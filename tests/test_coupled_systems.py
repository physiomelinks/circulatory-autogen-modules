"""
The coupled example system models (system_models/coupled, built by
tools/build_coupled_examples.py): 0D modules generated as C++ and coupled to the FEniCS model
tissue_diffusion_FEniCS (an external model: libcuflynx.coupling), checked against the same
problem with the tissue as a finite-volume grid of CellML cells, and timed.

* SN_NE_single_volume reproduces SN_NE_original: the varicosity NEexchange_v01 with one
  well-mixed volume of its own size is the original varicosity (CellML, myokit).
* FEniCS with fv_scheme 1 (DG0, two-point face fluxes) solves the same discrete equations as
  the CellML grid, so the coupled C++ + FEniCS run must match the all-CellML C++ run up to time
  stepping and coupling error. With Q1 finite elements (fv_scheme 0, refined) it differs by the
  spatial discretisation, which is reported.
* The run times of both are recorded in results/coupled_comparison.json of each FEniCS model.

Needs libcuflynx with libcuflynx.coupling (circulatory_autogen PR #520), CMake, a C++ compiler
and SUNDIALS; the FEniCS tests also need dolfinx (skipped without it). In a conda environment,
install cmake, cxx-compiler and sundials there too (the 0D model is loaded into its Python).
"""
import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pytest

from cam_testing.library import MODULES_DIR
from cam_testing.system import generate, load_system, simulate

pytestmark = pytest.mark.system_model

coupling = pytest.importorskip('libcuflynx.coupling', reason='needs libcuflynx with libcuflynx.coupling (#520)')


def _need_build_tools():
    if shutil.which('cmake') is None:
        pytest.skip('cmake not available')


def _need_dolfinx():
    pytest.importorskip('dolfinx')


def generate_cpp(name, work):
    """The system model generated as C++ (CVODE at its spec's tolerances); (folder, seconds)."""
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture
    s = load_system(name)
    spec = s.spec
    cfg = {'file_prefix': name, 'input_param_file': f'{name}_parameters.csv', 'model_type': 'cpp',
           'resources_dir': s.dir, 'generated_models_dir': str(work), 'DEBUG': False, 'couple_to_1d': False,
           'module_library_dirs': [MODULES_DIR], 'use_builtin_modules': False,
           'pre_time': spec.get('pre_time', 0.0), 'sim_time': spec['sim_time'], 'dt': spec['dt'],
           'solver_info': {'solver': 'CVODE', 'dt_solver': spec['dt'], 'MaximumNumberOfSteps': 100000,
                           'rtol': spec['solver_info']['rtol'], 'atol': spec['solver_info']['atol']}}
    t = time.perf_counter()
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        ok = generate_with_new_architecture(False, cfg)
    assert ok, log.getvalue()[-4000:]
    return os.path.join(str(work), name), time.perf_counter() - t


def _load(path):
    header = open(path).readline().lstrip('#').split(';')
    cols = {h.split(':', 1)[1].split('[')[0].strip(): i for i, h in enumerate(header) if ':' in h}
    return cols, np.loadtxt(path, ndmin=2)


def run_fv(name, work):
    """The all-CellML model as C++ through main0d: (t, {name: values}, timings)."""
    _need_build_tools()
    model_dir, t_gen = generate_cpp(name, work)
    build = os.path.join(model_dir, 'build')
    args = ['cmake', '-S', model_dir, '-B', build, '-DCMAKE_BUILD_TYPE=Release']
    if os.environ.get('SUNDIALS_DIR'):
        args.append(f"-DSUNDIALS_DIR={os.environ['SUNDIALS_DIR']}")
    if os.environ.get('CONDA_PREFIX'):
        args.append(f"-DCMAKE_PREFIX_PATH={os.environ['CONDA_PREFIX']}")
    t = time.perf_counter()
    for cmd in (args, ['cmake', '--build', build, '-j']):
        p = subprocess.run(cmd, capture_output=True, text=True)
        assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    t_build = time.perf_counter() - t
    out = os.path.join(str(work), 'out')
    t = time.perf_counter()
    p = subprocess.run([os.path.join(build, 'main0d'), '-outDir', out], capture_output=True, text=True)
    t_run = time.perf_counter() - t
    assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    cs, s = _load(os.path.join(out, 'sol0D_states.txt'))
    cv, v = _load(os.path.join(out, 'sol0D_variables.txt'))
    values = {k: s[:, i] for k, i in cs.items()}
    values.update({k: v[:, i] for k, i in cv.items()})
    return s[:, 0], values, {'generate': t_gen, 'build': t_build, 'run': t_run}


def run_fenics(name, work, settings=None, **params):
    """The FEniCS-coupled model (C++ 0D + FEniCS), with tissue_diffusion_FEniCS parameters
    overridden by ``params`` and coupling settings (e.g. subiterations) by ``settings``:
    (CoupledResult, timings)."""
    _need_build_tools()
    _need_dolfinx()
    model_dir, t_gen = generate_cpp(name, work)
    info_path = os.path.join(model_dir, 'external_models.json')
    info = json.load(open(info_path))
    info['external_models'][0]['parameters'].update(params)
    info['external_models'][0].update(settings or {})
    json.dump(info, open(info_path, 'w'), indent=2)
    result = coupling.run_coupled(model_dir, output_dir=os.path.join(str(work), 'out'), verbose=False)
    timings = dict(result.timings)
    timings['generate'] = t_gen
    return result, timings


def coupled_cells(name):
    """The grid cells of a CellML-grid system model that 0D modules exchange with, in the order
    of their sources (the order the FEniCS model numbers its regions)."""
    records = json.load(open(os.path.join(load_system(name).dir, f'{name}_vessel_array.json')))
    by_name = {r['name']: r for r in records}
    cells = {}
    for r in records:
        if r['module_type'] == 'tissue_diffusion_volume':
            for src in r['inp_instances']:
                if by_name[src]['module_type'] != 'tissue_diffusion_face':
                    cells[src] = r['name']
    order = [r['name'] for r in records if r['name'] in cells]
    return [cells[src] for src in order]


def _rel_diff(t_ref, ref, t, values):
    """Largest |values - ref| (values interpolated onto t_ref), relative to max |ref|."""
    return float(np.max(np.abs(np.interp(t_ref, t, values) - ref)) / np.max(np.abs(ref)))


def _record(name, entry):
    """Adds a comparison to system_models/coupled/<name>/results/coupled_comparison.json."""
    s = load_system(name)
    os.makedirs(s.results_dir, exist_ok=True)
    path = os.path.join(s.results_dir, 'coupled_comparison.json')
    data = json.load(open(path)) if os.path.isfile(path) else {}
    data.update(entry)
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)


# ---------------------------------------------------------------------------------------------
# NE around a varicosity
# ---------------------------------------------------------------------------------------------

def test_varicosity_NEexchange_with_one_volume_is_the_original(tmp_path):
    """NEexchange_v01 releasing into one well-mixed volume of the varicosity's size (CellML,
    myokit) reproduces the original varicosity, whose NE is an internal state."""
    def run(name, names):
        s = load_system(name)
        path = generate(s, str(tmp_path / name))
        return simulate(path, s.spec, s.spec['solver_info'], names=names)
    t0, a = run('SN_NE_original', ['var_SN/NE', 'var_SN/Cai'])
    t1, b = run('SN_NE_single_volume', ['ecs_0/C_P', 'var_SN/Cai'])
    assert np.max(a['var_SN/NE']) > 1e-4, 'no NE released: the stimulus does not fire the neuron'
    assert _rel_diff(t0, a['var_SN/NE'], t1, b['ecs_0/C_P']) < 1e-5
    assert _rel_diff(t0, a['var_SN/Cai'], t1, b['var_SN/Cai']) < 1e-5


@pytest.fixture(scope='module')
def ne_fv(tmp_path_factory):
    return run_fv('SN_NE_FV', tmp_path_factory.mktemp('ne_fv'))


@pytest.mark.slow
@pytest.mark.parametrize('scheme', ['fv', 'fe'])
def test_NE_FEniCS_matches_the_cellml_grid(scheme, ne_fv, tmp_path):
    t_ref, ref_values, ref_time = ne_fv
    (centre,) = coupled_cells('SN_NE_FV')       # the cell the varicosity exchanges with
    ref = ref_values[f'{centre}/C_P']
    # NE is released in ~1 ms pulses, a few coupling steps each: the DG0 comparison iterates
    # each step (second-order coupling) so what is left is the time stepping
    if scheme == 'fv':
        result, times = run_fenics('SN_NE_FEniCS', tmp_path, {'subiterations': 3, 'tol': 1e-10}, fv_scheme=1.0)
    else:
        result, times = run_fenics('SN_NE_FEniCS', tmp_path, fv_scheme=0.0, refine=2.0)
    ne = result.exchange['tissue/C_t'][:, 0]
    diff = _rel_diff(t_ref, ref, result.times, ne)
    peak_ratio = float(ne.max() / ref.max())
    _record('SN_NE_FEniCS', {f'NE_at_varicosity_{scheme}': {
        'max_rel_difference_to_cellml_grid': diff, 'peak_ratio_fenics_to_cellml_grid': peak_ratio,
        'peak_mM': {'cellml_grid': float(ref.max()), 'fenics': float(ne.max())},
        'seconds': {'cellml_grid_cpp': ref_time, 'fenics_coupled': times}}})
    print(f'\nNE ({scheme}): max difference {diff:.2%} of the peak; CellML grid C++ run {ref_time["run"]:.2f} s '
          f'(+ build {ref_time["build"]:.1f} s, generate {ref_time["generate"]:.1f} s); coupled run '
          f'{times["total"]:.2f} s (0D {times["zero_d"]:.2f} s, FEniCS {times["external"]:.2f} s)')
    assert np.max(ref) > 1e-6, 'no NE reached the extracellular space'
    if scheme == 'fv':
        # the same discrete equations: equal up to time stepping
        assert diff < 0.02 and abs(peak_ratio - 1.0) < 0.01, (diff, peak_ratio)
    else:
        # Q1 elements resolve the gradient around the point-like release that the 1 um CellML
        # cells can't (refining the elements moves the peak by ~2 %, the CellML grid is ~25 %
        # higher): the same magnitude, not the same numbers
        assert 0.6 < peak_ratio < 1.4, peak_ratio


# ---------------------------------------------------------------------------------------------
# tissue O2 around two capillaries
# ---------------------------------------------------------------------------------------------

@pytest.fixture(scope='module')
def o2_fv(tmp_path_factory):
    return run_fv('microvasc_O2_FV', tmp_path_factory.mktemp('o2_fv'))


@pytest.mark.slow
@pytest.mark.parametrize('scheme', ['fv', 'fe'])
def test_O2_FEniCS_matches_the_cellml_grid(scheme, o2_fv, tmp_path):
    t_ref, ref_values, ref_time = o2_fv
    cells = coupled_cells('microvasc_O2_FV')    # the cells the two capillaries exchange with
    params = {'fv_scheme': 1.0} if scheme == 'fv' else {'fv_scheme': 0.0, 'refine': 2.0}
    result, times = run_fenics('microvasc_O2_FEniCS', tmp_path, **params)
    C = result.exchange['tissue/C_t']
    diffs = [_rel_diff(t_ref, ref_values[f'{c}/C_P'], result.times, C[:, k]) for k, c in enumerate(cells)]
    _record('microvasc_O2_FEniCS', {f'tissue_O2_{scheme}': {
        'max_rel_difference_to_cellml_grid': max(diffs),
        'end_mM': {'cellml_grid': [float(ref_values[f'{c}/C_P'][-1]) for c in cells], 'fenics': C[-1].tolist()},
        'seconds': {'cellml_grid_cpp': ref_time, 'fenics_coupled': times}}})
    print(f'\nO2 ({scheme}): max difference {max(diffs):.2%}; CellML grid C++ run {ref_time["run"]:.2f} s '
          f'(+ build {ref_time["build"]:.1f} s, generate {ref_time["generate"]:.1f} s); coupled run '
          f'{times["total"]:.2f} s (0D {times["zero_d"]:.2f} s, FEniCS {times["external"]:.2f} s)')
    assert max(diffs) < (0.01 if scheme == 'fv' else 0.10), diffs


@pytest.mark.slow
def test_FEniCS_coupling_under_mpi_matches_serial(tmp_path):
    """The FEniCS mesh distributed over two ranks gives the same coupled run as one rank."""
    _need_build_tools()
    _need_dolfinx()
    if shutil.which('mpiexec') is None:
        pytest.skip('mpiexec not available')
    model_dir, _ = generate_cpp('microvasc_O2_FEniCS', tmp_path / 'gen')
    serial = coupling.run_coupled(model_dir, output_dir=str(tmp_path / 'out1'), verbose=False)
    script = tmp_path / 'run.py'
    script.write_text('import json, sys\nfrom libcuflynx.coupling import run_coupled\n'
                      f'r = run_coupled({model_dir!r}, output_dir={str(tmp_path / "out2")!r}, verbose=False)\n'
                      'from mpi4py import MPI\n'
                      'if MPI.COMM_WORLD.rank == 0:\n'
                      f'    json.dump(r.exchange["tissue/C_t"].tolist(), open({str(tmp_path / "mpi.json")!r}, "w"))\n')
    p = subprocess.run(['mpiexec', '-n', '2', sys.executable, str(script)], capture_output=True, text=True,
                       env=dict(os.environ))
    assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    parallel = np.array(json.load(open(tmp_path / 'mpi.json')))
    assert np.max(np.abs(parallel - serial.exchange['tissue/C_t'])) <= 1e-10 * np.max(np.abs(parallel))
