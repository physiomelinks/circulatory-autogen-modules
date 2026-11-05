"""
SN_simple (system_models/cellular/SN_simple), the heart-free test model for automatic smoothing: its
non-smooth sites are 16 abs(tau/1000) guards, the i_Ca_L GHK limit on |V| < 1, the Kv1.5 and Kv4.2
time-constant switches and the varicosity's release rectifier max(Cai - Cai_NE_base, 0).

- The smoothed model (rel 1e-4) simulates and stays within tol of the original on the range-based
  metric (checks.compare_smoothed) over 2 s.
- With the GHK limit left as written (its branches agree at |V| = 1: blending them evaluates the GHK
  form at V = 0, 0/0, and CVODES's sensitivities fail there), the Myokit CVODES FSA gradient of the final
  NE with respect to g_M_mod and rho_Kv4_2 matches central finite differences (relative step 1e-4, rtol
  1e-10) within 2%. (The original's FSA gradient misses the jump terms of the time-constant switches at
  each crossing: for these two parameters it has the wrong sign; not asserted here.)
- The CVODE steps and wall times of both models are printed, not asserted.
About 45 s.
"""
import contextlib
import io
import os
import time

import numpy as np
import pytest

from cam_testing import checks, module_array, smooth

SN = os.path.join(os.path.dirname(__file__), '..', 'system_models', 'cellular', 'SN_simple')
SPEC = {'sim_time': 2.0, 'pre_time': 0.0, 'dt': 1e-4, 'solver': 'CVODE_myokit', 'solver_info': {}}
OUTPUTS = ['soma_SN/V', 'soma_SN/Cai', 'var_SN/Cai', 'var_SN/NE']
PARAMS = ['parameters/g_M_mod_soma_SN', 'parameters/rho_Kv4_2_soma_SN']
GHK = [{'variable': 'soma_SN_module/i_Ca_L', 'reason': 'GHK limit: a removable-singularity guard'}]


@pytest.fixture(scope='module')
def sn_simple(tmp_path_factory):
    work = tmp_path_factory.mktemp('sn_simple')
    records = module_array.read_records(os.path.join(SN, 'SN_simple_module_array.json'))
    model = smooth._generate(str(work), 'SN_simple', records, os.path.join(SN, 'SN_simple_parameters.json'), [])
    return work, model


def _helper(path, solver_info, sim_time=2.0, dt=1e-3):
    from libcuflynx.solver_wrappers import get_simulation_helper
    with contextlib.redirect_stdout(io.StringIO()):
        return get_simulation_helper(model_path=path, solver='CVODE_myokit', model_type='cellml', dt=dt,
                                     sim_time=sim_time, solver_info=solver_info, pre_time=0.0)


def _final_ne(h, values):
    h.reset_states()
    h.set_param_vals(list(values), list(values.values()))
    assert h.run()
    return float(np.asarray(h.get_results([['var_SN/NE']], flatten=True)[0]).ravel()[-1])


def _fsa_and_fd(path):
    """{param: (FSA d(final NE)/dp at the default tolerance, central FD at rtol 1e-10)}."""
    h = _helper(path, {})
    nominal = {p: float(v[0] if isinstance(v, list) else v) for p, v in zip(PARAMS, h.get_init_param_vals(PARAMS))}
    h.enable_fsa(['var_SN/NE'], PARAMS)
    h.reset_states()
    ok = h.run()
    fsa = {p: float(np.asarray(s).ravel()[-1]) for p, s in h.get_sensitivities(['var_SN/NE'], PARAMS)['var_SN/NE'].items()} \
        if ok else {}
    tight = _helper(path, {'rtol': 1e-10, 'atol': 1e-12})
    out = {}
    for p in PARAMS:
        step = 1e-4 * abs(nominal[p])
        fd = (_final_ne(tight, {**nominal, p: nominal[p] + step}) - _final_ne(tight, {**nominal, p: nominal[p] - step})) / (2 * step)
        out[p] = (fsa.get(p, np.nan), fd)
    return out


def test_smoothed_sn_simple_simulates_within_tol_and_reports_its_speed(sn_simple):
    work, model = sn_simple
    status, message, metrics, _figs, _details = checks.compare_smoothed(
        model, SPEC, OUTPUTS, [{'label': 'nominal', 'params': {}}], {}, dict(checks.DEFAULT_SMOOTHING, rel=1e-4),
        str(work / 'smoothing'), str(work), 'SN_simple', floor_solver_info={'rtol': 1e-10, 'atol': 1e-12})
    assert status == checks.PASSED, message
    assert metrics['counts'] == {'abs': 17, 'piecewise': 4, 'maxmin': 1}
    for name, path in (('original', checks.flat_model_path(model)), ('smoothed', metrics['smoothed_model'])):
        h = _helper(path, {}, sim_time=20.0, dt=1e-3)
        t0 = time.perf_counter()
        assert h.run()
        print(f'{name}: {h.simulation.last_number_of_steps()} CVODE steps, '
              f'{h.simulation.last_number_of_evaluations()} RHS evaluations, {time.perf_counter() - t0:.2f} s for 20 s')


def test_fsa_gradient_of_the_smoothed_model_matches_finite_differences(sn_simple):
    from libcuflynx.utilities import smooth_cellml as sc
    work, model = sn_simple
    flat = checks.flat_model_path(model)
    probe = work / 'probe.cellml'
    probe.write_text(sc.probe(open(flat).read()))
    with contextlib.redirect_stdout(io.StringIO()):
        ranges = sc.reference_ranges(str(probe), [{'sim_time': 2.0, 'dt': 1e-4}])
    text, report = sc.smooth(open(flat).read(), ranges, rel=1e-3, exempt=GHK)
    assert report['counts'] == {'abs': 16, 'abs_exempt': 1, 'piecewise_exempt': 1, 'piecewise': 3, 'maxmin': 1}
    path = work / 'smoothed_ghk_exempt.cellml'
    path.write_text(text)
    for p, (fsa, fd) in _fsa_and_fd(str(path)).items():
        assert fsa == pytest.approx(fd, rel=2e-2), p
