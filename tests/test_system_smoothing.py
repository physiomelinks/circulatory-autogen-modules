"""
SN_simple (system_models/cellular/SN_simple), the heart-free test model for automatic smoothing. Every
non-smooth site is smoothed, with no exemption: 16 abs(tau/1000) guards (eps below each tau's smallest value),
the i_Ca_L GHK limit on |V| < 1 (abs(V) and its piecewise), the GHK form's division, whose denominator
1 - exp(2FV/RT) crosses 0 at V = 0 (guarded num den/(den^2 + eps_d^2)), the Kv1.5 and Kv4.2 time-constant
switches, and the varicosity's release rectifier max(Cai - Cai_NE_base, 0).

- The smoothed model (rel 1e-4, the division at rel 1e-6: at rel_div 1e-4 the guard moves spikes by 68 ms
  over 20 s) simulates and stays within tol of the original on the range-based metric over 2 s.
- The Myokit CVODES FSA gradient of the final NE with respect to g_M_mod and rho_Kv4_2 runs through V = 0
  (with the GHK branch blended but unguarded it stopped there) and matches central finite differences
  (relative step 1e-4, rtol 1e-10) within 2%.
- CasADi symbolic AD (slow, about 90 s, and needs casadi) matches them within 5%.
- The CVODE steps and wall times of both models are printed, not asserted.
About 2.5 minutes (generation and the reference runs take 40 s), and 4 with -m slow.
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


def _rel(rel):
    return {'abs': rel, 'maxmin': rel, 'piecewise': rel, 'division': 1e-6}


@pytest.fixture(scope='module')
def sn_simple(tmp_path_factory):
    from libcuflynx.utilities import smooth_cellml as sc
    work = tmp_path_factory.mktemp('sn_simple')
    records = module_array.read_records(os.path.join(SN, 'SN_simple_module_array.json'))
    model = smooth._generate(str(work), 'SN_simple', records, os.path.join(SN, 'SN_simple_parameters.json'), [])
    flat = checks.flat_model_path(model)
    probe = work / 'probe.cellml'
    probe.write_text(sc.probe(open(flat).read()))
    with contextlib.redirect_stdout(io.StringIO()):
        ranges = sc.reference_ranges(str(probe), [{'sim_time': 2.0, 'dt': 1e-4}])
    text, report = sc.smooth(open(flat).read(), ranges, rel=_rel(1e-3))
    smoothed = work / 'smoothed_rel1e-3.cellml'
    smoothed.write_text(text)
    return work, model, str(smoothed), report


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


def _fd(path, nominal):
    tight = _helper(path, {'rtol': 1e-10, 'atol': 1e-12})
    out = {}
    for p in PARAMS:
        step = 1e-4 * abs(nominal[p])
        out[p] = (_final_ne(tight, {**nominal, p: nominal[p] + step})
                  - _final_ne(tight, {**nominal, p: nominal[p] - step})) / (2 * step)
    return out


def test_smoothed_sn_simple_simulates_within_tol_and_reports_its_speed(sn_simple):
    work, model, _smoothed, _report = sn_simple
    status, message, metrics, _figs, _details = checks.compare_smoothed(
        model, SPEC, OUTPUTS, [{'label': 'nominal', 'params': {}}], {}, dict(checks.DEFAULT_SMOOTHING, rel=_rel(1e-4)),
        str(work / 'smoothing'), str(work), 'SN_simple', floor_solver_info={'rtol': 1e-10, 'atol': 1e-12})
    assert status == checks.PASSED, message
    assert metrics['counts'] == {'abs': 17, 'division': 1, 'piecewise': 4, 'maxmin': 1, 'division_skipped': 2}
    for name, path in (('original', checks.flat_model_path(model)), ('smoothed', metrics['smoothed_model'])):
        h = _helper(path, {}, sim_time=20.0, dt=1e-3)
        t0 = time.perf_counter()
        assert h.run()
        print(f'{name}: {h.simulation.last_number_of_steps()} CVODE steps, '
              f'{h.simulation.last_number_of_evaluations()} RHS evaluations, {time.perf_counter() - t0:.2f} s for 20 s')


def test_fsa_gradient_of_the_smoothed_model_matches_finite_differences(sn_simple):
    _work, _model, smoothed, report = sn_simple
    assert not [s for s in report['sites'] if s['status'] == 'exempt']
    h = _helper(smoothed, {})
    nominal = {p: float(v[0] if isinstance(v, list) else v) for p, v in zip(PARAMS, h.get_init_param_vals(PARAMS))}
    h.enable_fsa(['var_SN/NE'], PARAMS)
    h.reset_states()
    assert h.run(), 'the FSA run stopped (at V = 0, where an unguarded GHK form is 0/0)'
    fsa = {p: float(np.asarray(s).ravel()[-1]) for p, s in h.get_sensitivities(['var_SN/NE'], PARAMS)['var_SN/NE'].items()}
    fd = _fd(smoothed, nominal)
    for p in PARAMS:
        assert fsa[p] == pytest.approx(fd[p], rel=2e-2), p


@pytest.mark.slow
def test_casadi_gradient_of_the_smoothed_model_matches_finite_differences(sn_simple, tmp_path):
    ca = pytest.importorskip('casadi')
    from libcuflynx.generators.PythonGenerator import PythonGenerator
    from libcuflynx.solver_wrappers.casadi_python_solver_helper import SimulationHelper
    _work, _model, smoothed, _report = sn_simple
    with contextlib.redirect_stdout(io.StringIO()):
        py = PythonGenerator(smoothed, output_dir=str(tmp_path), module_name='SN_simple', casadi_compat=True).generate()
    h = SimulationHelper(py, 1e-3, 2.0, solver_info={'method': 'cvodes', 'rtol': 1e-8, 'atol': 1e-10,
                                                     'max_num_steps': 200000})
    h._create_param_subset(PARAMS)
    h.run()
    ne = h._extract_symb('var_SN/NE')
    f = ca.Function('f', [h.states_symb, h.variables_symb], [ca.jacobian(ne[-1], h.variables_symb_subset)])
    jac = np.array(f(ca.DM(h._x0_numeric()), ca.DM(h.variables))).ravel()
    nominal = dict(zip(PARAMS, h.variables_subset))
    fd = _fd(smoothed, nominal)
    for i, p in enumerate(PARAMS):
        assert jac[i] == pytest.approx(fd[p], rel=5e-2), p
