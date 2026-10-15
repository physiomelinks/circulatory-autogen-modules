"""
A version's tests driven by the protocol of an obs_data file (verification_config ``protocol``;
validation.<p>.baseline ``obs_data``): the heart valve heart/valve/pp driven by a prescribed inlet
pressure u_in(t) (a protocol_traces entry) that reverses the pressure difference across it, so it
opens and then closes. Results and plots go to a temporary directory.
"""
import json

import numpy as np
import pytest

from cam_testing import checks, harness
from cam_testing.library import load_version, protocol_problems

# u_out is 10666 J/m3 in the default parameterisation: u_in is 2000 above it, then ramps to 2000 below
TRACE_T = [0.0, 0.1, 0.2]
TRACE_U = [12666.0, 12666.0, 8666.0]
U_OUT = 10666.0


def _obs_data(path, protocol, data_items=()):
    path.write_text(json.dumps({'protocol_info': protocol, 'data_items': list(data_items), 'prediction_items': []}))
    return str(path)


REVERSAL = {'pre_times': [0.0], 'sim_times': [[0.2]], 'params_to_change': {'mod/u_in': [['dp']]},
            'protocol_traces': {'dp': {'t': TRACE_T, 'values': TRACE_U}}}


@pytest.fixture(scope='module')
def out_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp('protocol_results')
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(checks, 'result_path', lambda c, test, p=None: str(d / f'{getattr(p, "name", "v")}__{test}.json'))
        mp.setattr(checks, 'plot_path', lambda c, name, p=None: str(d / f'{name}.png'))
        yield d


@pytest.fixture(scope='module')
def cm(out_dir):
    version = load_version('valve', 'pp')
    version.spec = dict(
        version.spec, dt=0.005, protocol=_obs_data(out_dir / 'reversal_obs_data.json', REVERSAL),
        invariants=[
            {'expr': f'np.allclose(u_in, np.interp(t, {TRACE_T}, {TRACE_U}))', 'description': 'u_in follows the trace'},
            {'expr': 'np.max(zeta[t <= 0.1]) > 0.99', 'description': 'opens while u_in > u_out'},
            {'expr': 'np.all(np.diff(zeta[t >= 0.16]) < 0)', 'description': 'closes once u_in < u_out'},
        ],
        stability={'cvode': [{'rtol': 1e-6, 'atol': 1e-8}], 'solve_ivp': ['RK45'], 'fixed_step': {'schemes': ['rk4']},
                   'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]})
    return checks.ComponentModel(version, str(out_dir))


def test_the_protocol_sets_the_run_and_its_input_is_an_output(cm):
    assert cm.spec['sim_time'] == 0.2 and cm.spec['pre_time'] == 0.0
    assert cm.outputs()[-1] == 'u_in'
    t, out = cm.run()
    assert t[0] == 0.0 and t[-1] == pytest.approx(0.2)
    np.testing.assert_allclose(out['u_in'], np.interp(t, TRACE_T, TRACE_U))
    assert checks.run_test(cm).status == checks.PASSED


def test_invariants_see_the_trace_input(cm):
    r = checks.verification_test_invariants(cm)
    assert r.status == checks.PASSED, r.message


def test_the_driven_input_is_not_swept(cm):
    names, _ = checks.sweep_parameter_names(cm)
    assert 'u_in' not in names and 'u_out' in names


def test_fixed_step_tests_are_not_applicable(cm):
    r = checks.verification_test_timestep(cm)
    assert r.status == checks.NOT_APPLICABLE and 'protocol' in r.message
    r = checks.stability_test(cm)
    assert r.status == checks.PASSED, r.message
    stops = {f['family']: f['stop'] for f in r.metrics['families']}
    assert stops['solve_ivp (method=RK45)'] == stops['fixed_step_rk4'] == checks.PROTOCOL_CVODE_ONLY
    assert all(row['config']['solver'] == 'CVODE_myokit' for row in r.metrics['matrix'])


def test_sub_experiments_are_joined(cm, out_dir):
    # two sub-experiments with constant u_in: forward then reversed
    protocol = {'pre_times': [0.0], 'sim_times': [[0.1, 0.1]], 'params_to_change': {'mod/u_in': [[12666.0, 8666.0]]}}
    helper = harness.simulation_helper(cm.model_path, harness.with_protocol(harness.without_protocol(cm.spec), protocol))
    t, out, starts = harness.run_protocol(helper, ['zeta', 'u_in'], protocol)
    assert list(starts) == [0.0, 0.1] and len(t) == len(out['u_in']) == 41
    assert np.all(out['u_in'][t < 0.1] == 12666.0) and np.all(out['u_in'][t > 0.1] == 8666.0)
    assert out['zeta'][20] > 0.99 > 0.5 > out['zeta'][-1]


def test_structure_checks_the_protocol_files(tmp_path):
    _obs_data(tmp_path / 'one_obs_data.json', REVERSAL)
    _obs_data(tmp_path / 'two_obs_data.json', dict(REVERSAL, pre_times=[0.0, 0.0], sim_times=[[0.2], [0.2]]))
    assert protocol_problems(str(tmp_path), {'protocol': 'one_obs_data.json'}) == []
    assert protocol_problems(str(tmp_path), {'validation': {'p': {'baseline': {'obs_data': 'two_obs_data.json'}}}}) == []
    assert 'experiments' in protocol_problems(str(tmp_path), {'protocol': 'two_obs_data.json'})[0]
    assert 'not a file' in protocol_problems(str(tmp_path), {'validation': {'p': {'baseline': {'obs_data': 'x.json'}}}})[0]


def _baseline(cm, out_dir, data_items, protocol=REVERSAL):
    path = _obs_data(out_dir / 'baseline_obs_data.json', protocol, data_items)
    cm.component.spec['validation'] = {'default': {'baseline': {'status': 'active', 'source': 'synthetic',
                                                                 'obs_data': path, 'metric': 'nrmse', 'threshold': 0.01}}}
    return checks.validation_test_baseline(cm)


def _series(name, operands, values, obs_dt, operation=None):
    return {'data_item_name': name, 'data_type': 'series', 'operation': operation, 'operands': operands,
            'unit': 'J_per_m3', 'weight': 1.0, 'value': list(values), 'std': [1.0] * len(values), 'obs_dt': obs_dt,
            'experiment_idx': 0, 'subexperiment_idx': 0}


def _constant(name, operands, value, std, operation, sub=0):
    return {'data_item_name': name, 'data_type': 'constant', 'operation': operation, 'operands': operands,
            'unit': 'dimensionless', 'weight': 1.0, 'value': value, 'std': std, 'experiment_idx': 0,
            'subexperiment_idx': sub}


def test_baseline_scores_the_obs_data_against_the_protocol_run(cm, out_dir):
    ts = np.arange(0, 0.2001, 0.02)
    dp = np.interp(ts, TRACE_T, TRACE_U) - U_OUT
    items = [_series('dp', ['mod/u_in', 'mod/u_out'], dp, 0.02, operation='subtraction'),
             _constant('zeta_max', ['mod/zeta'], 1.0, 0.01, 'max')]
    r = _baseline(cm, out_dir, items)
    assert r.status == checks.PASSED, r.message
    assert r.metrics['nrmse']['dp'] < 1e-9
    assert [row['ok'] for row in r.metrics['targets']] == [True, True] and r.plots

    # data that does not match the run fails
    r = _baseline(cm, out_dir, [_series('dp', ['mod/u_in', 'mod/u_out'], -dp, 0.02, operation='subtraction'),
                                _constant('zeta_final', ['mod/zeta'], 1.0, 0.01, None)])
    assert r.status == checks.FAILED and '2 of 2 items missed' in r.message


def test_baseline_items_are_scored_on_their_sub_experiment(cm, out_dir):
    protocol = {'pre_times': [0.0], 'sim_times': [[0.1, 0.1]], 'params_to_change': {'mod/u_in': [[12666.0, 8666.0]]}}
    items = [_constant('open', ['mod/zeta'], 1.0, 0.01, None, sub=0), _constant('closed', ['mod/zeta'], 0.0, 0.2, None, sub=1),
             _constant('closing', ['mod/zeta'], 1.0, 0.01, 'max', sub=1)]
    r = _baseline(cm, out_dir, items, protocol)
    assert r.status == checks.PASSED, r.message
