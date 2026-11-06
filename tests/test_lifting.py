"""
cam_testing.lift: a version's smoothed model lifted to quadratic form (libcuflynx.utilities.lift_cellml).
Each test generates a small model (seconds to tens of seconds).
"""
import json
import os

import pytest

pytest.importorskip('sympy')
pytest.importorskip('qbee')

from cam_testing import checks, lift  # noqa: E402
from cam_testing.library import load_version  # noqa: E402


def test_fitzhugh_nagumo_lifts_with_v_squared_and_passes_over_a_short_run(tmp_path):
    """Lift & Learn eq. 59-60: one auxiliary, V^2. Over 2 s the lifted run stays on the manifold."""
    result, doc = lift.run('FitzHugh_Nagumo/FitzHugh1961_v01', output_dir=str(tmp_path),
                           settings_override={'sim_time': 2.0, 'time_budget': 10.0})
    assert result.status == checks.PASSED, result.message
    assert [a['expr'] for a in doc['auxiliaries']] == ['V**2']
    assert [q.split('.')[-1] for q in doc['states']] == ['V', 'R']
    saved = json.load(open(os.path.join(tmp_path, 'lifting_result.json')))
    assert saved['status'] == checks.PASSED and saved['metrics']['size']['n_w'] == 3


def test_a_long_run_that_leaves_the_manifold_is_reported_with_its_growth_rate(tmp_path):
    """w = V^2 is invariant but does not attract: over 20 s solver error grows off it."""
    result, _ = lift.run('FitzHugh_Nagumo/FitzHugh1961_v01', output_dir=str(tmp_path),
                         settings_override={'sim_time': 20.0, 'time_budget': 10.0})
    assert result.status == lift.UNSTABLE_MANIFOLD, result.message
    assert result.metrics['verification']['defect_growth']['rate'] > 0


def test_a_smoothed_site_of_an_input_alone_is_an_input_function(tmp_path):
    """cardiac_clock/controlled: ds/dt = 1/T_eff, T_eff = max(T + Delta_T, T_min) smoothed. Delta_T is an unconnected
    boundary condition (a held input), so the smoothed 1/T_eff is a function of the inputs: no auxiliary."""
    result, doc = lift.run('cardiac_clock/controlled', output_dir=str(tmp_path),
                           settings_override={'sim_time': 2.0, 'time_budget': 10.0})
    assert result.status == checks.PASSED, result.message
    kinds = {i['name']: i['kind'] for i in doc['inputs']}
    assert kinds.pop('Delta_T_mod') == 'variable'
    assert list(kinds.values()) == ['input_function']
    assert 'eps_sm' in [i for i in doc['inputs'] if i['kind'] == 'input_function'][0]['expr']
    assert doc['size']['n_aux'] == 0 and doc['operators']['present']['B']
    assert result.metrics['smoothing']['status'] == checks.PASSED


def test_held_inputs_make_the_boundary_condition_a_state(tmp_path):
    """cardiac_clock/controlled with held inputs: Delta_T is a state with dDelta_T/dt = 0, and the smoothed 1/T_eff,
    a function of it alone, one held variable: no inputs remain."""
    result, doc = lift.run('cardiac_clock/controlled', output_dir=str(tmp_path),
                           settings_override={'sim_time': 2.0, 'time_budget': 10.0, 'held_inputs': True})
    assert result.status == checks.PASSED, result.message
    assert doc['inputs'] == [] and doc['size']['n_held'] >= 1
    assert all(h['index'] >= doc['size']['n'] for h in doc['held_inputs'])


def test_a_lifting_that_cannot_finish_is_a_failed_result(tmp_path):
    result, doc = lift.run('FitzHugh_Nagumo/FitzHugh1961_v01', output_dir=str(tmp_path),
                           settings_override={'time_budget': 0.01, 'max_aux': 0})
    assert doc is None and result.status == checks.FAILED
    assert 'LiftError' in result.message


def test_spec_settings_override_the_defaults():
    s = lift.lifting_settings({'lifting': {'tol': 1e-3, 'rationale': 'x'}}, time_budget=5.0)
    assert s['tol'] == 1e-3 and s['time_budget'] == 5.0 and 'rationale' not in s
    assert s['max_aux'] == lift.DEFAULT_LIFTING['max_aux']


def test_unconnected_boundary_conditions_are_the_default_inputs(tmp_path):
    """A version alone: its boundary_condition variables are constants of the generated model."""
    version = load_version('FitzHugh_Nagumo', 'FitzHugh1961_v01')
    assert lift.default_inputs(version, checks.flat_model_path(checks.ComponentModel(version, str(tmp_path)).model_path)) == []
    clock = load_version('cardiac_clock', 'controlled')
    assert lift.default_inputs(clock, checks.flat_model_path(checks.ComponentModel(clock, str(tmp_path / 'c')).model_path)) \
        == ['mod_module/Delta_T']
