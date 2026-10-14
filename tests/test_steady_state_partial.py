"""
Unit tests of the periodic steady-state check, the time-step test's partial verdict, the report's
supermodule columns for a component and its status counts, and version listing. No model is
generated or simulated.
"""
from types import SimpleNamespace

import numpy as np

from cam_testing import checks, report
from cam_testing.library import _has_source


# ---- time step: failed / partial / passed -------------------------------------------------------

def _verdict(diffs, observed, cvode_diff=1e-4):
    return checks.timestep_problems(diffs, observed, cvode_diff, min_order=1.8, tol=1e-3, cvode_tol=1e-3,
                                    roundoff=1e-12)


def test_timestep_second_order_passes():
    assert _verdict([1e-2, 2.5e-3, 6.25e-4], [None, 2.0, 2.0]) == ([], [])


def test_timestep_first_order_but_convergent_is_partial():
    problems, order = _verdict([1e-2, 5e-3, 2.5e-3], [None, 1.0, 1.0])
    assert problems == []
    assert any('observed order 1.00' in p for p in order) and any('finest-step' in p for p in order)


def test_timestep_hard_problems_fail():
    assert _verdict([1e-3, 2e-3, 1e-3], [None, -1.0, 1.0])[0]                       # not decreasing
    assert _verdict([1e-2, 2.5e-3, 6.25e-4], [None, 2.0, 2.0], cvode_diff=1e-2)[0]   # disagrees with CVODE
    assert _verdict([np.nan, 1e-3], [None, None])[0]                                # non-finite


# ---- periodic steady state ----------------------------------------------------------------------

def test_period_changes_zero_at_steady_state_and_large_in_a_transient():
    t = np.linspace(0.0, 2.0, 2001)
    periodic = {'u': np.sin(2 * np.pi * t)}
    growing = {'u': (1 - np.exp(-t)) * np.sin(2 * np.pi * t)}
    assert checks.period_changes(t, periodic, 1.0)['u'] < 1e-3
    assert checks.period_changes(t, growing, 1.0)['u'] > 0.1


# ---- version listing ----------------------------------------------------------------------------

def test_a_directory_left_with_only_generated_output_is_not_a_version(tmp_path):
    left = tmp_path / 'vv'
    (left / 'results').mkdir(parents=True)
    (left / 'plots').mkdir()
    (left / 'chamber_vv.html').write_text('')
    assert not _has_source(left)
    (left / 'chamber_vv_modules_config.json').write_text('[]')
    assert _has_source(left)


# ---- a component's supermodule columns ----------------------------------------------------------

def _super(key, equivalent=True):
    name = key.split('/')[1]
    return SimpleNamespace(key=key, vessel_type=key.split('/')[0], name=name,
                           spec={'supermodule': {'equivalent': [{'model': 'm', 'reproduces': 'r'}]} if equivalent else {}})


def _with(monkeypatch, users, results):
    component = SimpleNamespace(vessel_type='chamber', name='v1')
    monkeypatch.setattr(report, '_supermodule_index', lambda: {('chamber', 'v1'): users} if users else {})
    monkeypatch.setattr(checks, 'load', lambda s, test, *a: results.get(s.key))
    monkeypatch.setattr(report, '_equivalence_result', lambda s: results.get(s.key))
    return component


def test_component_shows_its_supermodules_structure(monkeypatch):
    ok = checks.Result('supermodule_structure_test', checks.PASSED, 'fine')
    bad = checks.Result('supermodule_structure_test', checks.FAILED, 'broken')
    c = _with(monkeypatch, [_super('heart/a')], {'heart/a': ok})
    assert report._in_supermodules(c, 'supermodule_structure_test').status == checks.PASSED_IN_SUPER
    c = _with(monkeypatch, [_super('heart/a'), _super('heart/b')], {'heart/a': ok, 'heart/b': bad})
    r = report._in_supermodules(c, 'supermodule_structure_test')
    assert r.status == checks.FAILED_IN_SUPER and [x['key'] for x in r.metrics['supermodules']] == ['heart/a', 'heart/b']


def test_component_not_in_a_supermodule_or_not_run(monkeypatch):
    assert report._in_supermodules(_with(monkeypatch, [], {}), 'supermodule_structure_test') is None
    assert report._in_supermodules(_with(monkeypatch, [_super('heart/a')], {}), 'supermodule_structure_test') is None
    # reproduces: only supermodules with a system model to reproduce count
    r = report._in_supermodules(_with(monkeypatch, [_super('heart/a', equivalent=False)], {}),
                                'supermodule_equivalence_test')
    assert r.status == checks.NOT_APPLICABLE


# ---- status counts ------------------------------------------------------------------------------

def test_counts_include_partial_and_each_parameterisations_archive():
    comp = {'tests': [{'status': checks.PASSED}, {'status': checks.PARTIAL}, {'status': checks.PASSED_IN_SUPER}],
            'parameterisations': [
                {'tests': [{'status': checks.FAILED}], 'omex_status': checks.FAILED},
                {'tests': [{'status': checks.NOT_APPLICABLE}], 'omex_status': checks.PASSED},
                {'tests': [{'status': checks.NOT_APPLICABLE}], 'omex_status': None}]}
    counts = report._status_counts([comp])
    assert counts['passed'] == 3 and counts['partial'] == 1 and counts['failed'] == 2
    assert counts['not_applicable'] == 2
