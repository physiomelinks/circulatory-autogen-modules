"""
Unit tests of the boundary-condition sweep's machinery (cam_testing.checks): clipping sweep points
to valid ranges, classifying failures, and choosing recommended solver settings. No model is
generated; the simulation runs are faked.
"""
import math
import os
from types import SimpleNamespace

import numpy as np
import pytest

from cam_testing import checks, plots
from cam_testing.library import Parameter


def _param(name, units, value, kind='constant'):
    return Parameter(name, units, str(value), 'test', 'yes', kind)


def _component(tmp_path, params, spec):
    comp = SimpleNamespace(spec=spec, label='fake/v', dir=str(tmp_path), results_dir=str(tmp_path / 'results'),
                           plots_dir=str(tmp_path / 'plots'), id='fake_v')
    comp.parameters = lambda: params
    comp.boundary_conditions = lambda: [p for p in params if p.kind == 'boundary_condition']
    comp.config = {'module_format': 'cellml'}
    return comp


GATE_INVARIANT = {'expr': 'np.all((h >= -1e-9) & (h <= 1 + 1e-9))'}


# ---- clipping -----------------------------------------------------------------------------------

def test_unit_interval_states_needs_both_bounds():
    spec = {'invariants': [GATE_INVARIANT, 'np.all(x >= 0)', 'np.all(y <= 1)', 'np.all((z >= 0) & (z <= 1 + B))',
                           'np.all((w >= -0.5) & (w <= 1))', 'np.all((u >= 0) & (u <= 1.0))']}
    assert checks.unit_interval_states(spec) == {'h', 'u'}


def test_parameter_bounds_gate_heuristic_and_non_negative(tmp_path):
    params = [_param('h_init', 'dimensionless', 0.9), _param('n_init', 'dimensionless', 0.9),   # n not bounded
              _param('h_scale', 'dimensionless', 0.5), _param('g', 'microS', 1e-3),
              _param('V', 'milliV', -70), _param('C', 'millimolar', 0.1)]
    comp = _component(tmp_path, params, {'invariants': [GATE_INVARIANT]})
    b = checks.parameter_bounds(comp)
    assert b['h_init'][:2] == (0.0, 1.0)
    assert 'n_init' not in b and 'h_scale' not in b and 'V' not in b
    assert b['g'][:2] == (0.0, math.inf) and b['C'][:2] == (0.0, math.inf)


def test_explicit_bounds_win(tmp_path):
    params = [_param('h_init', 'dimensionless', 0.9), _param('frac', 'dimensionless', 0.2)]
    comp = _component(tmp_path, params, {'invariants': [GATE_INVARIANT],
                                         'bc_sweep': {'bounds': {'h_init': [0.1, 0.95], 'frac': [None, 1]}}})
    b = checks.parameter_bounds(comp)
    assert b['h_init'] == (0.1, 0.95, 'bc_sweep.bounds')
    assert b['frac'][:2] == (-math.inf, 1.0)


def test_point_skip_reason_bounds_and_constraints():
    bounds = {'h_init': (0.0, 1.0, 'gate')}
    env = {'B_init': 7.95, 'B_total': 8.0}
    assert checks.point_skip_reason('h_init', 0.99, bounds, [], env) is None
    assert 'outside valid range' in checks.point_skip_reason('h_init', 1.485, bounds, [], env)
    cons = ['B_init < B_total']
    assert checks.point_skip_reason('B_init', 7.95, bounds, cons, env) is None
    assert 'violates constraint' in checks.point_skip_reason('B_init', 11.9, bounds, cons, dict(env, B_init=11.9))
    assert 'violates constraint' in checks.point_skip_reason('B_total', 4.0, bounds, cons, dict(env, B_total=4.0))
    with pytest.raises(ValueError):
        checks.point_skip_reason('B_init', 1.0, bounds, ['nonexistent < 1'], env)


def test_invariant_excess_relative_to_tolerances():
    t = np.linspace(0, 1, 5)
    out = {'h': np.array([0.5, 0.9, 1.0 + 5e-9, 1.0, 0.9])}
    ex = checks.invariant_excess({'invariants': [GATE_INVARIANT]}, t, out, {}, rtol=1e-8, atol=1e-10)
    assert len(ex) == 1
    assert ex[0]['excess'] == pytest.approx(4e-9, rel=1e-3)
    assert ex[0]['relative'] == pytest.approx(4e-9 / (1e-10 + 1e-8 * (1 + 1e-9)), rel=1e-3)


# ---- candidate selection ------------------------------------------------------------------------

@pytest.mark.parametrize('n, first_pass', [(9, 0), (9, 4), (9, 8), (9, None), (1, 0), (1, None), (0, None)])
def test_cheapest_passing_binary_search(n, first_pass):
    calls = []

    def passes(i):
        calls.append(i)
        return first_pass is not None and i >= first_pass

    idx, seen = checks.cheapest_passing(n, passes)
    assert idx == first_pass
    assert len(calls) == len(set(calls)) <= max(1, math.ceil(math.log2(n + 1)))
    if idx is not None:
        assert seen[idx] is True


# ---- failure classification ---------------------------------------------------------------------

def _res(ok):
    return {'ok': ok, 'problems': [] if ok else ['h: violated'], 't': np.linspace(0, 1, 5),
            'out': {'h': np.full(5, 0.5)}, 'bad': [], 'seconds': 1e-3}


def test_classify_failures():
    failed = [{'pname': 'a', 'value': 2.0}, {'pname': 'b', 'value': 0.5}]
    ref = lambda pname, value: _res(pname == 'a')
    checks.classify_failures(failed, ref)
    assert [p['class'] for p in failed] == [checks.NUMERICAL_FAILURE, checks.PARAMETER_FAILURE]
    assert failed[1]['reference_problems'] == ['h: violated']


class _FakeCM:
    '''A ComponentModel stand-in: parameters a and b, one output h.'''
    def __init__(self, component):
        self.component, self.spec = component, component.spec
        self.nominal = {'a_mod': 1.0, 'b_mod': 1.0}
        self.helper = SimpleNamespace(model=SimpleNamespace(count_states=lambda: 1))
        self.model_path = 'unused'

    def outputs(self):
        return ['h']

    def units(self):
        return {'h': 'dimensionless'}

    def param_env(self, overrides=None):
        values = dict(self.nominal, **(overrides or {}))
        return {k[:-4]: v for k, v in values.items()}


def _fake_sweep_point(param_failure):
    '''Version settings (helper None) fail at a=2: a numerical failure that a MaximumStep fixes and the
    reference passes. With param_failure, b=0.5 fails with every setting.'''
    def run(cm, helper, pname, value):
        if param_failure and pname == 'b_mod' and value == 0.5:
            return _res(False)
        if pname == 'a_mod' and value == 2.0:
            return _res(helper is not None and bool(helper.get('MaximumStep')))
        return _res(True)
    return run


@pytest.mark.parametrize('param_failure', [False, True])
def test_sweep_classifies_and_recommends(tmp_path, monkeypatch, param_failure):
    params = [_param('a', 'per_s', 1.0), _param('b', 'per_s', 1.0)]
    spec = {'dt': 1e-3, 'sim_time': 1.0, 'solver': 'CVODE_myokit', 'solver_info': {'rtol': 1e-6, 'atol': 1e-8},
            'bc_sweep': {'factors': [0.5, 1.0, 2.0]}}
    cm = _FakeCM(_component(tmp_path, params, spec))
    monkeypatch.setattr(checks, '_sweep_point', _fake_sweep_point(param_failure))
    monkeypatch.setattr(checks, '_cvode_helper', lambda cm, cfg: dict(cfg))
    monkeypatch.setattr(plots, 'plot_sweep', lambda path, *a, **k: path)
    r = checks.verification_test_BC(cm)
    m = r.metrics
    assert m['n_numerical_failures'] == 1
    assert m['recommended'] is not None and m['recommended']['config'].get('MaximumStep')
    assert m['recommended']['failures'] == f'0/{m["n_runs"] - int(param_failure)}'
    # the recommendation is never looser than the version's settings
    assert m['recommended']['config']['rtol'] <= 1e-6 and m['recommended']['config']['atol'] <= 1e-8
    classes = {(f['parameter'], f['value']): f['class'] for f in m['failures']}
    assert classes[('a', 2.0)] == checks.NUMERICAL_FAILURE
    if param_failure:
        assert classes[('b', 0.5)] == checks.PARAMETER_FAILURE
        assert r.status == checks.FAILED and 'parameter failure' in r.message
    else:
        assert r.status == checks.PASSED and 'recommended settings' in r.message
    assert os.path.isfile(os.path.join(str(tmp_path), 'results', 'verification_test_BC.json'))


def test_sweep_skips_out_of_range(tmp_path, monkeypatch):
    params = [_param('h_init', 'dimensionless', 0.9), _param('b', 'per_s', 1.0)]
    spec = {'dt': 1e-3, 'sim_time': 1.0, 'solver_info': {'rtol': 1e-6, 'atol': 1e-8}, 'invariants': [GATE_INVARIANT],
            'bc_sweep': {'factors': [0.5, 1.0, 2.0], 'constraints': ['b > 0.6']}}
    cm = _FakeCM(_component(tmp_path, params, spec))
    cm.nominal = {'h_init_mod': 0.9, 'b_mod': 1.0}
    monkeypatch.setattr(checks, '_sweep_point', lambda cm, helper, pname, value: _res(True))
    monkeypatch.setattr(plots, 'plot_sweep', lambda path, *a, **k: path)
    r = checks.verification_test_BC(cm)
    assert r.status == checks.PASSED
    skipped = {(s['parameter'], s['value']) for s in r.metrics['skipped_out_of_range']}
    assert skipped == {('h_init', 1.8), ('b', 0.5)}
    assert r.metrics['n_runs'] == 4
