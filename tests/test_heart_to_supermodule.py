"""
tools/heart_to_supermodule.py replaces a system model's written-out heart (cardiac_clock, chamber and
valve records) with one heart supermodule record: the host neighbours go to per_submodule_inputs /
per_submodule_outputs, a flow_merge into a chamber and a flow_split out of one are folded into the
chamber's summing ports, and the submodules' names in params_for_id, obs_data and the spec become
heart/<submodule>/... . supermodule.equivalence compares wrapped and reciprocal outputs as a system
model's equivalence does.
"""
import importlib.util
import os
from types import SimpleNamespace

import numpy as np
import pytest

from cam_testing import supermodule

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location('heart_to_supermodule',
                                               os.path.join(ROOT, 'tools', 'heart_to_supermodule.py'))
h2s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h2s)

HEART = h2s.version_array(os.path.join(ROOT, 'modules'), 'Argus2026_v01')


def _r(name, mt, sub, inp=(), out=()):
    return {'name': name, 'module_type': mt, 'module_subtype': sub, 'parameterisation': 'default',
            'inp_instances': list(inp), 'out_instances': list(out)}


def _written_out(ra_inputs, la_outputs, chamber='Liang2009NoVisco_v01'):
    return [
        _r('heart_clock', 'cardiac_clock', 'Liang2009_v01', out=['ra', 'rv', 'la', 'lv']),
        _r('ra', 'chamber', chamber, inp=ra_inputs + ['heart_clock'], out=['trv', 'volume_sum']),
        _r('rv', 'chamber', chamber, inp=['trv', 'heart_clock'], out=['puv']),
        _r('la', 'chamber', chamber, inp=['pvn', 'heart_clock'], out=la_outputs),
        _r('lv', 'chamber', chamber, inp=['miv', 'heart_clock'], out=['aov']),
        _r('trv', 'valve', 'pp', inp=['ra'], out=['rv']),
        _r('puv', 'valve', 'pp', inp=['rv'], out=['par']),
        _r('miv', 'valve', 'pp', inp=[m for m in la_outputs if m == 'la_outflow'] or ['la'], out=['lv']),
        _r('aov', 'valve', 'pp', inp=['lv'], out=['aorta']),
    ]


def test_both_venae_cavae_go_straight_into_ra():
    records = [_r('ivc', 'venous', 'vp', out=['ra_inflow']), _r('svc', 'venous', 'vp', out=['ra_inflow']),
               _r('ra_inflow', 'flow_merge', 'vp', inp=['ivc', 'svc'], out=['ra']),
               _r('pvn', 'venous', 'vp', out=['la']), _r('par', 'arterial_simple', 'vp', inp=['puv']),
               _r('aorta', 'arterial', 'vv', inp=['aov']), _r('volume_sum', 'volume_sum', 'v1', inp=['ra', 'svc'])]
    records += _written_out(['ra_inflow'], ['miv'])
    out, local, folded = h2s.convert_array(records, HEART, 'Argus2026_v01')
    assert folded == {'ra_inflow'}
    assert local['heart_clock'] == 'clock' and local['ra'] == 'ra'
    heart = next(r for r in out if r['name'] == 'heart')
    assert heart['per_submodule_inputs'] == {'ra': ['ivc', 'svc'], 'la': ['pvn']}
    assert heart['per_submodule_outputs'] == {'ra': ['volume_sum'], 'puv': ['par'], 'aov': ['aorta']}
    by_name = {r['name']: r for r in out}
    assert set(by_name) == {'ivc', 'svc', 'pvn', 'par', 'aorta', 'volume_sum', 'heart'}
    assert by_name['ivc']['out_instances'] == ['heart'] and by_name['svc']['out_instances'] == ['heart']
    assert by_name['volume_sum']['inp_instances'] == ['heart', 'svc']


def test_an_asd_shunt_leaves_la_and_enters_ra_with_no_junction():
    records = [_r('svc', 'venous', 'vp', out=['ra_inflow']),
               _r('ra_inflow', 'flow_merge', 'vp', inp=['svc', 'asd'], out=['ra']),
               _r('la_outflow', 'flow_split', 'pv', inp=['la'], out=['miv', 'asd']),
               _r('asd', 'resistor', 'pp', inp=['la_outflow'], out=['ra_inflow']),
               _r('pvn', 'venous', 'vp', out=['la']), _r('par', 'arterial_simple', 'vp', inp=['puv']),
               _r('aorta', 'arterial', 'vv', inp=['aov'])]
    records += _written_out(['ra_inflow'], ['la_outflow'])
    out, _, folded = h2s.convert_array(records, HEART, 'Argus2026_v01')
    assert folded == {'ra_inflow', 'la_outflow'}
    heart = next(r for r in out if r['name'] == 'heart')
    assert heart['per_submodule_inputs']['ra'] == ['svc', 'asd']
    assert heart['per_submodule_outputs']['la'] == ['asd']
    assert 'miv' not in heart['per_submodule_inputs']
    asd = next(r for r in out if r['name'] == 'asd')
    assert asd['inp_instances'] == ['heart'] and asd['out_instances'] == ['heart']


def test_a_heart_wired_otherwise_is_refused():
    records = [_r('svc', 'venous', 'vp', out=['ra']), _r('pvn', 'venous', 'vp', out=['la']),
               _r('par', 'arterial_simple', 'vp', inp=['puv']), _r('aorta', 'arterial', 'vv', inp=['aov'])]
    records += _written_out(['svc'], ['miv'])
    records[-2]['out_instances'] = ['rv']        # miv into rv, not lv
    with pytest.raises(ValueError, match='wires'):
        h2s.convert_array(records, HEART, 'Argus2026_v01')


def test_names_become_paths_under_the_heart(tmp_path):
    local = {'heart_clock': 'clock', 'lv': 'lv', 'ra': 'ra'}
    assert h2s.rename_path('lv/E_A', local) == 'heart/lv/E_A'
    assert h2s.rename_path('heart_clock/T', local) == 'heart/clock/T'
    assert h2s.rename_path('aortic_root/C', local) == 'aortic_root/C'
    spec = tmp_path / 'm_system.yaml'
    spec.write_text('output_map: {heart/u_ra: ra/u, heart/s: heart_clock/s}\n'
                    'ignore: {ra/t: time, ra_inflow/t: time, heart/q: sum}\n')
    h2s.convert_spec(str(spec), local, {'ra_inflow'}, 'Argus2026_v01')
    assert spec.read_text() == ('output_map: {heart/u_ra: heart/ra/u, heart/s: heart/clock/s}\n'
                                'ignore: {heart/ra/t: time, heart/q: sum}\n')


def test_supermodule_equivalence_compares_wrapped_and_reciprocal(monkeypatch):
    from cam_testing import system as systems
    t = np.linspace(0, 1, 5)
    a_eff = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    ref = {'heart/mt': np.array([0.2, 0.6, 1.0 - 1e-9, 0.4, 0.8]), 'heart/L_aov': 1.0 / (a_eff + 1e-9)}
    # the phase wraps a step early; L = 1/(A_eff + eps) magnifies an A_eff difference of 1e-9
    new = {'h__clock/mt': np.array([0.2, 0.6, 0.0, 0.4, 0.8]), 'h__aov/L': 1.0 / (a_eff + 2e-9)}
    monkeypatch.setattr(systems, 'load_system', lambda key: SimpleNamespace(name=key, spec={}))
    monkeypatch.setattr(systems, 'generate', lambda system, work_dir: system.name)
    monkeypatch.setattr(systems, 'simulate', lambda path, spec, solver_info: (t, ref if path == 'ref' else new))
    version = type('V', (), {'vessel_type': 'h', 'submodule_names': ['clock', 'aov']})()
    entry = {'model': 'new', 'reproduces': 'ref', 'instance': 'h', 'tol': 1e-6,
             'output_map': {'heart/mt': 'h/clock/mt', 'heart/L_aov': 'h/aov/L'}}

    def ok(e):
        rows, missing = supermodule.equivalence(version, e, 'unused')
        assert not missing
        return {r['reference']: r['ok'] for r in rows}

    assert ok(entry) == {'heart/mt': False, 'heart/L_aov': False}
    assert ok(dict(entry, wrapped={'heart/mt': 1.0}, reciprocal=['heart/L_aov'])) == {'heart/mt': True,
                                                                                      'heart/L_aov': True}
