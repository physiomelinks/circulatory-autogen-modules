"""
Unit tests of what the reports show and of the supermodule-version machinery: every version (component
or supermodule) has the same test columns, a column that doesn't apply is N/A with a reason, a test
that didn't run is never shown as passed, the version-level calibration derived from supermodules
("Pass in super" / "Fail in super"), supermodule parameter naming, the supermodule BC-sweep default,
and the PhLynx "no supermodule support" failure. No model is generated.
"""
import json
import os
from types import SimpleNamespace

import pytest

from cam_testing import pytest_plugin
from cam_testing import checks, harness, phlynx, report
from cam_testing.library import Parameter, all_versions, load_version


# ---- uniform columns ----------------------------------------------------------------------------

EXPECTED_COLUMNS = ['run_test', 'verification_test_invariants', 'verification_test_BC', 'verification_test_timestep',
                    'stability_test', 'version_calibration', 'phlynx_export_test', 'cuflynx_simulate_test',
                    'phlynx_equivalence_test', 'supermodule_structure_test', 'supermodule_equivalence_test']


def test_report_columns_in_order():
    assert report.REPORT_TESTS == EXPECTED_COLUMNS
    assert [report.TEST_SHORT[k] for k in EXPECTED_COLUMNS] == [
        'Run', 'Invariants', 'BC sweep', 'Timestep', 'Stability', 'Calibration', 'PhLynx export', 'CUFLynx simulate',
        'PhLynx equivalence', 'Supermodule structure', 'Supermodule reproduces']


def test_every_version_has_the_same_test_keys():
    versions = all_versions()
    assert any(v.is_supermodule for v in versions) and any(not v.is_supermodule for v in versions)
    for v in versions:
        keys = [t['key'] for t in report.version_tests(v)]
        assert keys == EXPECTED_COLUMNS, v.key


def test_module_page_columns_are_the_report_columns():
    ctx = report.module_context('heart')        # component and supermodule versions
    assert ctx['test_keys'] == EXPECTED_COLUMNS
    for row in ctx['versions']:
        assert [t['key'] for t in row['tests']] == EXPECTED_COLUMNS, row['version']


def test_supermodule_columns_not_applicable_to_a_component():
    v = load_version('Lotka_Volterra', 'nn')
    tests = {t['key']: t for t in report.version_tests(v)}
    for k in report.SUPERMODULE_TESTS:
        assert tests[k]['status'] == checks.NOT_APPLICABLE
        assert tests[k]['message'] == report.NOT_A_SUPERMODULE


def test_supermodule_without_equivalent_is_not_applicable():
    v = load_version('soma', 'sympathetic')
    assert not (v.spec.get('supermodule') or {}).get('equivalent')
    assert report.not_applicable_reason(v, 'supermodule_equivalence_test') == report.NO_EQUIVALENT
    assert report.not_applicable_reason(v, 'supermodule_structure_test') is None
    assert report.not_applicable_reason(v, 'verification_test_BC') is None


def test_cpp_reasons():
    cpp = SimpleNamespace(is_supermodule=False, spec={}, config={'module_format': 'cpp'})
    assert report.not_applicable_reason(cpp, 'stability_test') == checks.CPP_NOT_APPLICABLE
    assert report.not_applicable_reason(cpp, 'phlynx_export_test') == phlynx.CPP_REASON
    assert report.not_applicable_reason(cpp, 'run_test') is None


def _fake_version(tmp_path, name='v', supermodule=False, spec=None, instances=(), submodules=()):
    vdir = tmp_path / name
    v = SimpleNamespace(name=name, vessel_type=f'mt_{name}', key=f'mt_{name}/{name}', id=f'mt_{name}__{name}',
                        dir=str(vdir), results_dir=str(vdir / 'results'), plots_dir=str(vdir / 'plots'),
                        is_supermodule=supermodule, spec=dict(spec or {}), config={'module_format': 'cellml'},
                        submodules=list(submodules))
    v.instances = lambda: [SimpleNamespace(name=i, has_obs_data=data, results_dir=str(vdir / 'results' / 'instances' / i))
                           for i, data in instances]
    return v


def test_a_test_that_did_not_run_is_not_run(tmp_path):
    v = _fake_version(tmp_path, instances=[('default', False)])
    tests = {t['key']: t for t in report.version_tests(v)}
    for k in ('run_test', 'verification_test_BC', 'stability_test', 'phlynx_export_test'):
        assert tests[k]['status'] is None and tests[k]['status_label'] == 'Not run'
    assert not any(t['status'] == checks.PASSED for t in tests.values())


# ---- version-level calibration from supermodules --------------------------------------------------

def _sub(v):
    return {'name': v.name, 'module_type': v.vessel_type, 'module_subtype': v.name}


def _save_calibration(v, instance, status):
    path = os.path.join(v.results_dir, 'instances', instance, 'validation_test_calibrate.json')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump({'test': 'validation_test_calibrate', 'status': status, 'message': status}, f)


def test_supermodule_index(tmp_path):
    a = _fake_version(tmp_path, 'a')
    b = _fake_version(tmp_path, 'b')
    s1 = _fake_version(tmp_path, 's1', True, submodules=[_sub(a), _sub(b)])
    s2 = _fake_version(tmp_path, 's2', True, submodules=[_sub(a)])
    index = checks.supermodule_index([a, b, s1, s2])
    assert [v.key for v in index[('mt_a', 'a')]] == [s1.key, s2.key]
    assert [v.key for v in index[('mt_b', 'b')]] == [s1.key]
    assert ('mt_s1', 's1') not in index


def test_supermodule_index_of_the_library():
    index = checks.supermodule_index()
    assert 'soma/sympathetic' in [v.key for v in index[('SN_membrane_soma', 'nn')]]
    assert 'neuron/sympathetic' in [v.key for v in index[('soma', 'sympathetic')]]


def test_calibration_in_super_is_transitive(tmp_path):
    chan = _fake_version(tmp_path, 'chan', instances=[('default', False)])
    soma = _fake_version(tmp_path, 'soma', True, instances=[('default', False)], submodules=[_sub(chan)])
    neuron = _fake_version(tmp_path, 'neuron', True, instances=[('default', True)], submodules=[_sub(soma)])
    index = checks.supermodule_index([chan, soma, neuron])

    # neuron has no result yet: pending, which is not a pass
    assert checks.version_calibration(neuron, index).status == checks.PENDING
    r = checks.version_calibration(chan, index)
    assert r.status == checks.FAILED_IN_SUPER and 'mt_soma/soma' in r.message

    _save_calibration(neuron, 'default', checks.PASSED)
    assert checks.version_calibration(neuron, index).status == checks.PASSED
    assert checks.version_calibration(soma, index).status == checks.PASSED_IN_SUPER
    r = checks.version_calibration(chan, index)
    assert r.status == checks.PASSED_IN_SUPER
    assert [row['key'] for row in r.metrics['supermodules']] == ['mt_soma/soma']
    assert r.metrics['supermodules'][0]['status'] == checks.PASSED_IN_SUPER


def test_calibration_in_super_any_passing_supermodule(tmp_path):
    chan = _fake_version(tmp_path, 'chan', instances=[('default', False)])
    s_nodata = _fake_version(tmp_path, 's_nodata', True, instances=[('default', False)], submodules=[_sub(chan)])
    s_pass = _fake_version(tmp_path, 's_pass', True, instances=[('default', True)], submodules=[_sub(chan)])
    index = checks.supermodule_index([chan, s_nodata, s_pass])
    # a supermodule without calibration data counts as failed
    assert checks.version_calibration(s_nodata, index).status == checks.FAILED
    assert checks.version_calibration(chan, index).status == checks.FAILED_IN_SUPER
    _save_calibration(s_pass, 'default', checks.PASSED)
    r = checks.version_calibration(chan, index)
    assert r.status == checks.PASSED_IN_SUPER and 'mt_s_pass/s_pass' in r.message


def test_own_calibration_data_wins(tmp_path):
    chan = _fake_version(tmp_path, 'chan', instances=[('default', True)])
    sup = _fake_version(tmp_path, 'sup', True, instances=[('default', True)], submodules=[_sub(chan)])
    index = checks.supermodule_index([chan, sup])
    _save_calibration(sup, 'default', checks.PASSED)
    _save_calibration(chan, 'default', checks.FAILED)
    assert checks.version_calibration(chan, index).status == checks.FAILED


def test_calibration_in_supermodule_opt_out(tmp_path):
    chan = _fake_version(tmp_path, 'chan', spec={'calibration_in_supermodule': False}, instances=[('default', False)])
    sup = _fake_version(tmp_path, 'sup', True, instances=[('default', True)], submodules=[_sub(chan)])
    index = checks.supermodule_index([chan, sup])
    _save_calibration(sup, 'default', checks.PASSED)
    r = checks.version_calibration(chan, index)
    assert r.status == checks.FAILED and r.message == checks.NO_INSTANCE_CALIBRATION


def test_not_a_submodule_keeps_the_plain_rule(tmp_path):
    v = _fake_version(tmp_path, 'v', instances=[('default', False)])
    assert checks.version_calibration(v, {}).status == checks.FAILED


def test_report_links_and_counts_for_calibration_in_super():
    v = load_version('SN_membrane_soma', 'nn')
    t = next(t for t in report.version_tests(v) if t['key'] == 'version_calibration')
    assert t['status'] in (checks.PASSED_IN_SUPER, checks.FAILED_IN_SUPER)
    assert t['status_label'] in ('Pass in super', 'Fail in super')
    link = next(l for l in t['links'] if l['key'] == 'soma/sympathetic')
    assert link['href'] == '../../../versions/sympathetic/soma_sympathetic.html'
    assert link['href_module'] == '../versions/sympathetic/soma_sympathetic.html'
    counts = report._status_counts([{'tests': [t], 'instances': []}])
    assert counts.get(checks.PASSED_IN_SUPER, 0) == 0 and counts.get(checks.FAILED_IN_SUPER, 0) == 0
    assert counts['passed'] + counts['failed'] == 1


# ---- supermodule versions in the test parametrisation ---------------------------------------------

def test_supermodule_versions_get_the_verification_tests():
    config = SimpleNamespace(getoption=lambda k: ['soma'] if k == '--module' else ([] if k == '--component' else True))
    ids = [p.id for p in pytest_plugin._selected_versions(config)]
    assert 'soma/sympathetic' in ids and 'soma/sympathetic_monolithic_v01' in ids
    inst = [p.id for p in pytest_plugin._selected_instances(config)]
    assert 'soma/sympathetic/default' in inst


# ---- supermodule parameter names ----------------------------------------------------------------

def test_supermodule_model_names():
    soma = load_version('soma', 'sympathetic')
    assert harness.supermodule_model_name(soma, 'T') == 'T'
    assert harness.supermodule_model_name(soma, 'Cm_membrane') == 'Cm_mod_membrane'
    assert harness.supermodule_model_name(soma, 'tau_z_BK_i_BK') == 'tau_z_BK_mod_i_BK'
    p = Parameter('Cm_membrane', 'picoF', '30', 'x')
    assert harness.model_parameter_name(soma, p) == 'Cm_mod_membrane'
    assert harness.parameter_name(Parameter('a', '', '1', '', model_name='a_mod_x')) == 'a_mod_x'


def test_split_model_name_nested():
    neuron = load_version('neuron', 'sympathetic')
    name, var, path, owner = harness.split_model_name(neuron, 'I_in_mod_soma_membrane')
    assert (name, var, path) == ('I_in_soma_membrane', 'I_in', 'soma_membrane')
    assert owner.key == 'SN_membrane_soma/nn'
    assert harness.split_model_name(neuron, 'C_mod_axon')[3].key == 'axon/sympathetic_monolithic_v01'
    assert harness.split_model_name(neuron, 'T') == ('T', 'T', None, None)


# ---- supermodule BC sweep default ---------------------------------------------------------------

def _fake_cm(supermodule, sweep_spec):
    params = [Parameter('T', 'kelvin', '310', '', kind='global_constant'),
              Parameter('I_in_membrane', 'nanoA', '0', '', kind='boundary_condition'),
              Parameter('Cm_membrane', 'picoF', '30', '', kind='constant')]
    comp = SimpleNamespace(is_supermodule=supermodule)
    return SimpleNamespace(component=comp, spec={'bc_sweep': sweep_spec}, parameters=lambda: params)


def test_supermodule_sweep_default_is_globals_and_bcs():
    names, kind = checks.sweep_parameter_names(_fake_cm(True, {}))
    assert names == ['T', 'I_in_membrane'] and 'global' in kind
    names, _ = checks.sweep_parameter_names(_fake_cm(True, {'extra_parameters': ['Cm_membrane']}))
    assert names == ['T', 'I_in_membrane', 'Cm_membrane']
    names, _ = checks.sweep_parameter_names(_fake_cm(True, {'sweep': 'all'}))
    assert set(names) == {'T', 'I_in_membrane', 'Cm_membrane'}


def test_component_sweep_default_is_all():
    names, _ = checks.sweep_parameter_names(_fake_cm(False, {}))
    assert names == ['I_in_membrane', 'T', 'Cm_membrane']
    with pytest.raises(ValueError):
        checks.sweep_parameter_names(_fake_cm(False, {'sweep': 'nonsense'}))


# ---- PhLynx -------------------------------------------------------------------------------------

def test_phlynx_export_fails_for_a_supermodule_as_a_known_issue():
    for v in all_versions():
        if not v.is_supermodule:
            continue
        r = phlynx.export_check(v, None)
        assert r.status == checks.FAILED and r.message == phlynx.NO_SUPERMODULE_SUPPORT
        assert phlynx.simulate_check(v, None).status == checks.SKIPPED
        assert (v.spec.get('expected_failures') or {}).get('phlynx_export_test') == phlynx.NO_SUPERMODULE_SUPPORT, v.key
