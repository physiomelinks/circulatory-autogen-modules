"""
verification_test_smoothing and the module-by-module smoothing pipeline (cam_testing.smooth).
Each test generates small models (seconds to tens of seconds).
"""
import json
import os
import textwrap

import pytest

from cam_testing import checks, smooth
from cam_testing.library import load_version, version_index

SYSTEMS = os.path.join(os.path.dirname(__file__), '..', 'system_models')


def test_a_version_with_a_site_passes_and_a_jumping_output_is_not_compared():
    """cardiac_clock/controlled: T_eff = max(T + Delta_T, T_min) is one site; mt = s - floor(s)
    jumps at each wrap, so it is left out of the comparison."""
    result = checks.verification_test_smoothing(checks.ComponentModel(load_version('cardiac_clock', 'controlled')))
    assert result.status == checks.PASSED, result.message
    m = result.metrics
    assert m['n_sites'] == 1 and m['counts'] == {'maxmin': 1, 'floor_not_smoothed': 1}
    assert m['not_compared'] == ['mt'] and m['n_points'] > 1
    assert m['fingerprint'] == checks.smoothing_fingerprint(load_version('cardiac_clock', 'controlled'),
                                                            checks.smoothing_settings({}))
    version = load_version('cardiac_clock', 'controlled')
    assert os.path.isfile(os.path.join(version.dir, m['smoothed_model']))
    assert checks.load(version, 'verification_test_smoothing').status == checks.PASSED


def test_a_smooth_version_records_no_sites():
    result = checks.verification_test_smoothing(checks.ComponentModel(load_version('Lotka_Volterra', 'Lotka1925_v01')))
    assert result.status == checks.PASSED
    assert result.message == checks.NO_SITES
    assert result.metrics['n_sites'] == 0


def test_a_tight_tol_fails_and_names_the_worst_site(monkeypatch):
    version = load_version('cardiac_clock', 'controlled')
    monkeypatch.setitem(version.spec, 'smoothing', {'rel': 0.3, 'tol': 1e-9, 'parameters': ['T']})
    result = checks.verification_test_smoothing(checks.ComponentModel(version))
    assert result.status == checks.FAILED
    assert result.metrics['attribution'][0]['variable'] == 'T_eff'
    assert 'worst sites: T_eff' in result.message


def test_the_library_versions_resolve_through_the_renames():
    index = version_index(include_libraries=True)
    assert smooth.resolve({'module_type': 'heart', 'module_subtype': 'vp'}, index).key == 'heart/vp'
    assert smooth.resolve({'module_type': 'heart_simple', 'module_subtype': 'vp'}, index).key == 'heart/vp_simple'
    assert smooth.resolve({'module_type': 'nowhere', 'module_subtype': 'v'}, index) is None


def test_the_pipeline_stops_before_the_system_when_a_module_fails(tmp_path, monkeypatch):
    def failing(version, settings_override=None, force=False):
        return checks.Result('verification_test_smoothing', checks.FAILED, 'over tol', {}), False

    def no_system(*args, **kwargs):
        raise AssertionError('the system step ran')

    monkeypatch.setattr(smooth, 'library_check', failing)
    monkeypatch.setattr(checks, 'compare_smoothed', no_system)
    lv = os.path.join(SYSTEMS, 'benchmarks', 'Lotka_Volterra')
    with pytest.raises(smooth.SmoothingStopped, match='modules failed: Lotka_Volterra'):
        smooth.run(os.path.join(lv, 'Lotka_Volterra_module_array.json'), os.path.join(lv, 'Lotka_Volterra_parameters.json'),
                   output_dir=str(tmp_path), work_dir=str(tmp_path / 'work'), spec={'sim_time': 1.0, 'pre_time': 0.0})
    report = json.loads((tmp_path / 'smoothing_report.json').read_text())
    assert report['system']['status'] == 'not run'
    assert report['modules']['Lotka_Volterra']['status'] == checks.FAILED


RELAXER_CELLML = textwrap.dedent('''\
    <?xml version='1.0' encoding='UTF-8'?>
    <model name="relaxer_v01" xmlns="http://www.cellml.org/cellml/1.1#" xmlns:cellml="http://www.cellml.org/cellml/1.1#">
        <component name="relaxer">
            <variable name="t" public_interface="in" units="second"/>
            <variable name="k" public_interface="in" units="dimensionless"/>
            <variable name="u" public_interface="in" units="dimensionless"/>
            <variable name="x_init" public_interface="in" units="dimensionless"/>
            <variable initial_value="x_init" name="x" public_interface="out" units="dimensionless"/>
            <math xmlns="http://www.w3.org/1998/Math/MathML">
                <apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>x</ci></apply>
                    <apply><divide/><apply><times/><ci>k</ci><apply><minus/><apply><abs/><ci>u</ci></apply><ci>x</ci></apply></apply>
                    <cn cellml:units="second">1</cn></apply></apply>
            </math>
        </component>
    </model>
    ''')


def _relaxer_entry():
    return {'module_type': 'relaxer', 'module_subtype': 'v01', 'module_format': 'cellml',
            'component_file': 'relaxer_v01_modules.cellml', 'component_type': 'relaxer',
            'entrance_ports': [], 'exit_ports': [], 'general_ports': [],
            'variables_and_units': [['k', 'dimensionless', 'access', 'constant'],
                                    ['u', 'dimensionless', 'access', 'boundary_condition'],
                                    ['x_init', 'dimensionless', 'access', 'constant'],
                                    ['x', 'dimensionless', 'access', 'variable']]}


def test_a_module_with_no_library_version_gets_the_ad_hoc_harness(tmp_path):
    ext = tmp_path / 'ext'
    ext.mkdir()
    (ext / 'relaxer_v01_modules.cellml').write_text(RELAXER_CELLML)
    (ext / 'relaxer_v01_modules_config.json').write_text(json.dumps([_relaxer_entry()]))
    array = tmp_path / 'relax_module_array.json'
    array.write_text(json.dumps([{'name': 'r', 'module_type': 'relaxer', 'module_subtype': 'v01',
                                  'inp_instances': [], 'out_instances': []}]))
    params = tmp_path / 'relax_parameters.json'
    params.write_text(json.dumps({'schema_version': 1, 'r': {
        'k': {'value': 2, 'units': 'dimensionless', 'data_reference': 'test'},
        'u': {'value': -0.5, 'units': 'dimensionless', 'data_reference': 'test'},
        'x_init': {'value': 0, 'units': 'dimensionless', 'data_reference': 'test'}}}))
    report = smooth.run(str(array), str(params), external_dirs=[str(ext)], output_dir=str(tmp_path / 'out'),
                        work_dir=str(tmp_path / 'work'), spec={'sim_time': 2.0, 'pre_time': 0.0, 'dt': 0.01})
    module = report['modules']['r']
    assert module['kind'] == smooth.AD_HOC and module['status'] == checks.PASSED, module['message']
    assert '1 sites (abs 1)' in module['message']
    assert '9 points' in module['message']        # nominal, and 4 factors for each of k and u (x_init is 0)
    assert report['system']['status'] == checks.PASSED
    assert report['system']['sites'][0]['module'] == {'module': 'relaxer/v01', 'kind': smooth.AD_HOC,
                                                      'status': checks.PASSED}


def test_spec_exemptions_leave_sites_as_written(monkeypatch):
    version = load_version('cardiac_clock', 'controlled')
    monkeypatch.setitem(version.spec, 'smoothing', {'exempt': [['T_eff', 'test: a switch meant to stay']],
                                                    'parameters': ['T']})
    result = checks.verification_test_smoothing(checks.ComponentModel(version))
    assert result.status == checks.PASSED, result.message
    assert result.metrics['counts'] == {'maxmin_exempt': 1, 'floor_not_smoothed': 1}
    report = json.loads(open(os.path.join(version.dir, result.metrics['report'])).read())
    assert report['sites'][0]['status'] == 'exempt'
    assert report['sites'][0]['reason'] == 'test: a switch meant to stay'


def test_the_system_exemptions_reach_the_ad_hoc_module_and_the_system(tmp_path):
    ext = tmp_path / 'ext'
    ext.mkdir()
    (ext / 'relaxer_v01_modules.cellml').write_text(RELAXER_CELLML)
    (ext / 'relaxer_v01_modules_config.json').write_text(json.dumps([_relaxer_entry()]))
    array = tmp_path / 'relax_module_array.json'
    array.write_text(json.dumps([{'name': 'r', 'module_type': 'relaxer', 'module_subtype': 'v01',
                                  'inp_instances': [], 'out_instances': []}]))
    params = tmp_path / 'relax_parameters.json'
    params.write_text(json.dumps({'schema_version': 1, 'r': {
        'k': {'value': 2, 'units': 'dimensionless', 'data_reference': 'test'},
        'u': {'value': -0.5, 'units': 'dimensionless', 'data_reference': 'test'},
        'x_init': {'value': 0, 'units': 'dimensionless', 'data_reference': 'test'}}}))
    report = smooth.run(str(array), str(params), external_dirs=[str(ext)], output_dir=str(tmp_path / 'out'),
                        work_dir=str(tmp_path / 'work'), spec={'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01},
                        exempt=[{'variable': 'r_module/x', 'reason': 'test'}])
    assert 'abs_exempt 1' in report['modules']['r']['message']
    assert report['system']['sites'][0]['status'] == 'exempt'
    assert report['exempt'] == [{'variable': 'r_module/x', 'reason': 'test'}]
