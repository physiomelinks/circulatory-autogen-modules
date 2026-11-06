"""
cam_testing.smooth: a module with no library version gets the ad hoc harness, and a failing module
stops the pipeline before the system step. Split so that each file runs in under about 6 minutes at a machine load of 60.
"""
import json
import os

import pytest

import smoothing_relaxer as relaxer
from cam_testing import checks, smooth

SYSTEMS = os.path.join(os.path.dirname(__file__), '..', 'system_models')


def test_a_module_with_no_library_version_gets_the_ad_hoc_harness(tmp_path):
    ext, array, params = relaxer.write_system(tmp_path, u=-0.5)
    # u is constant in the system: it is swept over its declared range, not factors x its value
    report = smooth.run(str(array), str(params), external_dirs=[str(ext)], output_dir=str(tmp_path / 'out'),
                        work_dir=str(tmp_path / 'work'), spec={'sim_time': 0.5, 'pre_time': 0.0, 'dt': 0.01},
                        input_ranges={'r/u': [-1, 0]})
    module = report['modules']['r']
    assert module['kind'] == smooth.AD_HOC and module['status'] == checks.PASSED, module['message']
    assert '1 sites (abs 1)' in module['message']
    assert '9 points' in module['message']        # nominal, 4 factors of k, and u at -1, -0.75, -0.25, 0
    assert report['system']['status'] == checks.PASSED
    assert report['system']['sites'][0]['module'] == {'module': 'relaxer/v01', 'kind': smooth.AD_HOC,
                                                      'status': checks.PASSED}


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
