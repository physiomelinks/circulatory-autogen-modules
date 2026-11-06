"""
cam_testing.smooth's system step: the system's exemptions reach the ad hoc module and the system, a
given system model is used instead of generating one, and the inputs' means are read through connections. Split so that each file runs in under about 6 minutes at a machine load of 60.
"""
import os

import pytest

import smoothing_relaxer as relaxer
from cam_testing import checks, smooth

SYSTEMS = os.path.join(os.path.dirname(__file__), '..', 'system_models')


def test_the_system_exemptions_reach_the_ad_hoc_module_and_the_system(tmp_path):
    ext, array, params = relaxer.write_system(tmp_path, u=-0.5)
    report = smooth.run(str(array), str(params), external_dirs=[str(ext)], output_dir=str(tmp_path / 'out'),
                        work_dir=str(tmp_path / 'work'), spec={'sim_time': 0.5, 'pre_time': 0.0, 'dt': 0.01},
                        exempt=[{'variable': 'r_module/x', 'reason': 'test'}])
    # the same system, given as a generated model instead of generated again
    again = smooth.run(str(array), str(params), external_dirs=[str(ext)], output_dir=str(tmp_path / 'out2'),
                       work_dir=str(tmp_path / 'work2'), spec={'sim_time': 0.5, 'pre_time': 0.0, 'dt': 0.01},
                       exempt=[{'variable': 'r_module/x', 'reason': 'test'}], modules_only=True,
                       system_model=str(tmp_path / 'work' / 'system' / 'generated_models' / 'relax' / 'relax.cellml'))
    assert not os.path.exists(tmp_path / 'work2' / 'system') and again['modules']['r']['status'] == checks.PASSED
    assert 'abs_exempt 1' in report['modules']['r']['message']
    assert report['system']['sites'][0]['status'] == 'exempt'
    assert report['exempt'] == [{'variable': 'r_module/x', 'reason': 'test'}]


def test_system_means_read_a_connected_input_from_its_source(tmp_path):
    """systemic_T/u_out is connected to venous_svc's u: the solver keeps one name for both."""
    d = os.path.join(SYSTEMS, 'closed_loop_cvs', '3compartment')
    records = smooth.module_array.read_records(os.path.join(d, '3compartment_module_array.json'))
    model = smooth._generate(str(tmp_path), '3compartment', records, os.path.join(d, '3compartment_parameters.json'), [])
    spec = dict(smooth.DEFAULT_SPEC, pre_time=0.0, sim_time=0.5, dt=0.01)
    means = smooth.system_means(model, spec, ['systemic_T/u_out', 'venous_svc/u'])
    assert set(means) == {'systemic_T/u_out', 'venous_svc/u'}
    assert means['systemic_T/u_out'] == pytest.approx(means['venous_svc/u'])
