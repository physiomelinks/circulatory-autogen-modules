"""
cam_testing.smooth: an ad hoc module's input that is constant in the system is swept over its
declared range, or reported as not exercised. Split so that each file runs in under about 6 minutes at a machine load of 60.
"""
import smoothing_relaxer as relaxer
from cam_testing import checks, smooth


def test_an_input_constant_in_the_system_is_swept_over_its_declared_range_or_reported(tmp_path):
    """u is 0 in the system: factors x 0 give no sweep, so abs(u) would be checked without ever moving.
    With no declared range the module is reported with u not exercised and the site is left as written;
    with one (a bound naming a system parameter) u is swept and the site smoothed."""
    ext, array, params = relaxer.write_system(tmp_path, u=0)
    spec = {'sim_time': 0.5, 'pre_time': 0.0, 'dt': 0.01}

    def module(name, **kwargs):
        report = smooth.run(str(array), str(params), external_dirs=[str(ext)], output_dir=str(tmp_path / name),
                            work_dir=str(tmp_path / f'work_{name}'), spec=spec, modules_only=True, **kwargs)
        return report['modules']['r']

    held = module('held')
    assert held['not_exercised_inputs'] == ['u'] and held['message'].startswith('input u not exercised')
    assert 'not exercised (x constant at every point, left as written): x' in held['message']
    assert 'u' in (tmp_path / 'held' / 'module_table.csv').read_text().splitlines()[1].split(',')
    # declared in the settings (run() also takes input_ranges, as --input-ranges gives it)
    swept = module('swept', settings={'input_ranges': {'r/u': {'min': -1, 'max': 'k', 'points': 3}}})
    assert swept['not_exercised_inputs'] == [] and swept['status'] == checks.PASSED, swept['message']
    assert '1 sites (abs 1)' in swept['message'] and '8 points' in swept['message']   # u at -1, 0.5, 2
