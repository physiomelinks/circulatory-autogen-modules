"""
verification_test_smoothing and its error metric. The pipeline (cam_testing.smooth) is tested in
test_smoothing_ad_hoc.py, test_smoothing_inputs.py and test_smoothing_system.py. Split so that each file runs in under about 6 minutes at a machine load of 60.
"""
import json
import os

import numpy as np
import pytest

from cam_testing import checks, smooth
from cam_testing.library import load_version, version_index


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


def test_spec_exemptions_leave_sites_as_written(monkeypatch):
    version = load_version('cardiac_clock', 'controlled')
    monkeypatch.setitem(version.spec, 'smoothing', {'exempt': [['T_eff', 'test: a switch meant to stay']],
                                                    'parameters': ['T']})
    result = checks.verification_test_smoothing(checks.ComponentModel(version))
    assert result.status == checks.PASSED, result.message
    assert result.metrics['counts'] == {'maxmin_exempt': 1, 'floor_not_smoothed': 1}
    report = json.loads(open(os.path.join(version.dir, result.metrics['report'])).read())
    site = next(x for x in report['sites'] if x['variable'] == 'T_eff')
    assert site['status'] == 'exempt' and site['reason'] == 'test: a switch meant to stay'


def test_an_output_is_judged_on_max_of_its_range_its_mean_and_atol():
    """RMS / max(R, M, atol), not / its own norm: a clamp that is 0 throughout is judged against atol, an
    output that is constant and nonzero (the afferent modules' outputs with their inputs held) against its
    value, a varying one against its range; abs_floor overrides per output."""
    t = np.linspace(0, 1, 101)
    originals = [{'clamp': np.zeros_like(t), 'const': np.full_like(t, 3.0), 'p': np.sin(2 * np.pi * t)},
                 {'clamp': np.zeros_like(t), 'const': np.full_like(t, 3.0), 'p': 2 * np.sin(2 * np.pi * t)}]
    smoothed = {'clamp': np.full_like(t, 2e-7), 'const': np.full_like(t, 3.0 + 3e-7), 'p': 2 * np.sin(2 * np.pi * t) + 1e-3}
    scales = checks.output_scales(originals, atol=1e-6)
    assert scales == {'clamp': 1e-6, 'const': pytest.approx(3.0), 'p': pytest.approx(4.0, rel=1e-3)}
    errors = checks._relative_errors(originals[1], smoothed, scales)
    assert errors['clamp'][0] == pytest.approx(0.2)
    assert errors['const'][0] == pytest.approx(1e-7)
    assert errors['p'][0] == pytest.approx(1e-3 / 4.0, rel=1e-3)
    stated = checks.output_scales(originals, {'clamp': 1e-3}, atol=1e-6)
    assert checks._relative_errors(originals[1], smoothed, stated)['clamp'][0] == pytest.approx(2e-4)
