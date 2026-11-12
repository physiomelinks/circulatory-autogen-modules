"""
verification_test_rhs: the generated right-hand side against reference values (verification_config
rhs_reference). The Van der Pol benchmark, dx/dt = y, dy/dt = mu ((1 - x^2) y - x), with reference
rates from those equations at points that also set mu. Results go to a temporary directory.
"""
import copy

import pytest

from cam_testing import checks, report
from cam_testing.library import load_version

POINTS = [(0.0, 2.0, 0.0, 100.0), (1.3, -0.7, 1.5, 0.5), (4.0, 0.2, -3.0, 7.0)]


def _reference():
    pts = [{'t': t, 'states': [x, y], 'parameters': [mu], 'rates': [y, mu * ((1 - x * x) * y - x)], 'observed': [x]}
           for t, x, y, mu in POINTS]
    return {'source': 'the Van der Pol equations', 'states': ['x', 'y'], 'parameters': ['mu'], 'observed': ['x'],
            'points': pts}


@pytest.fixture(scope='module')
def out_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp('rhs_results')
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(checks, 'result_path', lambda c, test, p=None: str(d / f'{test}.json'))
        mp.setattr(checks, 'plot_path', lambda c, name, p=None: str(d / f'{name}.png'))
        yield d


def _cm(out_dir, ref):
    version = load_version('VanDerPol', 'VanDerPol1926_v01')
    version.spec = dict(version.spec, rhs_reference=ref) if ref else \
        {k: v for k, v in version.spec.items() if k != 'rhs_reference'}
    return checks.ComponentModel(version, str(out_dir))


def test_matching_reference_passes(out_dir):
    r = checks.verification_test_rhs(_cm(out_dir, _reference()))
    assert r.status == checks.PASSED, r.message
    assert r.metrics['points'] == 3 and r.metrics['worst_relative_error'] <= 1.0


def test_a_wrong_rate_fails(out_dir):
    ref = _reference()
    ref['points'][1]['rates'][1] *= 1 + 1e-9
    r = checks.verification_test_rhs(_cm(out_dir, ref))
    assert r.status == checks.FAILED and '1 of 3 points' in r.message


def test_parameters_are_set_per_point(out_dir):
    # the same rates at the default mu (100) would be wrong at the second and third points
    ref = _reference()
    ref['parameters'] = []
    r = checks.verification_test_rhs(_cm(out_dir, ref))
    assert r.status == checks.FAILED and '2 of 3 points' in r.message


def test_a_wrong_observed_value_fails(out_dir):
    ref = copy.deepcopy(_reference())
    ref['points'][0]['observed'] = [2.5]
    assert checks.verification_test_rhs(_cm(out_dir, ref)).status == checks.FAILED


def test_states_must_be_the_models(out_dir):
    ref = _reference()
    ref['states'] = ['x', 'z']
    r = checks.verification_test_rhs(_cm(out_dir, ref))
    assert r.status == checks.FAILED and 'not the model' in r.message


def test_not_applicable_without_a_reference(out_dir):
    cm = _cm(out_dir, None)
    assert checks.verification_test_rhs(cm).status == checks.NOT_APPLICABLE
    assert report.not_applicable_reason(cm.component, 'verification_test_rhs') == checks.NO_RHS_REFERENCE
    assert report.not_applicable_reason(_cm(out_dir, _reference()).component, 'verification_test_rhs') is None

