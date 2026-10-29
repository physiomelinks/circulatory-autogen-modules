'''expected_failures keyed per parameterisation ("<test>[<parameterisation>]") or for every one ("<test>").'''
from types import SimpleNamespace

from cam_testing import checks


def _version(ef):
    return SimpleNamespace(spec={'expected_failures': ef})


def test_a_parameterisation_key_covers_only_that_parameterisation():
    v = _version({'validation_test_baseline[aov]': 'cross-model comparison'})
    assert checks.expected_failure(v, 'validation_test_baseline', 'aov') == 'cross-model comparison'
    assert checks.expected_failure(v, 'validation_test_baseline', SimpleNamespace(name='aov')) == 'cross-model comparison'
    assert checks.expected_failure(v, 'validation_test_baseline', 'miv') is None
    assert checks.expected_failure(v, 'validation_test_baseline') is None


def test_a_bare_key_covers_every_parameterisation_and_the_specific_one_wins():
    v = _version({'validation_test_baseline': 'all', 'validation_test_baseline[aov]': 'aov only'})
    assert checks.expected_failure(v, 'validation_test_baseline', 'miv') == 'all'
    assert checks.expected_failure(v, 'validation_test_baseline', 'aov') == 'aov only'
    assert checks.expected_failure(v, 'run_test') is None
