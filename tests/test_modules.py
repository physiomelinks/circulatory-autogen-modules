"""
The standard V&V tests, parametrised over every (module, component) that has a test spec.

    pytest                                  # every module
    pytest --module diffusion_volume        # one module (repeatable)
    pytest -m "not slow"                    # skip calibration

Each test writes its result JSON and plots into modules/<name>/results and plots/, which
the report generator (python -m cam_testing.report) turns into modules/<name>/<name>.html.
"""
import pytest

from cam_testing import checks

pytestmark = pytest.mark.module_vv


def _assert(result, component_model):
    if result.status == checks.PASSED:
        return
    if result.status in (checks.SKIPPED, checks.PENDING):
        pytest.skip(f'{result.status}: {result.message}')
    # A failure listed under expected_failures in the spec is a recorded known issue: it
    # stays 'failed' in the results and the report, but does not fail the run.
    known = (component_model.component.spec.get('expected_failures') or {}).get(result.test)
    detail = '\n'.join(str(d) for d in result.details[:20])
    if known:
        pytest.xfail(f'known issue: {known} ({result.message})')
    pytest.fail(f'{result.message}\n{detail}', pytrace=False)


def run_test(component_model):
    _assert(checks.run_test(component_model), component_model)


def verification_test_BC(component_model):
    _assert(checks.verification_test_BC(component_model), component_model)


def verification_test_timestep(component_model):
    _assert(checks.verification_test_timestep(component_model), component_model)


def stability_test(component_model):
    _assert(checks.stability_test(component_model), component_model)


def validation_test_baseline(component_model):
    _assert(checks.validation_test_baseline(component_model), component_model)


@pytest.mark.slow
def validation_test_calibrate(component_model):
    _assert(checks.validation_test_calibrate(component_model), component_model)
