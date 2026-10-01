"""
The standard V&V tests: per version (module_type/version) for verification, and per instance
(module_type/version/instance) for validation.

    pytest                                        # every reviewed version
    pytest --module Lotka_Volterra                # one module_type (or a category: --module cell)
    pytest --component Lotka_Volterra/nn          # one version
    pytest -m "not slow"                          # skip calibration

Each test writes its result JSON and plots into the version's results/ and plots/, which the
report generator (python -m cam_testing.report) turns into <module_type>_<version>.html and
<module_type>.html.
"""
import pytest

from cam_testing import checks, supermodule
from cam_testing.library import load_version

pytestmark = pytest.mark.module_vv


def _assert(result, version):
    if result.status == checks.PASSED:
        return
    if result.status in (checks.SKIPPED, checks.PENDING, checks.NOT_APPLICABLE):
        pytest.skip(f'{result.status}: {result.message}')
    # A failure listed under expected_failures in the spec is a recorded known issue: it
    # stays 'failed' in the results and the report, but does not fail the run.
    known = (version.spec.get('expected_failures') or {}).get(result.test)
    detail = '\n'.join(str(d) for d in result.details[:20])
    if known:
        pytest.xfail(f'known issue: {known} ({result.message})')
    pytest.fail(f'{result.message}\n{detail}', pytrace=False)


# ---- per version --------------------------------------------------------------------------------

def run_test(component_model):
    _assert(checks.run_test(component_model), component_model.component)


def verification_test_invariants(component_model):
    _assert(checks.verification_test_invariants(component_model), component_model.component)


def verification_test_BC(component_model):
    _assert(checks.verification_test_BC(component_model), component_model.component)


def verification_test_timestep(component_model):
    _assert(checks.verification_test_timestep(component_model), component_model.component)


def stability_test(component_model):
    _assert(checks.stability_test(component_model), component_model.component)


# ---- per instance -------------------------------------------------------------------------------

def validation_test_baseline(instance_model):
    _assert(checks.validation_test_baseline(instance_model), instance_model.component)


@pytest.mark.slow
def validation_test_calibrate(instance_model):
    _assert(checks.validation_test_calibrate(instance_model), instance_model.component)


# ---- supermodule versions -----------------------------------------------------------------------

def supermodule_structure_test(supermodule_key):
    version = load_version(*supermodule_key)
    _assert(supermodule.structure_check(version), version)


@pytest.mark.system_model
def supermodule_equivalence_test(equivalence_key, tmp_path):
    mt, v, model = equivalence_key
    version = load_version(mt, v)
    entry = next(e for e in version.spec['supermodule']['equivalent'] if e['model'] == model)
    _assert(supermodule.equivalence_check(version, entry, str(tmp_path)), version)
