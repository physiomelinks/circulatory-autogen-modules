"""
PhLynx -> CUFLynx pipeline, per version: the version's test network is built in PhLynx
(its real code, from a checkout) out of this library, exported as the .omex PhLynx sends to
CUFLynx, imported into a released CUFLynx and simulated there, and compared with libcuflynx's
model of the same network.

    pytest tests/test_phlynx.py --module coupling            # a category, or a module_type
    PHLYNX_DIR=../phlynx CUFLYNX_BIN=~/software/CUFLynx pytest tests/test_phlynx.py

Skips (not fails) when node, PhLynx or the CUFLynx binary isn't available. A supermodule version
has no PhLynx form (PhLynx has no supermodule support): its export test fails with that reason and
the other two are skipped, without running the pipeline. A failure listed
under a version's expected_failures (e.g. a module feature PhLynx doesn't support yet) is a
recorded known issue, as in the other tests.
"""
import pytest

from cam_testing import checks, phlynx
from cam_testing.library import load_module_type, load_version

pytestmark = pytest.mark.phlynx_pipeline

_runs = {}


@pytest.fixture
def pipeline(component_key, tmp_path_factory):
    # every version of a module_type in one PhLynx run and one CUFLynx run
    mt, version_name = component_key
    component = load_version(mt, version_name)
    if component.is_supermodule:
        # PhLynx has no supermodule form: the checks record that without running the pipeline
        return component, None, None
    if mt not in _runs:
        versions = [v for v in load_module_type(mt).versions() if not v.is_supermodule]
        work_dir = str(tmp_path_factory.mktemp(f'phlynx_{mt}'))
        try:
            _runs[mt] = (phlynx.run_versions(mt, versions, keep_dir=work_dir), work_dir)
        except phlynx.PipelineUnavailable as e:
            _runs[mt] = (e, work_dir)
    results, work_dir = _runs[mt]
    if isinstance(results, phlynx.PipelineUnavailable):
        pytest.skip(f'PhLynx/CUFLynx pipeline unavailable: {results}')
    return component, results.get(component.id), work_dir


def _record(component, result):
    checks.save(component, result)
    if result.status == checks.PASSED:
        return
    if result.status in (checks.SKIPPED, checks.NOT_APPLICABLE):
        pytest.skip(f'{result.status}: {result.message}')
    known = (component.spec.get('expected_failures') or {}).get(result.test)
    if known:
        pytest.xfail(f'known issue: {known} ({result.message})')
    pytest.fail(result.message + '\n' + '\n'.join(str(d) for d in result.details[:20]), pytrace=False)


def phlynx_export_test(pipeline):
    component, res, _ = pipeline
    _record(component, phlynx.export_check(component, res))


def cuflynx_simulate_test(pipeline):
    component, res, _ = pipeline
    _record(component, phlynx.simulate_check(component, res))


def phlynx_equivalence_test(pipeline, tmp_path):
    component, res, _ = pipeline
    _record(component, phlynx.equivalence_check(component, res, str(tmp_path)))
