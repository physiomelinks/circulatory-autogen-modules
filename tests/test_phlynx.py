"""
PhLynx -> CUFLynx pipeline, per component: the component's test network is built in PhLynx
(its real code, from a checkout) out of this library, exported as the .omex PhLynx sends to
CUFLynx, imported into a released CUFLynx and simulated there, and compared with libcuflynx's
model of the same network.

    pytest tests/test_phlynx.py --module coupling
    PHLYNX_DIR=../phlynx CUFLYNX_BIN=~/software/CUFLynx pytest tests/test_phlynx.py

Skips (not fails) when node, PhLynx or the CUFLynx binary isn't available. A failure listed
under a component's expected_failures (e.g. a module feature PhLynx doesn't support yet) is a
recorded known issue, as in the other tests.
"""
import pytest

from cam_testing import checks, phlynx
from cam_testing.library import load_module

pytestmark = pytest.mark.phlynx_pipeline

_runs = {}


@pytest.fixture
def pipeline(component_key, tmp_path_factory):
    module_name, component_id = component_key
    if module_name not in _runs:
        module = load_module(module_name)
        work_dir = str(tmp_path_factory.mktemp(f'phlynx_{module_name}'))
        try:
            _runs[module_name] = (phlynx.run_module(module, module.components(), keep_dir=work_dir), work_dir)
        except phlynx.PipelineUnavailable as e:
            _runs[module_name] = (e, work_dir)
    results, work_dir = _runs[module_name]
    if isinstance(results, phlynx.PipelineUnavailable):
        pytest.skip(f'PhLynx/CUFLynx pipeline unavailable: {results}')
    component = load_module(module_name).component(component_id)
    return component, results.get(component_id), work_dir


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
