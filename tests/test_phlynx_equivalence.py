# cam_testing.phlynx.equivalence_check, apart from PhLynx and CUFLynx.
from types import SimpleNamespace

from cam_testing import checks, harness, phlynx


def test_a_reference_libcuflynx_cannot_generate_is_a_failed_result(monkeypatch, tmp_path):
    """So the version's expected_failures can name it, rather than the exception ending the test."""
    def fail(component, work_dir):
        raise harness.GenerationFailed('generation returned False:\nVariable q has no units')

    monkeypatch.setattr(harness, 'generate', fail)
    monkeypatch.setattr(phlynx, '_cpp', lambda component: False)
    component = SimpleNamespace(is_supermodule=False, spec={'outputs': ['q']})
    r = phlynx.equivalence_check(component, {'cuflynx': {'ok': True}}, str(tmp_path))
    assert r.status == checks.FAILED
    assert r.message == 'libcuflynx could not generate the reference model: generation returned False:'
