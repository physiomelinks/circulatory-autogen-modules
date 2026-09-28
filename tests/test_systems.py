"""
System models (modules/system/<category>/<model>/): generated from this module library,
checked to run, and checked to reproduce circulatory_autogen's original model.

    pytest tests/test_systems.py                       # every system model
    pytest tests/test_systems.py -k 3compartment       # one model
"""
import numpy as np
import pytest

from cam_testing import checks, system as systems

pytestmark = pytest.mark.system_model

SYSTEMS = systems.system_models()


def _param(s):
    marks = []
    if s.spec.get('skip'):
        marks.append(pytest.mark.skip(reason=s.spec['skip']))
    return pytest.param(s, id=s.id, marks=marks)


PARAMS = [_param(s) for s in SYSTEMS]


def _known(s, test, message):
    reason = (s.spec.get('expected_failures') or {}).get(test)
    if reason:
        pytest.xfail(f'known issue: {reason} ({message[:300]})')
    pytest.fail(message, pytrace=False)


@pytest.mark.parametrize('system', PARAMS)
def system_run_test(system, tmp_path):
    try:
        t, out = systems.run_model(system, str(tmp_path))
    except Exception as e:
        systems.save_result(system, 'system_run_test', 'failed', f'{type(e).__name__}: {str(e)[:2000]}')
        _known(system, 'system_run_test', f'{type(e).__name__}: {str(e)[-2000:]}')
    bad = [k for k, v in out.items() if not np.all(np.isfinite(v))]
    status = 'failed' if bad else 'passed'
    msg = f'non-finite: {bad[:10]}' if bad else f'generated and ran {system.spec.get("sim_time")} s; {len(out)} outputs finite'
    systems.save_result(system, 'system_run_test', status, msg, {'n_outputs': len(out)})
    if bad:
        _known(system, 'system_run_test', msg)


@pytest.mark.parametrize('system', PARAMS)
def system_equivalence_test(system, tmp_path):
    eq = system.spec.get('equivalence') or {}
    if eq.get('status') in ('not_applicable', 'skipped'):
        systems.save_result(system, 'system_equivalence_test', eq['status'], eq.get('reason', ''))
        pytest.skip(eq.get('reason', ''))
    try:
        t, ref, new, rows, missing = systems.equivalence(system, str(tmp_path))
    except Exception as e:
        msg = f'{type(e).__name__}: {str(e)[-2000:]}'
        systems.save_result(system, 'system_equivalence_test', 'failed', msg)
        _known(system, 'system_equivalence_test', msg)
    worst = max(rows, key=lambda r: r['difference']) if rows else None
    bad = [r for r in rows if not r['ok']]
    metrics = {'compared': len(rows), 'not_in_model': missing, 'tol': eq.get('tol', 1e-6),
               'worst': worst, 'rows': sorted(rows, key=lambda r: -r['difference'])}
    if bad or missing or not rows:
        msg = (f'{len(bad)} of {len(rows)} outputs differ beyond {eq.get("tol", 1e-6)}; '
               f'{len(missing)} reference outputs not found in the model')
        systems.save_result(system, 'system_equivalence_test', 'failed', msg, metrics,
                            [f"{r['reference']} -> {r['model']}: {r['difference']:.3g}" for r in bad[:30]] + missing[:30])
        _known(system, 'system_equivalence_test', msg + '\n' + '\n'.join(f"{r['reference']}: {r['difference']:.3g}" for r in bad[:10]))
    systems.save_result(system, 'system_equivalence_test', 'passed',
                        f"all {len(rows)} outputs of circulatory_autogen's original reproduced; worst "
                        f"{worst['reference']} {worst['difference']:.2e} <= {eq.get('tol', 1e-6)}", metrics)


@pytest.mark.parametrize('system', PARAMS)
def system_invariants_test(system, tmp_path):
    invariants = system.spec.get('invariants') or []
    if not invariants:
        systems.save_result(system, 'system_invariants_test', 'skipped', 'no invariants in the spec')
        pytest.skip('no invariants in the spec')
    t, out = systems.run_model(system, str(tmp_path))
    env_out = {k.replace('/', '__'): v for k, v in out.items()}
    failures = checks._check_invariants(system.spec, t, env_out, {}, where='run')
    status = 'failed' if failures else 'passed'
    systems.save_result(system, 'system_invariants_test', status,
                        '; '.join(failures) if failures else f'all {len(invariants)} invariants hold')
    if failures:
        _known(system, 'system_invariants_test', '; '.join(failures))
