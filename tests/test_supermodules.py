"""
Supermodules (modules/supermodules/<name>/): each harness in the spec closes the supermodule
with host vessels, is flattened by cam_testing.supermodule.expand, generated and run.

    pytest tests/test_supermodules.py

  supermodule_structure_test    every internal vessel is a library (vessel_type, BC_type) and
                                every '@<port>' marker is a declared port
  supermodule_run_test          each harness generates and runs; outputs finite
  supermodule_equivalence_test  a harness with equivalent_to reproduces that system model
  supermodule_invariants_test   the harness invariants hold
"""
import os

import numpy as np
import pytest

from cam_testing import checks, library, supermodule as sm, system as systems

pytestmark = pytest.mark.system_model

SUPERMODULES = sm.supermodules()
HARNESSES = [pytest.param(s, h, id=f'{s.name}/{h["name"]}')
             for s in SUPERMODULES for h in ((s.spec.get('tests') or {}).get('harnesses') or [])]


def _harness_system(s, h, work_dir):
    res = os.path.join(work_dir, 'resources')
    sm.write_model(os.path.join(s.dir, h['vessel_array']), os.path.join(s.dir, h['parameters']), res, h['name'])
    spec = {'sim_time': h.get('sim_time', 2.0), 'pre_time': h.get('pre_time', 0.0), 'dt': h.get('dt', 0.01),
            'solver_info': h.get('solver_info') or {}, 'invariants': h.get('invariants') or []}
    return systems.System(h['name'], 'supermodules', res, spec)


def _run(s, h, work_dir):
    system = _harness_system(s, h, work_dir)
    path = systems.generate(system, work_dir)
    return systems.simulate(path, system.spec, system.spec['solver_info'])


@pytest.mark.parametrize('s', [pytest.param(s, id=s.name) for s in SUPERMODULES])
def supermodule_structure_test(s):
    known = {(e['vessel_type'], e['BC_type']) for n in library.module_names() for e in library.load_module(n).config}
    ports = s.spec.get('ports') or {}
    problems = []
    for r in sm.read_vessels(s.path('vessel_array.csv')):
        if (r['vessel_type'], r['BC_type']) not in known:
            problems.append(f'{r["name"]}: ({r["vessel_type"]}, {r["BC_type"]}) is not in the module library')
        for n in (r['inp_vessels'] + ' ' + r['out_vessels']).split():
            if n.startswith('@') and n[1:] not in ports:
                problems.append(f'{r["name"]}: {n} is not a declared port')
    names = {r['name'] for r in sm.read_vessels(s.path('vessel_array.csv'))}
    for p, d in ports.items():
        problems += [f'port {p}: {v} is not an internal vessel' for v in d.get('vessels', []) if v not in names]
    assert not problems, '\n'.join(problems)


@pytest.mark.parametrize('s, h', HARNESSES)
def supermodule_run_test(s, h, tmp_path):
    t, out = _run(s, h, str(tmp_path))
    bad = [k for k, v in out.items() if not np.all(np.isfinite(v))]
    assert out and not bad, f'non-finite outputs: {bad[:10]}'


@pytest.mark.parametrize('s, h', HARNESSES)
def supermodule_equivalence_test(s, h, tmp_path):
    if not h.get('equivalent_to'):
        pytest.skip('harness has no equivalent_to')
    target = systems.load_system(h['equivalent_to'])
    t, new = _run(s, h, str(tmp_path / 'harness'))
    path = systems.generate(target, str(tmp_path / 'target'))
    t_ref, ref = systems.simulate(path, target.spec, h.get('solver_info') or {})
    rows, missing = systems.compare(ref, new, {}, float(h.get('tol', 1e-9)))
    bad = [r for r in rows if not r['ok']]
    assert rows and not bad and not missing, (
        f'{len(bad)} of {len(rows)} outputs differ; missing {missing[:10]}\n'
        + '\n'.join(f"{r['reference']}: {r['difference']:.3g}" for r in bad[:10]))


@pytest.mark.parametrize('s, h', HARNESSES)
def supermodule_invariants_test(s, h, tmp_path):
    if not h.get('invariants'):
        pytest.skip('harness has no invariants')
    t, out = _run(s, h, str(tmp_path))
    failures = checks._check_invariants({'invariants': h['invariants']}, t,
                                        {k.replace('/', '__'): v for k, v in out.items()}, {}, where='run')
    assert not failures, '; '.join(failures)
