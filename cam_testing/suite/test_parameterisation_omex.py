"""
Per-parameterisation COMBINE archives (cam_testing/omex.py): each parameterisation's .omex is built
from the library's files, loaded into a released CUFLynx through its HTTP API (cam_testing/bridges/cuflynx)
and simulated, and must reproduce libcuflynx's model of the same parameterisation.

    pytest tests/test_parameterisation_omex.py --module Lotka_Volterra
    CUFLYNX_BIN=~/software/CUFLynx pytest tests/test_parameterisation_omex.py

cuflynx_parameterisation_omex_test checks that CUFLynx takes the archive's master model, its obs_data and
params_for_id (when the parameterisation has them), that the run is finite, and that the outputs match
libcuflynx within 1e-6 at the output times both have. Skips when the CUFLynx binary isn't there;
C++ versions are not applicable (CUFLynx runs CellML).
"""
import glob
import os

import numpy as np
import pytest

from cam_testing import checks, harness, omex, phlynx
from cam_testing import system as systems
from cam_testing.library import load_version

pytestmark = pytest.mark.phlynx_pipeline

#: parameterisations built and simulated in one CUFLynx run. Each run launches CUFLynx once (about
#: 10 s), so the session's parameterisations share runs rather than taking one per module_type; the
#: first test of a batch makes its run, so no test builds more than this many archives.
BATCH = 40

_batches = {}   # parameterisation key -> the keys of its batch
_runs = {}      # batch (tuple of keys) -> ((built, sims, models), work_dir), or PipelineUnavailable


def _session_keys(session):
    '''The parameterisation keys of this test that the session runs (not skipped), in order.'''
    keys = []
    for item in session.items:
        if getattr(item, 'originalname', None) != 'cuflynx_parameterisation_omex_test' \
                or item.get_closest_marker('skip') is not None:
            continue
        callspec = getattr(item, 'callspec', None)
        key = callspec.params.get('parameterisation_key') if callspec else None
        if key and key not in keys:
            keys.append(key)
    return keys


def _batch_of(key, session):
    '''The batch (tuple of keys) ``key`` is simulated in: the session's keys not yet in a batch,
    BATCH at a time.'''
    if key not in _batches:
        keys = [k for k in _session_keys(session) if k not in _batches]
        if key not in keys:
            keys = [key]
        for i in range(0, len(keys), BATCH):
            batch = tuple(keys[i:i + BATCH])
            _batches.update({k: batch for k in batch})
    return _batches[key]


def _run_batch(keys, work_dir):
    '''Builds the archives of the parameterisations ``keys`` and simulates them in one CUFLynx run.
    Returns (built, sims, models): models is each one's libcuflynx model, generated for its archive.'''
    from cam_testing.pytest_plugin import UNCHANGED   # skipped as unchanged (--changed-only): not built
    jobs, reports, built, models = [], {}, {}, {}
    for mt, version_name, inst_name in keys:
        v = load_version(mt, version_name)
        if (v.format != 'cellml' and not v.is_supermodule) or v.key in UNCHANGED:
            continue
        spec = v.spec
        inst = v.parameterisation(inst_name)
        key = f'{v.id}__{inst.name}'
        try:
            # at the parameterisation's own archive path (what make omex writes), so the report links
            # this archive and CI builds each one once
            path = omex.build(v, inst, work_dir=os.path.join(work_dir, key))
        except Exception as e:  # noqa: BLE001 - recorded per parameterisation
            built[key] = e
            continue
        built[key] = path
        # the model omex.build generated for the archive (harness.generate, as the comparison would)
        generated = [g for g in glob.glob(os.path.join(work_dir, key, 'generated_models', '*', '*.cellml'))
                     if not g.endswith('_flat.cellml')]
        if len(generated) == 1:
            models[key] = generated[0]
        reports[key] = {'omex': path}
        jobs.append({'id': key, 'sim_time': float(spec.get('sim_time', 1.0)), 'dt': float(spec.get('dt', 0.01)),
                     'extra_outputs': [harness.output_name(o) for o in (spec.get('outputs') or [])]})
    sims = phlynx.simulate(jobs, reports, work_dir) if jobs else {}
    return built, sims, models


def cuflynx_parameterisation_omex_test(parameterisation_key, tmp_path_factory, request):
    mt, version_name, inst_name = parameterisation_key
    v = load_version(mt, version_name)
    name = f'cuflynx_parameterisation_omex_test__{inst_name}'
    if v.format != 'cellml' and not v.is_supermodule:
        checks.save(v, checks.Result(name, checks.NOT_APPLICABLE, 'C++ version: CUFLynx runs CellML models'))
        pytest.skip('C++ version')
    batch = _batch_of(parameterisation_key, request.session)
    if batch not in _runs:
        work_dir = str(tmp_path_factory.mktemp('omex'))
        try:
            _runs[batch] = (_run_batch(batch, work_dir), work_dir)
        except phlynx.PipelineUnavailable as e:
            _runs[batch] = (e, work_dir)
    run, work_dir = _runs[batch]
    if isinstance(run, phlynx.PipelineUnavailable):
        pytest.skip(f'CUFLynx unavailable: {run}')
    built, sims, models = run
    inst = v.parameterisation(inst_name)
    key = f'{v.id}__{inst_name}'
    problems = []
    if isinstance(built.get(key), Exception):
        problems.append(f'archive not built: {type(built[key]).__name__}: {built[key]}')
    cf = sims.get(key) or {}
    if not problems and not cf.get('ok'):
        problems.append(f'CUFLynx failed at "{cf.get("stage")}": {(cf.get("error") or "")[:400]}')
    imported = cf.get('imported') or {}
    stem = f'{v.vessel_type}_{v.name}_{inst_name}'
    if cf.get('ok'):
        if imported.get('model_filename') not in (f'{stem}.cellml', f'./{stem}.cellml'):
            problems.append(f'CUFLynx took {imported.get("model_filename")!r} as the model, not {stem}.cellml')
        for role, path in (('obs_data', inst.obs_data_path), ('params_for_id', inst.params_for_id_path)):
            got = imported.get(role) or {}
            if path and os.path.isfile(path):
                if got.get('error') or os.path.basename(str(got.get('filename') or '')) != os.path.basename(path):
                    problems.append(f'CUFLynx {role}: took {got.get("filename")!r} (error {got.get("error")!r}), '
                                    f'not {os.path.basename(path)}')
        bad = [k for k, x in (cf.get('outputs') or {}).items() if not np.all(np.isfinite(np.asarray(x, dtype=float)))]
        if bad:
            problems.append(f'non-finite CUFLynx outputs: {bad[:10]}')
    metrics = {'imported': imported, 'n_points': cf.get('n_points'), 'cuflynx_version': cf.get('cuflynx_version')}
    if not problems:
        # the same parameterisation through libcuflynx, compared at the times both logged
        path = models.get(key) or harness.generate(v, str(tmp_path_factory.mktemp(f'lib_{key}')),
                                                   parameters=inst.parameters(), parameterisation=inst)
        lib_t, lib = systems.simulate(path, dict(v.spec, pre_time=0.0), phlynx.SOLVER_INFO)
        wanted = [harness.output_name(o) for o in (v.spec.get('outputs') or [])]
        if any(w not in lib for w in wanted):
            lib_t, extra = systems.simulate(path, dict(v.spec, pre_time=0.0), phlynx.SOLVER_INFO, names=wanted)
            lib = {**extra, **lib}
        rows, unmatched = phlynx.compare(cf, lib_t, lib, phlynx.EQUIVALENCE_TOL)
        metrics.update({'compared': len(rows), 'unmatched': unmatched,
                        'worst': max(rows, key=lambda r: r['difference']) if rows else None})
        if not rows:
            problems.append(f'no CUFLynx output matched libcuflynx\'s ({unmatched[:10]})')
        problems += [f"{r['cuflynx']} differs from libcuflynx by {r['difference']:.3g}" for r in rows if not r['ok']]
    if problems:
        checks.save(v, checks.Result(name, checks.FAILED, problems[0], metrics, details=problems[1:20]))
        ef = v.spec.get('expected_failures') or {}
        known = checks.expected_failure(v, 'cuflynx_parameterisation_omex_test', inst_name)
        if not known and isinstance(built.get(key), Exception) and ef.get('run_test'):
            # the version's model doesn't generate (a known run_test failure), so there is no archive
            known = f'run_test known failure: {ef["run_test"]}'
        if known:
            pytest.xfail(f'known issue: {known} ({problems[0]})')
        pytest.fail('\n'.join(problems[:20]), pytrace=False)
    checks.save(v, checks.Result(name, checks.PASSED,
                                 f'{stem}.omex loads in CUFLynx {cf.get("cuflynx_version")} and reproduces libcuflynx '
                                 f'({metrics.get("compared")} outputs, worst {metrics["worst"]["difference"]:.1e})', metrics))
