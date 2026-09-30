"""
The PhLynx -> CUFLynx pipeline for one module: can each component be built in PhLynx from this
library, exported as PhLynx's CUFLynx .omex, imported into CUFLynx and simulated there, and
does CUFLynx's run reproduce the model libcuflynx builds from the same test network?

Two bridges, each running the real application code from a checkout:
  tools/phlynx_bridge/export_omex.mjs    PhLynx's JS under node + jsdom (PHLYNX_DIR)
  tools/cuflynx_bridge/simulate_omex.py  a released CUFLynx binary (CUFLYNX_BIN, default
                                         ~/software/CUFLynx; CI downloads the latest Ubuntu
                                         release) driven over its HTTP API

A module's components are exported in one node run and simulated in one CUFLynx run
(run_module), and the per-component results are cached in the module's results dir for the
three tests (phlynx_export_test, cuflynx_simulate_test, phlynx_equivalence_test).
"""
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from cam_testing import harness, vessel_array
from cam_testing.library import MODULES_DIR

REPO_DIR = os.path.dirname(MODULES_DIR)
PHLYNX_BRIDGE = os.path.join(REPO_DIR, 'tools', 'phlynx_bridge', 'export_omex.mjs')
CUFLYNX_BRIDGE = os.path.join(REPO_DIR, 'tools', 'cuflynx_bridge', 'simulate_omex.py')
SIBLINGS = os.path.dirname(REPO_DIR)


class PipelineUnavailable(Exception):
    '''PhLynx or CUFLynx (or node) isn't available here: the tests skip, they don't fail.'''


def phlynx_dir():
    d = os.environ.get('PHLYNX_DIR') or os.path.join(SIBLINGS, 'phlynx')
    if not os.path.isfile(os.path.join(d, 'src', 'utils', 'cellml.js')):
        raise PipelineUnavailable(f'no PhLynx checkout at {d} (set PHLYNX_DIR)')
    if not os.path.isdir(os.path.join(d, 'node_modules', 'libcellml.js')):
        raise PipelineUnavailable(f'PhLynx at {d} has no node_modules (run yarn install there)')
    return d


def cuflynx_bin():
    b = os.environ.get('CUFLYNX_BIN') or os.path.expanduser('~/software/CUFLynx')
    if not os.access(b, os.X_OK):
        raise PipelineUnavailable(f'no CUFLynx binary at {b} (set CUFLYNX_BIN)')
    return b


def library_files():
    '''Every module's CellML, units and config: a harness can use neighbours from other modules.'''
    cellml, units, configs = [], [], []
    from cam_testing.library import module_dir, module_names
    for name in module_names():
        d = module_dir(name)
        cellml += sorted(glob.glob(os.path.join(d, f'{name}_modules.cellml')))
        units += sorted(glob.glob(os.path.join(d, f'{name}_units.cellml')))
        configs += sorted(glob.glob(os.path.join(d, f'{name}_modules_config.json')))
    return {'cellml': cellml, 'units': units, 'configs': configs}


def _read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def component_job(component, work_dir):
    '''The component's test network (harness) as a PhLynx instance array plus its parameters.'''
    prefix = f'{component.module.name}__{component.id}'
    res = os.path.join(work_dir, 'resources', component.id)
    harness._write_resources(component, res, prefix, {})
    records = vessel_array.read_records(vessel_array.find(res, prefix))
    # inp/out space-separated, as PhLynx's instance-array importer gives them
    instances = [{'name': r['name'], 'module_type': r['module_type'], 'module_subtype': r['module_subtype'],
                  'inp_instances': ' '.join(r.get('inp_instances', [])),
                  'out_instances': ' '.join(r.get('out_instances', []))} for r in records]
    params = [{k: (v or '') for k, v in r.items()} for r in _read_csv(os.path.join(res, f'{prefix}_parameters.csv'))]
    spec = component.spec
    extra = [harness.output_name(o) for o in (spec.get('outputs') or [])]
    return {'id': component.id, 'instances': instances, 'parameters': params, 'extra_outputs': extra,
            'sim_time': float(spec.get('sim_time', 1.0)), 'dt': float(spec.get('dt', 0.01)),
            'resources': res, 'prefix': prefix}


def _node():
    node = shutil.which('node')
    if not node:
        raise PipelineUnavailable('node is not installed (PhLynx bridge needs node >= 22.15)')
    return node


def export(jobs, work_dir, library=None):
    '''Runs the PhLynx bridge once per job; returns {id: report} and writes <work_dir>/<id>.omex.'''
    pdir = phlynx_dir()
    node = _node()
    reports = {}
    lib = library or library_files()
    for job in jobs:
        job_path = os.path.join(work_dir, f'{job["id"]}.phlynx_job.json')
        omex = os.path.join(work_dir, f'{job["id"]}.omex')
        with open(job_path, 'w') as f:
            json.dump(dict(lib, **{k: job[k] for k in ('instances', 'parameters', 'sim_time', 'dt')}), f)
        p = subprocess.run([node, PHLYNX_BRIDGE, job_path, omex], capture_output=True, text=True,
                           env=dict(os.environ, PHLYNX_DIR=pdir), timeout=600)
        if p.returncode == 3:
            raise PipelineUnavailable(p.stderr.strip()[:500])
        try:
            rep = json.loads(p.stdout.strip().splitlines()[-1])
        except (ValueError, IndexError):
            rep = {'ok': False, 'stage': 'bridge', 'error': (p.stderr or p.stdout)[-3000:]}
        rep['omex'] = omex if rep.get('ok') else None
        reports[job['id']] = rep
    return reports


def simulate(jobs, reports, work_dir, solver_info=None):
    '''Runs the CUFLynx bridge on every exported .omex; returns {id: result}.'''
    binary = cuflynx_bin()
    models = [{'id': j['id'], 'omex': reports[j['id']]['omex'], 'sim_time': j['sim_time'], 'dt': j['dt'],
               'extra_outputs': j.get('extra_outputs') or [],
               'out': os.path.join(work_dir, f'{j["id"]}.cuflynx.json')}
              for j in jobs if reports.get(j['id'], {}).get('omex')]
    if not models:
        return {}
    job_path = os.path.join(work_dir, 'cuflynx_job.json')
    with open(job_path, 'w') as f:
        json.dump({'models': models, 'solver_info': solver_info or SOLVER_INFO}, f)
    p = subprocess.run([sys.executable, CUFLYNX_BRIDGE, job_path], capture_output=True, text=True,
                       env=dict(os.environ, CUFLYNX_BIN=binary), timeout=3600, cwd=work_dir)
    if p.returncode == 3:
        raise PipelineUnavailable(p.stderr.strip()[:500])
    out = {}
    for m in models:
        try:
            with open(m['out']) as f:
                out[m['id']] = json.load(f)
        except (OSError, ValueError):
            out[m['id']] = {'ok': False, 'stage': 'bridge', 'error': (p.stderr or p.stdout)[-3000:]}
    return out


def match_outputs(cuflynx_outputs, lib_outputs):
    '''
    Pairs CUFLynx's outputs (the PhLynx model's names, 'component/variable') with libcuflynx's
    ('vessel/variable'): PhLynx names a node's component after the instance, so the pair is the
    same instance and variable when the component name contains the vessel name.
    '''
    pairs = {}
    for name in cuflynx_outputs:
        comp, _, var = name.rpartition('/')
        cands = [k for k in lib_outputs if k.split('/', 1)[-1] == var
                 and (comp == k.split('/', 1)[0] or comp.startswith(k.split('/', 1)[0] + '_')
                      or comp.endswith('_' + k.split('/', 1)[0]))]
        if len(cands) == 1:
            pairs[name] = cands[0]
    return pairs


def compare(cf, lib_t, lib_outputs, tol):
    '''
    Normalised max differences between CUFLynx's run and libcuflynx's on the shared outputs, at
    the output times both runs have (the two apps log at different intervals; interpolating a
    coarse log through a fast transient would invent differences).
    '''
    t = np.asarray(cf.get('time') or [], dtype=float)
    lib_t = np.asarray(lib_t, dtype=float)
    pairs = match_outputs(cf.get('outputs') or {}, lib_outputs)
    # indices of common times
    j = np.searchsorted(lib_t, t)
    j = np.clip(j, 0, max(len(lib_t) - 1, 0))
    common = np.abs(lib_t[j] - t) <= 1e-9 * max(1.0, float(np.max(np.abs(t))) if t.size else 1.0) if lib_t.size else np.zeros(t.shape, bool)
    ci_idx, li_idx = np.nonzero(common)[0], j[common]
    rows = []
    for c_name, l_name in pairs.items():
        c = np.asarray(cf['outputs'][c_name], dtype=float).ravel()
        if c.size == 1 and t.size > 1:   # CUFLynx returns a constant output as one value
            c = np.full(t.shape, c[0])
        l_vals = np.asarray(lib_outputs[l_name], dtype=float).ravel()
        if c.size != t.size or l_vals.size != lib_t.size or len(ci_idx) < 2:
            continue
        cv, lv = c[ci_idx], l_vals[li_idx]
        scale = max(float(np.max(np.abs(lv))), float(np.ptp(lv)), 1e-300)
        d = float(np.max(np.abs(cv - lv)) / scale)
        rows.append({'cuflynx': c_name, 'libcuflynx': l_name, 'difference': d, 'ok': d <= tol,
                     'common_times': int(len(ci_idx))})
    unmatched = sorted(set(cf.get('outputs') or {}) - set(pairs))
    return rows, unmatched


def cache_path(module):
    return os.path.join(module.dir, 'results', 'phlynx_pipeline.json')


def run_module(module, components, solver_info=None, keep_dir=None):
    '''Exports and simulates every component of a module; caches and returns {component id: result}.'''
    work_dir = keep_dir or tempfile.mkdtemp(prefix=f'cam_phlynx_{module.name}_')
    os.makedirs(work_dir, exist_ok=True)
    jobs, results = [], {}
    for c in components:
        if _cpp(c):
            continue
        try:
            jobs.append(component_job(c, work_dir))
        except Exception as e:  # noqa: BLE001 - e.g. parameters still TODO
            results[c.id] = {'export': {'ok': False, 'stage': 'job', 'error': f'{type(e).__name__}: {e}'}}
    reports = export(jobs, work_dir)
    sims = simulate(jobs, reports, work_dir, solver_info)
    for j in jobs:
        results[j['id']] = {'export': reports.get(j['id']), 'cuflynx': sims.get(j['id']), 'job': {
            k: j[k] for k in ('instances', 'sim_time', 'dt')}}
    os.makedirs(os.path.dirname(cache_path(module)), exist_ok=True)
    with open(cache_path(module), 'w') as f:
        json.dump(results, f, default=float)
    return results


# ---- the three pipeline checks, recorded like the other tests (results/<component>/<test>.json) ----

PIPELINE_TESTS = ['phlynx_export_test', 'cuflynx_simulate_test', 'phlynx_equivalence_test']
EQUIVALENCE_TOL = 1e-6
SOLVER_INFO = {'rtol': 1e-10, 'atol': 1e-12}


CPP_REASON = ('C++ 1D-solver component (module_format cpp): PhLynx builds CellML models only, so it has no '
              'PhLynx/CUFLynx form')


def _cpp(component):
    return component.config.get('module_format', 'cellml') == 'cpp'


def export_check(component, res):
    from cam_testing.checks import FAILED, NOT_APPLICABLE, PASSED, Result
    if _cpp(component):
        return Result('phlynx_export_test', NOT_APPLICABLE, CPP_REASON)
    rep = (res or {}).get('export') or {}
    problems = []
    if not rep.get('ok'):
        problems.append(f'PhLynx failed at "{rep.get("stage")}": {(rep.get("error") or "").splitlines()[0:3]}')
    if rep.get('missing_modules'):
        problems.append(f'modules PhLynx could not find: {rep["missing_modules"]}')
    if rep.get('stub_modules'):
        problems.append(f'modules with no math in PhLynx: {rep["stub_modules"]}')
    if rep.get('ok') and rep.get('edges') != rep.get('expected_edges'):
        problems.append(f'PhLynx made {rep.get("edges")} of the {rep.get("expected_edges")} connections')
    problems += [f'refused: {w}' for w in rep.get('refused_edges') or []]
    problems += [f'unsupported: {w}' for w in rep.get('multiport_warnings') or []]
    metrics = {k: rep.get(k) for k in ('nodes', 'edges', 'expected_edges', 'parameters_applied', 'configs_loaded',
                                        'omex_members', 'flattened_bytes')}
    if problems:
        return Result('phlynx_export_test', FAILED, problems[0], metrics, details=problems[1:] + (rep.get('warnings') or [])[:20])
    return Result('phlynx_export_test', PASSED,
                  f'built in PhLynx ({rep.get("nodes")} nodes, all {rep.get("edges")} connections) and exported as '
                  f'the .omex PhLynx sends to CUFLynx', metrics)


def simulate_check(component, res):
    from cam_testing.checks import FAILED, NOT_APPLICABLE, PASSED, SKIPPED, Result
    if _cpp(component):
        return Result('cuflynx_simulate_test', NOT_APPLICABLE, CPP_REASON)
    if not ((res or {}).get('export') or {}).get('ok'):
        return Result('cuflynx_simulate_test', SKIPPED, 'no .omex: the PhLynx export failed')
    cf = (res or {}).get('cuflynx') or {}
    metrics = {k: cf.get(k) for k in ('n_points', 'cuflynx_version', 'import_warnings')}
    metrics['outputs'] = sorted(cf.get('outputs') or {})
    if not cf.get('ok'):
        return Result('cuflynx_simulate_test', FAILED,
                      f'CUFLynx failed at "{cf.get("stage")}": {(cf.get("error") or "")[:500]}', metrics)
    bad = [k for k, v in cf['outputs'].items() if not np.all(np.isfinite(np.asarray(v, dtype=float)))]
    if not cf['outputs']:
        return Result('cuflynx_simulate_test', FAILED, 'CUFLynx ran the model but returned no outputs '
                      '(no states, and none of the component\'s outputs resolved)', metrics)
    if bad:
        return Result('cuflynx_simulate_test', FAILED, f'non-finite outputs in CUFLynx: {bad[:10]}', metrics)
    return Result('cuflynx_simulate_test', PASSED,
                  f'imported into CUFLynx {cf.get("cuflynx_version")} and simulated '
                  f'({cf.get("n_points")} points, {len(cf["outputs"])} outputs finite)', metrics)


def equivalence_check(component, res, work_dir):
    '''CUFLynx's run of the PhLynx-built model against libcuflynx's model of the same test network.'''
    from cam_testing import system as systems
    from cam_testing.checks import FAILED, NOT_APPLICABLE, PASSED, SKIPPED, Result
    if _cpp(component):
        return Result('phlynx_equivalence_test', NOT_APPLICABLE, CPP_REASON)
    cf = (res or {}).get('cuflynx') or {}
    if not cf.get('ok'):
        return Result('phlynx_equivalence_test', SKIPPED, 'no CUFLynx run to compare')
    path = harness.generate(component, work_dir)
    spec = dict(component.spec, pre_time=0.0)
    lib_t, lib = systems.simulate(path, spec, SOLVER_INFO)
    wanted = [harness.output_name(o) for o in (component.spec.get('outputs') or [])]
    if any(w not in lib for w in wanted):
        # the listing misses some variables (e.g. an algebraic-only model); ask for them by name
        lib_t, extra = systems.simulate(path, spec, SOLVER_INFO, names=wanted)
        lib = {**extra, **lib}
    rows, unmatched = compare(cf, lib_t, lib, EQUIVALENCE_TOL)
    bad = [r for r in rows if not r['ok']]
    worst = max(rows, key=lambda r: r['difference']) if rows else None
    metrics = {'compared': len(rows), 'unmatched': unmatched, 'tol': EQUIVALENCE_TOL, 'worst': worst,
               'rows': sorted(rows, key=lambda r: -r['difference'])}
    if not rows:
        return Result('phlynx_equivalence_test', FAILED,
                      f'no CUFLynx output could be matched to libcuflynx\'s ({unmatched[:10]})', metrics)
    if bad:
        return Result('phlynx_equivalence_test', FAILED,
                      f'{len(bad)} of {len(rows)} outputs differ from libcuflynx\'s model beyond {EQUIVALENCE_TOL}',
                      metrics, details=[f"{r['cuflynx']} vs {r['libcuflynx']}: {r['difference']:.3g}" for r in bad[:20]])
    return Result('phlynx_equivalence_test', PASSED,
                  f'all {len(rows)} compared outputs match libcuflynx\'s model (worst {worst["cuflynx"]} '
                  f'{worst["difference"]:.1e} <= {EQUIVALENCE_TOL})', metrics)
