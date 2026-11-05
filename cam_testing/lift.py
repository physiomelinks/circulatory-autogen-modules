"""
Lifts a version's model to quadratic form for Operator Inference (libcuflynx.utilities.lift_cellml).

    python -m cam_testing.lift <module_type>/<version> [--system-model MODEL.cellml] [--output-dir DIR]

1. The version is generated in its harness (its verification_config's network, else the version alone), or
   ``--system-model`` names a generated model to use instead (its flat CellML must sit next to it).
2. It is smoothed with the spec's smoothing settings and exemptions (verification_test_smoothing, whose current
   saved result is reused). A model with no non-smooth sites is lifted as it is.
3. It is lifted: the inputs are the version's ``lifting.inputs`` (verification_config), else its
   ``boundary_condition`` variables left unconnected in the generated model; time functions and exempt sites
   become inputs as well.
4. The lifted model is simulated against the smoothed one and checked (J f = lifted right-hand side, operator
   residual, distance from the manifold).

Writes results/lifting/ (the lifted CellML, lifting.json, lifting_functions.py) and results/lifting.json (the
result, as the other checks write theirs), or into --output-dir. Informational: not a test column.
"""
import argparse
import contextlib
import io
import json
import os
import sys
import time

from cam_testing import checks
from cam_testing.library import version_index

LIFTING = 'lifting'
#: lifting (verification_config): inputs, the variables kept as inputs u (default: the version's boundary_condition
#: variables); time_budget [s] for each QBee stage; max_aux; max_terms, the largest polynomial system to quadratize;
#: sim_time and dt of the verification run (default: the version's sim_time and dt); tol, the largest state error and
#: manifold defect of the lifted run for "passed"; rationale.
DEFAULT_LIFTING = {'time_budget': 60.0, 'max_aux': 400, 'max_terms': 20000, 'tol': 1e-4, 'algebraic_tol': 1e-8}
#: status of a lifting whose algebra checks out but whose lifted run leaves the manifold beyond tol
UNSTABLE_MANIFOLD = 'unstable_manifold'


def lifting_settings(spec, **overrides):
    settings = dict(DEFAULT_LIFTING)
    settings.update({k: v for k, v in (spec.get('lifting') or {}).items() if k != 'rationale'})
    settings.update({k: v for k, v in overrides.items() if v is not None})
    return settings


def resolve(key, index=None):
    """The Version of '<module_type>/<version>'."""
    index = version_index(include_libraries=True) if index is None else index
    module_type, _, version = key.partition('/')
    found = index.get((module_type, version))
    if found is None:
        raise KeyError(f'no version {key}')
    return found


def default_inputs(version, flat_path):
    """The version's boundary_condition variables that are constants of the generated model (no port connects
    them), as '<component>/<variable>' of the version's component."""
    import xml.etree.ElementTree as ET
    from libcuflynx.utilities import smooth_cellml as sc
    names = [v[0] for v in (version.config.get('variables_and_units') or []) if len(v) > 3 and v[3] == 'boundary_condition']
    if not names:
        return []
    root = ET.parse(flat_path).getroot()
    graph = sc._ModelGraph(root)
    out = []
    component = f'{_vessel()}_module'
    for name in names:
        if (component, name) in {(c.get('name'), v.get('name')) for c in root.findall(sc._C + 'component')
                                 for v in c.findall(sc._C + 'variable')} and graph.is_constant(component, name):
            out.append(f'{component}/{name}')
    return out


def _vessel():
    from cam_testing.harness import VESSEL
    return VESSEL


def smoothed_model(version, cm, force=False):
    """(path of the model to lift, smoothing report or None, how it was smoothed)."""
    from cam_testing.smooth import library_check
    result, reused = library_check(version, force=force)
    smoothed = result.metrics.get('smoothed_model')
    report = result.metrics.get('report')
    if smoothed:
        smoothed = smoothed if os.path.isabs(smoothed) else os.path.join(version.dir, smoothed)
        report = report if (report is None or os.path.isabs(report)) else os.path.join(version.dir, report)
        return smoothed, report, {'status': result.status, 'message': result.message, 'reused': reused}
    return checks.flat_model_path(cm.model_path), report if report and os.path.isfile(report) else None, \
        {'status': result.status, 'message': result.message, 'reused': reused, 'note': 'no non-smooth sites'}


def smooth_given_model(model_path, spec, settings, work_dir, out_dir):
    """Smooths a given generated model (``--system-model``) with the version's settings, compared at the nominal
    point only. Returns (path to lift, report path or None, how)."""
    from cam_testing import harness
    flat = checks.flat_model_path(model_path)
    with contextlib.redirect_stdout(io.StringIO()):
        helper = harness.simulation_helper(model_path, spec)
    outputs = harness.state_outputs(helper)
    status, message, metrics, _figs, _details = checks.compare_smoothed(
        model_path, spec, outputs, [{'label': 'nominal', 'params': {}}], {}, settings, out_dir, work_dir, 'system')
    if metrics.get('smoothed_model'):
        return metrics['smoothed_model'], metrics.get('report'), {'status': status, 'message': message}
    return flat, metrics.get('report') if os.path.isfile(metrics.get('report') or '') else None, \
        {'status': status, 'message': message, 'note': 'no non-smooth sites'}


def run(key, system_model=None, output_dir=None, settings_override=None, force=False, inputs=None):
    """The pipeline above; returns (Result, lifting.json dict)."""
    import tempfile
    from libcuflynx.utilities import lift_cellml as lc
    t0 = time.perf_counter()
    version = resolve(key)
    spec = version.spec
    settings = lifting_settings(spec, **(settings_override or {}))
    out_dir = output_dir or os.path.join(version.results_dir, LIFTING)
    os.makedirs(out_dir, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix=f'cam_lift_{version.id}_')
    smoothing = checks.smoothing_settings(spec)
    if system_model:
        model, report, how = smooth_given_model(system_model, spec, smoothing, work_dir,
                                                os.path.join(out_dir, 'smoothing'))
        flat = checks.flat_model_path(system_model)
    else:
        cm = checks.ComponentModel(version, work_dir=work_dir)
        model, report, how = smoothed_model(version, cm, force=force)
        flat = checks.flat_model_path(cm.model_path)
    t_smooth = time.perf_counter() - t0
    if inputs is None:
        inputs = settings.get('inputs')
    if inputs is None:
        inputs = default_inputs(version, flat)
    sim_time = float(settings.get('sim_time', spec.get('sim_time', 10.0)))
    dt = float(settings.get('dt', spec.get('dt', 0.01)))
    try:
        doc = lc.lift_cellml_file(model, out_dir, inputs=inputs, smoothing_report=report,
                                  time_budget=float(settings['time_budget']), max_aux=int(settings['max_aux']),
                                  sim_time=sim_time, dt=dt, max_terms=int(settings['max_terms']))
    except Exception as e:  # noqa: BLE001 -- a failed lifting is a result, not a crash
        result = checks.Result(LIFTING, checks.FAILED, f'{type(e).__name__}: {e}',
                               {'smoothing': how, 'inputs': inputs, 'time_s': time.perf_counter() - t0})
        if output_dir is None:
            checks.save(version, result)
        else:
            _write(os.path.join(out_dir, 'lifting_result.json'), result)
        return result, None
    v, c = doc['verification'], doc['check']
    algebraic = max(v['jf']['max'], c['max_relative_residual'])
    run_err = max(v['max_state_error'], v['max_manifold_defect'])
    issues = doc['libcellml']['validator'] + doc['libcellml']['analyser']
    present = ', '.join(k for k, on in doc['operators']['present'].items() if on)
    what = (f'{doc["size"]["n"]} -> {doc["size"]["n_w"]} states, {len(doc["inputs"])} inputs, operators {present}; '
            f'engine {doc["engine"]["polynomialize"]}/{doc["engine"]["quadratize"]}')
    if issues or not algebraic <= settings['algebraic_tol']:
        status = checks.FAILED
        message = (f'libCellML: {issues[0]}' if issues else
                   f'the lifted right-hand side differs from J f by {algebraic:.3g}') + f' [{what}]'
    elif run_err <= settings['tol']:
        status = checks.PASSED
        message = f'exact: J f and operators within {algebraic:.2g}; lifted run within {run_err:.2g} [{what}]'
    else:
        growth = v.get('defect_growth', {})
        leave = (growth.get('defect') or growth.get('state_error') or {}).get('1e-06')
        rate = growth.get('rate')
        status = UNSTABLE_MANIFOLD
        message = (f'exact (J f and operators within {algebraic:.2g}), but the lifted run leaves the manifold: '
                   f'error {run_err:.3g} over {sim_time:g}' + (f', past 1e-6 at t = {leave:g}' if leave is not None else '')
                   + (f', defect growth rate {rate:.3g}' if rate is not None else '') + f' [{what}]')
    metrics = {'size': doc['size'], 'engine': doc['engine'], 'inputs': [i['name'] for i in doc['inputs']],
               'operators': doc['operators']['present'], 'verification': {k: v[k] for k in (
                   'max_state_error', 'max_manifold_defect', 'defect_growth')}, 'jf': v['jf']['max'],
               'operator_residual': c['max_relative_residual'], 'smoothing': how, 'settings': settings,
               'timings': dict(doc['timings'], smooth=t_smooth, total=time.perf_counter() - t0),
               'lifting_json': os.path.join(out_dir, 'lifting.json'),
               'rationale': (spec.get('lifting') or {}).get('rationale', '')}
    result = checks.Result(LIFTING, status, message, metrics)
    if output_dir is None:
        metrics['lifting_json'] = os.path.relpath(metrics['lifting_json'], version.dir)
        checks.save(version, result)
    else:
        _write(os.path.join(out_dir, 'lifting_result.json'), result)
    return result, doc


def _write(path, result):
    from dataclasses import asdict
    with open(path, 'w') as f:
        json.dump(asdict(result), f, indent=2, default=checks._json_default)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('version', help='<module_type>/<version>')
    parser.add_argument('--system-model', help='a generated model to lift instead of the version in its harness')
    parser.add_argument('--output-dir', help='where to write (default: the version\'s results/lifting/)')
    parser.add_argument('--inputs', help='comma-separated variables kept as inputs ([component/]name)')
    parser.add_argument('--time-budget', type=float, help='s for each QBee stage')
    parser.add_argument('--max-aux', type=int)
    parser.add_argument('--max-terms', type=int, help='largest polynomial system (terms) to quadratize')
    parser.add_argument('--sim-time', type=float, help='verification run length')
    parser.add_argument('--force', action='store_true', help='rerun the smoothing check even when current')
    args = parser.parse_args(argv)
    override = {'time_budget': args.time_budget, 'max_aux': args.max_aux, 'max_terms': args.max_terms,
                'sim_time': args.sim_time}
    inputs = [s for s in args.inputs.split(',') if s] if args.inputs is not None else None
    result, _doc = run(args.version, args.system_model, args.output_dir, override, args.force, inputs)
    print(f'{args.version}: {result.status}: {result.message}')
    return 0 if result.status in (checks.PASSED, UNSTABLE_MANIFOLD) else 1


if __name__ == '__main__':
    sys.exit(main())
