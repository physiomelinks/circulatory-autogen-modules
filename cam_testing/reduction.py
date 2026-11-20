"""
Versions reduced from other versions: the annotations that say so, and the test that compares a reduced
version with every version it was reduced from.

A reduced version's verification config names what it was reduced from and how:

    "reduction_method": "symbolic_steady_state",
    "reduced_from": [
      {"module_type": "SGLT1", "version": "kinetic_Ai2025_v01",
       "outputs": {"v_SGLT1": "v_r6"},
       "points": [{"c_Glc_o": 1.0}, {"c_Glc_o": 10.0}],
       "rtol": 1e-6}
    ]

reduction_method   A name in METHODS, or {"name": ..., "description": ...} for a method that is not listed
                   (the description is then required).
reduced_from       One entry for each version it was reduced from, each with how the two are compared:
    module_type, version       the version it was reduced from (the parent)
    parameterisation           the parent's parameterisation the reduction was made from (default: its default)
    reduced_parameterisation   this version's parameterisation that corresponds to it (default: its default)
    compare                    "steady_state": the outputs at the end of each run, both models settled; or
                               "time_course": the outputs at every output time. Default: METHODS' for the
                               method, time_course for an unlisted one
    outputs                    what is compared: {this version's variable: the parent's}, or a list of
                               variables with the same name in both
    parameters                 {this version's parameter: the parent's}, for the parameters of `points` and
                               `sweep` whose names differ
    points                     operating points, [{parameter: value}], each added to the run_parameters of
                               this version (default: one point, the run_parameters alone)
    sweep                      {"ranges": {parameter: [min, max] | {min, max, points, scale} | {values}}}:
                               more points, one parameter at a time from the first point (checks.sweep_values)
    sim_time, pre_time, dt     run settings for both models (default: each version's own; a time_course
                               uses this version's for both)
    rtol, atol                 the two agree when |reduced - parent| <= atol + rtol * scale, scale being
                               |parent| (steady_state) or the largest |parent| of the run (time_course).
                               rtol is required: how close a reduction must be is the modeller's decision
    settled_rtol               steady_state: each output's change over the last tenth of the run, relative
                               to its final value, must be at most this (default 1e-6) plus atol

reduction_test fails when any output at any point of any parent differs by more than its tolerance, and when
a steady_state run has not settled. It is N/A for a version with no reduced_from.
"""
import os

import numpy as np

from cam_testing import checks, harness
from cam_testing.library import load_version

TEST = 'reduction_test'
REDUCTION_TESTS = [TEST]
NOT_REDUCED = 'not reduced from another version (no reduced_from in its verification config)'

STEADY_STATE, TIME_COURSE = 'steady_state', 'time_course'
# method -> (what it is, how a version reduced this way is compared with its parent by default)
METHODS = {
    'symbolic_steady_state': ('the parent\'s equations solved symbolically at steady state (some or all of its '
                              'states at equilibrium), so the reduced version has the parent\'s steady state',
                              STEADY_STATE),
    'operator_inference': ('a low-dimensional model fitted to the parent\'s simulated trajectories (operator '
                           'inference)', TIME_COURSE),
    'manual_term_removal': ('terms of the parent\'s equations removed by hand', TIME_COURSE),
    'manual_term_grouping': ('terms or parameters of the parent\'s equations grouped by hand into fewer',
                             TIME_COURSE),
}
ENTRY_KEYS = {'module_type', 'version', 'parameterisation', 'reduced_parameterisation', 'compare', 'outputs',
              'parameters', 'points', 'sweep', 'sim_time', 'pre_time', 'dt', 'rtol', 'atol', 'settled_rtol'}
DEFAULT_SETTLED_RTOL = 1e-6


def is_reduced(version):
    return bool(version.spec.get('reduced_from'))


def method_name(version):
    method = version.spec.get('reduction_method')
    return method.get('name') if isinstance(method, dict) else method


def method_description(version):
    method = version.spec.get('reduction_method')
    if isinstance(method, dict) and method.get('description'):
        return method['description']
    return (METHODS.get(method_name(version)) or ('',))[0]


def _compare(version, entry):
    return entry.get('compare') or (METHODS.get(method_name(version)) or (None, TIME_COURSE))[1]


def _outputs(entry):
    outputs = entry.get('outputs')
    return dict(outputs) if isinstance(outputs, dict) else {name: name for name in outputs or []}


def parents(version):
    '''The keys (<module_type>/<version>) of the versions the version says it was reduced from.'''
    return [f'{e.get("module_type")}/{e.get("version")}' for e in version.spec.get('reduced_from') or []
            if isinstance(e, dict)]


def annotation_problems(version, load=load_version):
    '''What is wrong with the version's reduction annotations (the structure check and the test report these).'''
    spec = version.spec
    method, entries = spec.get('reduction_method'), spec.get('reduced_from')
    if method is None and entries is None:
        return []
    problems = []
    if not entries or not isinstance(entries, list):
        problems.append('reduction_method without reduced_from: name each version it was reduced from')
        entries = []
    name = method_name(version)
    if not name:
        problems.append('reduced_from without reduction_method: say how it was reduced '
                        f'({", ".join(METHODS)}, or {{"name", "description"}} for another method)')
    elif name not in METHODS and not (isinstance(method, dict) and method.get('description')):
        problems.append(f'reduction_method {name!r} is not one of {", ".join(METHODS)}: give it as '
                        '{"name": ..., "description": ...}')
    for i, entry in enumerate(entries):
        where = f'reduced_from[{i}]'
        if not isinstance(entry, dict) or not entry.get('module_type') or not entry.get('version'):
            problems.append(f'{where}: needs module_type and version')
            continue
        key = f'{entry["module_type"]}/{entry["version"]}'
        where = f'reduced_from {key}'
        unknown = sorted(set(entry) - ENTRY_KEYS)
        if unknown:
            problems.append(f'{where}: unknown key(s) {", ".join(unknown)}')
        if key == version.key:
            problems.append(f'{where}: a version is not reduced from itself')
            continue
        try:
            parent = load(entry['module_type'], entry['version'])
        except Exception as e:
            problems.append(f'{where}: no such version ({type(e).__name__}: {e})')
            continue
        if version.key in _ancestors(parent, load):
            problems.append(f'{where}: it is itself reduced from {version.key} (a cycle)')
        for which, owner in (('parameterisation', parent), ('reduced_parameterisation', version)):
            if entry.get(which) and entry[which] not in owner.parameterisation_names():
                problems.append(f'{where}: {which} {entry[which]!r} is not a parameterisation of '
                                f'{owner.key}')
        if 'rtol' not in entry:
            problems.append(f'{where}: needs rtol, how close the reduced version must be to this one')
        if not _outputs(entry):
            problems.append(f'{where}: needs outputs, the variables compared')
        if _compare(version, entry) not in (STEADY_STATE, TIME_COURSE):
            problems.append(f'{where}: compare is {STEADY_STATE} or {TIME_COURSE}, not {entry.get("compare")!r}')
    return problems


def _ancestors(version, load, seen=()):
    '''Every version ``version`` is reduced from, directly or through others.'''
    found = set()
    for key in parents(version):
        if key in seen or key in found:
            continue
        found.add(key)
        try:
            found |= _ancestors(load(*key.split('/', 1)), load, seen + (version.key,))
        except Exception:      # a missing parent is reported where it is named
            pass
    return found


def _model(version, parameterisation, entry, work_dir, keys=('sim_time', 'pre_time', 'dt'), spec_from=None):
    '''The version's model at ``parameterisation``, with the entry's run settings (else ``spec_from``'s).'''
    cm = checks.ComponentModel(version, work_dir=work_dir,
                               parameterisation=version.parameterisation(parameterisation) if parameterisation else None)
    settings = {k: (spec_from or {}).get(k) for k in keys if (spec_from or {}).get(k) is not None}
    settings.update({k: entry[k] for k in keys if k in entry})
    cm.spec = dict(cm.spec, **settings)
    return cm


def _harness_names(cm, values):
    '''{variable name: value} as {the generated model's parameter name: value}.'''
    by_var = {p.variable_name: p for p in cm.parameters()}
    return {(harness.parameter_name(by_var[k]) if k in by_var else cm.resolve(k)): float(v) for k, v in values.items()}


def _points(cm, entry):
    '''The operating points [{parameter: value}] (this version's variable names), each on top of its
    run_parameters: `points`, then `sweep` one parameter at a time from the first.'''
    base = dict(cm.spec.get('run_parameters') or {})
    points = [dict(base, **p) for p in entry.get('points') or [{}]]
    sweep = entry.get('sweep') or {}
    nominal = cm.param_env(points[0])
    for name in (sweep.get('ranges') or {}):
        values, reason = checks.sweep_values(name, nominal.get(name, 0.0), sweep)
        if reason:
            raise ValueError(f'sweep of {name}: {reason}')
        points += [dict(points[0], **{name: v}) for v in values]
    return points


def _run(cm, point, outputs):
    t, out = cm.run(_harness_names(cm, point))
    missing = [name for name in outputs if harness.output_key(name) not in out]
    if missing:
        raise KeyError(f'{cm.component.key} has no output {", ".join(missing)}: add it to its verification '
                       f'config\'s outputs')
    return np.asarray(t, float), {name: np.asarray(out[harness.output_key(name)], float) for name in outputs}


def _unsettled(t, y, settled_rtol, atol):
    '''The trace's change over the last tenth of the run when it is more than a settled trace's, else None.'''
    late = y[t >= t[0] + 0.9 * (t[-1] - t[0])]
    change = float(np.max(np.abs(late - y[-1]))) if len(late) else 0.0
    return change if change > atol + settled_rtol * abs(float(y[-1])) else None


def compare_parent(version, entry, work_dir):
    '''Compares the version with one parent: (problems, metrics, the points' rows for the plot).'''
    parent = load_version(entry['module_type'], entry['version'])
    compare = _compare(version, entry)
    outputs = _outputs(entry)
    rtol, atol = float(entry['rtol']), float(entry.get('atol', 0.0))
    settled_rtol = float(entry.get('settled_rtol', DEFAULT_SETTLED_RTOL))
    own = _model(version, entry.get('reduced_parameterisation'), entry, os.path.join(work_dir, 'reduced'))
    # a time course needs the two on one time grid: the parent takes this version's run settings
    theirs = _model(parent, entry.get('parameterisation'), entry, os.path.join(work_dir, 'parent'),
                    spec_from=own.spec if compare == TIME_COURSE else None)
    names = entry.get('parameters') or {}
    problems, rows, worst = [], [], (0.0, None)
    for point in _points(own, entry):
        t_own, y_own = _run(own, point, list(outputs))
        t_par, y_par = _run(theirs, {names.get(k, k): v for k, v in point.items()}, list(outputs.values()))
        label = ', '.join(f'{k}={v:g}' for k, v in point.items()) or 'default parameters'
        row = {'point': point, 'outputs': {}}
        for name, parent_name in outputs.items():
            a, b = y_own[name], y_par[parent_name]
            if compare == STEADY_STATE:
                for who, t, y in ((f'{name}', t_own, a), (f'the parent\'s {parent_name}', t_par, b)):
                    change = _unsettled(t, y, settled_rtol, atol)
                    if change is not None:
                        problems.append(f'{label}: {who} has not settled (it still changes by {change:.3g} over '
                                        f'the last tenth of the run; final value {y[-1]:.6g}): run for longer')
                a, b = a[-1:], b[-1:]
            elif len(a) != len(b):
                problems.append(f'{label}: {name} has {len(a)} output times and the parent\'s {parent_name} '
                                f'{len(b)}: give sim_time and dt in the reduced_from entry')
                continue
            difference = float(np.max(np.abs(a - b)))
            scale = float(np.max(np.abs(b)))
            allowed = atol + rtol * scale
            relative = difference / scale if scale else (0.0 if difference == 0 else float('inf'))
            row['outputs'][name] = {'reduced': float(a[-1]), 'parent': float(b[-1]), 'difference': difference,
                                    'relative': relative}
            if relative > worst[0]:
                worst = (relative, f'{name} at {label}')
            if not difference <= allowed:
                problems.append(f'{label}: {name} differs from the parent\'s {parent_name} by {difference:.3g} '
                                f'({relative:.3g} of {scale:.3g}); allowed {allowed:.3g} (rtol {rtol:g}, atol {atol:g})')
        rows.append(row)
    metrics = {'parent': parent.key, 'compare': compare, 'rtol': rtol, 'atol': atol,
               'points': len(rows), 'worst_relative_difference': worst[0], 'worst_at': worst[1], 'rows': rows}
    return problems, metrics


def _plot(version, parent_metrics):
    '''One figure: each compared output of each parent against the point number, reduced and parent.'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    panels = [(m['parent'], name) for m in parent_metrics for name in (m['rows'][0]['outputs'] if m['rows'] else [])]
    if not panels:
        return []
    fig, axes = plt.subplots(len(panels), 1, figsize=(7, 2.6 * len(panels)), squeeze=False)
    for ax, (parent, name) in zip(axes[:, 0], panels):
        rows = next(m for m in parent_metrics if m['parent'] == parent)['rows']
        x = range(1, len(rows) + 1)
        ax.plot(x, [r['outputs'][name]['parent'] for r in rows if name in r['outputs']], 'o-', label=parent)
        ax.plot(x, [r['outputs'][name]['reduced'] for r in rows if name in r['outputs']], 'x--', label='this version')
        ax.set_ylabel(name)
        ax.set_xlabel('operating point')
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = checks.plot_path(version, 'reduction')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return [path]


def reduction_check(version, work_dir):
    '''reduction_test: the version against every version it was reduced from.'''
    if not is_reduced(version) and version.spec.get('reduction_method') is None:
        return checks.save(version, checks.Result(TEST, checks.NOT_APPLICABLE, NOT_REDUCED))

    def check():
        problems = annotation_problems(version)
        if problems:
            return checks.Result(TEST, checks.FAILED, f'reduction annotations: {problems[0]}', details=problems)
        all_metrics = []
        for i, entry in enumerate(version.spec['reduced_from']):
            found, metrics = compare_parent(version, entry, os.path.join(work_dir, f'parent_{i}'))
            problems += [f'{metrics["parent"]}: {p}' for p in found]
            all_metrics.append(metrics)
        metrics = {'method': method_name(version), 'method_description': method_description(version),
                   'parents': all_metrics}
        figures = _plot(version, all_metrics)
        summary = '; '.join(f'{m["parent"]} ({m["compare"].replace("_", " ")}, {m["points"]} point(s), largest '
                            f'relative difference {m["worst_relative_difference"]:.3g}, rtol {m["rtol"]:g})'
                            for m in all_metrics)
        if problems:
            return checks.Result(TEST, checks.FAILED, f'{len(problems)} difference(s) from what it was reduced '
                                 f'from: {problems[0]}', metrics, figures, problems)
        return checks.Result(TEST, checks.PASSED, f'reduced by {method_name(version)} from {summary}', metrics, figures)
    return checks._guard(version, TEST, check)
