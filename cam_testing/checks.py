"""
The standard verification & validation tests, run on one module component at a time.

Each check returns a Result and writes it to modules/<name>/results/<component>/<test>.json,
with any figures in modules/<name>/plots/. The pytest layer (tests/test_modules.py) turns a
Result into pass / fail / skip; the report generator reads the JSON files, so the HTML shows
exactly what the last test run found.
"""
import contextlib
import datetime
import io
import json
import math
import os
import tempfile
import traceback
from dataclasses import asdict, dataclass, field

import numpy as np

from cam_testing import fixed_step, harness, plots, ranges
from cam_testing.library import REPO_ROOT

PASSED, FAILED, SKIPPED, PENDING = 'passed', 'failed', 'skipped', 'pending'

TESTS = ['run_test', 'verification_test_BC', 'verification_test_timestep', 'stability_test',
         'validation_test_baseline', 'validation_test_calibrate']


@dataclass
class Result:
    test: str
    status: str
    message: str = ''
    metrics: dict = field(default_factory=dict)
    plots: list = field(default_factory=list)
    details: list = field(default_factory=list)
    timestamp: str = ''


def result_path(component, test):
    return os.path.join(component.module.results_dir, component.id, f'{test}.json')


def plot_path(component, name):
    return os.path.join(component.module.plots_dir, f'{component.id}__{name}.png')


def save(component, result):
    result.timestamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    # plots are stored relative to the module dir, which is where the HTML lives
    result.plots = [os.path.relpath(p, component.module.dir) for p in result.plots]
    path = result_path(component, result.test)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(asdict(result), f, indent=2, default=_json_default)
    return result


def load(component, test):
    path = result_path(component, test)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return Result(**json.load(f))


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


class ComponentModel(object):
    '''The component's generated model and simulation helper, built once and reused.'''

    def __init__(self, component, work_dir=None):
        self.component = component
        self.spec = component.spec
        self.work_dir = work_dir or tempfile.mkdtemp(prefix=f'cam_{component.module.name}_')
        self._model_path = None
        self._helper = None
        self.nominal = {harness.parameter_name(p): p.float_value
                        for p in component.parameters() if not p.is_todo}
        self._var_of = {harness.parameter_name(p): p.variable_name for p in component.parameters()}

    def param_env(self, overrides=None):
        '''Parameter values by variable name (e.g. alpha), for invariant expressions.'''
        values = dict(self.nominal)
        values.update(overrides or {})
        return {self._var_of[k]: v for k, v in values.items()}

    @property
    def model_path(self):
        if self._model_path is None:
            self._model_path = harness.generate(self.component, self.work_dir)
        return self._model_path

    @property
    def helper(self):
        if self._helper is None:
            self._helper = harness.simulation_helper(self.model_path, self.spec)
        return self._helper

    def outputs(self):
        wanted = self.spec.get('outputs') or [v[0] for v in self.component.variables()]
        return list(wanted)

    def units(self):
        return {v[0]: v[1] for v in self.component.config['variables_and_units']}

    def run(self, params=None):
        '''Runs at nominal parameters plus ``params``; always restores the nominal values.'''
        params = params or {}
        try:
            return harness.run(self.helper, self.outputs(), params={'parameters/' + k: v for k, v in params.items()})
        finally:
            if params:
                self.helper.set_param_vals(['parameters/' + k for k in params],
                                           [self.nominal[k] for k in params])


SAFE_BUILTINS = {'abs': abs, 'max': max, 'min': min, 'len': len, 'float': float}


def period(t, y, skip=1):
    '''
    Mean period of an oscillating trace from its upward crossings of its mid-range level,
    ignoring the first ``skip`` cycles (transient). NaN when there are too few cycles.
    '''
    t, y = np.asarray(t, float), np.asarray(y, float)
    level = 0.5 * (np.max(y) + np.min(y))
    up = np.where((y[:-1] < level) & (y[1:] >= level))[0]
    # linear interpolation of the crossing times
    tc = t[up] + (level - y[up]) * (t[up + 1] - t[up]) / (y[up + 1] - y[up])
    tc = tc[skip:]
    return float(np.mean(np.diff(tc))) if len(tc) >= 2 else float('nan')


INVARIANT_HELPERS = {'period': period}


def _check_invariants(spec, t, outputs, params=None, where='sweep'):
    '''
    Returns a list of failed invariant descriptions. An invariant is a numpy expression of
    t, the outputs and the component's parameters (by variable name), e.g. "q >= 0".
    '''
    failures = []
    env = {'np': np, 't': t}
    env.update(INVARIANT_HELPERS)
    env.update(params or {})
    env.update(outputs)
    for inv in spec.get('invariants') or []:
        # applies: all (default) | run -- 'run' invariants hold at the nominal parameters only
        if isinstance(inv, dict) and inv.get('applies', 'all') != 'all' and where not in inv['applies'].split(','):
            continue
        expr = inv['expr'] if isinstance(inv, dict) else inv
        try:
            ok = np.all(eval(expr, {'__builtins__': SAFE_BUILTINS}, env))  # noqa: S307 - expressions come from the repo's own yaml
        except Exception as e:
            failures.append(f'{expr}: could not evaluate ({e})')
            continue
        if not ok:
            failures.append(f'{expr}: violated')
    return failures


def _non_finite(outputs):
    return [name for name, y in outputs.items() if not np.all(np.isfinite(y))]


def _guard(component, test, fn):
    '''Runs a check; any unexpected exception is a failed Result, not a crash.'''
    plots.TIME_LABEL = component.spec.get('time_label', 'time [s]')
    try:
        return save(component, fn())
    except harness.MissingParameters as e:
        return save(component, Result(test, SKIPPED, str(e)))
    except Exception as e:
        return save(component, Result(test, FAILED, f'{type(e).__name__}: {e}',
                                      details=[traceback.format_exc()[-3000:]]))


# ----------------------------------------------------------------------------------------------
# run_test
# ----------------------------------------------------------------------------------------------

def run_test(cm):
    def check():
        t, outputs = cm.run()
        fig = plots.plot_outputs(plot_path(cm.component, 'run'), t, outputs, cm.units(),
                                 f'{cm.component.label}: nominal parameters')
        bad = _non_finite(outputs)
        inv = _check_invariants(cm.spec, t, outputs, cm.param_env(), where='run')
        metrics = {'n_time_points': len(t), 'sim_time': float(cm.spec['sim_time']),
                   'final_values': {k: float(v[-1]) for k, v in outputs.items()}}
        if bad or inv:
            msg = '; '.join(([f'non-finite output: {", ".join(bad)}'] if bad else []) + inv)
            return Result('run_test', FAILED, msg, metrics, [fig], inv)
        return Result('run_test', PASSED, f'generated and simulated {cm.spec["sim_time"]} s; '
                      f'all {len(outputs)} outputs finite', metrics, [fig])
    return _guard(cm.component, 'run_test', check)


# ----------------------------------------------------------------------------------------------
# verification_test_BC
# ----------------------------------------------------------------------------------------------

def sweep_values(param_name, nominal, sweep_spec):
    """Values to sweep: bc_sweep.ranges[param] as [min, max] or {min, max, points, scale: log},
    else bc_sweep.factors times the nominal value."""
    ranges = sweep_spec.get('ranges') or {}
    points = int(sweep_spec.get('points', 5))
    if param_name in ranges:
        r = ranges[param_name]
        if isinstance(r, dict):
            lo, hi = float(r['min']), float(r['max'])
            n = int(r.get('points', points))
            if r.get('scale') == 'log':
                return list(np.geomspace(lo, hi, n)), None
            return list(np.linspace(lo, hi, n)), None
        lo, hi = r
        return list(np.linspace(float(lo), float(hi), points)), None
    if nominal == 0:
        return [], 'nominal value is 0; set bc_sweep.ranges to sweep it'
    return [nominal * float(f) for f in sweep_spec.get('factors', [0.5, 1.0, 2.0])], None


def verification_test_BC(cm):
    def check():
        component, spec = cm.component, cm.spec
        sweep_spec = spec.get('bc_sweep') or {}
        # 'all' (default): boundary conditions and every constant, so each parameter gets a
        # verified range; 'bcs': boundary conditions only.
        mode = sweep_spec.get('sweep', 'all')
        bcs = [p.variable_name for p in component.boundary_conditions()]
        if mode == 'bcs':
            names, swept_kind = bcs, 'boundary conditions'
        else:
            names = bcs + [p.variable_name for p in component.parameters()
                           if p.kind != 'boundary_condition' and p.variable_name not in bcs]
            swept_kind = 'boundary conditions and constants' if bcs else 'constants (no boundary conditions)'
        names = [n for n in names if n not in (sweep_spec.get('exclude') or [])]
        names += [n for n in sweep_spec.get('extra_parameters', []) if n not in names]
        if not names:
            return Result('verification_test_BC', SKIPPED, 'component has no boundary conditions or constants to sweep')
        by_var = {p.variable_name: p for p in component.parameters()}
        plot_output = sweep_spec.get('plot_output') or cm.outputs()[0]

        sweeps, failures, notes, n_runs, verified = {}, [], [], 0, {}
        for var in names:
            pname = harness.parameter_name(by_var[var])
            values, why_not = sweep_values(var, cm.nominal[pname], sweep_spec)
            if why_not:
                notes.append(f'{var}: {why_not}')
                continue
            runs, passed = [], []
            for value in values:
                n_runs += 1
                try:
                    t, outputs = cm.run({pname: value})
                except Exception as e:
                    failures.append(f'{var}={value:.4g}: simulation failed ({type(e).__name__}: {e})')
                    runs.append((value, None, None))
                    passed.append((value, False))
                    continue
                bad = _non_finite(outputs)
                inv = _check_invariants(spec, t, outputs, cm.param_env({pname: value}))
                if bad:
                    failures.append(f'{var}={value:.4g}: non-finite {", ".join(bad)}')
                failures += [f'{var}={value:.4g}: {m}' for m in inv]
                runs.append((value, t, outputs[plot_output] if not bad else None))
                passed.append((value, not bad and not inv))
            sweeps[var] = runs
            verified[var] = ranges.verified_ranges_from_sweep(passed, cm.nominal[pname])

        # The ranges recorded in <name>_parameters.csv must be backed by this run.
        unsupported = []
        for p in component.parameters():
            rec = ranges.recorded_range(p, 'verified')
            if rec is None:
                continue
            got = verified.get(p.variable_name)
            tol = 1e-9 * max(abs(rec[0]), abs(rec[1]), 1.0)
            if got is None or rec[0] < got[0] - tol or rec[1] > got[1] + tol:
                unsupported.append(f'{p.variable_name}: recorded verified range [{rec[0]:.4g}, {rec[1]:.4g}] '
                                   f'not supported by this sweep ({got})')

        figs = []
        if sweeps:
            figs.append(plots.plot_sweep(plot_path(component, 'bc_sweep'), sweeps, plot_output, cm.units(),
                                         f'{component.label}: boundary-condition sweep'))
        metrics = {'swept_kind': swept_kind, 'swept': list(sweeps), 'n_runs': n_runs, 'n_failures': len(failures),
                   'verified_ranges': verified, 'unsupported_recorded_ranges': unsupported,
                   'sweep_values': {k: [r[0] for r in v] for k, v in sweeps.items()}}
        if unsupported:
            return Result('verification_test_BC', FAILED,
                          f'{len(unsupported)} recorded verified range(s) not supported by the latest sweep '
                          f'(re-run `make ranges` after checking why)', metrics, figs, unsupported + failures + notes)
        if failures:
            # Failures outside the verified span are information, not a defect: the verified
            # range records where the module works. Only a failure at the nominal value fails.
            nominal_failures = [v for v, r in verified.items() if r is None]
            if nominal_failures:
                return Result('verification_test_BC', FAILED,
                              f'fails at the nominal value of: {", ".join(nominal_failures)}',
                              metrics, figs, failures + notes)
            return Result('verification_test_BC', PASSED,
                          f'{n_runs} runs over {len(sweeps)} {swept_kind}; {len(failures)} run(s) outside the '
                          f'verified ranges failed (see ranges)', metrics, figs, failures + notes)
        if not sweeps:
            return Result('verification_test_BC', SKIPPED, 'no boundary condition could be swept', metrics,
                          figs, notes)
        return Result('verification_test_BC', PASSED,
                      f'{n_runs} runs over {len(sweeps)} {swept_kind}: all finite, invariants hold',
                      metrics, figs, notes)
    return _guard(cm.component, 'verification_test_BC', check)


# ----------------------------------------------------------------------------------------------
# verification_test_timestep
# ----------------------------------------------------------------------------------------------

def verification_test_timestep(cm):
    def check():
        component = cm.component
        ts = cm.spec.get('timestep') or {}
        scheme = ts.get('scheme', 'rk4')
        order = fixed_step.SCHEMES[scheme][0]
        dts = [float(d) for d in ts.get('dts', [4e-3, 2e-3, 1e-3, 5e-4])]
        t_end = float(ts.get('t_end', cm.spec['sim_time']))
        min_order = float(ts.get('min_order', 0.8))
        tol = float(ts.get('tol', 1e-4))
        cvode_tol = float(ts.get('cvode_tol', 1e-3))
        roundoff = float(ts.get('roundoff', 1e-11))

        rhs, y0, state_qnames, observe = fixed_step.compile_rhs(cm.helper.model)
        if not state_qnames:
            return Result('verification_test_timestep', SKIPPED, 'component has no state variables (algebraic only)')

        outputs = cm.outputs()
        qnames = [cm.helper._resolve_name(harness.output_name(o))[1] for o in outputs]
        stride_base = dts[0]
        runs = []
        for dt in dts:
            t, Y = fixed_step.integrate(rhs, y0, t_end, dt, scheme)
            stride = int(round(stride_base / dt))
            idx = np.arange(0, len(t), stride)
            vals = np.array([observe(t[i], Y[i], qnames) for i in idx], dtype=float)
            runs.append((dt, t[idx], vals))

        finest = runs[-1][2]
        scale = np.maximum(np.max(np.abs(finest), axis=0), np.ptp(finest, axis=0))
        scale[scale == 0] = 1.0
        diffs = [float(np.max(np.abs(runs[k][2] - runs[k + 1][2]) / scale)) for k in range(len(runs) - 1)]
        observed = []
        for k in range(len(diffs) - 1):
            if diffs[k + 1] > roundoff and diffs[k] > roundoff:
                observed.append(math.log(diffs[k] / diffs[k + 1]) / math.log(dts[k] / dts[k + 1]))
            else:
                observed.append(None)

        # cross-check against libcuflynx's CVODE run with tight tolerances
        with contextlib.redirect_stdout(io.StringIO()):
            tight = harness.simulation_helper(cm.model_path, dict(cm.spec, dt=stride_base, sim_time=t_end, pre_time=0.0),
                                              solver_info={'rtol': 1e-10, 'atol': 1e-12})
        t_cv, cv = harness.run(tight, outputs)
        n = min(len(t_cv), finest.shape[0])
        cvode_diff = float(np.max(np.abs(np.column_stack([cv[o] for o in outputs])[:n] - finest[:n]) / scale))

        figs = [plots.plot_convergence(plot_path(component, 'timestep_convergence'), dts[:-1], diffs, order,
                                       next((o for o in reversed(observed) if o is not None), None),
                                       f'{component.label}: {scheme}'),
                plots.plot_timestep_traces(plot_path(component, 'timestep_traces'),
                                           [(dt, t, v[:, 0]) for dt, t, v in runs], outputs[0], cm.units(),
                                           f'{component.label}: {outputs[0]} vs step size')]
        metrics = {'scheme': scheme, 'dts': dts, 'successive_differences': diffs, 'observed_orders': observed,
                   'expected_order': order, 'cvode_tight_vs_finest': cvode_diff, 't_end': t_end}

        problems = []
        if not all(np.isfinite(diffs)):
            problems.append('non-finite solution at some step size (explicit scheme unstable? use smaller dts)')
        else:
            if any(diffs[k + 1] > diffs[k] and diffs[k + 1] > roundoff for k in range(len(diffs) - 1)):
                problems.append(f'differences do not decrease with dt: {["%.2e" % d for d in diffs]}')
            last = next((o for o in reversed(observed) if o is not None), None)
            if last is not None and last < min_order:
                problems.append(f'observed order {last:.2f} < {min_order}')
            if diffs[-1] > tol:
                problems.append(f'finest-step difference {diffs[-1]:.2e} > tol {tol:.1e}')
            if cvode_diff > cvode_tol:
                problems.append(f'fixed-step and CVODE solutions differ by {cvode_diff:.2e} > {cvode_tol:.1e}')
        last = next((o for o in reversed(observed) if o is not None), None)
        summary = (f'{scheme}: differences {diffs[0]:.1e} -> {diffs[-1]:.1e}'
                   + (f', observed order {last:.2f} (expected {order})' if last is not None else ', at round-off')
                   + f'; agrees with CVODE to {cvode_diff:.1e}')
        if problems:
            return Result('verification_test_timestep', FAILED, '; '.join(problems), metrics, figs, problems)
        return Result('verification_test_timestep', PASSED, summary, metrics, figs)
    return _guard(cm.component, 'verification_test_timestep', check)


# ----------------------------------------------------------------------------------------------
# stability_test
# ----------------------------------------------------------------------------------------------

DEFAULT_STABILITY = {
    'cvode': [
        {'rtol': 1e-4, 'atol': 1e-6},
        {'rtol': 1e-6, 'atol': 1e-8},
        {'rtol': 1e-8, 'atol': 1e-10},
    ],
    'solve_ivp': ['RK45', 'BDF', 'LSODA'],
    'fixed_step': {'schemes': ['euler', 'heun', 'rk4'], 'dt_start': 1e-2},   # halved until it works
    'max_step_start': 1e-2,             # CVODE / solve_ivp: MaximumStep tried after the unconstrained run, /10 each
    'time_budget': None,                # seconds per run before refinement stops (default 20 s per state)
    'tol': 1e-2,                        # max normalised difference from the reference to count as working
    # configurations the module is declared to work with; each must pass
    'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}],
}


def config_label(cfg):
    extra = ', '.join(f'{k}={v:g}' if isinstance(v, float) else f'{k}={v}'
                      for k, v in cfg.items() if k not in ('solver',))
    return f"{cfg['solver']}" + (f' ({extra})' if extra else '')


def _matches(declared, cfg):
    return all(cfg.get(k) == v or (isinstance(v, float) and isinstance(cfg.get(k), float)
                                   and math.isclose(cfg[k], v, rel_tol=1e-9))
               for k, v in declared.items())


def stability_test(cm):
    """
    Runs the component with a matrix of solvers and settings and records which work:
    a run "works" when it finishes, every output is finite, and it stays within ``tol``
    (normalised max difference) of a tight-tolerance CVODE reference. The configurations
    listed under ``stability.supported`` in the spec must work; the full matrix is reported
    so the HTML shows which timesteps, tolerances and solvers the module can be run with.
    """
    def check():
        from cam_testing.library import deep_merge
        component = cm.component
        st = deep_merge(DEFAULT_STABILITY, cm.spec.get('stability') or {})
        for key in ('cvode', 'solve_ivp', 'supported'):   # lists replace, not merge
            if key in (cm.spec.get('stability') or {}):
                st[key] = cm.spec['stability'][key]
        tol = float(st['tol'])
        outputs = cm.outputs()
        t_end = float(st.get('t_end', cm.spec['sim_time']))
        base = dict(cm.spec, pre_time=0.0, sim_time=t_end)

        def compare(t, out, t_ref, ref, scale):
            y = np.column_stack([np.interp(t_ref, t, out[o]) for o in outputs])
            if not np.all(np.isfinite(y)):
                return math.inf
            return float(np.max(np.abs(y - ref) / scale))

        ref_helper = harness.simulation_helper(cm.model_path, base, solver_info={'rtol': 1e-10, 'atol': 1e-12})
        t_ref, ref_out = harness.run(ref_helper, outputs)
        ref = np.column_stack([ref_out[o] for o in outputs])
        scale = np.maximum(np.max(np.abs(ref), axis=0), np.ptp(ref, axis=0))
        scale[scale == 0] = 1.0

        rows = []

        def record(cfg, fn):
            import time
            start = time.perf_counter()
            try:
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    t, out = fn()
                err = compare(np.asarray(t, float), out, t_ref, ref, scale)
                works = bool(np.isfinite(err) and err <= tol)
                note = '' if works else ('non-finite output' if not np.isfinite(err) else f'difference {err:.1e} > tol')
            except Exception as e:
                err, works, note = None, False, f'{type(e).__name__}: {str(e)[:200]}'
            row = {'config': cfg, 'label': config_label(cfg), 'works': works, 'difference': err,
                   'seconds': round(time.perf_counter() - start, 3), 'note': note}
            rows.append(row)
            return row

        # Each solver family is refined (smaller step) until it works, or until the next run
        # is predicted to take longer than the budget: 20 s per state variable by default.
        rhs, y0, state_qnames, observe = fixed_step.compile_rhs(cm.helper.model)
        n_states = max(len(state_qnames), 1)
        budget = float(st.get('time_budget') or 20.0 * n_states)
        min_step = float(st.get('min_step', 1e-7))
        families = []

        def refine(family, make_cfg, run_fn, first_step, factor, try_unconstrained):
            """Runs make_cfg(step) for step = first_step, first_step/factor, ... (after an
            unconstrained try when try_unconstrained) until one works or the budget stops it."""
            steps = ([None] if try_unconstrained else []) + [first_step]
            last = None
            stop = ''
            while True:
                step = steps.pop(0) if steps else last['step'] / factor
                if step is not None and step < min_step:
                    stop = f'reached min_step {min_step:g}'
                    break
                if last is not None and last['row']['seconds'] * (factor if last['step'] else 2.0) > budget:
                    stop = (f'stopped: next run (step {step:g}) predicted over the '
                            f'{budget:.0f} s budget ({n_states} states x 20 s)')
                    break
                row = record(make_cfg(step), lambda step=step: run_fn(step))
                last = {'step': step, 'row': row}
                if row['works']:
                    break
                if row['seconds'] > budget:
                    stop = f'stopped: run took {row["seconds"]:.0f} s > {budget:.0f} s budget'
                    break
            works_at = last['row']['label'] if last and last['row']['works'] else None
            families.append({'family': family, 'works_at': works_at, 'stop': '' if works_at else stop})

        cv_start = float(st.get('max_step_start', 1e-2))
        for info in st['cvode']:
            info = dict(info)
            given_step = info.pop('MaximumStep', None)

            def cv_cfg(step, info=info):
                cfg = dict({'solver': 'CVODE_myokit'}, **info)
                if step is not None:
                    cfg['MaximumStep'] = step
                return cfg

            def cv_run(step, info=info):
                si = dict(info, **({'MaximumStep': step} if step is not None else {}))
                return harness.run(harness.simulation_helper(cm.model_path, base, solver_info=si), outputs)
            family = config_label(cv_cfg(None))
            refine(family, cv_cfg, cv_run, given_step or cv_start, 10.0, given_step is None)

        if st['solve_ivp']:
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    py_path = harness.generate(component, cm.work_dir, model_type='python')
                py_error = None
            except Exception as e:
                py_path, py_error = None, f'python generation failed: {str(e).splitlines()[0][:200]}'
            for method in st['solve_ivp']:
                def py_cfg(step, method=method):
                    return dict({'solver': 'solve_ivp', 'method': method}, **({'max_step': step} if step else {}))
                if py_path is None:
                    rows.append({'config': py_cfg(None), 'label': config_label(py_cfg(None)), 'works': False,
                                 'difference': None, 'seconds': 0.0, 'note': py_error})
                    families.append({'family': config_label(py_cfg(None)), 'works_at': None, 'stop': py_error})
                    continue

                def py_run(step, method=method):
                    si = dict({'method': method}, **({'max_step': step} if step else {}))
                    return harness.run(harness.simulation_helper(py_path, base, solver='solve_ivp', model_type='python',
                                                                 solver_info=si), outputs)
                refine(config_label(py_cfg(None)), py_cfg, py_run, cv_start, 10.0, True)

        fs = st.get('fixed_step') or {}
        if fs.get('schemes') and state_qnames:
            qnames = [cm.helper._resolve_name(harness.output_name(o))[1] for o in outputs]
            out_dt = float(cm.spec['dt'])
            for scheme in fs['schemes']:
                def fx_run(dt, scheme=scheme):
                    with np.errstate(all='ignore'):
                        t, Y = fixed_step.integrate(rhs, y0, t_end, dt, scheme)
                        stride = max(1, int(round(out_dt / dt)))
                        idx = np.arange(0, len(t), stride)
                        vals = np.array([observe(t[i], Y[i], qnames) for i in idx], dtype=float)
                    return t[idx], {o: vals[:, k] for k, o in enumerate(outputs)}
                refine(f'fixed_step_{scheme}', lambda dt, scheme=scheme: {'solver': f'fixed_step_{scheme}', 'dt': dt},
                       fx_run, float(fs.get('dt_start', (fs.get('dts') or [1e-2])[0])), 2.0, False)

        # declared-supported configurations are always run, even if refinement stopped earlier
        for declared in st.get('supported') or []:
            if any(_matches(declared, r['config']) for r in rows):
                continue
            d = dict(declared)
            solver = d.pop('solver')
            if solver == 'CVODE_myokit':
                record(dict(declared), lambda d=d: harness.run(harness.simulation_helper(cm.model_path, base, solver_info=d), outputs))

        supported = st.get('supported') or []
        broken = []
        for declared in supported:
            matching = [r for r in rows if _matches(declared, r['config'])]
            if not matching:
                broken.append(f'{config_label(declared)}: declared supported but not in the tested matrix')
            broken += [f"{r['label']}: declared supported but {r['note']}" for r in matching if not r['works']]

        working = [r['label'] for r in rows if r['works']]
        metrics = {'matrix': rows, 'families': families, 'tol': tol, 'budget_seconds': budget, 't_end': t_end,
                   'supported': [config_label(d) for d in supported],
                   'n_working': len(working), 'n_configs': len(rows)}
        if broken:
            return Result('stability_test', FAILED, '; '.join(broken), metrics, [], broken)
        n_fam_ok = sum(1 for f in families if f['works_at'])
        return Result('stability_test', PASSED,
                      f'{n_fam_ok} of {len(families)} solver families work at some step '
                      f'({len(working)} of {len(rows)} configurations tried); all declared-supported ones do',
                      metrics)
    return _guard(cm.component, 'stability_test', check)


# ----------------------------------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------------------------------

def _validation_status(component, kind, test):
    v = (component.spec.get('validation') or {}).get(kind) or {}
    status = v.get('status', PENDING)
    if status in (PENDING, SKIPPED):
        return Result(test, status, v.get('reason', 'no validation data chosen yet')), v
    return None, v


def validation_test_baseline(cm):
    '''
    Compare the model with baseline data. Spec (validation.baseline):
        status: active
        source: citation / URL
        data: validation/<file>.csv        time column + one column per compared output
        time_column: t        time_offset: 0
        variables: {model_output: csv_column}
        parameters: {variable_name: value}  the parameter set the data are compared at
        sim_time / dt: optional (defaults: data span, spec dt)
        metric: nrmse | log_rmse        threshold: 0.1
    '''
    def check():
        from cam_testing.calibrate import load_data, score
        component = cm.component
        early, v = _validation_status(component, 'baseline', 'validation_test_baseline')
        if early:
            return early
        t_data, data = load_data(component.module.dir, v)
        by_var = {p.variable_name: p for p in component.parameters()}
        params = {k: float(val) for k, val in (v.get('parameters') or {}).items()}
        overrides = {'parameters/' + harness.parameter_name(by_var[k]): val for k, val in params.items()}
        run_spec = dict(cm.spec, pre_time=0.0, sim_time=float(v.get('sim_time', np.max(t_data))),
                        dt=float(v.get('dt', cm.spec['dt'])))
        helper = harness.simulation_helper(cm.model_path, run_spec, solver_info={'rtol': 1e-8, 'atol': 1e-10})
        t, out = harness.run(helper, list(v['variables']), params=overrides)
        metric = v.get('metric', 'nrmse')
        threshold = float(v.get('threshold', 0.1))
        errors = {var: score(metric, np.interp(t_data, t, out[var]), data[var]) for var in v['variables']}
        fig = plots.plot_model_vs_data(plot_path(component, 'validation_baseline'), t, out,
                                       {k: t_data for k in data}, data, cm.units(),
                                       f'{component.label}: model vs {v.get("source_short", "data")}')
        metrics = {metric: errors, 'threshold': threshold, 'source': v.get('source', ''),
                   'validated_values': params}
        worst = max(errors.values())
        if worst > threshold:
            return Result('validation_test_baseline', FAILED, f'worst {metric} {worst:.3f} > {threshold}', metrics, [fig])
        return Result('validation_test_baseline', PASSED, f'worst {metric} {worst:.3f} <= {threshold}', metrics, [fig])
    return _guard(cm.component, 'validation_test_baseline', check)


def validation_test_calibrate(cm):
    '''
    Calibrate to one data set with libcuflynx parameter identification, then validate on
    held-out data. Implemented per module when its validation data is chosen; until then the
    spec's status (pending / skipped) and reason are reported.
    '''
    def check():
        early, v = _validation_status(cm.component, 'calibrate', 'validation_test_calibrate')
        if early:
            return early
        from cam_testing import calibrate
        return calibrate.run(cm, v, plot_path)
    return _guard(cm.component, 'validation_test_calibrate', check)


CHECKS = {
    'run_test': run_test,
    'verification_test_BC': verification_test_BC,
    'verification_test_timestep': verification_test_timestep,
    'stability_test': stability_test,
    'validation_test_baseline': validation_test_baseline,
    'validation_test_calibrate': validation_test_calibrate,
}
