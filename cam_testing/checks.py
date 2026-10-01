"""
The standard verification & validation tests, run on one module version at a time (the
validation tests on one of its instances).

Each check returns a Result and writes it to <version dir>/results/<test>.json (instance tests:
results/instances/<instance>/<test>.json), with any figures in <version dir>/plots/. The pytest
layer (tests/test_modules.py) turns a Result into pass / fail / skip; the report generator reads
the JSON files, so the HTML shows exactly what the last test run found.
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

from cam_testing import fixed_step, harness, plots
from cam_testing.library import REPO_ROOT

PASSED, FAILED, SKIPPED, PENDING, NOT_APPLICABLE = 'passed', 'failed', 'skipped', 'pending', 'not_applicable'
# validation.<kind>.status: 'active' (confirmed data) or 'proposed' (run, but the data awaits confirmation)
PROPOSED = 'proposed'

# per version
VERSION_TESTS = ['run_test', 'verification_test_invariants', 'verification_test_BC', 'verification_test_timestep',
                 'stability_test']
# per instance
INSTANCE_TESTS = ['validation_test_baseline', 'validation_test_calibrate']
TESTS = VERSION_TESTS + INSTANCE_TESTS
NO_CALIBRATION_DATA = 'no calibration data in this instance'
NO_BASELINE_DATA = 'no baseline data in this instance'


@dataclass
class Result:
    test: str
    status: str
    message: str = ''
    metrics: dict = field(default_factory=dict)
    plots: list = field(default_factory=list)
    details: list = field(default_factory=list)
    timestamp: str = ''


def result_path(component, test, instance=None):
    if instance is not None:
        return os.path.join(instance.results_dir, f'{test}.json')
    return os.path.join(component.results_dir, f'{test}.json')


def plot_path(component, name, instance=None):
    prefix = f'{instance.name}__' if instance is not None else ''
    return os.path.join(component.plots_dir, f'{prefix}{name}.png')


def save(component, result, instance=None):
    result.timestamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    # plots are stored relative to the version dir, which is where the HTML lives
    result.plots = [os.path.relpath(p, component.dir) for p in result.plots]
    path = result_path(component, result.test, instance)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(asdict(result), f, indent=2, default=_json_default)
    return result


def load(component, test, instance=None):
    path = result_path(component, test, instance)
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
    '''A version's generated model and simulation helper, built once and reused. The parameter
    values are an instance's (default: the version's default instance).'''

    def __init__(self, component, work_dir=None, instance=None):
        self.component = component
        self.instance = instance or component.default_instance
        self.spec = component.spec
        self.work_dir = work_dir or tempfile.mkdtemp(prefix=f'cam_{component.id}_')
        self._model_path = None
        self._helper = None
        self._params = component.parameters() if self.instance.is_default else self.instance.parameters()
        self.nominal = {harness.parameter_name(p): p.float_value
                        for p in self._params if not p.is_todo}
        self._var_of = {harness.parameter_name(p): p.variable_name for p in self._params}
        # the test network's own parameters (e.g. an outlet's flow v_vout), by their full names, so a
        # sweep can vary them (bc_sweep.extra_parameters) and invariants can use them
        for name, _units, value, *_ref in (self.spec.get('harness') or {}).get('parameters') or []:
            if name not in self.nominal:
                self.nominal[name] = float(value)
                self._var_of[name] = name

    def param_env(self, overrides=None):
        '''Parameter values by variable name (e.g. alpha), for invariant expressions.'''
        values = dict(self.nominal)
        values.update(overrides or {})
        return {self._var_of[k]: v for k, v in values.items()}

    @property
    def model_path(self):
        if self._model_path is None:
            if getattr(self, '_generation_error', None) is not None:
                raise self._generation_error          # don't regenerate a model that failed to generate
            try:
                self._model_path = harness.generate(self.component, self.work_dir, parameters=self._instance_params())
            except harness.GenerationFailed as e:
                self._generation_error = e
                raise
        return self._model_path

    @property
    def helper(self):
        if self._helper is None:
            self._helper = harness.simulation_helper(self.model_path, self.spec)
        return self._helper

    def parameters(self):
        '''The parameters the model was generated with (the instance's).'''
        return self._params

    def _instance_params(self):
        return None if self.instance.is_default else self._params

    def generate(self, work_dir, **kw):
        '''Another model of this version at this instance's parameters (e.g. a python model).'''
        return harness.generate(self.component, work_dir, parameters=self._instance_params(), **kw)

    def run_point_params(self):
        '''The spec's run_parameters as simulation-helper overrides ({'parameters/<name>': value}).'''
        return {'parameters/' + k: v for k, v in _run_point(self).items()}

    def point_model(self):
        '''A copy of the Myokit model with the run_parameters applied (for fixed-step integration).'''
        model = self.helper.model.clone()
        for name, value in _run_point(self).items():
            for comp in ('parameters', 'parameters_global'):
                if model.has_variable(f'{comp}.{name}'):
                    model.get(f'{comp}.{name}').set_rhs(value)
                    break
            else:
                raise KeyError(f'run parameter {name} not found in the model')
        return model

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


def cumint(t, y):
    '''Cumulative trapezoidal integral of y over t, starting at 0 (for balance invariants).'''
    t, y = np.asarray(t, float), np.asarray(y, float)
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])


INVARIANT_HELPERS = {'period': period, 'cumint': cumint}


def _check_invariants(spec, t, outputs, params=None, where='sweep'):
    '''
    Returns a list of failed invariant descriptions. An invariant is a numpy expression of
    t, the outputs and the component's parameters (by variable name), e.g. "q >= 0".
    '''
    failures = []
    env = {'np': np, 't': t, '__builtins__': SAFE_BUILTINS}
    env.update(INVARIANT_HELPERS)
    env.update(params or {})
    env.update(outputs)
    env['param'] = dict(params or {})     # param['x']: the parameter even when an output is also named x
    for inv in spec.get('invariants') or []:
        # applies: all (default) | run -- 'run' invariants hold at the nominal parameters only
        if isinstance(inv, dict) and inv.get('applies', 'all') != 'all' and where not in inv['applies'].split(','):
            continue
        expr = inv['expr'] if isinstance(inv, dict) else inv
        try:
            # one namespace, so comprehensions/lambdas in an expression can see the parameters and np
            ok = np.all(eval(expr, env))  # noqa: S307 - expressions come from the repo's own yaml
        except Exception as e:
            failures.append(f'{expr}: could not evaluate ({e})')
            continue
        if not ok:
            failures.append(f'{expr}: violated')
    return failures


def _non_finite(outputs):
    return [name for name, y in outputs.items() if not np.all(np.isfinite(y))]


def _guard(component, test, fn, instance=None):
    '''Runs a check; any unexpected exception is a failed Result, not a crash.'''
    plots.TIME_LABEL = component.spec.get('time_label', 'time [s]')
    try:
        return save(component, fn(), instance)
    except harness.MissingParameters as e:
        return save(component, Result(test, SKIPPED, str(e)), instance)
    except Exception as e:
        return save(component, Result(test, FAILED, f'{type(e).__name__}: {e}',
                                      details=[traceback.format_exc()[-3000:]]), instance)


# ----------------------------------------------------------------------------------------------
# run_test
# ----------------------------------------------------------------------------------------------

def reference_tolerances(cm):
    '''Tolerances of the tight CVODE reference: rtol 1e-10 and atol min(1e-12, the component's
    atol), so tiny state scales (e.g. microvessel volumes) keep error control without asking
    CVODE for an absolute accuracy beyond the component's own.'''
    atol = float((cm.spec.get('solver_info') or {}).get('atol', 1e-8))
    ref = cm.spec.get('reference_solver_info') or {}
    return {'rtol': float(ref.get('rtol', 1e-10)), 'atol': float(ref.get('atol', min(1e-12, atol)))}


def _run_point(cm):
    """run_parameters from the spec: the operating point for the run and invariant tests
    (e.g. off-equilibrium boundary conditions so the dynamics are exercised); the library
    defaults stay as they are."""
    by_var = {p.variable_name: p for p in cm.parameters()}
    return {(harness.parameter_name(by_var[k]) if k in by_var else k): float(v) for k, v in (cm.spec.get('run_parameters') or {}).items()}


CPP_NOT_APPLICABLE = ('C++ 1D-solver component (module_format cpp): libcuflynx generates the 0D model coupled to it '
                      '(checked by run_test), but the 1D solver it couples to is not runnable here, so there is no '
                      'simulation to check')


def is_cpp(component):
    return component.config.get('module_format', 'cellml') == 'cpp'


def _cpp_not_applicable(cm, test):
    return save(cm.component, Result(test, NOT_APPLICABLE, CPP_NOT_APPLICABLE))


def cpp_generation_check(cm):
    '''
    run_test for a C++ 1D-solver component: its test network (harness) is generated through
    libcuflynx's C++ 0D-1D path (model_type cpp, couple_to_1d) from this library, and the generated
    model0d.cc must compile (g++ -fsyntax-only against the SUNDIALS headers).
    '''
    import shutil
    import subprocess
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture

    def check():
        component = cm.component
        prefix = component.id
        work = os.path.join(cm.work_dir, 'cpp')
        res = os.path.join(work, 'resources')
        harness._write_resources(component, res, prefix, {})
        cfg = {'file_prefix': prefix, 'input_param_file': f'{prefix}_parameters.csv', 'model_type': 'cpp',
               'solver': 'RK4', 'couple_to_1d': True, 'resources_dir': res,
               'generated_models_dir': os.path.join(work, 'generated_models'),
               'cpp_generated_models_dir': os.path.join(work, 'cpp_out'), 'cpp_1d_model_config_path': None,
               'module_library_dirs': [harness.MODULES_DIR], 'use_builtin_modules': False, 'DEBUG': False,
               'dt': float(cm.spec.get('dt', 1e-3)),
               'solver_info': {'dt_solver': 1e-4, 'MaximumNumberOfSteps': 5000, 'solver': 'RK4'}}
        log = io.StringIO()
        try:
            with contextlib.redirect_stdout(log):
                ok = generate_with_new_architecture(False, cfg)
        except SystemExit as e:
            ok = False
            log.write(f'\nlibcuflynx exited ({e.code})')
        tail = [l for l in log.getvalue().splitlines() if l.strip()][-6:]
        if not ok:
            return Result('run_test', FAILED, 'C++ 0D-1D generation failed: ' + (tail[-1] if tail else ''), details=tail)
        src = os.path.join(cfg['cpp_generated_models_dir'], 'model0d.cc')
        coupler = [f for f in os.listdir(cfg['cpp_generated_models_dir']) if f.endswith('_coupler1d0d.json')]
        metrics = {'generated': sorted(os.listdir(cfg['cpp_generated_models_dir'])), 'coupler': coupler}
        gxx = shutil.which('g++')
        if not gxx:
            return Result('run_test', PASSED, 'generated the C++ 0D model and 0D-1D coupler (g++ not available, '
                          'compilation not checked)', metrics)
        p = subprocess.run([gxx, '-std=c++17', '-fsyntax-only', src], capture_output=True, text=True)
        if p.returncode != 0:
            return Result('run_test', FAILED, 'generated model0d.cc does not compile: '
                          + (p.stderr.strip().splitlines() or [''])[0][:300], metrics, details=p.stderr.splitlines()[:20])
        return Result('run_test', PASSED, 'generated the C++ 0D model coupled to the 1D vessel (model0d.cc, '
                      f'{", ".join(coupler)}) and it compiles; the 1D solver itself is not runnable here', metrics)
    return _guard(cm.component, 'run_test', check)


def run_test(cm):
    """Generates the component with libcuflynx, simulates it, and checks every output is finite."""
    if is_cpp(cm.component):
        return cpp_generation_check(cm)

    def check():
        point = _run_point(cm)
        t, outputs = cm.run(point)
        fig = plots.plot_outputs(plot_path(cm.component, 'run'), t, outputs, cm.units(),
                                 f'{cm.component.label}: ' + ('run parameters' if point else 'nominal parameters'))
        bad = _non_finite(outputs)
        metrics = {'n_time_points': len(t), 'sim_time': float(cm.spec['sim_time']),
                   'final_values': {k: float(v[-1]) for k, v in outputs.items()}}
        if bad:
            return Result('run_test', FAILED, f'non-finite output: {", ".join(bad)}', metrics, [fig])
        # A component whose outputs never change may be a silent no-op (e.g. parameters that don't
        # reach the equations); flag it so the invariants/review look at it. Not a failure: an
        # equilibrium run point is legitimate.
        flat = [k for k, v in outputs.items() if np.ptp(v) == 0]
        zero = [k for k in flat if np.all(outputs[k] == 0)]
        metrics['constant_outputs'] = flat
        warning = ''
        if flat and len(flat) == len(outputs):
            warning = f'; WARNING: every output is constant{" (all zero)" if len(zero) == len(outputs) else ""}'
        elif zero:
            warning = f'; note: identically zero: {", ".join(zero)}'
        return Result('run_test', PASSED, f'generated and simulated {cm.spec["sim_time"]} s; '
                      f'all {len(outputs)} outputs finite' + warning, metrics, [fig],
                      [f'constant outputs: {", ".join(flat)}'] if flat else [])
    return _guard(cm.component, 'run_test', check)


# ----------------------------------------------------------------------------------------------
# verification_test_invariants
# ----------------------------------------------------------------------------------------------

def verification_test_invariants(cm):
    if is_cpp(cm.component):
        return _cpp_not_applicable(cm, 'verification_test_invariants')
    """
    Checks the simulation against what the component is supposed to do: the spec's
    invariants (exact solutions, conservation laws, bounds, delays, ...) evaluated on the run
    at the run parameters. The BC sweep checks the same invariants across parameter ranges.
    """
    def check():
        invariants = cm.spec.get('invariants') or []
        if not invariants:
            return Result('verification_test_invariants', SKIPPED, 'no invariants defined in the spec for this component')
        point = _run_point(cm)
        t, outputs = cm.run(point)
        with np.errstate(all='ignore'):
            failures = _check_invariants(cm.spec, t, outputs, cm.param_env(point), where='run')
        checked = [inv for inv in invariants
                   if not (isinstance(inv, dict) and inv.get('applies', 'all') != 'all' and 'run' not in inv['applies'].split(','))]
        results = []
        for inv in checked:
            expr = inv['expr'] if isinstance(inv, dict) else inv
            desc = (inv.get('description') if isinstance(inv, dict) else None) or expr
            failed = [f for f in failures if f.startswith(expr)]
            results.append({'invariant': desc, 'expr': expr, 'holds': not failed, 'note': failed[0][len(expr):].strip(': ') if failed else ''})
        metrics = {'invariants': results}
        n_fail = sum(not r['holds'] for r in results)
        if n_fail:
            names = '; '.join(r['invariant'] for r in results if not r['holds'])
            return Result('verification_test_invariants', FAILED, f'{n_fail} of {len(results)} invariants violated: {names}',
                          metrics, [], failures)
        return Result('verification_test_invariants', PASSED, f'all {len(results)} invariants hold', metrics)
    return _guard(cm.component, 'verification_test_invariants', check)


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
        if isinstance(r, dict) and 'values' in r:
            return [float(v) for v in r['values']], None      # discrete (e.g. 0/1 flags)
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
    if is_cpp(cm.component):
        return _cpp_not_applicable(cm, 'verification_test_BC')
    def check():
        component, spec = cm.component, cm.spec
        sweep_spec = spec.get('bc_sweep') or {}
        # 'all' (default): boundary conditions and every constant; 'bcs': boundary conditions only.
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

        sweeps, failures, notes, n_runs, failed_runs = {}, [], [], 0, 0
        for var in names:
            # a component parameter, or (extra_parameters) one of its test network's, by full name
            pname = harness.parameter_name(by_var[var]) if var in by_var else var
            if pname not in cm.nominal:
                notes.append(f'{var}: not a parameter of the component or its test network')
                continue
            values, why_not = sweep_values(var, cm.nominal[pname], sweep_spec)
            if why_not:
                notes.append(f'{var}: {why_not}')
                continue
            runs = []
            for value in values:
                n_runs += 1
                try:
                    t, outputs = cm.run({pname: value})
                except Exception as e:
                    failed_runs += 1
                    failures.append(f'{var}={value:.4g}: simulation failed ({type(e).__name__}: {e})')
                    runs.append((value, None, None))
                    continue
                bad = _non_finite(outputs)
                inv = _check_invariants(spec, t, outputs, cm.param_env({pname: value}))
                if bad or inv:
                    failed_runs += 1
                if bad:
                    failures.append(f'{var}={value:.4g}: non-finite {", ".join(bad)}')
                failures += [f'{var}={value:.4g}: {m}' for m in inv]
                runs.append((value, t, outputs[harness.output_key(plot_output)] if not bad else None))
            sweeps[var] = runs

        figs = []
        if sweeps:
            figs.append(plots.plot_sweep(plot_path(component, 'bc_sweep'), sweeps, plot_output, cm.units(),
                                         f'{component.label}: boundary-condition sweep'))
        metrics = {'swept_kind': swept_kind, 'swept': list(sweeps), 'n_runs': n_runs, 'n_failures': len(failures),
                   'sweep_values': {k: [r[0] for r in v] for k, v in sweeps.items()}}
        if failures:
            return Result('verification_test_BC', FAILED, f'{failed_runs} of {n_runs} runs failed',
                          metrics, figs, failures + notes)
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
    if is_cpp(cm.component):
        return _cpp_not_applicable(cm, 'verification_test_timestep')
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

        rhs, y0, state_qnames, observe = fixed_step.compile_rhs(cm.point_model())
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
        # sawtooth outputs (timestep.wrapped: {output: period}, e.g. a cardiac phase) are compared modulo their
        # period, so a wrap one step earlier at one step size isn't a difference of a whole period
        periods = np.array([float((ts.get('wrapped') or {}).get(o, 0.0)) for o in outputs])

        def _diff(a, b):
            d = np.abs(a - b)
            w = periods > 0
            if w.any():
                d[..., w] = d[..., w] % periods[w]
                d[..., w] = np.minimum(d[..., w], periods[w] - d[..., w])
            return d

        diffs = [float(np.max(_diff(runs[k][2], runs[k + 1][2]) / scale)) for k in range(len(runs) - 1)]
        observed = []
        for k in range(len(diffs) - 1):
            if diffs[k + 1] > roundoff and diffs[k] > roundoff:
                observed.append(math.log(diffs[k] / diffs[k + 1]) / math.log(dts[k] / dts[k + 1]))
            else:
                observed.append(None)

        # cross-check against libcuflynx's CVODE run with tight tolerances
        with contextlib.redirect_stdout(io.StringIO()):
            tight = harness.simulation_helper(cm.model_path, dict(cm.spec, dt=stride_base, sim_time=t_end, pre_time=0.0),
                                              solver_info=reference_tolerances(cm))
        t_cv, cv = harness.run(tight, outputs, params=cm.run_point_params())
        n = min(len(t_cv), finest.shape[0])
        cvode_diff = float(np.max(_diff(np.column_stack([cv[harness.output_key(o)] for o in outputs])[:n], finest[:n]) / scale))

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
    if is_cpp(cm.component):
        return _cpp_not_applicable(cm, 'stability_test')
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
            y = np.column_stack([np.interp(t_ref, t, out[harness.output_key(o)]) for o in outputs])
            if not np.all(np.isfinite(y)):
                return math.inf
            return float(np.max(np.abs(y - ref) / scale))

        ref_helper = harness.simulation_helper(cm.model_path, base, solver_info=reference_tolerances(cm))
        point = cm.run_point_params()         # integrate at the run point, like the run/invariant tests
        t_ref, ref_out = harness.run(ref_helper, outputs, params=point)
        ref = np.column_stack([ref_out[harness.output_key(o)] for o in outputs])
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
        rhs, y0, state_qnames, observe = fixed_step.compile_rhs(cm.point_model())
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
                return harness.run(harness.simulation_helper(cm.model_path, base, solver_info=si), outputs, params=point)
            family = config_label(cv_cfg(None))
            refine(family, cv_cfg, cv_run, given_step or cv_start, 10.0, given_step is None)

        if st['solve_ivp']:
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    py_path = cm.generate(cm.work_dir, model_type='python')
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
                                                                 solver_info=si), outputs, params=point)
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
                    return t[idx], {harness.output_key(o): vals[:, k] for k, o in enumerate(outputs)}
                refine(f'fixed_step_{scheme}', lambda dt, scheme=scheme: {'solver': f'fixed_step_{scheme}', 'dt': dt},
                       fx_run, float(fs.get('dt_start', (fs.get('dts') or [1e-2])[0])), 2.0, False)

        # declared-supported configurations are always run, even if refinement stopped earlier
        for declared in st.get('supported') or []:
            if any(_matches(declared, r['config']) for r in rows):
                continue
            d = dict(declared)
            solver = d.pop('solver')
            if solver == 'CVODE_myokit':
                record(dict(declared), lambda d=d: harness.run(harness.simulation_helper(cm.model_path, base, solver_info=d),
                                                               outputs, params=point))

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

def _validation_status(cm, kind, test):
    '''(early Result or None, spec block) for validation.<instance>.<kind> of the model's instance.'''
    inst = cm.instance
    v = inst.validation.get(kind) or {}
    if kind == 'calibrate' and not inst.has_obs_data:
        # an instance without obs_data has nothing to calibrate to: recorded as a failure (the
        # instance stays usable; pytest marks it xfail so CI isn't blocked)
        why = v.get('note') or v.get('reason')
        return Result(test, FAILED, NO_CALIBRATION_DATA, details=[why] if why else []), v
    if not v:
        return Result(test, NOT_APPLICABLE, NO_BASELINE_DATA if kind == 'baseline' else NO_CALIBRATION_DATA), v
    status = v.get('status', PENDING)
    if status in (PENDING, SKIPPED, NOT_APPLICABLE):
        return Result(test, status, v.get('reason', 'no validation data chosen yet')), v
    return None, v


def validation_test_baseline(cm):
    '''
    Compare the model, at the instance's parameters, with baseline data. Spec
    (validation.<instance>.baseline in the version's verification_config.json):
        status: active
        source: citation / URL
        data: instances/<instance>/<file>.csv   time column + one column per compared output
        time_column: t        time_offset: 0
        variables: {model_output: csv_column}
        parameters: {variable_name: value}  the parameter set the data are compared at
        sim_time / dt: optional (defaults: data span, spec dt)
        metric: nrmse | log_rmse        threshold: 0.1
    '''
    def check():
        from cam_testing.calibrate import load_data, score
        component = cm.component
        early, v = _validation_status(cm, 'baseline', 'validation_test_baseline')
        if early:
            return early
        if v.get('targets'):
            return _baseline_targets(cm, v)
        # data paths are relative to the version directory (instances/<instance>/<file>)
        t_data, data = load_data(cm.component.dir, v)
        by_var = {p.variable_name: p for p in cm.parameters()}
        params = {k: float(val) for k, val in (v.get('parameters') or {}).items()}
        overrides = {'parameters/' + harness.parameter_name(by_var[k]): val for k, val in params.items()}
        run_spec = dict(cm.spec, pre_time=0.0, sim_time=float(v.get('sim_time', np.max(t_data))),
                        dt=float(v.get('dt', cm.spec['dt'])))
        helper = harness.simulation_helper(cm.model_path, run_spec, solver_info={'rtol': 1e-8, 'atol': 1e-10})
        t, out = harness.run(helper, list(v['variables']), params=overrides)
        metric = v.get('metric', 'nrmse')
        threshold = float(v.get('threshold', 0.1))
        errors = {var: score(metric, np.interp(t_data, t, out[var]), data[var]) for var in v['variables']}
        fig = plots.plot_model_vs_data(plot_path(component, 'validation_baseline', cm.instance), t, out,
                                       {k: t_data for k in data}, data, cm.units(),
                                       f'{component.label} [{cm.instance.name}]: model vs {v.get("source_short", "data")}')
        metrics = {metric: errors, 'threshold': threshold, 'source': v.get('source', ''),
                   'validated_values': params}
        if v.get('kind') == 'fit_check':
            # the published values were fitted to this same data: goodness of fit, not validation
            metrics['kind'] = 'fit_check'
            metrics['note'] = v.get('note', '')
        worst = max(errors.values())
        if worst > threshold:
            return Result('validation_test_baseline', FAILED, f'worst {metric} {worst:.3f} > {threshold}', metrics, [fig])
        return Result('validation_test_baseline', PASSED, f'worst {metric} {worst:.3f} <= {threshold}', metrics, [fig])
    return _guard(cm.component, 'validation_test_baseline', check, cm.instance)


FEATURE_OPS = {
    'max': np.max, 'min': np.min, 'mean': np.mean, 'final': lambda y: y[-1],
    'ptp': np.ptp,
}


def _baseline_targets(cm, v):
    '''
    Scalar validation targets, e.g. steady-state clinical values, evaluated on the logged run
    (after pre_time). Spec (validation.<instance>.baseline):
        pre_time: 20            sim_time: 5            parameters: {var: value}   (optional)
        targets:
          - {name: LV end-diastolic volume, expr: 'np.max(q_lv)*1e6', value: 142, std: 21, units: ml}
          - {name: aortic systolic pressure, expr: 'np.max(vessel__u)/133.322', range: [100, 140], units: mmHg}
    Each expr is a numpy expression of t and the outputs (as in invariants). A target passes when
    |model - value| <= z_threshold * std, or when the model lies within range.
    '''
    by_var = {p.variable_name: p for p in cm.parameters()}
    params = {k: float(val) for k, val in (v.get('parameters') or {}).items()}
    overrides = {'parameters/' + harness.parameter_name(by_var[k]): val for k, val in params.items()}
    run_spec = dict(cm.spec, **{k: v[k] for k in ('sim_time', 'pre_time', 'dt') if k in v})
    helper = harness.simulation_helper(cm.model_path, run_spec)
    t, out = harness.run(helper, cm.outputs(), params=overrides)
    env = {'np': np, 't': t, '__builtins__': SAFE_BUILTINS}
    env.update(INVARIANT_HELPERS)
    env.update(cm.param_env({harness.parameter_name(by_var[k]): val for k, val in params.items()}))
    env.update(out)
    z = float(v.get('z_threshold', 2.0))
    rows, problems = [], []
    for tg in v['targets']:
        try:
            model = float(eval(tg['expr'], env))  # noqa: S307 - expressions come from the repo's own yaml
        except Exception as e:
            problems.append(f"{tg['name']}: could not evaluate ({e})")
            continue
        if 'range' in tg:
            lo, hi = map(float, tg['range'])
            ok = lo <= model <= hi
            target = f'[{lo:g}, {hi:g}]'
        else:
            ok = abs(model - float(tg['value'])) <= z * float(tg['std'])
            target = f"{float(tg['value']):g} ± {float(tg['std']):g}"
        rows.append({'name': tg['name'], 'model': model, 'target': target, 'units': tg.get('units', ''), 'ok': ok})
        if not ok:
            problems.append(f"{tg['name']}: model {model:.4g} vs {target} {tg.get('units', '')}")
    metrics = {'targets': rows, 'z_threshold': z, 'source': v.get('source', ''), 'validated_values': params}
    if problems:
        return Result('validation_test_baseline', FAILED, f'{len(problems)} of {len(v["targets"])} targets missed',
                      metrics, [], problems)
    return Result('validation_test_baseline', PASSED, f'all {len(rows)} targets met', metrics)


def validation_test_calibrate(cm):
    '''
    Calibrate to the instance's obs_data with libcuflynx parameter identification, then validate
    on its held-out data; writes the instance's calibrated parameters. An instance without
    obs_data is recorded as failed ("no calibration data in this instance").
    '''
    def check():
        early, v = _validation_status(cm, 'calibrate', 'validation_test_calibrate')
        if early:
            return early
        from cam_testing import calibrate
        return calibrate.run(cm, v, plot_path)
    return _guard(cm.component, 'validation_test_calibrate', check, cm.instance)


CHECKS = {
    'run_test': run_test,
    'verification_test_invariants': verification_test_invariants,
    'verification_test_BC': verification_test_BC,
    'verification_test_timestep': verification_test_timestep,
    'stability_test': stability_test,
    'validation_test_baseline': validation_test_baseline,
    'validation_test_calibrate': validation_test_calibrate,
}
