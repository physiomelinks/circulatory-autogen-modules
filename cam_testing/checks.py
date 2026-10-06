"""
The standard verification & validation tests, run on one module version at a time (the
validation tests on one of its instances).

Each check returns a Result and writes it to <version dir>/results/<test>.json (instance tests:
results/instances/<instance>/<test>.json), with any figures in <version dir>/plots/. The pytest
layer (tests/test_modules.py) turns a Result into pass / fail / skip; the report generator reads
the JSON files, so the HTML shows exactly what the last test run found.
"""
import ast
import contextlib
import datetime
import io
import json
import math
import os
import re
import tempfile
import time
import traceback
from dataclasses import asdict, dataclass, field

import numpy as np

from cam_testing import fixed_step, harness, plots
from cam_testing.library import library_dirs

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
NO_INSTANCE_CALIBRATION = 'no instance has calibration data'


# version_calibration of a version calibrated only as part of supermodules (below): pass-styled /
# fail-styled, but distinct from a plain pass / fail so it's clear the version wasn't calibrated alone
PASSED_IN_SUPER, FAILED_IN_SUPER = 'passed_in_super', 'failed_in_super'


def supermodule_index(versions=None):
    '''{(module_type, module_subtype): [supermodule versions whose submodules use it]} over ``versions``
    (default: the whole library).'''
    if versions is None:
        from cam_testing.library import all_versions
        versions = all_versions()
    index = {}
    for v in versions:
        if not v.is_supermodule:
            continue
        for sub in v.submodules:
            users = index.setdefault((sub['module_type'], sub['module_subtype']), [])
            if all(u.key != v.key for u in users):
                users.append(v)
    return index


def version_calibration(version, index=None, _seen=()):
    '''
    The version's calibration, from its instances: failed when no instance has calibration data
    (obs_data); otherwise failed if any instance's calibration failed, passed if they ran and
    passed, pending if not run yet.

    A version with no calibration data of its own that is a submodule of supermodule versions
    (``index``: supermodule_index(), built from the library when not given) is calibrated as part of
    them instead: "passed_in_super" when any of those supermodules' version_calibration passes
    (plainly or in its own supermodules: the rule applies transitively), "failed_in_super" when none
    does (a supermodule with no calibration data counts as failed). tests.yaml
    ``calibration_in_supermodule: false`` opts a version out (the plain rule).
    '''
    with_data = [i for i in version.instances() if i.has_obs_data]
    if not with_data:
        return _calibration_in_supermodules(version, index, _seen) \
            or Result('version_calibration', FAILED, NO_INSTANCE_CALIBRATION)
    results = {i.name: load(version, 'validation_test_calibrate', i) for i in with_data}
    failed = [n for n, r in results.items() if r is not None and r.status == FAILED]
    passed = [n for n, r in results.items() if r is not None and r.status == PASSED]
    if failed:
        return Result('version_calibration', FAILED, f'calibration failed: {", ".join(failed)}')
    if passed and len(passed) == len(results):
        return Result('version_calibration', PASSED, f'calibrated: {", ".join(passed)}')
    missing = [n for n, r in results.items() if r is None or r.status not in (PASSED, FAILED)]
    return Result('version_calibration', PENDING, f'calibration not run yet: {", ".join(missing)}')


def _calibration_in_supermodules(version, index, seen):
    if version.spec.get('calibration_in_supermodule', True) is False:
        return None
    if index is None:
        index = supermodule_index()
    seen = tuple(seen) + (version.key,)
    users = [s for s in index.get((version.vessel_type, version.name)) or [] if s.key not in seen]
    if not users:
        return None
    rows = []
    for s in users:
        r = version_calibration(s, index, seen)
        rows.append({'key': s.key, 'module_type': s.vessel_type, 'version': s.name, 'status': r.status,
                     'message': r.message})
    passing = [r for r in rows if r['status'] in (PASSED, PASSED_IN_SUPER)]
    listed = '; '.join(f'{r["key"]} {r["status"].replace("_", " ")} ({r["message"]})' for r in rows)
    metrics = {'supermodules': rows}
    if passing:
        return Result('version_calibration', PASSED_IN_SUPER,
                      f'no calibration data of its own; calibrated in supermodule {", ".join(r["key"] for r in passing)}: '
                      + listed, metrics)
    return Result('version_calibration', FAILED_IN_SUPER,
                  f'no calibration data of its own, and no supermodule using it passes calibration: {listed}', metrics)


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
        # the instance's own rows: what the model is generated with
        self._raw_params = component.parameters() if self.instance.is_default else self.instance.parameters()
        self._params = None
        if not component.is_supermodule:
            self._load_parameters()
        # A supermodule's parameters are those of its generated (flattened) model, read once it is
        # generated: its own instance's, its submodules' instances' and the boundary conditions no
        # internal connection closes (harness.supermodule_parameters).

    def _load_parameters(self):
        if self._params is not None:
            return
        if self.component.is_supermodule:
            params = harness.supermodule_parameters(self.component, self.model_path)
        else:
            params = self._raw_params
        nominal, var_of = {}, {}
        for p in params:
            var_of[harness.parameter_name(p)] = p.variable_name
            if p.is_todo:
                continue
            try:
                nominal[harness.parameter_name(p)] = p.float_value
            except ValueError:
                pass
        # the test network's own parameters (e.g. an outlet's flow v_vout), by their full names, so a
        # sweep can vary them (bc_sweep.extra_parameters) and invariants can use them
        for name, _units, value, *_ref in (self.spec.get('harness') or {}).get('parameters') or []:
            if name not in nominal:
                nominal[name] = float(value)
                var_of[name] = name
        self._params, self._nominal, self._var_of_map = params, nominal, var_of

    @property
    def nominal(self):
        '''Parameter values by their names in the generated model.'''
        self._load_parameters()
        return self._nominal

    @property
    def _var_of(self):
        self._load_parameters()
        return self._var_of_map

    def param_env(self, overrides=None):
        '''Parameter values by variable name (e.g. alpha; a supermodule's: <var>_<submodule>), for invariant expressions.'''
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
        '''The model's parameters: the instance's (a supermodule's: every parameter of its flattened
        model, named <var>_<submodule path> or as the global it is).'''
        self._load_parameters()
        return self._params

    def _instance_params(self):
        return None if self.instance.is_default else self._raw_params

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
        '''The spec's outputs; by default the version's variables (a supermodule's: every state of its model).'''
        wanted = self.spec.get('outputs')
        if not wanted:
            wanted = harness.state_outputs(self.helper) if self.component.is_supermodule \
                else [v[0] for v in self.component.variables()]
        return list(wanted)

    def units(self):
        if self.component.is_supermodule:
            return harness.supermodule_units(self.component, self.outputs())
        return {v[0]: v[1] for v in self.component.config.get('variables_and_units') or []}

    def run(self, params=None, helper=None):
        '''Runs at nominal parameters plus ``params``; always restores the nominal values.
        ``helper``: another simulation helper of this model (e.g. other solver settings).'''
        params = params or {}
        helper = helper or self.helper
        try:
            return harness.run(helper, self.outputs(), params={'parameters/' + k: v for k, v in params.items()})
        finally:
            if params:
                helper.set_param_vals(['parameters/' + k for k in params],
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
               'module_library_dirs': library_dirs(), 'use_builtin_modules': False, 'DEBUG': False,
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
        figs, details = [fig], [f'constant outputs: {", ".join(flat)}'] if flat else []
        rest = rest_check(cm)
        if rest:
            warning += f'; {rest["text"]}'
            metrics['rest_check'] = {k: v for k, v in rest.items() if k != 'fig'}
            details.append(rest['text'])
            figs += [rest['fig']] if rest.get('fig') else []
        return Result('run_test', PASSED, f'generated and simulated {cm.spec["sim_time"]} s; '
                      f'all {len(outputs)} outputs finite' + warning, metrics, figs, details)
    return _guard(cm.component, 'run_test', check)


def rest_check(cm):
    '''
    spec ``rest_check: {sim_time, window, voltage, parameters, threshold (mV, default 0), dt}``: a long
    run (e.g. at 0 pA injected) reported by run_test, information rather than a pass/fail gate: the
    spikes (upward crossings of the threshold by ``voltage``) and their rate in the last ``window`` s.
    '''
    rc = cm.spec.get('rest_check')
    if not rc:
        return None
    sim_time, window, thr = float(rc['sim_time']), float(rc.get('window', 10.0)), float(rc.get('threshold', 0.0))
    label = rc.get('label', f'rest ({sim_time:g} s)')
    try:
        spec = dict(cm.spec, sim_time=sim_time, dt=float(rc.get('dt', cm.spec['dt'])))
        helper = harness.simulation_helper(cm.model_path, spec)
        by_var = {p.variable_name: p for p in cm.parameters()}
        params = {'parameters/' + (harness.parameter_name(by_var[k]) if k in by_var else k): float(v)
                  for k, v in (rc.get('parameters') or {}).items()}
        t, out = harness.run(helper, [rc['voltage']], params)
    except Exception as e:
        return {'text': f'{label}: did not run ({type(e).__name__}: {str(e)[:200]})'}
    v = out[harness.output_key(rc['voltage'])]
    up = np.flatnonzero((v[:-1] < thr) & (v[1:] >= thr)) + 1
    n = int(np.sum(t[up] >= t[-1] - window))
    text = (f'{label}: fires at {n / window:.3g} Hz in the last {window:g} s ({n} spikes)' if n
            else f'{label}: silent in the last {window:g} s (no spikes; V ends at {v[-1]:.1f} mV)')
    fig = plots.plot_outputs(plot_path(cm.component, 'rest'), t, {rc['voltage']: v}, {rc['voltage']: 'mV'},
                             f'{cm.component.label}: {label}')
    return {'text': text, 'spikes_in_window': n, 'rate_Hz': n / window, 'window_s': window, 'sim_time_s': sim_time,
            'spikes_total': int(up.size), 'fig': fig}


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


# Valid ranges. A sweep value outside its parameter's valid range is not run; it is reported as
# skipped. The range is bc_sweep.bounds[param] when given, else a generic rule (below);
# bc_sweep.constraints are relations between parameters (e.g. "B_Ca_init < B_Ca_total").

# units of physically non-negative quantities: concentrations, amounts, volumes, conductances,
# capacitances, resistances, times (time constants) and absolute temperature
NON_NEGATIVE_UNITS = {
    'millimolar', 'micromolar', 'nanomolar', 'molar', 'mol_per_m3', 'mole', 'fmol',
    'm3', 'litre', 'microm3',
    'microS', 'nanoS', 'milliS', 'siemens',
    'picoF', 'nanoF', 'microF', 'farad',
    'megaOhm', 'kiloOhm', 'ohm',
    'second', 'millis',
    'kelvin',
}
GATE_SLACK = 1e-3      # an invariant's "x >= -eps" / "x <= 1 + eps" counts as [0, 1] for eps up to this
_NUM = r'\d+(?:\.\d*)?(?:[eE][-+]?\d+)?'
_LOWER0 = re.compile(rf'\b([A-Za-z_]\w*)\s*>=?\s*(-?\s*{_NUM})(?![\w.])(?!\s*[-+*/(])')
_UPPER1 = re.compile(rf'\b([A-Za-z_]\w*)\s*<=?\s*1(?:\.0*)?(?:\s*\+\s*({_NUM})(?![\w.]))?(?![\w.])(?!\s*[-+*/(])')


def unit_interval_states(spec):
    '''Variables the spec's invariants keep in [0, 1]: an invariant has both "x >= 0" (or >= -eps)
    and "x <= 1" (or <= 1 + eps), with eps <= GATE_SLACK. Used for the gate heuristic.'''
    lower, upper = set(), set()
    for inv in spec.get('invariants') or []:
        expr = inv['expr'] if isinstance(inv, dict) else inv
        for name, val in _LOWER0.findall(expr):
            if abs(float(val.replace(' ', ''))) <= GATE_SLACK:
                lower.add(name)
        for name, eps in _UPPER1.findall(expr):
            if not eps or float(eps) <= GATE_SLACK:
                upper.add(name)
    return lower & upper


def parameter_bounds(component, params=None):
    '''
    {variable: (lo, hi, why)}: the valid range of each parameter, for clipping the sweep.
      1. bc_sweep.bounds {variable: [min, max]} (null for no limit) when given;
      2. else a gate's initial value: a dimensionless "<x>_init" whose state x an invariant keeps
         in [0, 1] (unit_interval_states) lies in [0, 1];
      3. else a quantity in NON_NEGATIVE_UNITS with a non-negative value stays >= 0.
    ``params``: the parameters to bound (default: the version's default instance; a supermodule's
    model gives its flattened parameters, ComponentModel.parameters()).
    '''
    spec = component.spec
    gates = unit_interval_states(spec)
    out = {}
    for p in (component.parameters() if params is None else params):
        name, units = p.variable_name, (p.units or '').strip()
        try:
            value = None if p.is_todo else p.float_value
        except ValueError:
            value = None
        if units == 'dimensionless' and name.endswith('_init') and name[:-len('_init')] in gates:
            out[name] = (0.0, 1.0, f'gate initial value: an invariant keeps {name[:-len("_init")]} in [0, 1]')
        elif units in NON_NEGATIVE_UNITS and value is not None and value >= 0:
            out[name] = (0.0, math.inf, f'non-negative quantity ({units})')
    for name, (lo, hi) in ((spec.get('bc_sweep') or {}).get('bounds') or {}).items():
        out[name] = (-math.inf if lo is None else float(lo), math.inf if hi is None else float(hi), 'bc_sweep.bounds')
    return out


def point_skip_reason(var, value, bounds, constraints, env):
    '''Why a sweep point is not run (outside its valid range, or violating a constraint), or None.
    ``env``: the point's parameter values by variable name, for the constraints.'''
    if var in bounds:
        lo, hi, why = bounds[var]
        if not lo <= value <= hi:
            return f'outside valid range [{lo:g}, {hi:g}] ({why})'
    for expr in constraints:
        names = {'np': np, '__builtins__': SAFE_BUILTINS}
        names.update(env)
        try:
            ok = bool(np.all(eval(expr, names)))  # noqa: S307 - expressions come from the repo's own spec
        except Exception as e:
            raise ValueError(f'bc_sweep.constraints: could not evaluate {expr!r} ({type(e).__name__}: {e})')
        if not ok:
            return f'violates constraint {expr}'
    return None


def invariant_excess(spec, t, outputs, params, rtol, atol, only=None):
    '''
    For each invariant (in ``only``, when given), how far its comparisons miss: the largest
    violation over its "a <= b" / "a >= b" sub-expressions, as {expr, excess, bound, relative}
    with relative = excess / (atol + rtol*|bound|). A relative size of order 1-100 hints at solver
    error; a much larger one at the model itself. Expressions without comparisons are left out.
    '''
    env = {'np': np, 't': t, '__builtins__': SAFE_BUILTINS}
    env.update(INVARIANT_HELPERS)
    env.update(params or {})
    env.update(outputs)
    env['param'] = dict(params or {})
    out = []
    for inv in spec.get('invariants') or []:
        expr = inv['expr'] if isinstance(inv, dict) else inv
        if only is not None and expr not in only:
            continue
        try:
            tree = ast.parse(expr, mode='eval')
        except SyntaxError:
            continue
        worst = None
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Compare) and len(node.ops) == 1):
                continue
            op = node.ops[0]
            if not isinstance(op, (ast.Lt, ast.LtE, ast.Gt, ast.GtE)):
                continue
            try:
                with np.errstate(all='ignore'):
                    a = np.asarray(eval(compile(ast.Expression(node.left), '<inv>', 'eval'), env), float)  # noqa: S307
                    b = np.asarray(eval(compile(ast.Expression(node.comparators[0]), '<inv>', 'eval'), env), float)  # noqa: S307
                    miss = (a - b) if isinstance(op, (ast.Lt, ast.LtE)) else (b - a)
                    k = np.nanargmax(miss) if np.size(miss) and np.any(np.isfinite(miss)) else None
            except Exception:
                continue
            if k is None:
                continue
            m = float(np.ravel(miss)[k])
            if not m > 0:
                continue
            bound = float(np.max(np.abs(b))) if np.size(b) else 0.0
            rel = m / (atol + rtol * bound)
            if worst is None or rel > worst['relative']:
                worst = {'excess': m, 'bound': bound, 'relative': rel}
        if worst:
            out.append(dict(expr=expr, **worst))
    return out


def cheapest_passing(n, passes):
    '''
    Binary search of a cost-ordered list of n candidates for the cheapest with passes(i) true,
    assuming that once a candidate passes the more expensive (more accurate) ones do too.
    Returns (index or None, {index: result} of the candidates evaluated).
    '''
    lo, hi, seen = 0, n, {}
    while lo < hi:
        mid = (lo + hi) // 2
        if mid not in seen:
            seen[mid] = bool(passes(mid))
        if seen[mid]:
            hi = mid
        else:
            lo = mid + 1
    return (lo if lo < n else None), seen


PARAMETER_FAILURE, NUMERICAL_FAILURE = 'parameter', 'numerical'
HARD_POINTS_CAP = 10           # points a candidate solver setting is evaluated on before the full confirmation
HARD_FAILURES_CAP = 8          # of which at most this many are numerical failures (leaving room for nominal/extremes)


def bc_reference_config(cm):
    '''The tight reference a failed sweep point is rerun with: the stability test's reference
    tolerances (reference_tolerances) and a MaximumStep of a tenth of the output interval dt.'''
    step = float(cm.spec['dt']) / 10.0
    given = (cm.spec.get('solver_info') or {}).get('MaximumStep')
    if given:
        step = min(step, float(given))
    return dict({'solver': 'CVODE_myokit'}, **reference_tolerances(cm), MaximumStep=step)


def solver_candidates(cm):
    '''
    CVODE settings to try where the version's own fail numerically: the CVODE configurations the
    stability test found working (results/stability_test.json, if run), then a default ladder,
    rtol 1e-6 / 1e-8 / 1e-10 (atol rtol/100, or the version's atol if smaller) x MaximumStep
    none / 1e-4 / 1e-5. Duplicates removed; they are ordered by measured cost later.
    '''
    cands = []
    r = load(cm.component, 'stability_test')
    for row in ((r.metrics or {}).get('matrix') or []) if r is not None else []:
        cfg = row.get('config') or {}
        if row.get('works') and cfg.get('solver') == 'CVODE_myokit':
            cands.append(dict(cfg))
    atol_v = float((cm.spec.get('solver_info') or {}).get('atol', 1e-8))
    for rtol in (1e-6, 1e-8, 1e-10):
        for step in (None, 1e-4, 1e-5):
            cfg = {'solver': 'CVODE_myokit', 'rtol': rtol, 'atol': min(rtol * 1e-2, atol_v)}
            if step:
                cfg['MaximumStep'] = step
            cands.append(cfg)
    out, seen = [], set()
    for c in cands:
        if config_label(c) not in seen:
            seen.add(config_label(c))
            out.append(c)
    return out


def _stricter(a, b):
    '''Config a asks at least as much of CVODE as b (tighter or equal rtol and MaximumStep).'''
    inf = math.inf
    return (float(a.get('rtol', inf)) <= float(b.get('rtol', inf))
            and float(a.get('atol', inf)) <= float(b.get('atol', inf))
            and float(a.get('MaximumStep') or inf) <= float(b.get('MaximumStep') or inf))


def classify_failures(failed, reference):
    '''Classifies each failed point by one rerun at the reference settings (``reference(pname,
    value)`` -> a _sweep_point result): "parameter" when the reference fails too (the model breaks
    there), "numerical" when it passes (the version's solver settings are inadequate there).'''
    for p in failed:
        ref = reference(p['pname'], p['value'])
        p['class'] = NUMERICAL_FAILURE if ref['ok'] else PARAMETER_FAILURE
        if not ref['ok']:
            p['reference_problems'] = ref['problems']
    return failed


def _cvode_helper(cm, cfg):
    info = {k: v for k, v in cfg.items() if k != 'solver'}
    return harness.simulation_helper(cm.model_path, cm.spec, solver_info=info, solver='CVODE_myokit')


def _sweep_point(cm, helper, pname, value):
    '''Runs one sweep point (pname None: the nominal point) with a helper (None: the version's
    settings). Returns {ok, problems, t, out, bad, seconds}.'''
    params = {} if pname is None else {pname: value}
    start = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            t, out = cm.run(params, helper=helper)
    except Exception as e:
        return {'ok': False, 'problems': [f'simulation failed ({type(e).__name__}: {str(e)[:300]})'], 't': None,
                'out': None, 'bad': [], 'seconds': time.perf_counter() - start}
    seconds = time.perf_counter() - start
    bad = _non_finite(out)
    with np.errstate(all='ignore'):
        inv = _check_invariants(cm.spec, t, out, cm.param_env(params))
    problems = ([f'non-finite {", ".join(bad)}'] if bad else []) + inv
    return {'ok': not problems, 'problems': problems, 't': t, 'out': out, 'bad': bad, 'seconds': seconds}


def _accuracy(res, ref):
    '''Max over outputs and time of |run - reference| / the reference's magnitude (as the stability test).'''
    if res['out'] is None or ref['out'] is None:
        return math.inf
    keys = list(ref['out'])
    R = np.column_stack([ref['out'][k] for k in keys])
    scale = np.maximum(np.max(np.abs(R), axis=0), np.ptp(R, axis=0))
    scale[scale == 0] = 1.0
    with np.errstate(all='ignore'):
        Y = np.column_stack([np.interp(ref['t'], res['t'], res['out'][k]) for k in keys])
    if not np.all(np.isfinite(Y)):
        return math.inf
    return float(np.max(np.abs(Y - R) / scale))


def _recommend_settings(cm, points, numerical, reference, tol, version_cfg):
    '''
    The cheapest CVODE setting that fixes the numerical failures: candidates (solver_candidates)
    at least as strict as the version's settings (no looser rtol, atol or MaximumStep: a looser
    setting that happens to pass is no fix) and accurate at the nominal point are ordered by
    measured run time; a binary search finds the
    cheapest that passes the hard points (the numerical failures, nominal, and each parameter's
    extreme values) within ``tol`` of the reference; the winner is then confirmed over every sweep
    point. Runs predicted over the stability budget (20 s per state) are not made.
    Returns (recommended or None, [candidate rows], [notes]).
    '''
    n_states = max(cm.helper.model.count_states(), 1)
    budget = float((cm.spec.get('stability') or {}).get('time_budget') or 20.0 * n_states)
    usable = [p for p in points if p.get('class') != PARAMETER_FAILURE]
    hard = []

    def add(pname, value, cap=HARD_POINTS_CAP):
        if (pname, value) not in hard and len(hard) < cap:
            hard.append((pname, value))

    for p in numerical:
        add(p['pname'], p['value'], HARD_FAILURES_CAP)
    add(None, None)
    by_param = {}
    for p in usable:
        by_param.setdefault(p['pname'], []).append(p['value'])
    for pname, values in by_param.items():
        add(pname, min(values))
        add(pname, max(values))

    rows, notes, over = [], [], []
    measured = []
    ref_nominal = reference(None, None)
    for cfg in solver_candidates(cm):
        if not _stricter(cfg, version_cfg) or config_label(cfg) == config_label(version_cfg):
            continue
        row = {'config': cfg, 'label': config_label(cfg), 'seconds': None, 'passed_hard': None, 'note': ''}
        rows.append(row)
        if any(_stricter(cfg, o) for o in over):
            row['note'] = f'not run: predicted over the {budget:.0f} s budget ({n_states} states x 20 s)'
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                helper = _cvode_helper(cm, cfg)
        except Exception as e:
            row['note'] = f'{type(e).__name__}: {str(e)[:200]}'
            continue
        res = _sweep_point(cm, helper, None, None)
        if res['seconds'] > budget:
            over.append(cfg)
            row['note'] = f'run took {res["seconds"]:.0f} s > {budget:.0f} s budget'
            continue
        res2 = _sweep_point(cm, helper, None, None) if res['ok'] else res     # second timing: less noise
        row['seconds'] = min(res['seconds'], res2['seconds'])
        acc = _accuracy(res, ref_nominal) if ref_nominal['ok'] else math.nan
        if not res['ok'] or not acc <= tol:
            row['note'] = 'fails at the nominal point' if not res['ok'] else f'inaccurate at the nominal point ({acc:.1e} > tol)'
            continue
        measured.append({'row': row, 'helper': helper, 'cfg': cfg})
    # cost order; ties (equal times) go to the less strict setting first
    strictness = lambda c: (-math.log10(float(c['cfg'].get('rtol', 1))), -math.log10(float(c['cfg'].get('MaximumStep') or 1e3)))
    measured.sort(key=lambda c: (round(c['row']['seconds'], 3), strictness(c)))

    def passes(i):
        c = measured[i]
        worst = 0.0
        for k, (pname, value) in enumerate(hard):
            res = _sweep_point(cm, c['helper'], pname, value)
            if res['seconds'] > budget:
                c['row']['note'] = f'run took {res["seconds"]:.0f} s > {budget:.0f} s budget'
                c['row']['passed_hard'] = False
                return False
            ref = reference(pname, value)
            acc = _accuracy(res, ref) if ref['ok'] else 0.0
            if not res['ok'] or not acc <= tol:
                c['row']['passed_hard'] = False
                c['row']['note'] = (f'fails at hard point {k + 1}/{len(hard)}'
                                    + ('' if res['ok'] else f' ({res["problems"][0][:120]})')
                                    + ('' if not res['ok'] else f' (difference {acc:.1e} > tol)'))
                return False
            worst = max(worst, acc)
        c['row']['passed_hard'] = True
        c['row']['accuracy'] = worst
        c['row']['note'] = f'passes all {len(hard)} hard points'
        return True

    start_at = 0
    for _attempt in range(3):
        idx, _seen = cheapest_passing(len(measured) - start_at, lambda i: passes(start_at + i))
        if idx is None:
            notes.append('no candidate setting passes the hard points')
            return None, rows, notes
        c = measured[start_at + idx]
        fails, secs = [], []
        for p in usable:
            res = _sweep_point(cm, c['helper'], p['pname'], p['value'])
            secs.append(res['seconds'])
            if not res['ok']:
                fails.append(p)
        c['row']['confirmation'] = f'{len(fails)}/{len(usable)} failures'
        if not fails:
            return {'config': c['cfg'], 'label': c['row']['label'],
                    'mean_seconds': float(np.mean(secs)) if secs else None,
                    'n_failures': 0, 'n_runs': len(usable), 'failures': f'0/{len(usable)}',
                    'accuracy': c['row'].get('accuracy'), 'tol': tol, 'n_hard_points': len(hard),
                    'budget_seconds': budget}, rows, notes
        notes.append(f'{c["row"]["label"]} passed the hard points but failed {len(fails)} of {len(usable)} '
                     f'sweep points; those points were added to the hard points')
        for p in fails:
            add(p['pname'], p['value'], cap=len(hard) + len(fails))
        start_at += idx + 1
        if start_at >= len(measured):
            break
    notes.append('no candidate setting passes every sweep point')
    return None, rows, notes


# bc_sweep.sweep: what the boundary-condition sweep varies
#   all             every boundary condition and constant (a component version's default)
#   bcs             the boundary conditions only
#   globals_and_bcs the global constants and the boundary conditions (a supermodule version's
#                   default: each submodule's own parameters are swept in that submodule version's
#                   own tests, so the supermodule sweeps what it adds -- its globals and the boundary
#                   conditions no internal connection closes)
# bc_sweep.extra_parameters are swept as well, whatever the mode.
SWEEP_MODES = ('all', 'bcs', 'globals_and_bcs')


def default_sweep_mode(component):
    return 'globals_and_bcs' if component.is_supermodule else 'all'


def sweep_parameter_names(cm):
    '''(variable names to sweep, what they are) for the version's bc_sweep spec.'''
    component, sweep_spec = cm.component, cm.spec.get('bc_sweep') or {}
    params = cm.parameters()
    mode = sweep_spec.get('sweep', default_sweep_mode(component))
    if mode not in SWEEP_MODES:
        raise ValueError(f'bc_sweep.sweep {mode!r}: expected one of {", ".join(SWEEP_MODES)}')
    bcs = [p.variable_name for p in params if p.kind == 'boundary_condition']
    if mode == 'bcs':
        names, swept_kind = bcs, 'boundary conditions'
    elif mode == 'globals_and_bcs':
        globals_ = [p.variable_name for p in params if p.kind == 'global_constant' and p.variable_name not in bcs]
        names = globals_ + bcs
        swept_kind = 'global constants and boundary conditions (submodule parameters are swept in their own versions)'
    else:
        names = bcs + [p.variable_name for p in params if p.kind != 'boundary_condition' and p.variable_name not in bcs]
        swept_kind = 'boundary conditions and constants' if bcs else 'constants (no boundary conditions)'
    names = [n for n in names if n not in (sweep_spec.get('exclude') or [])]
    names += [n for n in sweep_spec.get('extra_parameters', []) if n not in names]
    return names, swept_kind


def verification_test_BC(cm):
    if is_cpp(cm.component):
        return _cpp_not_applicable(cm, 'verification_test_BC')
    def check():
        component, spec = cm.component, cm.spec
        sweep_spec = spec.get('bc_sweep') or {}
        cm.helper           # generate once: a model that doesn't generate fails here, not at every point
        names, swept_kind = sweep_parameter_names(cm)
        if not names:
            return Result('verification_test_BC', SKIPPED, 'component has no boundary conditions or constants to sweep')
        by_var = {p.variable_name: p for p in cm.parameters()}
        plot_output = sweep_spec.get('plot_output') or cm.outputs()[0]

        bounds = parameter_bounds(component, cm.parameters())
        constraints = list(sweep_spec.get('constraints') or [])
        si = spec.get('solver_info') or {}
        version_cfg = dict({'solver': spec.get('solver', 'CVODE_myokit')}, **si)
        rtol, atol = float(si.get('rtol', 1e-6)), float(si.get('atol', 1e-8))

        sweeps, points, skipped, notes = {}, [], [], []
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
                value = float(value)
                why = point_skip_reason(var, value, bounds, constraints, cm.param_env({pname: value}))
                if why:
                    skipped.append({'parameter': var, 'value': value, 'reason': why})
                    continue
                res = _sweep_point(cm, None, pname, value)
                point = {'var': var, 'pname': pname, 'value': value, 'ok': res['ok'], 'problems': res['problems']}
                if not res['ok'] and res['out'] is not None:
                    failed_inv = {m.rsplit(': ', 1)[0] for m in res['problems']}
                    point['violation'] = invariant_excess(spec, res['t'], res['out'], cm.param_env({pname: value}),
                                                          rtol, atol, only=failed_inv)
                points.append(point)
                runs.append((value, res['t'], res['out'][harness.output_key(plot_output)]
                             if res['out'] is not None and not res['bad'] else None))
            if runs:
                sweeps[var] = runs

        # classify each failure with one run at the reference settings
        ref_cfg = bc_reference_config(cm)
        ref_cache, ref_helper = {}, []

        def reference(pname, value):
            if (pname, value) not in ref_cache:
                if not ref_helper:
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        ref_helper.append(_cvode_helper(cm, ref_cfg))
                ref_cache[(pname, value)] = _sweep_point(cm, ref_helper[0], pname, value)
            return ref_cache[(pname, value)]

        failed = classify_failures([p for p in points if not p['ok']], reference)
        numerical = [p for p in failed if p['class'] == NUMERICAL_FAILURE]
        parameter = [p for p in failed if p['class'] == PARAMETER_FAILURE]

        recommended, candidates, search_notes = None, [], []
        if numerical:
            tol = float((spec.get('stability') or {}).get('tol', DEFAULT_STABILITY['tol']))
            recommended, candidates, search_notes = _recommend_settings(cm, points, numerical, reference, tol,
                                                                         version_cfg)
        fixed = recommended is not None

        figs = []
        if sweeps:
            figs.append(plots.plot_sweep(plot_path(component, 'bc_sweep'), sweeps, plot_output, cm.units(),
                                         f'{component.label}: boundary-condition sweep'))

        def fmt_point(p):
            return f'{p["var"]}={p["value"]:.4g}'

        failure_rows = []
        for p in failed:
            failure_rows.append({'parameter': p['var'], 'value': p['value'], 'class': p['class'],
                                 'problems': p['problems'], 'violation': p.get('violation') or [],
                                 'reference_problems': p.get('reference_problems') or [],
                                 'fixed_by_recommended': bool(fixed and p['class'] == NUMERICAL_FAILURE)})
        details = []
        for p in failed:
            hint = ''
            if p.get('violation'):
                v = max(p['violation'], key=lambda x: x['relative'])
                hint = f' [violation {v["excess"]:.2g} = {v["relative"]:.2g} x (atol + rtol*|bound|)]'
            tag = p['class'] + (', fixed by the recommended settings' if fixed and p['class'] == NUMERICAL_FAILURE else '')
            details += [f'[{tag}] {fmt_point(p)}: {m}{hint if m.endswith("violated") else ""}' for m in p['problems']]
        details += [f'[skipped] {s["parameter"]}={s["value"]:.4g}: {s["reason"]}' for s in skipped]
        details += search_notes + notes

        n_runs = len(points)
        metrics = {'swept_kind': swept_kind, 'swept': list(sweeps), 'n_runs': n_runs, 'n_failures': len(failed),
                   'sweep_mode': sweep_spec.get('sweep', default_sweep_mode(component)),
                   'sweep_values': {k: [r[0] for r in v] for k, v in sweeps.items()},
                   'skipped_out_of_range': skipped,
                   'bounds': {k: [None if math.isinf(lo) else lo, None if math.isinf(hi) else hi, why]
                              for k, (lo, hi, why) in bounds.items() if k in names},
                   'constraints': constraints,
                   'failures': failure_rows,
                   'n_parameter_failures': len(parameter), 'n_numerical_failures': len(numerical),
                   'version_settings': config_label(version_cfg), 'reference_settings': config_label(ref_cfg),
                   'recommended': recommended, 'candidates': candidates}

        skip_note = (f'; {len(skipped)} point{"" if len(skipped) == 1 else "s"} skipped (outside the valid range)'
                     if skipped else '')
        if not sweeps:
            return Result('verification_test_BC', SKIPPED, 'no boundary condition could be swept' + skip_note,
                          metrics, figs, details)
        if parameter or (numerical and not fixed):
            parts = []
            if parameter:
                parts.append(f'{len(parameter)} parameter failure{"" if len(parameter) == 1 else "s"} '
                             f'(the reference solver fails there too: the model breaks)')
            if numerical:
                parts.append(f'{len(numerical)} numerical failure{"" if len(numerical) == 1 else "s"} '
                             f'(the reference solver passes there) '
                             + (f'that the recommended settings {recommended["label"]} fix' if fixed
                                else 'that no candidate solver setting fixes'))
            return Result('verification_test_BC', FAILED, f'{len(failed)} of {n_runs} runs failed: ' + '; '.join(parts)
                          + skip_note, metrics, figs, details)
        msg = f'{n_runs} runs over {len(sweeps)} {swept_kind}: '
        if numerical:
            one = len(numerical) == 1
            msg += (f'{len(numerical)} run{"" if one else "s"} failed at the version\'s solver settings '
                    f'({config_label(version_cfg)}) but pass{"es" if one else ""} with the reference solver '
                    f'(numerical failure{"" if one else "s"}), and every point passes with the recommended '
                    f'settings {recommended["label"]}: a recommendation to adopt into the version\'s solver_info, '
                    f'not yet applied')
        else:
            msg += 'all finite, invariants hold' if spec.get('invariants') else 'all finite (no invariants defined in the spec)'
        return Result('verification_test_BC', PASSED, msg + skip_note, metrics, figs, details)
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
        # an instance without obs_data has nothing to calibrate to: not applicable for that
        # instance. A version with no calibration data in any instance fails the version-level
        # calibration (version_calibration, in the report's test overview).
        why = v.get('note') or v.get('reason')
        return Result(test, NOT_APPLICABLE, NO_CALIBRATION_DATA, details=[why] if why else []), v
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
    on its held-out data; writes the instance's calibrated parameters. Not applicable to an
    instance without obs_data ("no calibration data in this instance"); see version_calibration.
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


# ----------------------------------------------------------------------------------------------
# unit consistency (libcellml): informational, not a test column; no simulation
# ----------------------------------------------------------------------------------------------

UNIT_CONSISTENCY = 'unit_consistency'
UNITS_CONSISTENT, UNITS_INCONSISTENT, UNITS_NOT_CHECKED = 'consistent', 'inconsistent', 'not_checked'
# libcellml 0.6 ANALYSER_UNITS messages: "The units in '<expr>'[ in '<parent>' ...][ in equation '<eq>'] in component
# '<c>' are not equivalent. '<a>' is in '<u>' (i.e. '<base>') while '<b>' is '<u>'." and "The unit of '<expr>' ... is not
# dimensionless. '<expr>' is in '<u>'."
_UNITS_ISSUE = re.compile(r"^The units? (?:in|of) '(?P<expression>[^']*)'(?P<chain>.*?) in component '(?P<component>[^']*)' "
                          r"(?:are|is) (?P<what>[^.]*)\.\s*(?P<detail>.*)$", re.S)
_UNITS_EQUATION = re.compile(r" in equation '(?P<equation>[^']*)'")
_UNITS_SIDE = re.compile(r"'(?P<expression>[^']*)' is (?:in )?'(?P<units>[^']*)'(?: \(i\.e\. '(?P<base>[^']*)'\))?")


def unit_consistency_path(version):
    return os.path.join(version.results_dir, f'{UNIT_CONSISTENCY}.json')


def _units_model(version):
    '''The version's CellML component as a stand-alone libcellml model: its units file's units added,
    every input made a local constant (initial value 1) unless an equation defines it or it is the
    variable of integration, so the analyser can check the equations' units on their own.'''
    import xml.etree.ElementTree as ET
    import libcellml
    from cam_testing.mathml import CELLML_NS, MATHML_NS
    parser = libcellml.Parser(False)        # CellML 1.1 -> 2.0
    with open(version.cellml_path) as f:
        model = parser.parseModel(f.read())
    if os.path.isfile(version.units_path):
        with open(version.units_path) as f:
            units_model = parser.parseModel(f.read())
        for i in range(units_model.unitsCount()):
            u = units_model.units(i)
            if not model.hasUnits(u.name()):
                model.addUnits(u.clone())
    model.linkUnits()
    if not re.match(r'[A-Za-z_]', model.name() or ''):
        model.setName('m_' + (model.name() or 'model'))   # a CellML 2.0 identifier can't start with a digit
    root = ET.parse(version.cellml_path).getroot()
    defined = {}
    from cam_testing.mathml import equation_target
    for comp in root.iter(f'{{{CELLML_NS}}}component'):
        bvars = {(ci.text or '').strip() for b in comp.iter(f'{{{MATHML_NS}}}bvar') for ci in b.iter(f'{{{MATHML_NS}}}ci')}
        targets = {equation_target(eq)[0] for m in comp.iter(f'{{{MATHML_NS}}}math') for eq in m}
        defined[comp.get('name')] = bvars | targets
    for ci in range(model.componentCount()):
        c = model.component(ci)
        keep = defined.get(c.name(), set())
        for vi in range(c.variableCount()):
            var = c.variable(vi)
            var.removeInterfaceType()
            init = var.initialValue()
            if var.name() not in keep and not init:
                var.setInitialValue('1')
            elif init in keep:
                # a state initialised by a computed variable (e.g. a gate's steady state): the analyser
                # stops at that before checking any units, and an initial value doesn't enter the equations
                var.setInitialValue('1')
    return model


def unit_consistency(version, save=True):
    '''
    libcellml's unit checks on the version's CellML component (no simulation): the Validator's units
    issues (undefined or invalid units) and the Analyser's ANALYSER_UNITS issues (the two sides of an
    equation, or the terms of a sum, in units that are not equivalent or differ by a scaling factor).
    Each analyser issue is mapped to the equation it is about (equation_index, in the order of
    mathml.component_equations, and the variable that equation defines). Written to
    results/unit_consistency.json; informational, not a test column.
    '''
    import libcellml
    from cam_testing.mathml import component_equation_targets
    out = {'status': UNITS_NOT_CHECKED, 'message': '', 'n_failures': 0, 'n_equations': 0, 'failures': [],
           'definition_issues': [], 'other_issues': [], 'libcellml_version': libcellml.versionString()}
    if version.is_supermodule:
        out.update(status=NOT_APPLICABLE, message='a supermodule has no CellML of its own: each part is checked in its '
                                                  'own version')
    elif version.format != 'cellml' or not os.path.isfile(version.cellml_path):
        out.update(status=NOT_APPLICABLE, message=f'no CellML component to check (module_format {version.format})')
    else:
        try:
            out.update(_unit_issues(version, component_equation_targets(version.cellml_path, version.module_type)))
        except Exception as e:  # noqa: BLE001
            out.update(status=UNITS_NOT_CHECKED, message=f'libcellml could not check the units: {type(e).__name__}: {e}')
    out['timestamp'] = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    if save:
        path = unit_consistency_path(version)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            json.dump(out, f, indent=2)
    return out


def _unit_issues(version, targets):
    import libcellml
    model = _units_model(version)
    rule = libcellml.Issue.ReferenceRule
    error = libcellml.Issue.Level.ERROR
    validator = libcellml.Validator()
    validator.validateModel(model)
    definition = []
    for i in range(validator.issueCount()):
        it = validator.issue(i)
        name = next((n for n in dir(rule) if getattr(rule, n) == it.referenceRule() and n.isupper()), '')
        if it.level() == error and 'UNIT' in name:
            definition.append({'rule': name, 'message': it.description()})
    analyser = libcellml.Analyser()
    analyser.analyseModel(model)
    failures, other = [], []
    for i in range(analyser.issueCount()):
        it = analyser.issue(i)
        if it.referenceRule() == rule.ANALYSER_UNITS:
            failures.append(_unit_failure(it.description(), targets))
        elif it.level() == error:
            other.append(it.description())
    n_equations = len({f['equation'] or f['message'] for f in failures})
    if failures or definition:
        status = UNITS_INCONSISTENT
    elif other:
        # libcellml's analyser checks units only once the model is otherwise valid
        return {'status': UNITS_NOT_CHECKED, 'n_failures': 0, 'n_equations': 0, 'failures': [], 'definition_issues': [],
                'other_issues': other,
                'message': 'libcellml could not check the equations\' units: ' + '; '.join(other[:3])
                           + (f' (and {len(other) - 3} more)' if len(other) > 3 else '')}
    else:
        return {'status': UNITS_CONSISTENT, 'message': 'every equation is unit-consistent', 'n_failures': 0,
                'n_equations': 0, 'failures': [], 'definition_issues': [], 'other_issues': []}
    parts = []
    if failures:
        parts.append(f'{n_equations} equation{"" if n_equations == 1 else "s"} with inconsistent units'
                     + (f' ({len(failures)} issues)' if len(failures) != n_equations else ''))
    if definition:
        parts.append(f'{len(definition)} undefined or invalid units' + (
            ', so libcellml could not check the equations' if not failures and other else ''))
    return {'status': status, 'message': '; '.join(parts), 'n_failures': n_equations + len(definition),
            'n_equations': n_equations, 'failures': failures, 'definition_issues': definition, 'other_issues': other}


def _unit_failure(description, targets):
    '''One ANALYSER_UNITS issue -> {equation, variable, equation_index, kind, sides, message}.'''
    m = _UNITS_ISSUE.match(description.strip())
    row = {'message': description.strip(), 'equation': None, 'expression': None, 'variable': None,
           'equation_index': None, 'kind': 'not equivalent', 'sides': []}
    if not m:
        return row
    # a sub-expression's issue names its equation: "The units in '<sub>' in '<...>' in equation '<eq>' ..."
    e = _UNITS_EQUATION.search(m['chain'])
    eq = e['equation'] if e else m['expression']
    row['equation'] = eq
    row['expression'] = m['expression'] if e else None
    what = m['what']
    row['kind'] = ('not dimensionless' if 'dimensionless' in what else
                   'scaling factor' if 'scaling' in what or 'multiplication' in m['detail'] else 'not equivalent')
    row['sides'] = [{'expression': x['expression'], 'units': x['units'], 'base': x['base']}
                    for x in _UNITS_SIDE.finditer(m['detail'])]
    lhs = eq.split(' = ', 1)[0].strip()
    r = re.fullmatch(r'd(\w+)/d\w+', lhs)
    var = r[1] if r else (lhs if re.fullmatch(r'\w+', lhs) else None)
    row['variable'] = var
    if var is not None:
        row['equation_index'] = next((k for k, (v, _) in enumerate(targets) if v == var), None)
    return row


def load_unit_consistency(version):
    '''The recorded unit check (results/unit_consistency.json), recomputed when missing or older than the
    version's CellML, units or config file.'''
    path = unit_consistency_path(version)
    sources = [p for p in (version.cellml_path, version.units_path, version.config_path) if os.path.isfile(p)]
    if os.path.isfile(path) and all(os.path.getmtime(path) >= os.path.getmtime(p) for p in sources):
        with open(path) as f:
            return json.load(f)
    return unit_consistency(version)
