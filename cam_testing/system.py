"""
System models: system_models/<category>/<model>/ (at the repo root; not modules, not supermodules).

Each system dir holds the model built from this module library (<model>_module_array.json, whose
records (instances) name a module_type, its version (module_subtype) and the parameterisation
it uses), its parameter sets under parameterisations/ (values that win over the module
parameterisations'), the reference it must reproduce under reference/, and a spec <model>_system.yaml.

Parameter sets:

  parameterisations/default/default_parameters.json
      the model's values: what generation, the tests and users get unless the spec's
      "parameterisation" names another set. For a model with a calibration workflow it is the
      workflow's result (python -m cam_testing.system write-default), and nothing else writes it;
      default_calibration.json records the inputs it was calibrated from.
  parameterisations/uncalibrated/
      only for a model with a calibration workflow: the workflow's input, uncalibrated_parameters.json,
      with <model>_calibration_workflow.json, its obs_data, params_for_id, protocols/ and SOURCES.md.
      A module recalibration (python -m cam_testing.calibrate apply) writes here, never into default,
      and leaves the workflow stale (workflow_stale) until it is run again. The model is generated with every module library
(cam_testing.paths: this repo's modules/, then the extra ones), so its records may name versions of
either. The reference is one of:

  - circulatory_autogen's original: its module array (and parameters) under reference/, and the model
    generated from them with libcuflynx's built-in modules, committed as reference/<model>.cellml
    (python -m cam_testing.system write-reference, while libcuflynx still has its built-ins), and
    simulated as a generated model;
  - a ready CellML model, simulated as it is (e.g. the flattened model PhLynx exports):
    equivalence.reference_cellml (a path relative to the system dir), or reference/<model>.cellml
    when reference/ holds no module array. Its outputs are "<component>/<variable>" for every
    component but those in equivalence.ignore_components (default: environment,
    global_parameters, instance_parameters), mapped through output_map as for an original.

The tests:

  system_run_test          the model generates with this library and runs; outputs finite
  system_equivalence_test  the model reproduces its reference (circulatory_autogen's original, as
                           committed CellML, or a ready CellML model): every
                           variable the reference logs,
                           mapped through equivalence.output_map (e.g. heart/u_lv -> lv/u),
                           agrees to equivalence.tol (normalised max difference) with both
                           simulated at equivalence.solver_info
  system_invariants_test   the spec's invariants (as for modules)
  system_workflow_current_test
                           a model with a calibration workflow: its default set is the workflow's
                           result for the inputs in uncalibrated/ now (workflow_stale)
"""
import contextlib
import glob
import io
import json
import os
import tempfile
from dataclasses import dataclass

import numpy as np
import yaml

from cam_testing import paths
from cam_testing.library import library_dirs


def system_dir():
    '''system_models/ of the repo under test (cam_testing.paths).'''
    return paths.roots().system_models_dir


def legacy_system_dir():
    '''System models not moved yet (another session's modules/system/poiseuille): read where they are.'''
    return os.path.join(paths.roots().modules_dir, 'system')


def __getattr__(name):
    # the former constants
    if name == 'SYSTEM_DIR':
        return system_dir()
    if name == 'LEGACY_SYSTEM_DIR':
        return legacy_system_dir()
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


@dataclass
class System:
    name: str
    category: str
    dir: str
    spec: dict

    @property
    def id(self):
        return f'{self.category}/{self.name}'

    @property
    def results_dir(self):
        return os.path.join(self.dir, 'results')

    @property
    def plots_dir(self):
        return os.path.join(self.dir, 'plots')

    @property
    def parameterisation(self):
        '''The parameter set the model is generated with: the spec's "parameterisation", else default.'''
        return self.spec.get('parameterisation') or DEFAULT

    @property
    def parameters_path(self):
        return parameters_path(self.dir, self.parameterisation)

    @property
    def workflow_path(self):
        '''uncalibrated/<model>_calibration_workflow.json, None when the model has no workflow.'''
        path = os.path.join(self.dir, PARAMETERISATIONS, UNCALIBRATED, f'{self.name}_calibration_workflow.json')
        return path if os.path.isfile(path) else None

    @property
    def calibration_dir(self):
        '''Where the model's calibration data are (obs_data, SOURCES.md): uncalibrated/ when the model has a
        workflow, else the model's directory.'''
        return os.path.dirname(self.workflow_path) if self.workflow_path else self.dir


PARAMETERISATIONS, DEFAULT, UNCALIBRATED = 'parameterisations', 'default', 'uncalibrated'
# beside default_parameters.json: the sha256 of each workflow input it was calibrated from
CALIBRATION_RECORD = 'default_calibration.json'


def parameters_path(model_dir, parameterisation=DEFAULT):
    '''parameterisations/<p>/<p>_parameters.json of the system model in model_dir.'''
    return os.path.join(model_dir, PARAMETERISATIONS, parameterisation, f'{parameterisation}_parameters.json')


def _workflow_inputs(system):
    '''{path relative to uncalibrated/: sha256} of every file there but notes (*.md).'''
    from libcuflynx.calibration_workflow.runner import sha256
    root = os.path.dirname(system.workflow_path)
    out = {}
    for d, dirs, files in os.walk(root):
        dirs.sort()
        for f in sorted(files):
            if not f.endswith('.md') and not f.startswith('.'):
                path = os.path.join(d, f)
                out[os.path.relpath(path, root).replace(os.sep, '/')] = sha256(path)
    return out


def workflow_stale(system):
    '''Why the model's default set is not its workflow's result for the inputs in uncalibrated/ now
    ([] when it is; None when the model has no workflow): no record of a calibration, or inputs
    changed, added or removed since (e.g. a module recalibration wrote into uncalibrated/).'''
    if not system.workflow_path:
        return None
    record_path = os.path.join(system.dir, PARAMETERISATIONS, DEFAULT, CALIBRATION_RECORD)
    if not os.path.isfile(record_path):
        return [f'default has no {CALIBRATION_RECORD}: it was not written by the workflow']
    with open(record_path) as f:
        then = json.load(f).get('inputs') or {}
    now = _workflow_inputs(system)
    return ([f'{p} changed' for p in sorted(now) if p in then and then[p] != now[p]]
            + [f'{p} added' for p in sorted(set(now) - set(then))]
            + [f'{p} removed' for p in sorted(set(then) - set(now))])


def write_default(system, run_dir):
    '''Writes the workflow run in run_dir (complete and current) as the default set, with its record;
    returns the default parameters path.'''
    import datetime
    from libcuflynx.calibration_workflow.runner import write_calibrated_parameters_to
    if not system.workflow_path:
        raise ValueError(f'{system.id} has no calibration workflow ({UNCALIBRATED}/{system.name}_calibration_workflow.json)')
    path = parameters_path(system.dir)
    write_calibrated_parameters_to(system.workflow_path, run_dir, path, module_library_dirs=library_dirs())
    record = {'date': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d'),
              'inputs': _workflow_inputs(system)}
    with open(os.path.join(os.path.dirname(path), CALIBRATION_RECORD), 'w') as f:
        json.dump(record, f, indent=1)
        f.write('\n')
    return path


def system_models():
    out = []
    sdir, legacy = system_dir(), legacy_system_dir()
    specs = sorted(glob.glob(os.path.join(sdir, '*', '*', '*_system.yaml')))
    moved = {os.path.relpath(os.path.dirname(p), sdir) for p in specs}
    specs += sorted(p for p in glob.glob(os.path.join(legacy, '*', '*', '*_system.yaml'))
                    if os.path.relpath(os.path.dirname(p), legacy) not in moved)
    for spec_path in specs:
        d = os.path.dirname(spec_path)
        with open(spec_path) as f:
            spec = yaml.safe_load(f) or {}
        out.append(System(os.path.basename(d), os.path.basename(os.path.dirname(d)), d, spec))
    return out


def load_system(system_id):
    for s in system_models():
        if s.id == system_id or s.name == system_id:
            return s
    raise KeyError(system_id)


def generate(system, work_dir, reference=False, quiet=True):
    '''Generates the model, or with reference=True circulatory_autogen's original from reference/ with
    libcuflynx's built-in modules (write_reference); returns the .cellml path.'''
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture
    from cam_testing.harness import GenerationFailed

    resources = os.path.join(system.dir, 'reference') if reference else system.dir
    out_dir = os.path.join(work_dir, 'reference' if reference else 'library')
    params = f'{system.name}_parameters.json' if reference else os.path.relpath(system.parameters_path, system.dir)
    config = {
        'file_prefix': system.name,
        'input_param_file': params,
        'model_type': 'cellml', 'solver': 'CVODE_myokit',
        'resources_dir': resources, 'generated_models_dir': out_dir, 'DEBUG': False,
    }
    if not reference:
        config.update({'module_library_dirs': library_dirs(), 'use_builtin_modules': False})
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(log) if quiet else contextlib.nullcontext():
            ok = generate_with_new_architecture(False, config)
    except SystemExit as e:
        raise GenerationFailed(f'libcuflynx exited during generation ({e.code}):\n{log.getvalue()[-4000:]}')
    except Exception as e:
        raise GenerationFailed(f'{type(e).__name__}: {e}\n{log.getvalue()[-4000:]}') from e
    if not ok:
        raise GenerationFailed(f'generation returned False:\n{log.getvalue()[-4000:]}')
    return os.path.join(out_dir, system.name, f'{system.name}.cellml')


# components of a ready CellML reference that are not vessels (PhLynx's export)
REFERENCE_CELLML_IGNORED_COMPONENTS = ('environment', 'global_parameters', 'instance_parameters')


# How a spec's expected_failures.system_equivalence_test starts when circulatory_autogen cannot generate
# the original, so it has no reference/<model>.cellml: '<prefix> (<circulatory_autogen's error>)'
CANNOT_GENERATE = 'circulatory_autogen cannot generate the original'


def original_cannot_generate(system):
    '''Whether the spec says circulatory_autogen cannot generate the original (an expected failure of
    system_equivalence_test whose reason says so), so reference/ has no CellML of it.'''
    reason = (system.spec.get('expected_failures') or {}).get('system_equivalence_test') or ''
    return CANNOT_GENERATE in reason


def has_original(system):
    '''Whether reference/ holds circulatory_autogen's original (a module array).'''
    ref = os.path.join(system.dir, 'reference')
    return any(glob.glob(os.path.join(ref, f'*_module_array.{ext}')) for ext in ('json', 'csv'))


def original_cellml_path(system):
    '''Where circulatory_autogen's original is committed as generated CellML.'''
    return os.path.join(system.dir, 'reference', f'{system.name}.cellml')


def reference_cellml(system):
    '''The CellML model the system model must reproduce, or None (it has no reference):
    equivalence.reference_cellml, else reference/<model>.cellml (circulatory_autogen's original as
    generated, when reference/ has its module array; else a ready CellML model).'''
    eq = system.spec.get('equivalence') or {}
    if eq.get('reference_cellml'):
        path = os.path.join(system.dir, eq['reference_cellml'])
        if not os.path.isfile(path):
            raise FileNotFoundError(f'equivalence.reference_cellml: no file {path}')
        return path
    candidate = original_cellml_path(system)
    if os.path.isfile(candidate):
        return candidate
    if has_original(system):
        raise FileNotFoundError(
            f'{system.id}: reference/ has circulatory_autogen\'s module array but not the model generated from it, '
            f'{candidate}: write it with python -m cam_testing.system write-reference {system.id} '
            f'(needs a libcuflynx that still has its built-in modules)')
    return None


def write_reference(system, work_dir=None):
    '''Generates circulatory_autogen's original (reference/'s module array and parameters) with
    libcuflynx's built-in modules and writes it, flattened into one CellML file, to
    reference/<model>.cellml; returns that path.'''
    import shutil
    if not has_original(system):
        raise ValueError(f'{system.id}: reference/ has no module array to generate')
    work_dir = work_dir or tempfile.mkdtemp(prefix=f'cam_reference_{system.name}_')
    path = generate(system, work_dir, reference=True)
    flat = path[:-len('.cellml')] + '_flat.cellml'
    dst = original_cellml_path(system)
    shutil.copy(flat, dst)
    return dst


def all_component_names(helper, ignore_components=REFERENCE_CELLML_IGNORED_COMPONENTS):
    '''"<component>/<variable>" for every variable of a ready CellML model, but of ignore_components.'''
    out = []
    for n in helper.get_all_variable_names():
        comp, _, var = n.partition('.')
        if var and comp not in ignore_components:
            out.append(f'{comp}/{var}')
    return out


def time_variable_names(helper):
    '''The names ('vessel/var' and 'component/var') of the model's free variable, which simulate()
    returns as t. Myokit binds one variable to time; a helper without a Myokit model has none here.'''
    model = getattr(helper, 'model', None)
    time = model.time() if model is not None and hasattr(model, 'time') else None
    if time is None:
        return set()
    qname = time.qname()
    return {qname.replace('_module.', '/', 1), qname.replace('.', '/', 1)}


def simulate(model_path, spec, solver_info, names=None, ignore_components=None):
    '''Runs a model; returns (t, {qualified name 'vessel/var': array}). names=None: a generated
    model's vessel outputs; with ignore_components (a ready CellML reference): every component's.'''
    from libcuflynx.solver_wrappers import get_simulation_helper
    with contextlib.redirect_stdout(io.StringIO()):
        helper = get_simulation_helper(model_path=model_path, solver='CVODE_myokit', model_type='cellml',
                                       dt=spec.get('dt', 0.01), sim_time=spec.get('sim_time', 2.0),
                                       solver_info=solver_info, pre_time=spec.get('pre_time', 0.0))
        helper.run()
    if names is None and ignore_components is not None:
        names = all_component_names(helper, ignore_components)
    if names is None:
        # vessel outputs, plus the multi_port "sum" totals (components multiport_sum_<module>_<variable>)
        names = [n.replace('_module.', '/', 1) if '_module.' in n else n.replace('.', '/', 1)
                 for n in helper.get_all_variable_names() if '_module.' in n or n.startswith('multiport_sum_')]
    # time is returned as t, not as an output: the solver merges every module's t into one variable
    # and names it after whichever component it picks, which differs between runs
    time_names = time_variable_names(helper)
    names = [n for n in names if n not in time_names]
    t = np.asarray(helper.get_time(), dtype=float)
    out = {}
    for name in names:
        try:
            values = np.asarray(helper.get_results([[name]], flatten=True)[0], dtype=float).ravel()
        except Exception:
            continue
        out[name] = np.full(t.shape, values[0]) if values.size == 1 and t.size > 1 else values
    return t, out


def compare(ref, new, output_map, tol, ignore=None, wrapped=None, reciprocal=None):
    '''
    Normalised max differences for every reference variable that exists in the new model.
    ``ignore``: {name: reason} not compared. ``wrapped``: {name: period} for sawtooth outputs
    (e.g. the cardiac phase mt), compared modulo the period so a wrap a few ULPs earlier in
    one model isn't a difference of a whole period. ``reciprocal``: names compared as 1/x,
    for outputs like a valve's L = rho l/(A_eff + eps) that spike near closure and amplify
    tiny differences in A_eff by up to 1/eps.
    '''
    from cam_testing import harness
    # names as the models name them: a nested vessel may be written by its path (SN/soma/membrane/V)
    output_map = harness.model_output_map(output_map or {})
    ignore, wrapped = harness.model_output_keys(ignore or {}), harness.model_output_keys(wrapped or {})
    reciprocal = set(harness.model_output_keys(reciprocal or []))
    rows, missing = [], []
    for name, r in ref.items():
        if name in ignore:
            continue
        target = output_map.get(name, name)
        if target not in new:
            missing.append(name)
            continue
        n = new[target]
        if name in reciprocal:
            r, n = 1.0 / r, 1.0 / n
        scale = max(float(np.max(np.abs(r))), float(np.ptp(r)), 1e-300)
        d = r - n if r.shape == n.shape else np.array([np.inf])
        if name in wrapped:
            period = float(wrapped[name])
            d = np.abs(d) % period
            d = np.minimum(d, period - d)
        diff = float(np.max(np.abs(d)) / scale)
        rows.append({'reference': name, 'model': target, 'difference': diff, 'ok': diff <= tol})
    return rows, missing


def equivalence(system, work_dir=None):
    work_dir = work_dir or tempfile.mkdtemp(prefix=f'cam_system_{system.name}_')
    eq = system.spec.get('equivalence') or {}
    solver_info = eq.get('solver_info') or {'rtol': 1e-10, 'atol': 1e-12}
    ref_path = reference_cellml(system)
    if ref_path is None:
        raise ValueError(f'{system.id} has no reference')
    if not eq.get('reference_cellml') and has_original(system):
        # circulatory_autogen's original, as generated: outputs named as the library model's
        t_ref, ref = simulate(ref_path, system.spec, solver_info)
    else:
        ignored = eq.get('ignore_components', REFERENCE_CELLML_IGNORED_COMPONENTS)
        t_ref, ref = simulate(ref_path, system.spec, solver_info, ignore_components=tuple(ignored))
    new_path = generate(system, work_dir, reference=False)
    from cam_testing import harness
    output_map = harness.model_output_map(eq.get('output_map') or {})
    wanted = sorted({output_map.get(n, n) for n in ref})
    t_new, new = simulate(new_path, system.spec, solver_info, names=wanted)
    rows, missing = compare(ref, new, output_map, float(eq.get('tol', 1e-6)), eq.get('ignore'), eq.get('wrapped'),
                           eq.get('reciprocal'))
    return t_ref, ref, new, rows, missing


def save_result(system, test, status, message, metrics=None, details=None):
    import datetime
    os.makedirs(system.results_dir, exist_ok=True)
    result = {'test': test, 'status': status, 'message': message, 'metrics': metrics or {},
              'details': details or [],
              'timestamp': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
    with open(os.path.join(system.results_dir, f'{test}.json'), 'w') as f:
        json.dump(result, f, indent=1, default=float)
    return result


def run_model(system, work_dir):
    path = generate(system, work_dir)
    return simulate(path, system.spec, system.spec.get('solver_info') or {})


def main(argv=None):
    import argparse
    import sys
    ap = argparse.ArgumentParser(prog='python -m cam_testing.system')
    sub = ap.add_subparsers(dest='cmd', required=True)
    w = sub.add_parser('write-reference', help="write reference/<model>.cellml from circulatory_autogen's original "
                                               "(with libcuflynx's built-in modules)")
    w.add_argument('systems', nargs='*', help='category/model or model; default: every system model with an original')
    d = sub.add_parser('write-default', help="write a calibration workflow's run as the model's default parameter set")
    d.add_argument('system', help='category/model or model')
    d.add_argument('run_dir', help="the workflow's output directory (a complete, current run)")
    args = ap.parse_args(argv)
    if args.cmd == 'write-default':
        print(f'wrote {write_default(load_system(args.system), args.run_dir)}')
        return 0
    chosen = [load_system(s) for s in args.systems] or [s for s in system_models() if has_original(s)]
    failed = 0
    for s in chosen:
        try:
            path = write_reference(s)
            print(f'{s.id}: wrote {path} ({os.path.getsize(path) / 1e3:.0f} kB)')
        except Exception as e:
            failed += 1
            print(f'{s.id}: FAILED {type(e).__name__}: {str(e)[:500]}', file=sys.stderr)
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
