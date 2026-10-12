"""
System models: system_models/<category>/<model>/ (at the repo root; not modules, not supermodules).

Each system dir holds the model built from this module library (<model>_module_array.json, whose
records (instances) name a module_type, its version (module_subtype) and the parameterisation
it uses, and <model>_parameters.csv, whose values win over the parameterisations'), the reference it must reproduce under
reference/, and a spec <model>_system.yaml. The model is generated with every module library
(cam_testing.paths: this repo's modules/, then the extra ones), so its records may name versions of
either. The reference is one of:

  - circulatory_autogen's original: a vessel array (and parameters) under reference/, generated
    with libcuflynx's bundled modules;
  - a ready CellML model, simulated as it is (e.g. the flattened model PhLynx exports):
    equivalence.reference_cellml (a path relative to the system dir), or reference/<model>.cellml
    when reference/ holds no vessel array. Its outputs are "<component>/<variable>" for every
    component but those in equivalence.ignore_components (default: environment,
    global_parameters, instance_parameters), mapped through output_map as for an original.

The tests:

  system_run_test          the model generates with this library and runs; outputs finite
  system_equivalence_test  the model reproduces its reference (circulatory_autogen's original, built
                           with libcuflynx's bundled modules, or a ready CellML model): every
                           variable the reference logs,
                           mapped through equivalence.output_map (e.g. heart/u_lv -> lv/u),
                           agrees to equivalence.tol (normalised max difference) with both
                           simulated at equivalence.solver_info
  system_invariants_test   the spec's invariants (as for modules)
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
    '''Generates the model (or circulatory_autogen's original); returns the .cellml path.'''
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture
    from cam_testing.harness import GenerationFailed

    resources = os.path.join(system.dir, 'reference') if reference else system.dir
    out_dir = os.path.join(work_dir, 'reference' if reference else 'library')
    config = {
        'file_prefix': system.name,
        'input_param_file': f'{system.name}_parameters.csv',
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


def reference_cellml(system):
    '''The ready CellML model the system model must reproduce, or None (its reference is a vessel array
    to generate): equivalence.reference_cellml, else reference/<model>.cellml when reference/ has no
    vessel array.'''
    eq = system.spec.get('equivalence') or {}
    if eq.get('reference_cellml'):
        path = os.path.join(system.dir, eq['reference_cellml'])
        if not os.path.isfile(path):
            raise FileNotFoundError(f'equivalence.reference_cellml: no file {path}')
        return path
    ref = os.path.join(system.dir, 'reference')
    candidate = os.path.join(ref, f'{system.name}.cellml')
    has_array = any(glob.glob(os.path.join(ref, f'*_{n}.{ext}')) for n in ('module_array', 'vessel_array')
                    for ext in ('json', 'csv'))
    return candidate if os.path.isfile(candidate) and not has_array else None


def all_component_names(helper, ignore_components=REFERENCE_CELLML_IGNORED_COMPONENTS):
    '''"<component>/<variable>" for every variable of a ready CellML model, but of ignore_components.'''
    out = []
    for n in helper.get_all_variable_names():
        comp, _, var = n.partition('.')
        if var and comp not in ignore_components:
            out.append(f'{comp}/{var}')
    return out


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
        # vessel outputs, plus the multi_port "sum" totals (component sum_blood_volume)
        names = [n.replace('_module.', '/', 1).replace('sum_blood_volume.', 'sum_blood_volume/', 1)
                 for n in helper.get_all_variable_names() if '_module.' in n or n.startswith('sum_blood_volume.')]
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
    ignore, wrapped, reciprocal = ignore or {}, wrapped or {}, set(reciprocal or [])
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
        ref_path = generate(system, work_dir, reference=True)
        t_ref, ref = simulate(ref_path, system.spec, solver_info)
    else:
        ignored = eq.get('ignore_components', REFERENCE_CELLML_IGNORED_COMPONENTS)
        t_ref, ref = simulate(ref_path, system.spec, solver_info, ignore_components=tuple(ignored))
    new_path = generate(system, work_dir, reference=False)
    output_map = eq.get('output_map') or {}
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
