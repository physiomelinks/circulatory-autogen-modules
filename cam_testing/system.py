"""
System models: system_models/<category>/<model>/ (at the repo root; not modules, not supermodules).

Each system dir holds the model built from this module library (<model>_module_array.json, whose
records name a module_type, its version (module_subtype) and an instance, and
<model>_parameters.csv, whose values win over the instances'), circulatory_autogen's original under reference/, and a spec
<model>_system.yaml. The tests:

  system_run_test          the model generates with this library and runs; outputs finite
  system_equivalence_test  the model reproduces circulatory_autogen's original (built with
                           libcuflynx's bundled modules): every variable the original logs,
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

from cam_testing.library import MODULES_DIR, SYSTEM_MODELS_DIR

SYSTEM_DIR = SYSTEM_MODELS_DIR
# system models not moved yet (another session's modules/system/poiseuille): read where they are
LEGACY_SYSTEM_DIR = os.path.join(MODULES_DIR, 'system')


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
    paths = sorted(glob.glob(os.path.join(SYSTEM_DIR, '*', '*', '*_system.yaml')))
    moved = {os.path.relpath(os.path.dirname(p), SYSTEM_DIR) for p in paths}
    paths += sorted(p for p in glob.glob(os.path.join(LEGACY_SYSTEM_DIR, '*', '*', '*_system.yaml'))
                    if os.path.relpath(os.path.dirname(p), LEGACY_SYSTEM_DIR) not in moved)
    for spec_path in paths:
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
        config.update({'module_library_dirs': [MODULES_DIR], 'use_builtin_modules': False})
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


def simulate(model_path, spec, solver_info, names=None):
    '''Runs a generated model; returns (t, {qualified name 'vessel/var': array}).'''
    from libcuflynx.solver_wrappers import get_simulation_helper
    with contextlib.redirect_stdout(io.StringIO()):
        helper = get_simulation_helper(model_path=model_path, solver='CVODE_myokit', model_type='cellml',
                                       dt=spec.get('dt', 0.01), sim_time=spec.get('sim_time', 2.0),
                                       solver_info=solver_info, pre_time=spec.get('pre_time', 0.0))
        helper.run()
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
    ref_path = generate(system, work_dir, reference=True)
    new_path = generate(system, work_dir, reference=False)
    t_ref, ref = simulate(ref_path, system.spec, solver_info)
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
