"""
Builds a single module component into a runnable model with libcuflynx, and runs it.

A component is instantiated alone, as one vessel with no connections, so every
boundary_condition variable becomes a parameter the tests can set. Models are generated
from this repo's modules only (``use_builtin_modules: false``), never libcuflynx's
bundled copies.
"""
import contextlib
import csv
import io
import os

import numpy as np

from cam_testing.library import MODULES_DIR

# The single vessel's name. Not 'heart': libcuflynx special-cases a vessel called heart.
VESSEL = 'mod'


class MissingParameters(Exception):
    pass


class GenerationFailed(Exception):
    pass


def parameter_name(param):
    '''The name libcuflynx gives a parameter in the generated model.'''
    if param.kind == 'global_constant' or param.vessel_type == 'global':
        return param.variable_name
    return f'{param.variable_name}_{VESSEL}'


def output_name(variable):
    return f'{VESSEL}/{variable}'


def _write_resources(component, resources_dir, prefix, overrides):
    os.makedirs(resources_dir, exist_ok=True)
    with open(os.path.join(resources_dir, f'{prefix}_vessel_array.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'BC_type', 'vessel_type', 'inp_vessels', 'out_vessels'])
        writer.writerow([VESSEL, component.BC_type, component.vessel_type, '', ''])

    params = component.parameters()
    todo = [p.variable_name for p in params if p.is_todo and parameter_name(p) not in overrides]
    if todo:
        raise MissingParameters(f'parameters still TODO in {component.module.name}_parameters.csv: {", ".join(todo)}')
    with open(os.path.join(resources_dir, f'{prefix}_parameters.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variable_name', 'units', 'value', 'data_reference'])
        for p in params:
            name = parameter_name(p)
            value = overrides.get(name, p.value)
            writer.writerow([name, p.units, value, p.data_reference or 'cam_testing'])


def generate(component, work_dir, overrides=None, quiet=True, model_type='cellml'):
    '''
    Generates the component's model in ``work_dir``; returns the path of the .cellml file,
    or of the .py file for ``model_type='python'``.
    '''
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture

    overrides = overrides or {}
    prefix = f'{component.module.name}__{component.id}'
    if model_type != 'cellml':
        work_dir = os.path.join(work_dir, model_type)
    resources_dir = os.path.join(work_dir, 'resources')
    generated_dir = os.path.join(work_dir, 'generated_models')
    _write_resources(component, resources_dir, prefix, overrides)
    config = {
        'file_prefix': prefix,
        'input_param_file': f'{prefix}_parameters.csv',
        'model_type': model_type,
        'solver': component.spec.get('solver', 'CVODE_myokit') if model_type == 'cellml' else 'solve_ivp',
        'resources_dir': resources_dir,
        'generated_models_dir': generated_dir,
        'module_library_dirs': [MODULES_DIR],
        'use_builtin_modules': False,
        'DEBUG': False,
    }
    log = io.StringIO()
    redirect = contextlib.redirect_stdout(log) if quiet else contextlib.nullcontext()
    try:
        with redirect:
            ok = generate_with_new_architecture(False, config)
    except SystemExit as e:
        raise GenerationFailed(f'libcuflynx exited during generation ({e.code}):\n{log.getvalue()[-4000:]}')
    except Exception as e:
        raise GenerationFailed(f'{type(e).__name__}: {e}\n{log.getvalue()[-4000:]}') from e
    if not ok:
        raise GenerationFailed(f'generation returned False:\n{log.getvalue()[-4000:]}')
    suffix = '.cellml' if model_type == 'cellml' else '.py'
    return os.path.join(generated_dir, prefix, f'{prefix}{suffix}')


def simulation_helper(model_path, spec, solver_info=None, solver=None, model_type='cellml'):
    from libcuflynx.solver_wrappers import get_simulation_helper

    info = dict(spec.get('solver_info') or {}) if solver is None else {}
    info.update(solver_info or {})
    with contextlib.redirect_stdout(io.StringIO()):
        return get_simulation_helper(model_path=model_path, solver=solver or spec.get('solver', 'CVODE_myokit'),
                                     model_type=model_type, dt=spec['dt'], sim_time=spec['sim_time'],
                                     solver_info=info, pre_time=spec.get('pre_time', 0.0))


def run(helper, outputs, params=None):
    '''
    Runs the helper from the model's initial state (optionally with parameter overrides);
    returns (t, {var: array}). Without the reset, a helper continues from the previous run's
    final state.
    '''
    if hasattr(helper, 'reset_states'):
        helper.reset_states()
    if params:
        names = list(params)
        helper.set_param_vals(names, [params[n] for n in names])
    with contextlib.redirect_stdout(io.StringIO()):
        ok = helper.run()
    if ok is False:
        raise RuntimeError('simulation failed')
    t = np.asarray(helper.get_time(), dtype=float)
    results = {}
    for var in outputs:
        values = np.asarray(helper.get_results([[output_name(var)]], flatten=True)[0], dtype=float).ravel()
        if values.size == 1 and t.size > 1:
            # an output that depends only on constants is logged once: broadcast it over time
            values = np.full(t.shape, values[0])
        results[var] = values
    return t, results
