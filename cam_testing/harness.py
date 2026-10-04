"""
Builds a single module version into a runnable model with libcuflynx, and runs it.

A version is instantiated alone, as one vessel with no connections, so every
boundary_condition variable becomes a parameter the tests can set. Its parameter values are
those of an instance (by default the version's default instance), written to the model's
parameters file. Models are generated from this repo's modules only
(``use_builtin_modules: false``), never libcuflynx's bundled copies.
"""
import contextlib
import csv
import io
import os

import numpy as np

from cam_testing import vessel_array
from cam_testing.library import MODULES_DIR

# The single vessel's name. Not 'heart': libcuflynx special-cases a vessel called heart.
VESSEL = 'mod'


class MissingParameters(Exception):
    pass


class GenerationFailed(Exception):
    pass


def parameter_name(param):
    '''The name libcuflynx gives a parameter in the generated model.'''
    if getattr(param, 'model_name', ''):
        return param.model_name            # known from the generated model (a supermodule's flattened parameters)
    if param.is_global:
        return param.variable_name
    return f'{param.variable_name}_{VESSEL}'


# ---- supermodule versions ----------------------------------------------------------------------
#
# A supermodule generated alone is the vessel VESSEL ("mod"); libcuflynx expands it into the vessels
# mod_<submodule> (nested supermodules: mod_<submodule>_<subsubmodule> ...), and names a submodule's
# parameter <var>_mod_<submodule>. Its own instance names that parameter <var>_<submodule> (the
# supermodule-level name used in the spec: run_parameters, bc_sweep, validation); globals keep
# their names in both.

def supermodule_model_name(version, name, vessel=VESSEL):
    '''A supermodule instance's parameter name (<var>_<submodule>, or a global) -> its name in a model in
    which the supermodule is the vessel ``vessel`` (<var>_<vessel>_<submodule>).'''
    if name in version.supermodule_globals:
        return name
    if name in version.shared_parameter_names:
        return f'{name}_{vessel}'
    for sub in sorted(version.submodule_names, key=len, reverse=True):
        if name.endswith('_' + sub) and len(name) > len(sub) + 1:
            return f'{name[:-len(sub)]}{vessel}_{sub}'
    return name


def model_parameter_name(component, param):
    '''The name an instance's parameter has in the version's generated model.'''
    if component.is_supermodule and not getattr(param, 'model_name', ''):
        return supermodule_model_name(component, param.variable_name)
    return parameter_name(param)


def _submodule_version(sub):
    from cam_testing.library import load_version
    return load_version(sub['module_type'], sub['module_subtype'])


def resolve_submodule_path(version, path):
    '''"ra" or (nested) "soma_Ca" -> the component version at the end of that submodule path, or None.'''
    for sub in sorted(version.submodules, key=lambda s: len(s['name']), reverse=True):
        name = sub['name']
        if path == name or path.startswith(name + '_'):
            try:
                target = _submodule_version(sub)
            except Exception:  # noqa: BLE001 - a broken reference is the structure test's to report
                return None
            if path == name:
                return target
            if target.is_supermodule:
                found = resolve_submodule_path(target, path[len(name) + 1:])
                if found is not None:
                    return found
    return None


def split_model_name(version, full, vessel=VESSEL):
    '''A generated model's parameter name -> (supermodule-level name, variable, submodule path or None,
    owning component version or None): "V_in_mod_membrane" -> ("V_in_membrane", "V_in", "membrane",
    <its version>); a global (no "_mod_<submodule>" suffix) -> (name, name, None, None).'''
    for shared in version.shared_parameter_names:
        if full == f'{shared}_{vessel}':
            return shared, shared, '', None
    marker = f'_{vessel}_'
    start = 0
    while True:
        i = full.find(marker, start)
        if i <= 0:
            return full, full, None, None
        path = full[i + len(marker):]
        owner = resolve_submodule_path(version, path)
        if owner is not None:
            return f'{full[:i]}_{path}', full[:i], path, owner
        start = i + 1


def generated_parameters_path(model_path):
    '''The parameters file libcuflynx writes next to a generated model.'''
    d = os.path.dirname(model_path)
    return os.path.join(d, f'{os.path.basename(d)}_parameters.csv')


def supermodule_parameters(version, model_path, vessel=VESSEL):
    '''
    Every parameter of a supermodule's generated model (the flattened model: the supermodule's own
    instance, its submodules' instances, and the boundary conditions no internal connection closes),
    as Parameters named at the supermodule level (<var>_<submodule path>, globals as they are), with
    model_name the generated model's name, and kind from the owning submodule's config
    (global_constant for the globals).
    '''
    from cam_testing.library import Parameter
    own = {p.variable_name: p for p in version.default_instance.parameters()}
    out = []
    with open(generated_parameters_path(model_path), newline='') as f:
        for row in csv.DictReader(f):
            full = (row.get('variable_name') or '').strip()
            if not full:
                continue
            name, var, path, owner = split_model_name(version, full, vessel)
            kind = ('constant' if path == '' else 'global_constant') if owner is None else owner.kinds.get(var, '')
            mine = own.get(name)
            out.append(Parameter(name, (row.get('units') or '').strip(), (row.get('value') or '').strip(),
                                 (row.get('data_reference') or '').strip(), mine.sourced if mine else '',
                                 kind, model_name=full))
    return out


def state_outputs(helper):
    '''Every state of a generated model as an output name ("<vessel>/<var>"): a supermodule's default outputs.'''
    out = []
    for state in helper.model.states():
        comp, var = state.parent().name(), state.name()
        if comp.endswith('_module'):
            comp = comp[:-len('_module')]
        out.append(f'{comp}/{var}')
    return out


def supermodule_units(version, outputs, vessel=VESSEL):
    '''{output_key: units} for "<vessel>_<submodule path>/<var>" outputs, from the owning submodule's config.'''
    units = {}
    for o in outputs:
        comp, _, var = o.rpartition('/')
        if not comp.startswith(vessel + '_'):
            continue
        owner = resolve_submodule_path(version, comp[len(vessel) + 1:])
        if owner is not None:
            u = {v[0]: v[1] for v in owner.config.get('variables_and_units') or []}.get(var)
            if u:
                units[output_key(o)] = u
    return units


def output_name(variable):
    '''"var" is a variable of the component under test; "vessel/var" one of a neighbour in its test network.'''
    return variable if '/' in variable else f'{VESSEL}/{variable}'


def output_key(variable):
    '''The name an output has in invariant expressions ("vessel/var" -> "vessel__var").'''
    return variable.replace('/', '__')


def _write_resources(component, resources_dir, prefix, overrides, parameters=None, instances=True):
    '''
    The test network: by default the version alone. A version that only works with
    neighbours (its inputs are variables another vessel supplies) gives a small network in its
    spec, in which the version under test is the vessel named "mod":

        harness:
          vessel_array:    # [name, module_subtype (version), module_type, inp, out, instance]
            - [pressure_in, nn_constant, inlet_pressure, '', mod, default]
            - [mod, pv_0D_1D, coupler, pressure_in, constant_1D, default]
          parameters:              # values for the neighbours' parameters
            - [P_pressure_in, J_per_m3, 2000, source]

    ``parameters``: the instance parameters to write (default: the version's default instance).
    ``instances=False`` leaves "instance" out of the records (each record then gets its version's
    default_instance). libcuflynx before circulatory_autogen #535's baad9e73 needed this for the C++
    0D-1D split, which appended 5-column rows.
    '''
    os.makedirs(resources_dir, exist_ok=True)
    network = (component.spec.get('harness') or {})
    rows = network.get('vessel_array') or [[VESSEL, component.BC_type, component.vessel_type, '', '', 'default']]
    if not any(r[0] == VESSEL for r in rows):
        raise ValueError(f'harness.vessel_array must contain the component under test as vessel "{VESSEL}"')
    records = vessel_array.from_rows(rows)
    if not instances:
        records = [{k: v for k, v in r.items() if k != 'instance'} for r in records]
    vessel_array.write_records(os.path.join(resources_dir, f'{prefix}_vessel_array.json'), records)
    if not vessel_array.libcuflynx_reads_json():
        # a libcuflynx without JSON vessel-array support reads the CSV
        with open(os.path.join(resources_dir, f'{prefix}_vessel_array.csv'), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['name', 'BC_type', 'vessel_type', 'inp_vessels', 'out_vessels'])
            writer.writerows([str(x) for x in r[:5]] for r in rows)

    params = component.parameters() if parameters is None else parameters
    todo = [p.variable_name for p in params if p.is_todo and model_parameter_name(component, p) not in overrides]
    if todo:
        raise MissingParameters(f'parameters still TODO in the instance parameters of {component.key}: {", ".join(todo)}')
    with open(os.path.join(resources_dir, f'{prefix}_parameters.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variable_name', 'units', 'value', 'data_reference'])
        written = set()
        for p in params:
            name = model_parameter_name(component, p)
            value = overrides.get(name, p.value)
            writer.writerow([name, p.units, value, p.data_reference or 'cam_testing'])
            written.add(name)
        for name, units, value, *ref in network.get('parameters') or []:
            if name not in written:
                writer.writerow([name, units, overrides.get(name, value), (ref[0] if ref else 'cam_testing harness')])


def generate(component, work_dir, overrides=None, quiet=True, model_type='cellml', parameters=None):
    '''
    Generates the version's model in ``work_dir``; returns the path of the .cellml file,
    or of the .py file for ``model_type='python'``. ``parameters``: an instance's parameters
    (default: the version's default instance).
    '''
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture

    overrides = overrides or {}
    prefix = component.id
    if model_type != 'cellml':
        work_dir = os.path.join(work_dir, model_type)
    resources_dir = os.path.join(work_dir, 'resources')
    generated_dir = os.path.join(work_dir, 'generated_models')
    _write_resources(component, resources_dir, prefix, overrides, parameters)
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
        results[output_key(var)] = values
    return t, results
