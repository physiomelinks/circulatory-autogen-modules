"""
Builds a single module version into a runnable model with libcuflynx, and runs it.

A version is instantiated alone, as one vessel with no connections, so every
boundary_condition variable becomes a parameter the tests can set. Its parameter values are
those of a parameterisation (by default the version's default parameterisation), written to the model's
parameters file. Models are generated from this repo's modules (and the extra libraries of
cam_testing.paths, which a harness may use) only (``use_builtin_modules: false``), never
libcuflynx's bundled copies.
"""
import contextlib
import csv
import io
import os

import numpy as np

from cam_testing import module_array
from cam_testing.library import library_dirs

# The single vessel's name. Not 'heart': libcuflynx special-cases a vessel called heart.
VESSEL = 'mod'


class MissingParameters(Exception):
    pass


class GenerationFailed(Exception):
    pass


def libcuflynx_model_name(written):
    '''The model's name for a parameter written in a parameters file (<vessel>/<var>, or a global).'''
    from libcuflynx.utilities.parameter_names import model_name
    return model_name(written)


def model_output(written):
    '''The model's name of an output written <vessel path>/<var> (SN/soma/membrane/V ->
    SN__soma__membrane/V); a <vessel>/<var> or flat name is unchanged.'''
    from libcuflynx.utilities.parameter_names import output_model_name
    return output_model_name(written)


def model_output_map(output_map):
    '''An equivalence output_map {reference output: model output} with its names as the models name them.'''
    return {model_output(k): model_output(v) for k, v in (output_map or {}).items()}


def model_output_keys(names):
    '''An ignore/wrapped dict (names to a reason or period) or a list of names, with the names as the models name them.'''
    if isinstance(names, dict):
        return {model_output(k): v for k, v in names.items()}
    return [model_output(n) for n in names or []]


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
# mod__<submodule> (nested supermodules: mod__<submodule>__<subsubmodule> ...: LEVEL_SEPARATOR), and
# names a submodule's parameter <var>_mod__<submodule>. Its own parameterisation names that parameter
# <submodule>/<var> (nested: <submodule>/<subsubmodule>/<var>), the supermodule-level name the spec
# uses too (run_parameters, bc_sweep, validation); globals keep their names in both.

LEVEL_SEPARATOR = '__'


def supermodule_model_name(version, name, vessel=VESSEL):
    '''A supermodule parameterisation's parameter name (<submodule path>/<var>, or a global) -> its
    name in a model in which the supermodule is the vessel ``vessel`` (<var>_<vessel>__<submodule path>).'''
    if name.startswith('global/'):
        return name[len('global/'):]
    path, _, var = name.rpartition('/')
    if not path:
        return name
    return f'{var}_{LEVEL_SEPARATOR.join([vessel] + path.split("/"))}'


def model_parameter_name(component, param):
    '''The name a parameterisation's parameter has in the version's generated model.'''
    if component.is_supermodule and not getattr(param, 'model_name', ''):
        return supermodule_model_name(component, param.variable_name)
    return parameter_name(param)


def written_parameter_name(component, param):
    '''The name a parameterisation's parameter is written with in the model's parameters file:
    <vessel>/<var> (a supermodule's: mod/<submodule path>/<var>), a global as it is.'''
    if component.is_supermodule:
        name = param.variable_name
        return f'{VESSEL}/{name}' if '/' in name and not name.startswith('global/') else name
    return param.variable_name if param.is_global else f'{VESSEL}/{param.variable_name}'


def _submodule_version(sub):
    from cam_testing.library import load_version
    return load_version(sub['module_type'], sub['module_subtype'])


def resolve_submodule_path(version, path):
    '''"ra" or (nested) "soma/Ca" -> the component version at the end of that submodule path, or None.'''
    head, _, rest = path.partition('/')
    sub = next((s for s in version.submodules if s['name'] == head), None)
    if sub is None:
        return None
    try:
        target = _submodule_version(sub)
    except Exception:  # noqa: BLE001 - a broken reference is the structure test's to report
        return None
    if not rest:
        return target
    return resolve_submodule_path(target, rest) if target.is_supermodule else None


def split_model_name(version, full, vessel=VESSEL):
    '''A generated model's parameter name -> (supermodule-level name, variable, submodule path or None,
    owning component version or None): "V_in_mod__membrane" -> ("membrane/V_in", "V_in", "membrane",
    <its version>); a global (no "_mod__<submodule>" suffix) -> (name, name, None, None).'''
    marker = f'_{vessel}{LEVEL_SEPARATOR}'
    i = full.rfind(marker)
    if i <= 0:
        return full, full, None, None
    var, path = full[:i], full[i + len(marker):].replace(LEVEL_SEPARATOR, '/')
    owner = resolve_submodule_path(version, path)
    if owner is None:
        return full, full, None, None
    return f'{path}/{var}', var, path, owner


def generated_parameters_path(model_path):
    '''The parameters file libcuflynx writes next to a generated model.'''
    d = os.path.dirname(model_path)
    return os.path.join(d, f'{os.path.basename(d)}_parameters.csv')


def supermodule_parameters(version, model_path, vessel=VESSEL):
    '''
    Every parameter of a supermodule's generated model (the flattened model: the supermodule's own
    parameterisation, its submodules' default parameterisations, and the boundary conditions no internal connection closes),
    as Parameters named at the supermodule level (<var>_<submodule path>, globals as they are), with
    model_name the generated model's name, and kind from the owning submodule's config
    (global_constant for the globals).
    '''
    from cam_testing.library import Parameter
    own = {p.variable_name: p for p in version.default_parameterisation.parameters()}
    out = []
    with open(generated_parameters_path(model_path), newline='') as f:
        for row in csv.DictReader(f):
            full = (row.get('variable_name') or '').strip()
            if not full:
                continue
            name, var, _path, owner = split_model_name(version, full, vessel)
            kind = 'global_constant' if owner is None else owner.kinds.get(var, '')
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
    '''{output_key: units} for "<vessel>__<submodule path>/<var>" outputs, from the owning submodule's config.'''
    units = {}
    prefix = vessel + LEVEL_SEPARATOR
    for o in outputs:
        comp, _, var = o.rpartition('/')
        if not comp.startswith(prefix):
            continue
        owner = resolve_submodule_path(version, comp[len(prefix):].replace(LEVEL_SEPARATOR, '/'))
        if owner is not None:
            u = {v[0]: v[1] for v in owner.config.get('variables_and_units') or []}.get(var)
            if u:
                units[output_key(o)] = u
    return units


def output_name(variable):
    '''"var" is a variable of the component under test; "vessel/var" one of a neighbour in its test
    network, or of a submodule by its path ("sub/var", "a/b/var").'''
    return model_output(variable) if '/' in variable else f'{VESSEL}/{variable}'


def output_key(variable):
    '''The name an output has in invariant expressions ("vessel/var" -> "vessel__var").'''
    return variable.replace('/', '__')


def network(component, parameterisation=None):
    '''The test network of ``component`` at ``parameterisation`` (a Parameterisation or its name): the
    parameterisation's own (validation.<parameterisation>.harness in the verification config, e.g. an experiment's set-up for a
    calibration) when it has one, else the version's (harness; empty: the version alone).'''
    name = getattr(parameterisation, 'name', parameterisation)
    own = (((component.spec.get('validation') or {}).get(name) or {}).get('harness')) if name else None
    return own if own is not None else (component.spec.get('harness') or {})


def _write_resources(component, resources_dir, prefix, overrides, parameters=None, parameterisation=None):
    '''
    The test network: by default the version alone. A version that only works with
    neighbours (its inputs are variables another vessel supplies) gives a small network in its
    spec, in which the version under test is the vessel named "mod":

        harness:
          module_array:    # [name, module_subtype (version), module_type, inp, out, parameterisation]
            - [pressure_in, constant, inlet_pressure, '', mod, default]
            - [mod, pv_0D_1D, coupler, pressure_in, constant_1D, default]
          parameters:              # values for the neighbours' parameters, as <vessel>/<var>
            - [pressure_in/P, J_per_m3, 2000, source]

    ``parameters``: the parameterisation's parameters to write (default: the version's default
    parameterisation). A parameterisation may have its own network instead
    (validation.<parameterisation>.harness, see ``network``); ``parameterisation`` selects it.
    '''
    os.makedirs(resources_dir, exist_ok=True)
    network = globals()['network'](component, parameterisation)
    rows = module_array.harness_rows(network) or [[VESSEL, component.BC_type, component.vessel_type, '', '', 'default']]
    if not any(r[0] == VESSEL for r in rows):
        raise ValueError(f'harness.module_array must contain the component under test as vessel "{VESSEL}"')
    records = module_array.from_rows(rows)
    module_array.write_records(os.path.join(resources_dir, f'{prefix}_module_array.json'), records)

    params = component.parameters() if parameters is None else parameters
    todo = [p.variable_name for p in params if p.is_todo and model_parameter_name(component, p) not in overrides]
    if todo:
        raise MissingParameters(f'parameters still TODO in the parameterisation parameters of {component.key}: {", ".join(todo)}')
    with open(os.path.join(resources_dir, f'{prefix}_parameters.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variable_name', 'units', 'value', 'data_reference'])
        written = set()
        for p in params:
            name = model_parameter_name(component, p)
            value = overrides.get(name, p.value)
            writer.writerow([written_parameter_name(component, p), p.units, value, p.data_reference or 'cam_testing'])
            written.add(name)
        for written_name, units, value, *ref in network.get('parameters') or []:
            name = libcuflynx_model_name(written_name)
            if name not in written:
                writer.writerow([written_name, units, overrides.get(name, value), (ref[0] if ref else 'cam_testing harness')])


def generate(component, work_dir, overrides=None, quiet=True, model_type='cellml', parameters=None, parameterisation=None):
    '''
    Generates the version's model in ``work_dir``; returns the path of the .cellml file,
    or of the .py file for ``model_type='python'``. ``parameters``: a parameterisation's
    parameters (default: the version's default parameterisation). ``parameterisation``: the
    parameterisation whose own test network (validation.<parameterisation>.harness) is used, if it has
    one (``instance``: its former name, still accepted).
    '''
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture

    overrides = overrides or {}
    prefix = component.id
    if model_type != 'cellml':
        work_dir = os.path.join(work_dir, model_type)
    resources_dir = os.path.join(work_dir, 'resources')
    generated_dir = os.path.join(work_dir, 'generated_models')
    _write_resources(component, resources_dir, prefix, overrides, parameters, parameterisation=parameterisation)
    config = {
        'file_prefix': prefix,
        'input_param_file': f'{prefix}_parameters.csv',
        'model_type': model_type,
        'solver': component.spec.get('solver', 'CVODE_myokit') if model_type == 'cellml' else 'solve_ivp',
        'resources_dir': resources_dir,
        'generated_models_dir': generated_dir,
        'module_library_dirs': library_dirs(),     # this repo's modules/, then the extra libraries
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
