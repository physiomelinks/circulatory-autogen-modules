"""Shared writer for module versions in this library's layout (CellML, units, config, instances, spec, references),
used by the tools/build_*_modules.py and tools/build_adaptive_cell.py builders.

Equations are infix text converted by tools/cellml_writer.py. Bib entries are registered in BIB by the builders; every
key must follow the library's convention (cam_testing.bib.convention_key).
"""
import contextlib
import csv
import io
import json
import os
import re

import yaml

from cellml_writer import component_model, units_model

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIB = {}
BUILTIN_UNITS = {'second', 'kelvin', 'watt', 'dimensionless', 'volt', 'ampere', 'coulomb', 'farad', 'siemens', 'metre',
                 'mole', 'joule', 'kilogram', 'newton', 'pascal'}

UNITS = {
    # electrophysiology, as in the library's sympathetic-neuron modules
    'milliV': [('volt', 'milli', 1, 1)],
    'nanoA': [('ampere', 'nano', 1, 1)],
    'picoF': [('farad', 'pico', 1, 1)],
    'picoS': [('siemens', 'pico', 1, 1)],
    'microS': [('siemens', 'micro', 1, 1)],
    'nanoC': [('coulomb', 'nano', 1, 1)],
    'um': [('metre', 'micro', 1, 1)],
    'um2': [('metre', 'micro', 2, 1)],
    'per_um2': [('metre', 'micro', -2, 1)],
    'nm2': [('metre', 'nano', 2, 1)],
    'm_per_s': [('metre', None, 1, 1), ('second', None, -1, 1)],
    'pascal_second': [('pascal', None, 1, 1), ('second', None, 1, 1)],
    'newton_s_per_m': [('newton', None, 1, 1), ('second', None, 1, 1), ('metre', None, -1, 1)],
    'm_per_um': [('metre', None, 1, 1), ('metre', 'micro', -1, 1)],
    'um2_per_nm2': [('metre', 'micro', 2, 1), ('metre', 'nano', -2, 1)],
    'uF_per_cm2': [('farad', 'micro', 1, 1), ('metre', 'centi', -2, 1)],
    'picoF_cm2_per_uF_um2': [('farad', 'pico', 1, 1), ('metre', 'centi', 2, 1), ('farad', 'micro', -1, 1),
                             ('metre', 'micro', -2, 1)],
    'm3_per_um3': [('metre', None, 3, 1), ('metre', 'micro', -3, 1)],
    'milliV_per_kiloV': [('volt', 'milli', 1, 1), ('volt', 'kilo', -1, 1)],
    'milliV_per_V': [('milliV', None, 1, 1), ('volt', None, -1, 1)],
    'V_per_milliV': [('volt', None, 1, 1), ('milliV', None, -1, 1)],
    'A_per_nanoA': [('ampere', 'nano', -1, 1), ('ampere', None, 1, 1)],
    'nanoA_per_A': [('ampere', 'nano', 1, 1), ('ampere', None, -1, 1)],
    'microS_per_picoS': [('siemens', 'micro', 1, 1), ('siemens', 'pico', -1, 1)],
    'C_per_mol': [('coulomb', None, 1, 1), ('mole', None, -1, 1)],
    'J_per_mol_K': [('joule', None, 1, 1), ('mole', None, -1, 1), ('kelvin', None, -1, 1)],
    'per_s': [('second', None, -1, 1)],
    'gray_per_s': [('joule', None, 1, 1), ('kilogram', None, -1, 1), ('second', None, -1, 1)],
    'per_gray': [('kilogram', None, 1, 1), ('joule', None, -1, 1)],
    'per_millimolar_s': [('millimolar', None, -1, 1), ('second', None, -1, 1)],
    'mol_per_s_per_um2': [('mole', None, 1, 1), ('second', None, -1, 1), ('metre', 'micro', -2, 1)],
    'per_mol': [('mole', None, -1, 1)],
    'microS_per_um2': [('siemens', 'micro', 1, 1), ('metre', 'micro', -2, 1)],
    'millimolar': [('mole', 'milli', 1, 1), ('litre', None, -1, 1)],
    'millimolar_per_s': [('millimolar', None, 1, 1), ('second', None, -1, 1)],
    'm2': [('metre', None, 2, 1)],
    'm3': [('metre', None, 3, 1)],
    'mol_per_s': [('mole', None, 1, 1), ('second', None, -1, 1)],
    'J_per_mol': [('joule', None, 1, 1), ('mole', None, -1, 1)],
    'kg_per_m3': [('kilogram', None, 1, 1), ('metre', None, -3, 1)],
    'J_per_kg_K': [('joule', None, 1, 1), ('kilogram', None, -1, 1), ('kelvin', None, -1, 1)],
    'W_per_m_K': [('watt', None, 1, 1), ('metre', None, -1, 1), ('kelvin', None, -1, 1)],
    'W_per_K': [('watt', None, 1, 1), ('kelvin', None, -1, 1)],
}



def register_bib(entries):
    from cam_testing.bib import convention_key
    for key, fields in entries.items():
        expected = convention_key(fields)
        assert expected == key, f'{key}: convention key is {expected}'
        BIB[key] = fields


def kept_record_keys(config_path):
    """The "licence" and "creator" of the config at config_path (given in review), so a rebuild keeps them; the
    defaults (CC0-1.0, no creator) when the file is missing."""
    keys = {'licence': 'CC0-1.0', 'creator': []}
    if os.path.isfile(config_path):
        entry = json.load(open(config_path))[0]
        keys.update({k: entry[k] for k in keys if k in entry})
    return keys


def parameters_json(csv_path):
    """Converts a parameters CSV just written at csv_path to the library's nested JSON beside it (libcuflynx's
    cuflynx-upgrade parameters, which reads the version's config for its globals), removing the CSV. The
    version's config must be written first. Returns the JSON's path."""
    from libcuflynx.upgrade import upgrade_parameters
    json_path = os.path.splitext(csv_path)[0] + '.json'
    if os.path.exists(json_path):
        os.remove(json_path)
    with contextlib.redirect_stdout(io.StringIO()):
        upgrade_parameters(csv_path)
    return json_path


def supermodule_layout(config_path):
    """Converts a supermodule version just written with "submodules" in its config (and its default parameterisation
    written) to libcuflynx's layout (cuflynx-upgrade module-config): its instances in <stem>_module_array.json with
    submodule_versions, and their values nested by instance in its default parameterisation."""
    from libcuflynx.upgrade import upgrade_module_config
    array_path = config_path.replace('_modules_config.json', '_module_array.json')
    if os.path.exists(array_path):
        os.remove(array_path)
    with contextlib.redirect_stdout(io.StringIO()):
        upgrade_module_config(config_path, module_library_dirs=(os.path.join(REPO, 'modules'),))


def parameter_rows(json_path):
    """A parameters JSON's rows as dicts (variable_name, units, value, data_reference, sourced), a global by its
    plain name, values as text."""
    from cam_testing.library import read_parameters
    return [{'variable_name': p.variable_name, 'units': p.units, 'value': p.value,
             'data_reference': p.data_reference, 'sourced': p.sourced} for p in read_parameters(json_path)]


def bib_text(keys):
    out = []
    for k in keys:
        f = BIB[k]
        fields = [(n, v) for n, v in f.items() if n != 'type']
        width = max(len(n) for n, _ in fields)
        body = ',\n'.join(f'  {n.ljust(width)} = {{{v}}}' for n, v in fields)
        out.append(f'@{f["type"]}{{{k},\n{body}\n}}')
    return '\n\n'.join(out) + '\n'


# ------------------------------------------------------------------ writing helpers
def write_version(category, module_type, component, variables, equations, config_extra, params, verification, tests,
                  bib_keys, version, extra_instances=None):
    d = os.path.join(REPO, 'modules', category, module_type, 'versions', version)
    os.makedirs(os.path.join(d, 'parameterisations', 'default'), exist_ok=True)
    stem = f'{module_type}_{version}'
    with open(os.path.join(d, f'{stem}_modules.cellml'), 'w') as f:
        f.write(component_model(stem, component, variables, equations))
    literal_units = {m for _, rhs in equations for m in re.findall(r"\bu\([^,]+,\s*(\w+)\)", rhs)}
    used_units = sorted(({v[1] for v in variables} | literal_units) - BUILTIN_UNITS)
    with open(os.path.join(d, f'{stem}_units.cellml'), 'w') as f:
        f.write(units_model({u: UNITS[u] for u in UNITS if u in used_units or any(
            p[0] == u for x in used_units for p in UNITS.get(x, []))}))
    roles = {name: role for name, role in config_extra.pop('roles').items()}
    # a multi_port "sum" is set by the modules connected through its port, or by its parameter when none is, so it is
    # a boundary_condition (libcuflynx refuses a "variable" sum with nothing connected)
    for key in ('entrance_ports', 'exit_ports', 'general_ports'):
        for port in config_extra.get(key) or []:
            multi = port.get('multi_port')
            kinds = multi if isinstance(multi, list) else [multi] * len(port['variables'])
            for var, kind in zip(port['variables'], kinds):
                if kind == 'sum':
                    roles[var] = 'boundary_condition'
    var_units = {v[0]: v[1] for v in variables}
    config = {
        'module_type': module_type, 'module_subtype': version, 'module_format': 'cellml',
        'component_file': f'{stem}_modules.cellml', 'component_type': component, 'default_parameterisation': 'default',
        **kept_record_keys(os.path.join(d, f'{stem}_modules_config.json')), **config_extra,
        'variables_and_units': [[n, var_units[n], 'access', roles[n]] for n in roles],
    }
    with open(os.path.join(d, f'{stem}_modules_config.json'), 'w') as f:
        json.dump([config], f, indent=2)
        f.write('\n')
    for inst, rows in [('default', params)] + list((extra_instances or {}).items()):
        os.makedirs(os.path.join(d, 'parameterisations', inst), exist_ok=True)
        csv_path = os.path.join(d, 'parameterisations', inst, f'{inst}_parameters.csv')
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f, lineterminator='\n')
            w.writerow(['variable_name', 'units', 'value', 'data_reference', 'sourced'])
            for name, value, ref, sourced in rows:
                w.writerow([name, var_units[name], repr(value) if isinstance(value, float) else value, ref,
                            'yes' if sourced else 'no'])
        parameters_json(csv_path)
    with open(os.path.join(d, f'{stem}_verification_config.json'), 'w') as f:
        json.dump({'module_type': module_type, 'version': version, **verification}, f, indent=2)
        f.write('\n')
    tests = dict(tests)
    described = tests.pop('parameterisation_notes', None)
    if described:
        # each parameterisation in a sentence of the notes
        tests['notes'] = (tests.get('notes', '') + ' Parameterisations: '
                          + '; '.join(f'{k}: {v}' for k, v in described.items()) + '.').strip()
    with open(os.path.join(d, f'{stem}_tests.yaml'), 'w') as f:
        yaml.safe_dump({'module_type': module_type, 'version': version, **tests}, f, sort_keys=False, width=120)
    with open(os.path.join(d, f'{stem}_references.bib'), 'w') as f:
        f.write(bib_text(bib_keys))
    return d


SOLVER = {'solver': 'CVODE_myokit', 'solver_info': {'rtol': 1e-10, 'atol': 1e-14, 'MaximumStep': 0.1}}
