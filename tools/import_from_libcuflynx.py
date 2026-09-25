"""
Import (or re-sync) modules from a circulatory_autogen / libcuflynx checkout into the
per-module layout of this repo:

    modules/<name>/<name>_modules.cellml
    modules/<name>/<name>_modules_config.json
    modules/<name>/<name>_units.cellml        units used by the module, from CA's units files
    modules/<name>/<name>_parameters.csv      seeded from CA resources/*_parameters.csv
    modules/<name>/<name>_tests.yaml          per-component test spec (created once, never overwritten)

Usage:
    python tools/import_from_libcuflynx.py --ca-dir ../circulatory_autogen [--only heart BG]

The CellML, config and units files are overwritten from CA. The parameters file is only
re-seeded with --reseed-parameters, and the tests spec is never overwritten, because both
are edited by hand once a module has been reviewed.
"""
import argparse
import csv
import glob
import json
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES_DIR = os.path.join(REPO_ROOT, 'modules')

CELLML_NS = 'http://www.cellml.org/cellml/1.1#'

# Built-in (package data) modules, then module_config_user modules: name -> (cellml, config).
BUILTIN_SKIP = {'test'}
USER_SKIP = {'tutorial', 'Simple_ODE_Benchmark',
             # a tutorial copy of two components of the built-in vessel_properties module
             'vessel_properties'}

# Units every CellML 1.1 processor knows; they are never written to a module's units file.
STANDARD_UNITS = {
    'ampere', 'becquerel', 'candela', 'celsius', 'coulomb', 'dimensionless', 'farad', 'gram',
    'gray', 'henry', 'hertz', 'joule', 'katal', 'kelvin', 'kilogram', 'liter', 'litre', 'lumen',
    'lux', 'meter', 'metre', 'mole', 'newton', 'ohm', 'pascal', 'radian', 'second', 'siemens',
    'sievert', 'steradian', 'tesla', 'volt', 'watt', 'weber',
}

UNITS_BLOCK_RE = re.compile(r'<units\b[^>]*?\bname="([^"]+)"[^>]*?(?:/>|>.*?</units>)', re.S)


def find_modules(ca_dir):
    builtin_dir = os.path.join(ca_dir, 'src', 'libcuflynx', 'generators', 'resources')
    user_dir = os.path.join(ca_dir, 'module_config_user')
    found = {}
    for src_dir, skip in ((builtin_dir, BUILTIN_SKIP), (user_dir, USER_SKIP)):
        for cellml in sorted(glob.glob(os.path.join(src_dir, '*_modules.cellml'))):
            name = os.path.basename(cellml)[:-len('_modules.cellml')]
            if name in skip:
                continue
            configs = [os.path.join(src_dir, f'{name}{suffix}') for suffix in
                       ('_modules_config.json', '_module_config.json', '_config.json')]
            configs = [c for c in configs if os.path.isfile(c)]
            if not configs:
                print(f'WARNING: no config for {cellml}, skipped')
                continue
            if name in found:
                raise RuntimeError(f'module {name} found in more than one CA directory')
            found[name] = (cellml, configs[0])
    return found


def read_units_blocks(paths):
    blocks = {}
    for path in paths:
        with open(path) as f:
            text = f.read()
        for match in UNITS_BLOCK_RE.finditer(text):
            blocks.setdefault(match.group(1), match.group(0))
    return blocks


def units_used(cellml_path, config):
    with open(cellml_path) as f:
        text = f.read()
    used = set(re.findall(r'\bunits="([^"]+)"', text))
    for entry in config:
        used.update(v[1] for v in entry['variables_and_units'])
    return used


def units_closure(names, blocks):
    '''names plus every unit their definitions depend on, minus the standard units.'''
    todo, needed, missing = list(names), set(), set()
    while todo:
        name = todo.pop()
        if name in needed or name in STANDARD_UNITS:
            continue
        if name not in blocks:
            missing.add(name)
            continue
        needed.add(name)
        todo.extend(re.findall(r'\bunits="([^"]+)"', blocks[name]))
    return needed, missing


def write_units_file(path, needed, blocks):
    # Keep CA's order so re-syncs give small diffs.
    with open(path, 'w') as f:
        f.write("<?xml version='1.0' encoding='UTF-8'?>\n")
        f.write(f'<model name="Units" xmlns="{CELLML_NS}" xmlns:cellml="{CELLML_NS}">\n')
        for name, block in blocks.items():
            if name in needed:
                f.write(f'    {block}\n')
        f.write('</model>\n')


def load_ca_parameter_sources(ca_dir):
    '''
    Every (vessel_array, parameters) pair in CA's resources/, as a list of
    (prefix, vessels: [(name, vessel_type, BC_type)], params: {variable_name: row}).
    '''
    resources = os.path.join(ca_dir, 'resources')
    sources = []
    for vessel_path in sorted(glob.glob(os.path.join(resources, '*_vessel_array.csv'))):
        prefix = os.path.basename(vessel_path)[:-len('_vessel_array.csv')]
        param_path = os.path.join(resources, f'{prefix}_parameters.csv')
        if not os.path.isfile(param_path):
            continue
        try:
            with open(vessel_path) as f:
                reader = csv.DictReader(f, skipinitialspace=True)
                vessels = [(r['name'].strip(), r['vessel_type'].strip(), r['BC_type'].strip())
                           for r in reader if r.get('name')]
            with open(param_path) as f:
                reader = csv.DictReader(f, skipinitialspace=True)
                params = {r['variable_name'].strip(): r for r in reader if r.get('variable_name')}
        except (KeyError, csv.Error, UnicodeDecodeError) as e:
            print(f'WARNING: could not read {prefix}: {e}')
            continue
        sources.append((prefix, vessels, params))
    return sources


def seed_parameters(config, sources):
    '''
    One row per constant / boundary_condition of each (vessel_type, BC_type), and per
    global_constant, with the first value found in a CA example model that uses it.
    '''
    rows = []
    globals_seen = {}
    for entry in config:
        vessel_type, bc_type = entry['vessel_type'], entry['BC_type']
        instances = [(prefix, vessel, params) for prefix, vessels, params in sources
                     for vessel, vt, bc in vessels if vt == vessel_type and bc == bc_type]
        for var, units, _access, kind in entry['variables_and_units']:
            kind = kind.strip()
            if kind == 'global_constant':
                if var not in globals_seen:
                    globals_seen[var] = (units, [(p, params) for p, _v, params in instances])
                else:
                    globals_seen[var][1].extend((p, params) for p, _v, params in instances)
                continue
            if kind not in ('constant', 'boundary_condition'):
                continue
            value, ref = 'TODO', 'TODO'
            for prefix, vessel, params in instances:
                row = params.get(f'{var}_{vessel}')
                if row is not None:
                    value = row['value'].strip()
                    ref = f"{prefix}: {(row.get('data_reference') or '').strip()}"
                    break
            rows.append({'vessel_type': vessel_type, 'BC_type': bc_type, 'variable_name': var,
                         'units': units, 'value': value, 'kind': kind, 'data_reference': ref})
    for var, (units, instances) in globals_seen.items():
        value, ref = 'TODO', 'TODO'
        for prefix, params in instances:
            row = params.get(var)
            if row is not None:
                value = row['value'].strip()
                ref = f"{prefix}: {(row.get('data_reference') or '').strip()}"
                break
        rows.append({'vessel_type': 'global', 'BC_type': '', 'variable_name': var, 'units': units,
                     'value': value, 'kind': 'global_constant', 'data_reference': ref})
    return rows


# Defects found in the CA modules at import (see tests/test_structure.py); reviewed per module.
KNOWN_ISSUES = {
    'BVC_Kidney': ["variable q_vc_W has kind 'variable '"],
    'cell': ['undefined units: C_per_M', 'undefined units: J_per_MK', 'undefined units: M_per_m3',
             'undefined units: m3_per_M_millis'],
    'cvs_mvp': ['module_type heart_simplesimple_type_OLD is not a component in the library'],
    'input_stimulation': ['undefined units: V'],
    'ion_channel': ['undefined units: C_per_M', 'undefined units: J_per_MK'],
    'parasympathetic': ['undefined units: milliL'],
}

PARAMETER_COLUMNS = ['vessel_type', 'BC_type', 'variable_name', 'units', 'value', 'kind', 'data_reference',
                     'verified_min', 'verified_max', 'validated_min', 'validated_max']


def write_parameters(path, rows):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=PARAMETER_COLUMNS, restval='')
        writer.writeheader()
        writer.writerows(rows)


def default_tests_spec(name, config, known_issues=None):
    return {
        'module': name,
        # Set to true once the module has been walked through: until then its V&V tests
        # are collected but skipped, so CI stays meaningful while modules are reviewed.
        'reviewed': False,
        'known_issues': known_issues or [],
        'defaults': {
            'sim_time': 2.0,
            'pre_time': 0.0,
            'dt': 0.01,
            'solver': 'CVODE_myokit',
            'bc_sweep': {'factors': [0.5, 0.75, 1.0, 1.5, 2.0], 'ranges': {}},
            'timestep': {'scheme': 'rk4', 'dts': [4e-3, 2e-3, 1e-3, 5e-4], 'min_order': 0.8, 'tol': 1e-4,
                         'cvode_tol': 1e-3},
            'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]},
        },
        'components': [
            {
                'vessel_type': e['vessel_type'],
                'BC_type': e['BC_type'],
                'module_type': e['module_type'],
                'outputs': [v[0] for v in e['variables_and_units'] if v[3].strip() == 'variable'][:6],
                'invariants': [],
                **({'skip': f"module_format {e['module_format']}: not a CellML component, cannot be tested in isolation"}
                   if e.get('module_format', 'cellml') != 'cellml' else {}),
                'validation': {
                    'baseline': {'status': 'pending', 'reason': 'validation data not chosen yet'},
                    'calibrate': {'status': 'pending', 'reason': 'validation data not chosen yet'},
                },
            }
            for e in config
        ],
    }


def import_module(name, cellml_src, config_src, unit_blocks, sources, reseed):
    dest = os.path.join(MODULES_DIR, name)
    os.makedirs(dest, exist_ok=True)
    cellml_dest = os.path.join(dest, f'{name}_modules.cellml')
    shutil.copyfile(cellml_src, cellml_dest)

    with open(config_src) as f:
        config_text = f.read()
    config = json.loads(config_text)
    # The module_file field must name the file as it is in this repo.
    config_text = re.sub(r'"module_file"\s*:\s*"[^"]*"', f'"module_file": "{name}_modules.cellml"', config_text)
    with open(os.path.join(dest, f'{name}_modules_config.json'), 'w') as f:
        f.write(config_text)

    needed, missing = units_closure(units_used(cellml_src, config), unit_blocks)
    if missing:
        print(f'WARNING: {name}: units not defined in CA units files: {sorted(missing)}')
    write_units_file(os.path.join(dest, f'{name}_units.cellml'), needed, unit_blocks)

    params_path = os.path.join(dest, f'{name}_parameters.csv')
    if reseed or not os.path.isfile(params_path):
        rows = seed_parameters(config, sources)
        write_parameters(params_path, rows)
        todo = sum(r['value'] == 'TODO' for r in rows)
        print(f'{name}: {len(config)} components, {len(rows)} parameters ({todo} TODO)')

    tests_path = os.path.join(dest, f'{name}_tests.yaml')
    if not os.path.isfile(tests_path):
        with open(tests_path, 'w') as f:
            yaml.safe_dump(default_tests_spec(name, config, KNOWN_ISSUES.get(name)), f, sort_keys=False, width=100, default_flow_style=None)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ca-dir', default=os.path.join(REPO_ROOT, '..', 'circulatory_autogen'))
    parser.add_argument('--only', nargs='*', help='only import these module names')
    parser.add_argument('--reseed-parameters', action='store_true',
                        help='overwrite existing <name>_parameters.csv files')
    args = parser.parse_args(argv)

    ca_dir = os.path.abspath(args.ca_dir)
    modules = find_modules(ca_dir)
    if args.only:
        unknown = set(args.only) - set(modules)
        if unknown:
            sys.exit(f'unknown modules: {sorted(unknown)}')
        modules = {k: v for k, v in modules.items() if k in args.only}

    unit_blocks = read_units_blocks([
        os.path.join(ca_dir, 'src', 'libcuflynx', 'generators', 'resources', 'units.cellml'),
        os.path.join(ca_dir, 'module_config_user', 'user_units.cellml'),
    ])
    sources = load_ca_parameter_sources(ca_dir)
    for name, (cellml_src, config_src) in sorted(modules.items()):
        import_module(name, cellml_src, config_src, unit_blocks, sources, args.reseed_parameters)


if __name__ == '__main__':
    main()
