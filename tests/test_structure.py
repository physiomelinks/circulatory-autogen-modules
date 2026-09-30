"""
Fast static checks on the module library (no libcuflynx needed).

A problem listed (as a substring) under ``known_issues`` in a module's spec is reported as
xfail rather than a failure, so known defects stay visible without blocking CI; once fixed,
remove the entry.
"""
import functools
import glob
import json
import os
import re
import xml.etree.ElementTree as ET
from collections import defaultdict

import pytest

from cam_testing.library import MODULES_DIR, REPO_ROOT, load_module, module_names
from cam_testing import bib
from cam_testing.mathml import component_names

CELLML_NS = 'http://www.cellml.org/cellml/1.1#'
STANDARD_UNITS = {
    'ampere', 'becquerel', 'candela', 'celsius', 'coulomb', 'dimensionless', 'farad', 'gram',
    'gray', 'henry', 'hertz', 'joule', 'katal', 'kelvin', 'kilogram', 'liter', 'litre', 'lumen',
    'lux', 'meter', 'metre', 'mole', 'newton', 'ohm', 'pascal', 'radian', 'second', 'siemens',
    'sievert', 'steradian', 'tesla', 'volt', 'watt', 'weber',
}
VARIABLE_KINDS = {'variable', 'constant', 'global_constant', 'boundary_condition'}


# unit names that denote the same unit (numerically identical definitions)
EQUIVALENT_UNITS = [{'Hz', 'per_s', 'per_second'}, {'mol_per_m3', 'millimolar', 'mM'},
                    {'J_per_m3', 'Pa', 'pascal'}]


def same_units(a, b):
    return a == b or any(a in group and b in group for group in EQUIVALENT_UNITS)


def _units_defined(path):
    root = ET.parse(path).getroot()
    return {u.get('name'): u for u in root.iter(f'{{{CELLML_NS}}}units')}


def _canonical(units_el):
    return tuple(sorted(tuple(sorted(c.attrib.items())) for c in units_el))


@functools.lru_cache(maxsize=1)
def _all_component_names():
    return {c for n in module_names() for c in component_names(load_module(n).cellml_path)}


def module_problems(name, library_components=None):
    """
    Returns (errors, warnings). Errors break model generation or are plainly wrong;
    warnings are worth a look but can be legitimate (e.g. a component defined in another
    module, which works because libcuflynx merges the whole library).
    """
    module = load_module(name)
    errors, warnings = [], []
    try:
        components = set(component_names(module.cellml_path))
    except ET.ParseError as e:
        return [f'{name}_modules.cellml does not parse: {e}'], []
    if library_components is None:
        library_components = _all_component_names()
    cellml_text = open(module.cellml_path).read()
    defined = set(_units_defined(module.units_path)) if os.path.isfile(module.units_path) else set()
    local_units = set(re.findall(r'<units\b[^>]*\bname="([^"]+)"', cellml_text))

    used_units = set(re.findall(r'\bunits="([^"]+)"', cellml_text))
    for entry in module.config:
        key = f"{entry['vessel_type']}/{entry['BC_type']}"
        if entry.get('module_format', 'cellml') != 'cellml':
            continue  # cpp / external modules are not CellML components
        if entry['module_file'] != f'{name}_modules.cellml':
            errors.append(f"{key}: module_file is {entry['module_file']}")
        if entry['module_type'] not in components:
            if entry['module_type'] in library_components:
                warnings.append(f"{key}: uses component {entry['module_type']} from another module")
            else:
                errors.append(f"{key}: module_type {entry['module_type']} is not a component in the library")
        names = {v[0] for v in entry['variables_and_units']}
        for v in entry['variables_and_units']:
            used_units.add(v[1])
            if v[3] not in VARIABLE_KINDS:
                errors.append(f"{key}: variable {v[0]} has kind {v[3]!r}")
        for port in entry.get('entrance_ports', []) + entry.get('exit_ports', []) + entry.get('general_ports', []):
            for var in port['variables']:
                if var not in names:
                    # libcuflynx stops ("the port variable ... is not a variable") as soon as such
                    # a port is connected, so this is an error even if no model connects it yet
                    errors.append(f"{key}: port variable {var} not in variables_and_units")
    for unit in sorted(used_units - defined - STANDARD_UNITS - local_units):
        errors.append(f'undefined units: {unit}')

    # a reviewed module's sourced parameters must cite an entry of its references.bib
    if module.reviewed:
        keys = set(bib.read(bib.bib_path(module)))
        for p in module.parameters:
            if not p.is_sourced or p.data_reference.lower().startswith('definitional'):
                continue
            key = bib.reference_key(p.data_reference)
            if key not in keys:
                errors.append(f'{p.vessel_type}/{p.variable_name}: sourced, but its reference '
                              f'"{p.data_reference[:40]}" is not a key in {module.name}_references.bib')

    # every identifier used in a component's equations must be one of its declared variables
    root = ET.parse(module.cellml_path).getroot()
    for comp in root.iter(f'{{{CELLML_NS}}}component'):
        declared = {v.get('name') for v in comp.findall(f'{{{CELLML_NS}}}variable')}
        used = {ci.text.strip() for ci in comp.iter('{http://www.w3.org/1998/Math/MathML}ci') if ci.text}
        for name in sorted(used - declared):
            warnings.append(f"component {comp.get('name')}: equations use undeclared variable {name}")

    # config variables must exist in the CellML component with the same units and a public interface
    comps = {c.get('name'): c for c in root.iter(f'{{{CELLML_NS}}}component')}
    for entry in module.config:
        comp = comps.get(entry['module_type'])
        if comp is None or entry.get('module_format', 'cellml') != 'cellml':
            continue
        cvars = {v.get('name'): v for v in comp.findall(f'{{{CELLML_NS}}}variable')}
        key = f"{entry['vessel_type']}/{entry['BC_type']}"
        declared_cfg = {v[0] for v in entry['variables_and_units']}
        for name, v in cvars.items():
            if v.get('public_interface') == 'in' and name not in declared_cfg and name not in ('t', 'time'):
                warnings.append(f'{key}: CellML input {name} of {entry["module_type"]} is not declared in the config (left unset)')
        seen_names = [v[0] for v in entry['variables_and_units']]
        for name in sorted({n for n in seen_names if seen_names.count(n) > 1}):
            warnings.append(f'{key}: variable {name} listed more than once in variables_and_units')
        for name, units, *_ in entry['variables_and_units']:
            v = cvars.get(name)
            if v is None:
                warnings.append(f'{key}: config variable {name} is not in CellML component {entry["module_type"]}')
                continue
            if not same_units(v.get('units'), units):
                warnings.append(f'{key}: {name} has units {units} in the config but {v.get("units")} in the CellML')
            if not v.get('public_interface'):
                warnings.append(f'{key}: config variable {name} has no public_interface in the CellML')

    spec_keys = {(c['vessel_type'], c['BC_type']) for c in module.spec.get('components', [])}
    for entry in module.config:
        if (entry['vessel_type'], entry['BC_type']) not in spec_keys:
            errors.append(f"{entry['vessel_type']}/{entry['BC_type']}: not listed in {name}_tests.yaml")
    return errors, warnings


def write_structure_results(name, errors, warnings, known):
    """Saved for the HTML report, next to the other test results."""
    module = load_module(name)
    os.makedirs(module.results_dir, exist_ok=True)
    with open(os.path.join(module.results_dir, 'structure.json'), 'w') as f:
        json.dump({'errors': errors, 'warnings': warnings, 'known_issues': known}, f, indent=2)


@pytest.mark.parametrize('name', module_names())
def test_module_structure(name):
    known = load_module(name).spec.get('known_issues') or []
    errors, warnings = module_problems(name)
    write_structure_results(name, errors, warnings, known)
    unknown = [p for p in errors if not any(k in p for k in known)]
    assert not unknown, '\n'.join(unknown)
    if errors:
        pytest.xfail('known issues: ' + '; '.join(errors))


def test_library_wide_uniqueness():
    '''libcuflynx merges every module: components, (vessel_type, BC_type) and units must not clash.'''
    component_owner, pair_owner, units_seen = {}, {}, {}
    problems = []
    for name in module_names():
        module = load_module(name)
        for comp in component_names(module.cellml_path):
            if comp in component_owner:
                problems.append(f'component {comp} is in both {component_owner[comp]} and {name}')
            component_owner[comp] = name
        for entry in module.config:
            pair = (entry['vessel_type'], entry['BC_type'])
            if pair in pair_owner:
                problems.append(f'{pair} is in both {pair_owner[pair]} and {name}')
            pair_owner[pair] = name
        if os.path.isfile(module.units_path):
            for unit, el in _units_defined(module.units_path).items():
                if unit in units_seen and units_seen[unit][1] != _canonical(el):
                    problems.append(f'units {unit} differ between {units_seen[unit][0]} and {name}')
                units_seen.setdefault(unit, (name, _canonical(el)))
    assert not problems, '\n'.join(problems)


@pytest.mark.parametrize('manifest', sorted(glob.glob(os.path.join(REPO_ROOT, 'manifests', '*.json'))),
                         ids=os.path.basename)
def test_manifest_paths_exist(manifest):
    with open(manifest) as f:
        data = json.load(f)
    missing = [e['path'] for entries in data['collections'].values() for e in entries
               if not os.path.isfile(os.path.join(REPO_ROOT, e['path']))]
    assert not missing, f'missing: {missing}'


# not modules: system models (tests/test_systems.py) and supermodules (tests/test_supermodules.py)
NON_MODULE_DIRS = {'template', 'system', 'supermodules'}


def test_every_module_dir_has_spec():
    # a module directory is any directory with a <name>_modules.cellml (modules may sit in category dirs)
    cellml = [p for p in glob.glob(os.path.join(MODULES_DIR, '**', '*_modules.cellml'), recursive=True)
              if os.path.relpath(p, MODULES_DIR).split(os.sep)[0] not in NON_MODULE_DIRS]
    missing = []
    for p in cellml:
        d = os.path.dirname(p)
        name = os.path.basename(d)
        if os.path.basename(p) != f'{name}_modules.cellml' or not os.path.isfile(os.path.join(d, f'{name}_tests.yaml')):
            missing.append(os.path.relpath(p, MODULES_DIR))
    assert not missing, f'module files not in a <name>/ directory with <name>_tests.yaml: {missing}'


# ----------------------------------------------------------------------------------------------
# JSON files against libcuflynx's schemas (src/libcuflynx/schemas/*.schema.json)
# ----------------------------------------------------------------------------------------------

def _libcuflynx_schema(name):
    try:
        import importlib.resources as ir
        path = ir.files('libcuflynx').joinpath('schemas', name)
        return json.loads(path.read_text()) if path.is_file() else None
    except (ImportError, FileNotFoundError, ModuleNotFoundError):
        return None


VESSEL_ARRAYS = sorted(p for p in glob.glob(os.path.join(MODULES_DIR, '**', '*_vessel_array.json'), recursive=True))
MODULE_CONFIGS = sorted(glob.glob(os.path.join(MODULES_DIR, '**', '*_modules_config.json'), recursive=True))


def _validate(path, schema_name):
    jsonschema = pytest.importorskip('jsonschema')
    schema = _libcuflynx_schema(schema_name)
    if schema is None:
        pytest.skip(f'the installed libcuflynx has no schemas/{schema_name}')
    with open(path) as f:
        data = json.load(f)
    errors = sorted(jsonschema.Draft202012Validator(schema).iter_errors(data), key=lambda e: list(e.path))
    assert not errors, '\n'.join(f'{list(e.path)}: {e.message}' for e in errors[:20])


@pytest.mark.parametrize('path', VESSEL_ARRAYS, ids=lambda p: os.path.relpath(p, MODULES_DIR))
def test_vessel_array_matches_libcuflynx_schema(path):
    _validate(path, 'vessel_array.schema.json')


@pytest.mark.parametrize('path', MODULE_CONFIGS, ids=lambda p: os.path.relpath(p, MODULES_DIR))
def test_module_config_matches_libcuflynx_schema(path):
    _validate(path, 'module_config.schema.json')


def test_no_csv_vessel_arrays_left():
    left = [os.path.relpath(p, MODULES_DIR) for p in glob.glob(os.path.join(MODULES_DIR, '**', '*_vessel_array.csv'), recursive=True)
            if 'poiseuille' not in p]
    assert not left, f'CSV vessel arrays left (run tools/convert_vessel_arrays.py): {left[:10]}'
