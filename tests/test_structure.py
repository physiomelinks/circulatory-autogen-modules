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


def test_every_module_dir_has_spec():
    dirs = [d for d in os.listdir(MODULES_DIR) if os.path.isdir(os.path.join(MODULES_DIR, d)) and d != 'template']
    missing = [d for d in dirs if not os.path.isfile(os.path.join(MODULES_DIR, d, f'{d}_tests.yaml'))]
    assert not missing
