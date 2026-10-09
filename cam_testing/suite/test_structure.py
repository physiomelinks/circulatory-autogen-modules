"""
Fast static checks on the module library (no libcuflynx generation needed).

  - the directory layout of modules/ and system_models/ against modules/directory_schema.json
    and its rules (version == module_subtype, default_parameterisation, obs_data_name, unique names,
    nested module_types used only within their parent, system-model records resolve);
  - each version's CellML, config and units (test_version_structure);
  - library-wide uniqueness, manifests, and the JSON files against libcuflynx's schemas.

The repo checked is cam_testing.paths' (its modules/, system_models/ and manifests/); with extra
module libraries, uniqueness is checked over every library (libcuflynx merges them) and system-model
records and harnesses may name their versions. The directory schema is the repo's
modules/directory_schema.json, or cam_testing's copy when it has none.

    pytest tests/test_structure.py                  # in this repo
    pytest --pyargs cam_testing.suite.test_structure   # in any repo (see README.md)

A problem listed (as a substring) under ``known_issues`` in a version's spec is reported as xfail
rather than a failure, so known defects stay visible without blocking CI; once fixed, remove the
entry.
"""
import csv
import functools
import glob
import json
import os
import re
import xml.etree.ElementTree as ET

import pytest

from cam_testing import bib, library, module_array, paths
from cam_testing.library import (IDENTITY_KEYS, INSTANCE_COLUMNS, LICENCES, TESTS_KEYS, VERIFICATION_KEYS, VERSIONS,
                                 all_versions, ancestors_of, load_version, misplaced_spec_keys, module_relpath,
                                 module_type_names, parent_of, read_spec_files)
from cam_testing.mathml import component_names

CELLML_NS = 'http://www.cellml.org/cellml/1.1#'
STANDARD_UNITS = {
    'ampere', 'becquerel', 'candela', 'celsius', 'coulomb', 'dimensionless', 'farad', 'gram',
    'gray', 'henry', 'hertz', 'joule', 'katal', 'kelvin', 'kilogram', 'liter', 'litre', 'lumen',
    'lux', 'meter', 'metre', 'mole', 'newton', 'ohm', 'pascal', 'radian', 'second', 'siemens',
    'sievert', 'steradian', 'tesla', 'volt', 'watt', 'weber',
}
VARIABLE_KINDS = {'variable', 'constant', 'global_constant', 'boundary_condition'}
ROOTS = paths.roots()
REPO_ROOT, MODULES_DIR, SYSTEM_MODELS_DIR = ROOTS.repo_root, ROOTS.modules_dir, ROOTS.system_models_dir
with open(ROOTS.directory_schema_path) as _f:
    SCHEMA = json.load(_f)
LEVELS = SCHEMA['levels']
if ROOTS.directory_schema_path == paths.PACKAGED_DIRECTORY_SCHEMA:
    # a repo without its own schema is checked against cam_testing's, and needn't hold a copy
    _root_files = LEVELS['modules_root']['files']
    _root_files['required'] = [f for f in _root_files['required'] if f != 'directory_schema.json']
    _root_files['optional'] = [*_root_files.get('optional', []), 'directory_schema.json']
EXCLUDED = {os.path.join(REPO_ROOT, p) for p in SCHEMA['excluded']}

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
def _versions():
    '''The versions of the repo's own modules/ (the ones checked).'''
    return all_versions()


@functools.lru_cache(maxsize=1)
def _library_versions():
    '''The versions of every library (the repo's and the extra ones): what libcuflynx merges.'''
    return all_versions(include_libraries=True)


VERSION_KEYS = [v.key for v in _versions()]


@functools.lru_cache(maxsize=1)
def _all_component_names():
    return {c for v in _library_versions() if os.path.isfile(v.cellml_path) for c in component_names(v.cellml_path)}


# ----------------------------------------------------------------------------------------------
# the directory layout
# ----------------------------------------------------------------------------------------------

def _fill(pattern, names):
    for k, v in names.items():
        pattern = pattern.replace('{' + k + '}', v)
    return pattern


def _check_files(level, d, names, problems, supermodule=False):
    spec = LEVELS[level].get('files') or {}
    files = sorted(f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f)) and not f.startswith('.'))
    required = [_fill(f, names) for f in spec.get('required', [])]
    if not supermodule:
        required += [_fill(f, names) for f in spec.get('required_unless_supermodule', [])]
    allowed = set(required) | {_fill(f, names) for f in spec.get('optional', [])} \
        | {_fill(f, names) for f in spec.get('required_unless_supermodule', [])}
    rel = os.path.relpath(d, REPO_ROOT)
    for f in required:
        if f not in files:
            problems.append(f'{rel}: missing {f}')
    extra_pat = spec.get('data_file_pattern')
    ignored = spec.get('ignored_pattern')
    for f in files:
        if f in allowed or (ignored and re.match(ignored, f)):
            continue
        if extra_pat and re.match(extra_pat, f):
            continue
        problems.append(f'{rel}: unexpected file {f} (not allowed at the {level} level)')


def _subdirs(d):
    return sorted(x for x in os.listdir(d) if os.path.isdir(os.path.join(d, x)) and not x.startswith(('.', '__')))


def _name_ok(level, name, d, problems):
    pat = LEVELS[level].get('name_pattern')
    if pat and not re.match(pat, name):
        problems.append(f'{os.path.relpath(d, REPO_ROOT)}: name {name!r} does not match the {level} pattern {pat}')


def walk_modules():
    '''(problems, categories {path: name}, module_types {name: [paths]}, parents {module_type path: the
    path of the module_type it is nested in, or None}) for modules/.'''
    problems, cats, mts, parents = [], {}, {}, {}
    _check_files('modules_root', MODULES_DIR, {}, problems)

    def category(d):
        name = os.path.basename(d)
        _name_ok('category', name, d, problems)
        _check_files('category', d, {}, problems)
        cats[d] = name
        children = [c for c in _subdirs(d) if os.path.join(d, c) not in EXCLUDED]
        if not children:
            problems.append(f'{os.path.relpath(d, REPO_ROOT)}: an empty category')
        for c in children:
            p = os.path.join(d, c)
            if os.path.isdir(os.path.join(p, VERSIONS)):
                module_type(p)
            elif c in LEVELS or c in SCHEMA['generic_directory_names']:
                problems.append(f'{os.path.relpath(p, REPO_ROOT)}: a generic directory name where a category belongs')
            else:
                category(p)

    def module_type(d, parent=None):
        name = os.path.basename(d)
        mts.setdefault(name, []).append(d)
        parents[d] = parent
        _name_ok('module_type', name, d, problems)
        _check_files('module_type', d, {'module_type': name}, problems)
        # its other subdirectories are nested module_types (each with versions/); anything else is stray
        for c in _subdirs(d):
            p = os.path.join(d, c)
            if c == VERSIONS:
                continue
            if os.path.isdir(os.path.join(p, VERSIONS)):
                module_type(p, d)
            else:
                problems.append(f'{os.path.relpath(p, REPO_ROOT)}: unexpected directory in a module_type '
                                f'(only versions/ and nested module_types, which have versions/)')
        vroot = os.path.join(d, VERSIONS)
        if not _subdirs(vroot):
            problems.append(f'{os.path.relpath(vroot, REPO_ROOT)}: no versions')
        _check_files('versions', vroot, {}, problems)
        for v in _subdirs(vroot):
            version(os.path.join(vroot, v), name)

    def version(d, mt):
        v = os.path.basename(d)
        names = {'module_type': mt, 'version': v}
        _name_ok('version', v, d, problems)
        rel = os.path.relpath(d, REPO_ROOT)
        cfg = os.path.join(d, f'{mt}_{v}_modules_config.json')
        supermodule = False
        if os.path.isfile(cfg):
            entries = json.load(open(cfg))
            if len(entries) != 1:
                problems.append(f'{rel}: the config has {len(entries)} entries (a version has one)')
            for e in entries[:1]:
                supermodule = e.get('module_format') == 'supermodule'
                if e.get('module_type') != mt or e.get('module_subtype') != v:
                    problems.append(f'{rel}: config entry is ({e.get("module_type")}, {e.get("module_subtype")}), '
                                    f'not ({mt}, {v}): the version is its module_subtype')
                di = e.get('default_parameterisation')
                if 'default_instance' in e:
                    problems.append(f'{rel}: config entry has "default_instance" (the former name): rename it '
                                    f'"default_parameterisation"')
                if not di:
                    problems.append(f'{rel}: config entry has no default_parameterisation')
                elif not os.path.isfile(os.path.join(d, 'parameterisations', di, f'{di}_parameters.csv')):
                    problems.append(f'{rel}: default_parameterisation {di} has no parameterisations/{di}/{di}_parameters.csv')
        _check_files('version', d, names, problems, supermodule)
        # the spec's two files: every key in its own file, identity keys matching the directories
        tests, verification = read_spec_files(d, f'{mt}_{v}')
        problems.extend(f'{rel}: {p}' for p in misplaced_spec_keys(tests, verification))
        review = tests.get('review')
        if isinstance(review, str) and not (review.startswith('reviews/') and os.path.isfile(os.path.join(REPO_ROOT, review))):
            problems.append(f'{rel}: review: {review} is not a file in reviews/')
        for fname, content in (('tests.yaml', tests), ('verification_config.json', verification)):
            if content and (content.get('module_type'), content.get('version')) != (mt, v):
                problems.append(f'{rel}: {fname} names ({content.get("module_type")}, {content.get("version")}), '
                                f'not ({mt}, {v})')
        allowed = set(LEVELS['version']['directories'])
        for c in _subdirs(d):
            if c not in allowed:
                problems.append(f'{rel}: unexpected directory {c}')
        if not os.path.isdir(os.path.join(d, 'parameterisations')):
            problems.append(f'{rel}: missing parameterisations/'
                            + (' (instances/ is the former name: rename it)' if os.path.isdir(os.path.join(d, 'instances')) else ''))
            return
        if os.path.isdir(os.path.join(d, 'risk')):
            _check_files('risk', os.path.join(d, 'risk'), names, problems)
        iroot = os.path.join(d, 'parameterisations')
        _check_files('parameterisations', iroot, {}, problems)
        if not _subdirs(iroot):
            problems.append(f'{rel}: no parameterisations')
        for i in _subdirs(iroot):
            parameterisation(os.path.join(iroot, i))

    def parameterisation(d):
        i = os.path.basename(d)
        _name_ok('parameterisation', i, d, problems)
        _check_files('parameterisation', d, {'parameterisation': i}, problems)
        rel = os.path.relpath(d, REPO_ROOT)
        allowed = LEVELS['parameterisation'].get('subdirectories') or {}
        for c in _subdirs(d):
            if c not in allowed:
                problems.append(f'{rel}: unexpected directory {c}')
                continue
            sub = allowed[c]
            for f in sorted(os.listdir(os.path.join(d, c))):
                if f not in sub.get('required', []) and not re.match(sub['file_pattern'], f):
                    problems.append(f'{rel}/{c}: unexpected file {f}')
            for f in sub.get('required', []):
                if not os.path.isfile(os.path.join(d, c, f)):
                    problems.append(f'{rel}/{c}: missing {f}')
        pp = os.path.join(d, f'{i}_parameters.csv')
        if os.path.isfile(pp):
            with open(pp, newline='') as f:
                header = next(csv.reader(f), [])
            if [h.strip() for h in header] != LEVELS['parameterisation']['parameters_columns']:
                problems.append(f'{rel}: {i}_parameters.csv columns are {header}, not {INSTANCE_COLUMNS}')
        p = os.path.join(d, f'{i}_obs_data.json')
        if os.path.isfile(p):
            name = json.load(open(p)).get('obs_data_name')
            if name != i:
                problems.append(f'{rel}: {i}_obs_data.json has obs_data_name {name!r}; the parameterisation is named by it ({i!r})')

    for c in _subdirs(MODULES_DIR):
        p = os.path.join(MODULES_DIR, c)
        if p in EXCLUDED:
            continue
        if os.path.isdir(os.path.join(p, VERSIONS)):
            module_type(p)     # a top-level module_type (heart), with the module_types nested in it
        else:
            category(p)
    return problems, cats, mts, parents


def test_modules_directory_layout():
    problems, _, _, _ = walk_modules()
    assert not problems, '\n'.join(problems[:60]) + (f'\n... {len(problems)} problems' if len(problems) > 60 else '')


def test_no_directory_name_twice():
    '''No category or module_type name (nested module_types included) appears twice in modules/, and no
    category is named as a module_type.'''
    _, cats, mts, _ = walk_modules()
    problems = []
    seen = {}
    for path, name in cats.items():
        seen.setdefault(name, []).append(os.path.relpath(path, REPO_ROOT))
    for name, paths in mts.items():
        seen.setdefault(name, []).extend(os.path.relpath(p, REPO_ROOT) for p in paths)
    generic = set(SCHEMA['generic_directory_names'])
    for name, paths in sorted(seen.items()):
        if name in generic:
            problems.append(f'{name}: a generic directory name used as a category or module_type ({paths})')
        elif len(paths) > 1:
            problems.append(f'{name} appears {len(paths)} times: {paths}')
    assert not problems, '\n'.join(problems)


def test_no_nn_versions():
    '''No version name starts with nn (directory_schema.json rules.version_names): a former nn_<x> is <x>,
    a former plain nn is named by its source. Only the pairs libcuflynx writes itself are kept
    (nn_versions_allowed).'''
    allowed = set(SCHEMA.get('nn_versions_allowed') or {})
    bad = [f'{module_relpath(v.vessel_type)}/{v.name}' for v in _versions() if v.name.startswith('nn')]
    assert not [b for b in bad if b not in allowed], f'nn versions: {[b for b in bad if b not in allowed]}'
    if ROOTS.directory_schema_path != paths.PACKAGED_DIRECTORY_SCHEMA:   # the list names this library's versions
        stale = sorted(allowed - set(bad))
        assert not stale, f'nn_versions_allowed lists versions that are gone: {stale}'


def test_type_names_carry_no_year():
    '''A module_type is named by its mechanism, never by its author or year (rules.placement).'''
    bad = sorted(n for n in module_type_names() if re.search(r'(19|20)\d\d', n) or n.endswith(('_OLD', '_Gee', '_Ursino')))
    assert not bad, f'module_type names with a year or source: {bad}'


def test_spec_key_lists_match_the_schema():
    '''The key lists cam_testing splits the spec by are the ones directory_schema.json documents.'''
    keys = SCHEMA['spec_keys']
    assert (keys['both'], keys['verification_config'], keys['tests']) == \
        (list(IDENTITY_KEYS), list(VERIFICATION_KEYS), list(TESTS_KEYS))


def test_categories_are_not_module_types():
    _, cats, mts, _ = walk_modules()
    clash = sorted(set(cats.values()) & set(mts))
    assert not clash, f'category names that are also module_type names: {clash}'


def test_nested_module_types_found_by_the_library():
    '''cam_testing.library finds the same module_types, and the same nesting, as the walk of modules/.'''
    _, _, mts, parents = walk_modules()
    walked = {name: os.path.relpath(paths[0], MODULES_DIR).replace(os.sep, '/') for name, paths in mts.items()}
    assert walked == {n: module_relpath(n) for n in module_type_names()}
    walked_parents = {os.path.basename(d): (os.path.basename(p) if p else None) for d, p in parents.items()}
    assert walked_parents == {n: parent_of(n) for n in module_type_names()}


def _uses_of_module_types(version):
    '''(module_type, where) for every module_type a version's supermodule submodules and harness name.'''
    out = [(s.get('module_type') or s.get('vessel_type'), f'submodule {s.get("name")}') for s in version.submodules]
    # the version's harness, and the parameterisations' own (validation.<parameterisation>.harness)
    networks = [('harness', version.spec.get('harness'))]
    networks += [(f'validation.{inst}.harness', (kinds or {}).get('harness'))
                 for inst, kinds in (version.spec.get('validation') or {}).items() if isinstance(kinds, dict)]
    for where, network in networks:
        for row in (module_array.harness_rows(network) or []):
            if isinstance(row, (list, tuple)) and len(row) > 2:
                out.append((row[2], f'{where} record {row[0]}'))
            elif isinstance(row, dict):
                out.append((row.get('module_type') or row.get('vessel_type'), f'{where} record {row.get("name")}'))
    return out


def test_nested_module_types_used_only_within_their_parent():
    '''A nested module_type is used only within its parent: every supermodule submodule and harness record
    naming it belongs to a version of the parent, of a module_type nested (at any depth) in the parent, or
    of the nested module_type itself. System models are exempt (they may wire a parent's parts explicitly).'''
    names = set(module_type_names(include_libraries=True))
    problems = []
    for version in _versions():
        user = version.vessel_type
        for used, where in _uses_of_module_types(version):
            parent = parent_of(used) if used in names else None
            if parent is None or used == user:
                continue
            if user == parent or parent in ancestors_of(user):
                continue
            problems.append(f'{version.key}: {where} uses {used}, which is nested in {parent} '
                            f'({module_relpath(used)}); {user} ({module_relpath(user)}) is not inside {parent}')
    assert not problems, '\n'.join(problems)


def test_system_models_directory_layout():
    if not os.path.isdir(SYSTEM_MODELS_DIR):
        pytest.skip(f'no system_models/ in {REPO_ROOT}')
    problems = []
    _check_files('system_models_root', SYSTEM_MODELS_DIR, {}, problems)
    for cat in _subdirs(SYSTEM_MODELS_DIR):
        cdir = os.path.join(SYSTEM_MODELS_DIR, cat)
        _name_ok('system_category', cat, cdir, problems)
        _check_files('system_category', cdir, {}, problems)
        for model in _subdirs(cdir):
            mdir = os.path.join(cdir, model)
            _name_ok('system_model', model, mdir, problems)
            _check_files('system_model', mdir, {'model': model}, problems)
            for c in _subdirs(mdir):
                if c not in LEVELS['system_model']['directories']:
                    problems.append(f'{os.path.relpath(mdir, REPO_ROOT)}: unexpected directory {c}')
    assert not problems, '\n'.join(problems)


SYSTEM_ARRAYS = sorted(glob.glob(os.path.join(SYSTEM_MODELS_DIR, '*', '*', '*_module_array.json')))


@pytest.mark.parametrize('path', SYSTEM_ARRAYS, ids=lambda p: os.path.relpath(os.path.dirname(p), SYSTEM_MODELS_DIR))
def test_system_model_records_resolve(path):
    '''Every record (an instance) names a (module_type, version) of the library and a
    parameterisation of it. A model
    listing modules that are not in this library (known.not_in_library in its spec) is exempt for those.'''
    import yaml
    model_dir = os.path.dirname(path)
    spec_path = glob.glob(os.path.join(model_dir, '*_system.yaml'))
    spec = yaml.safe_load(open(spec_path[0])) if spec_path else {}
    index = {}
    for v in _library_versions():     # a record may name a version of an extra library
        index.setdefault((v.vessel_type, v.name), v)
    problems, outside = [], []
    for r in json.load(open(path)):
        key = (r.get('module_type') or r.get('vessel_type'), r.get('module_subtype') or r.get('BC_type'))
        v = index.get(key)
        if v is None:
            outside.append(f'{r["name"]}: ({key[0]}, {key[1]}) is not a version in the library')
            continue
        inst = module_array.parameterisation_of(r)
        if inst is None:
            problems.append(f'{r["name"]}: names no parameterisation')
        elif inst not in v.parameterisation_names():
            problems.append(f'{r["name"]}: {v.key} has no parameterisation {inst} (it has {v.parameterisation_names()})')
    if outside and not (spec.get('expected_failures') or spec.get('skip')):
        problems += outside
    assert not problems, '\n'.join(problems)


# ----------------------------------------------------------------------------------------------
# each version
# ----------------------------------------------------------------------------------------------

def version_problems(version, library_components=None):
    """
    Returns (errors, warnings) for one version. Errors break model generation or are plainly
    wrong; warnings are worth a look but can be legitimate.
    """
    errors, warnings = [], []
    entry = version.config
    key = version.key
    if version.is_supermodule:
        return errors, warnings
    try:
        components = set(component_names(version.cellml_path))
    except ET.ParseError as e:
        return [f'{os.path.basename(version.cellml_path)} does not parse: {e}'], []
    if library_components is None:
        library_components = _all_component_names()
    cellml_text = open(version.cellml_path).read()
    defined = set(_units_defined(version.units_path)) if os.path.isfile(version.units_path) else set()
    local_units = set(re.findall(r'<units\b[^>]*\bname="([^"]+)"', cellml_text))
    used_units = set(re.findall(r'\bunits="([^"]+)"', cellml_text))
    if entry.get('module_format', 'cellml') == 'cellml':
        if entry['module_file'] != os.path.basename(version.cellml_path):
            errors.append(f"{key}: module_file is {entry['module_file']}")
        if entry['module_type'] not in components:
            if entry['module_type'] in library_components:
                warnings.append(f"{key}: uses component {entry['module_type']} from another version")
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
        extra = sorted(components - {entry['module_type']})
        if extra:
            errors.append(f'{key}: components other than its own in the CellML file: {extra}')
    for unit in sorted(used_units - defined - STANDARD_UNITS - local_units):
        errors.append(f'undefined units: {unit}')

    # a reviewed version's sourced parameters must cite an entry of its references.bib
    if version.reviewed:
        keys = set(bib.read(bib.bib_path(version)))
        for inst in version.parameterisations():
            for p in inst.parameters():
                if not p.is_sourced or p.data_reference.lower().startswith('definitional'):
                    continue
                k = bib.reference_key(p.data_reference)
                if k not in keys:
                    errors.append(f'{inst.name}/{p.variable_name}: sourced, but its reference '
                                  f'"{p.data_reference[:40]}" is not a key in {os.path.basename(bib.bib_path(version))}')

    # every parameterisation's parameters are variables of the version
    kinds = version.kinds
    for inst in version.parameterisations():
        for p in inst.parameters():
            if p.variable_name not in kinds:
                warnings.append(f'parameterisation {inst.name}: parameter {p.variable_name} is not a variable of the config entry')

    root = ET.parse(version.cellml_path).getroot()
    for comp in root.iter(f'{{{CELLML_NS}}}component'):
        declared = {v.get('name') for v in comp.findall(f'{{{CELLML_NS}}}variable')}
        used = {ci.text.strip() for ci in comp.iter('{http://www.w3.org/1998/Math/MathML}ci') if ci.text}
        for name in sorted(used - declared):
            warnings.append(f"component {comp.get('name')}: equations use undeclared variable {name}")

    # config variables must exist in the CellML component with the same units and a public interface
    comps = {c.get('name'): c for c in root.iter(f'{{{CELLML_NS}}}component')}
    comp = comps.get(entry.get('module_type'))
    if comp is not None and entry.get('module_format', 'cellml') == 'cellml':
        cvars = {v.get('name'): v for v in comp.findall(f'{{{CELLML_NS}}}variable')}
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
    return errors, warnings


def write_structure_results(version, errors, warnings, known):
    """Saved for the HTML report, next to the other test results."""
    os.makedirs(version.results_dir, exist_ok=True)
    with open(os.path.join(version.results_dir, 'structure.json'), 'w') as f:
        json.dump({'errors': errors, 'warnings': warnings, 'known_issues': known}, f, indent=2)


@pytest.mark.parametrize('key', VERSION_KEYS)
def test_version_structure(key):
    version = load_version(*key.split('/', 1))
    known = version.spec.get('known_issues') or []
    errors, warnings = version_problems(version)
    write_structure_results(version, errors, warnings, known)
    unknown = [p for p in errors if not any(k in p for k in known)]
    assert not unknown, '\n'.join(unknown)
    if errors:
        pytest.xfail('known issues: ' + '; '.join(errors))


def test_library_wide_uniqueness():
    '''libcuflynx merges every version: components, (module_type, version) pairs and units must not clash,
    across every library it is given (this repo's and the extra ones). Only clashes involving this
    repo's versions are reported: an extra library checks its own.'''
    component_owner, pair_owner, units_seen = {}, {}, {}
    clashes = []      # (the two versions' directories, problem)

    def where(v):
        return v.key if v.mtype.is_primary else f'{v.key} ({v.mtype.library})'

    for version in _library_versions():
        if os.path.isfile(version.cellml_path):
            for comp in component_names(version.cellml_path):
                if comp in component_owner:
                    other = component_owner[comp]
                    clashes.append(((other.dir, version.dir), f'component {comp} is in both {where(other)} and {where(version)}'))
                component_owner[comp] = version
        pair = (version.vessel_type, version.name)
        if pair in pair_owner:
            other = pair_owner[pair]
            clashes.append(((other.dir, version.dir), f'{pair} is in both {where(other)} and {where(version)}'))
        pair_owner[pair] = version
        if os.path.isfile(version.units_path):
            for unit, el in _units_defined(version.units_path).items():
                if unit in units_seen and units_seen[unit][1] != _canonical(el):
                    other = units_seen[unit][0]
                    clashes.append(((other.dir, version.dir), f'units {unit} differ between {where(other)} and {where(version)}'))
                units_seen.setdefault(unit, (version, _canonical(el)))
    own = {v.dir for v in _versions()}
    problems = [p for dirs, p in clashes if own & set(dirs)]
    assert not problems, '\n'.join(problems)


def test_every_cellml_file_is_in_a_version():
    '''Every *_modules.cellml under modules/ is <module_type>_<version>_modules.cellml in its version directory.'''
    stray = []
    for p in glob.glob(os.path.join(MODULES_DIR, '**', '*_modules.cellml'), recursive=True):
        if any(p.startswith(e + os.sep) for e in EXCLUDED):
            continue
        d = os.path.dirname(p)
        v, mt = os.path.basename(d), os.path.basename(os.path.dirname(os.path.dirname(d)))
        if os.path.basename(os.path.dirname(d)) != VERSIONS or os.path.basename(p) != f'{mt}_{v}_modules.cellml':
            stray.append(os.path.relpath(p, MODULES_DIR))
    assert not stray, f'module files outside a version directory: {stray}'


@pytest.mark.parametrize('manifest', sorted(glob.glob(os.path.join(REPO_ROOT, 'manifests', '*.json'))),
                         ids=os.path.basename)
def test_manifest_paths_exist(manifest):
    with open(manifest) as f:
        data = json.load(f)
    missing = [e['path'] for entries in data['collections'].values() for e in entries
               if not os.path.isfile(os.path.join(REPO_ROOT, e['path']))]
    assert not missing, f'missing: {missing}'


# ----------------------------------------------------------------------------------------------
# JSON files against libcuflynx's schemas (libcuflynx/schemas/*.schema.json)
# ----------------------------------------------------------------------------------------------

def _libcuflynx_schema(name):
    try:
        import importlib.resources as ir
        path = ir.files('libcuflynx').joinpath('schemas', name)
        return json.loads(path.read_text()) if path.is_file() else None
    except (ImportError, FileNotFoundError, ModuleNotFoundError):
        return None


def _not_excluded(p):
    return not any(p.startswith(e + os.sep) for e in EXCLUDED)


MODULE_ARRAYS = SYSTEM_ARRAYS
MODULE_CONFIGS = sorted(p for p in glob.glob(os.path.join(MODULES_DIR, '**', '*_modules_config.json'), recursive=True)
                        if _not_excluded(p))
OBS_DATA = sorted(p for p in glob.glob(os.path.join(MODULES_DIR, '**', 'parameterisations', '*', '*obs_data.json'), recursive=True)
                  if _not_excluded(p))


def _validate(path, schema_name, *older_names):
    '''older_names: the schema's name in older libcuflynx releases (vessel_array before #549).'''
    jsonschema = pytest.importorskip('jsonschema')
    schema = next((s for s in map(_libcuflynx_schema, (schema_name,) + older_names) if s is not None), None)
    if schema is None:
        pytest.skip(f'the installed libcuflynx has no schemas/{schema_name}')
    with open(path) as f:
        data = json.load(f)
    errors = sorted(jsonschema.Draft202012Validator(schema).iter_errors(data), key=lambda e: list(e.path))
    assert not errors, '\n'.join(f'{list(e.path)}: {e.message}' for e in errors[:20])


@pytest.mark.parametrize('path', MODULE_ARRAYS, ids=lambda p: os.path.relpath(p, REPO_ROOT))
def test_module_array_matches_libcuflynx_schema(path):
    _validate(path, 'module_array.schema.json', 'vessel_array.schema.json')


@pytest.mark.parametrize('path', MODULE_CONFIGS, ids=lambda p: os.path.relpath(p, MODULES_DIR))
def test_module_config_matches_libcuflynx_schema(path):
    _validate(path, 'module_config.schema.json')


def licence_and_creator_problems(entry):
    '''The config entry's "licence" (an SPDX id from library.LICENCES) and "creator" (a list of names).'''
    problems = []
    if 'licence' not in entry:
        problems.append('no "licence"')
    elif entry['licence'] not in LICENCES:
        problems.append(f'licence {entry["licence"]!r} is not one of {sorted(LICENCES)}')
    if 'creator' not in entry:
        problems.append('no "creator"')
    elif not (isinstance(entry['creator'], list) and all(isinstance(c, str) and c.strip() for c in entry['creator'])):
        problems.append(f'creator {entry["creator"]!r} is not a list of names')
    return problems


def test_every_config_entry_has_licence_and_creator():
    '''Every version's config entry has "licence" (CC0-1.0, CC-BY-4.0, MIT, Apache-2.0 or 0BSD) and
    "creator" (a list of names, empty until given in review).'''
    problems = []
    for path in MODULE_CONFIGS:
        with open(path) as f:
            for e in json.load(f):
                problems += [f'{os.path.relpath(path, MODULES_DIR)}: {p}' for p in licence_and_creator_problems(e)]
    assert not problems, '\n'.join(problems[:30])


def test_licence_and_creator_rules():
    assert licence_and_creator_problems({'licence': 'CC0-1.0', 'creator': []}) == []
    assert licence_and_creator_problems({'licence': 'MIT', 'creator': ['A. Author']}) == []
    assert licence_and_creator_problems({'creator': []}) == ['no "licence"']
    assert 'is not one of' in licence_and_creator_problems({'licence': 'GPL-3.0', 'creator': []})[0]
    assert 'not a list of names' in licence_and_creator_problems({'licence': 'MIT', 'creator': 'A. Author'})[0]
    assert licence_and_creator_problems({'licence': 'MIT'}) == ['no "creator"']


# ----------------------------------------------------------------------------------------------
# required_citations and the bib key convention (2026-10-06; modules/README.md, "Citations")
# ----------------------------------------------------------------------------------------------

def required_citations_problems(entry, bib_keys):
    '''The config entry's "required_citations" (a non-empty list of keys in the version's references.bib)
    and "required_citations_uncertain" (optional: a non-empty subset of them).'''
    problems = []
    rc = entry.get('required_citations')
    if not rc:
        problems.append('no "required_citations"' if rc is None else '"required_citations" is empty')
    elif not (isinstance(rc, list) and all(isinstance(k, str) and k for k in rc)):
        problems.append(f'required_citations {rc!r} is not a list of bib keys')
    else:
        if len(set(rc)) != len(rc):
            problems.append(f'required_citations repeats a key: {rc}')
        missing = [k for k in rc if k not in bib_keys]
        if missing:
            problems.append(f'required_citations not in the references.bib: {missing}')
    if 'required_citations_uncertain' in entry:
        unc = entry['required_citations_uncertain']
        if not (isinstance(unc, list) and unc and all(isinstance(k, str) for k in unc)):
            problems.append(f'required_citations_uncertain {unc!r} is not a non-empty list of keys (leave it out when none is)')
        elif not set(unc) <= set(rc or []):
            problems.append(f'required_citations_uncertain not in required_citations: {sorted(set(unc) - set(rc or []))}')
    return problems


def test_every_version_has_required_citations():
    '''Every version's config entry names the papers its model is built on: "required_citations", a
    non-empty list of keys that are in the version's references.bib ("required_citations_uncertain"
    marks the best guesses among them).'''
    problems = []
    for v in _versions():
        with open(v.path('modules_config.json')) as f:
            entries = json.load(f)
        keys = set(bib.read(bib.bib_path(v)))
        for e in entries:
            problems += [f'{v.key}: {p}' for p in required_citations_problems(e, keys)]
    assert not problems, '\n'.join(problems[:40]) + (f'\n... {len(problems)} in all' if len(problems) > 40 else '')


def test_required_citations_rules():
    keys = {'paci2013computational', 'tao2011model'}
    assert required_citations_problems({'required_citations': ['paci2013computational']}, keys) == []
    assert required_citations_problems({'required_citations': ['tao2011model', 'paci2013computational'],
                                        'required_citations_uncertain': ['tao2011model']}, keys) == []
    assert required_citations_problems({}, keys) == ['no "required_citations"']
    assert required_citations_problems({'required_citations': []}, keys) == ['"required_citations" is empty']
    assert 'not in the references.bib' in required_citations_problems({'required_citations': ['Paci2013']}, keys)[0]
    assert 'not in required_citations' in required_citations_problems(
        {'required_citations': ['tao2011model'], 'required_citations_uncertain': ['paci2013computational']}, keys)[0]
    assert 'non-empty' in required_citations_problems(
        {'required_citations': ['tao2011model'], 'required_citations_uncertain': []}, keys)[0]


def test_config_key_lists_match_the_schema():
    '''The record keys cam_testing checks are the ones directory_schema.json documents.'''
    assert list(SCHEMA['config_entry_record_keys']) == list(library.RECORD_KEYS) + list(library.CITATION_KEYS)


# Keys allowed to break the <surname><year><word> shape, because their entries have no year (web pages
# without a publication date); each still equals cam_testing.bib.convention_key of its entry. Only add a
# key here when the source really has no date.
UNDATED_BIB_KEYS = {
    'ranjanchannelpedia',                 # Channelpedia Kv1.5 model 21 (no date on the model page)
    'channelpedia',                       # Channelpedia Kv4.2 model 40 (no author or date on the model page)
    'usnationallibraryofmedicinerapid',   # MedlinePlus encyclopedia page "Rapid shallow breathing"
}


def _bib_files():
    return sorted(p for p in glob.glob(os.path.join(MODULES_DIR, '**', '*_references*.bib'), recursive=True)
                  if _not_excluded(p))


def bib_key_problems(key, raw_fields, n_words_max=4):
    '''A key follows the convention: lowercase <first author's surname><year><first significant title
    word> (cam_testing.bib.convention_key), or with further title words when two papers would share it.'''
    allowed = [bib.convention_key(raw_fields, n) for n in range(1, n_words_max + 1)]
    if key not in allowed:
        return [f'{key}: not the convention key of its entry ({allowed[0]})']
    if not bib.KEY_RE.match(key) and key not in UNDATED_BIB_KEYS:
        return [f'{key}: not <surname><year><word> (no year in the entry? add one, or list it in UNDATED_BIB_KEYS)']
    return []


def test_bib_keys_follow_the_convention():
    '''Every key in every references.bib / references_proposed.bib is the lowercase
    <surname><year><first title word> of its entry (tools/rekey_bib.py rekeyed the library on 2026-10-06).'''
    problems = []
    for p in _bib_files():
        with open(p, encoding='utf-8') as f:
            text = f.read()
        problems += [f'{os.path.relpath(p, MODULES_DIR)}: {m}'
                     for _, key, _, _, raw in bib.entries(text) for m in bib_key_problems(key, raw)]
    assert not problems, '\n'.join(problems[:40]) + (f'\n... {len(problems)} in all' if len(problems) > 40 else '')


def test_one_key_per_paper():
    '''No two keys cite the same paper (same DOI) anywhere in the library.'''
    by_doi = {}
    for p in _bib_files():
        for key, fields in bib.read(p).items():
            doi = (fields.get('doi') or '').strip().lower()
            if doi:
                by_doi.setdefault(doi, set()).add(key)
    dup = {d: sorted(k) for d, k in by_doi.items() if len(k) > 1}
    assert not dup, f'papers cited under two keys (merge them): {dup}'


def test_bib_key_rules():
    raw = {'author': 'Belluzzi, O. and Sacchi, O.', 'year': '1986',
           'title': 'A quantitative description of the sodium current in the rat sympathetic neurone'}
    assert bib.convention_key(raw) == 'belluzzi1986quantitative'
    assert bib.convention_key({'author': 'van der Pol, Balthasar', 'year': '1926', 'title': 'On relaxation-oscillations'}) \
        == 'vanderpol1926relaxation'
    assert bib.convention_key({'author': 'Hern{\\\'a}ndez-Cruz, A.', 'year': '1997', 'title': 'Ca2+ release'}) == 'hernandezcruz1997ca'
    assert bib_key_problems('belluzzi1986quantitative', raw) == []
    assert bib_key_problems('belluzzi1986quantitativedescription', raw) == []      # disambiguated
    assert bib_key_problems('BelluzziSacchi1986', raw)


@pytest.mark.parametrize('path', OBS_DATA, ids=lambda p: os.path.relpath(p, MODULES_DIR))
def test_obs_data_matches_libcuflynx_schema(path):
    _validate(path, 'obs_data.schema.json')


def test_no_csv_module_arrays_left():
    left = [os.path.relpath(p, REPO_ROOT) for name in module_array.NAMES
            for p in glob.glob(os.path.join(SYSTEM_MODELS_DIR, '**', f'*_{name}.csv'), recursive=True)
            if os.sep + 'reference' + os.sep not in p]
    assert not left, f'CSV module arrays left (run tools/convert_module_arrays.py): {left[:10]}'


def test_module_arrays_use_the_new_name():
    '''Module arrays were called vessel arrays: the library's own files and verification configs use the new
    name. reference/ keeps circulatory_autogen's originals as they were, and the readers take both.'''
    old_files = [os.path.relpath(p, REPO_ROOT)
                 for p in glob.glob(os.path.join(SYSTEM_MODELS_DIR, '**', '*_vessel_array.*'), recursive=True)
                 if os.sep + 'reference' + os.sep not in p]
    old_keys = [os.path.relpath(p, REPO_ROOT)
                for p in glob.glob(os.path.join(MODULES_DIR, '**', '*_verification_config.json'), recursive=True)
                if '"vessel_array"' in open(p).read()]
    assert not old_files and not old_keys, (f'rename to module_array: files {old_files[:10]}, '
                                            f'harness keys in {old_keys[:10]}')


def test_old_module_array_names_are_still_read(tmp_path):
    (tmp_path / 'm_vessel_array.json').write_text('[]')
    assert module_array.find(str(tmp_path), 'm').endswith('m_vessel_array.json')
    (tmp_path / 'm_module_array.json').write_text('[]')
    assert module_array.find(str(tmp_path), 'm').endswith('m_module_array.json')
    rows = [['mod', 'nn', 'x', '', '']]
    assert module_array.harness_rows({'vessel_array': rows}) == rows
    assert module_array.harness_rows({'module_array': rows}) == rows


def test_module_type_names_listed():
    assert module_type_names(), 'no module_types found under modules/'


# Parameterisations whose data were extracted from a publication but have no source screenshot yet.
# Remove an entry when its parameterisations/<p>/source_figures/ is added; never add new ones.
MISSING_SOURCE_FIGURES = {
    'capillary/pp_micro::default', 'heart/vp::default', 'heart/vp_Ca::default',
    'heart/vp_new_valve::default', 'heart/vp_wCont::default', 'heart/vp_wCont_nonstiff::default',
    'inlet_flow/adan::boileau2015_adan56_inflow', 'inlet_flow/adan_2::boileau2015_adan56_inflow',
    'inlet_flow/aorticbif::boileau2015_ibif_inflow', 'Lotka_Volterra/Lotka1925_v01::carpenter2018',
    'Lotka_Volterra/Lotka1925_v01::hudson_bay_lynx_hare', 'pulmonary_GE/Albanese2016_v01::pulmonary_GE_normal_blood_gases',
    # wanaverbecq2003plasma: no PDF on disk (not in ~/Zotero/storage or ~/Documents/papers; 2026-10-07); values from the full text
    'PMCA/Colegrove2000_v01::wanaverbecq2003_scg',
}


def _publication_parameterisations():
    return [(f'{v.key}::{i.name}', i) for v in library.all_versions() for i in v.parameterisations() if i.needs_source_figures()]


@pytest.mark.parametrize('key,inst', _publication_parameterisations(), ids=lambda x: x if isinstance(x, str) else '')
def test_publication_data_has_source_figures(key, inst):
    '''Validation or calibration data extracted from a paper or book carries a screenshot of the
    figure/table it came from (parameterisations/<p>/source_figures/, listed in source_figures.json), shown
    in the report beside the validation plots (modules/README.md, "Source figures").'''
    if key in MISSING_SOURCE_FIGURES:
        pytest.xfail('source screenshot not added yet')
    figures = inst.source_figures()
    assert figures, f'{key}: data from a publication but no {library.SOURCE_FIGURES}/{library.SOURCE_FIGURES_INDEX}'
    for f in figures:
        assert f.get('source') and f.get('file'), f'{key}: each source figure needs "file" and "source"'
        assert os.path.isfile(os.path.join(inst.version.dir, f['file'])), f'{key}: {f["file"]} missing'


def test_missing_source_figures_list_is_current():
    '''The known-gap list only shrinks: every entry still needs figures and still lacks them. (Entries
    are this library's versions; another repo's run ignores them.)'''
    have = {k for k, i in _publication_parameterisations() if not i.source_figures()}
    if os.path.realpath(REPO_ROOT) != os.path.realpath(paths.PACKAGE_CHECKOUT):
        pytest.skip('the known-gap list is circulatory-autogen-modules\' own')
    stale = MISSING_SOURCE_FIGURES - have
    assert not stale, f'remove from MISSING_SOURCE_FIGURES (now has figures or no longer needs them): {sorted(stale)}'


def _ca_resources():
    '''circulatory_autogen's resources/ directory, from the installed libcuflynx (or the sibling clone).'''
    try:
        import libcuflynx
        d = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(libcuflynx.__file__))), 'resources')
        if os.path.isdir(d):
            return d
    except Exception:  # noqa: BLE001
        pass
    d = os.path.join(REPO_ROOT, '..', 'circulatory_autogen', 'resources')
    return d if os.path.isdir(d) else None


def _reference_copies():
    return sorted(p for p in glob.glob(os.path.join(SYSTEM_MODELS_DIR, '*', '*', 'reference', '*'))
                  if p.endswith(('_parameters.csv', '_obs_data.json', '_params_for_id.csv')))


@pytest.mark.parametrize('path', _reference_copies(), ids=lambda p: os.path.relpath(p, SYSTEM_MODELS_DIR))
def test_reference_copies_are_circulatory_autogen_originals_in_library_keys(path):
    '''system_models/**/reference/ holds circulatory_autogen's original files with only the bib keys
    changed to this library's (tools/bib_rekey_map.json): each equals the original after the same rekey
    (tools/rekey_bib.rekey_string), so the copies stay exact and re-importing reproduces them.'''
    import importlib.util
    res = _ca_resources()
    original = os.path.join(res, os.path.basename(path)) if res else None
    if not original or not os.path.isfile(original):
        pytest.skip('no circulatory_autogen original to compare with')
    spec = importlib.util.spec_from_file_location('rekey_bib', os.path.join(REPO_ROOT, 'tools', 'rekey_bib.py'))
    rekey = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rekey)
    with open(original, encoding='utf-8') as f:
        expected = rekey.rekey_string(path, f.read())
    with open(path, encoding='utf-8') as f:
        assert f.read() == expected, f'{path} differs from {original} (after rekeying) beyond its bib keys'
