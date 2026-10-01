"""
The module library on disk: module_types, their versions and each version's instances.

Layout (modules/README.md and modules/directory_schema.json describe it in full):

    modules/<category path>/<module_type>/
        <module_type>.html                          report across the versions (generated)
        versions/<version>/                         version == the config entry's module_subtype
            <module_type>_<version>_modules.cellml
            <module_type>_<version>_modules_config.json   one entry; "default_instance"
            <module_type>_<version>_units.cellml
            <module_type>_<version>_verification_config.json   what the checks run
            <module_type>_<version>_tests.yaml            review and record-keeping
            <module_type>_<version>_references.bib        (+ _references_proposed.bib)
            <module_type>_<version>.html                  version report (generated)
            risk/  plots/  results/                       (plots, results generated)
            instances/<instance>/                         instance == obs_data_name
                <instance>_parameters.csv                 variable_name,units,value,data_reference,sourced
                <instance>_obs_data.json                  optional: calibration data (data_items) and
                                                          held-out data (prediction_items with a value)
                <instance>_params_for_id.csv              optional
                <instance>_calibrated_parameters.csv      written by calibration (committed)
                <instance>_calibration.json               written by calibration (committed)

A module_type directory is any directory under modules/ with a versions/ subdirectory; every
other directory above one is a category. System models live in system_models/.
"""
import copy
import csv
import functools
import json
import os
import re
from dataclasses import dataclass, field

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES_DIR = os.path.join(REPO_ROOT, 'modules')
SYSTEM_MODELS_DIR = os.path.join(REPO_ROOT, 'system_models')
# reviews shared by several versions (a whole pre-versions module's review): tests.yaml's
# review: reviews/<name>_review.yaml points to one copy
REVIEWS_DIR = os.path.join(REPO_ROOT, 'reviews')

VERSIONS, INSTANCES, DEFAULT_INSTANCE = 'versions', 'instances', 'default'
# directories under modules/ that are not categories or module_types (another session's work in
# the old layout, until it is moved)
UNMIGRATED = ('poiseuille_transport', 'system')
GENERATED_DIRS = ('plots', 'results')
# a supermodule's globals when its spec doesn't list supermodule.globals
SUPERMODULE_DEFAULT_GLOBALS = ('T', 'rho', 'l_eff')


def deep_merge(base, override):
    out = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def safe_id(text):
    '''A filesystem/URL-safe id, e.g. for work directories and HTML anchors.'''
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', text)


# ----------------------------------------------------------------------------------------------
# discovery
# ----------------------------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _module_type_dirs():
    '''module_type name -> directory, for every modules/**/<name>/ that has versions/.'''
    found = {}
    if not os.path.isdir(MODULES_DIR):
        return found
    for root, dirs, files in os.walk(MODULES_DIR):
        rel = os.path.relpath(root, MODULES_DIR)
        if rel != '.' and rel.split(os.sep)[0] in UNMIGRATED:
            dirs[:] = []
            continue
        if VERSIONS in dirs:
            name = os.path.basename(root)
            if name in found:
                raise ValueError(f'two module_types are called {name}: {found[name]} and {root}')
            found[name] = root
            dirs[:] = []
            continue
        dirs[:] = sorted(d for d in dirs if not d.startswith('.') and d not in GENERATED_DIRS)
    return found


def module_type_names():
    '''Every module_type, in case-insensitive order.'''
    return sorted(_module_type_dirs(), key=str.lower)


# the old name: a "module" is now a module_type
module_names = module_type_names


def module_type_dir(name):
    return _module_type_dirs().get(name) or os.path.join(MODULES_DIR, name)


module_dir = module_type_dir


def category_of(name):
    '''The category path of a module_type, e.g. "cell/neurons".'''
    return os.path.relpath(os.path.dirname(module_type_dir(name)), MODULES_DIR).replace(os.sep, '/')


def module_relpath(name):
    '''The module_type directory relative to modules/, e.g. "cell/neurons/soma".'''
    return os.path.relpath(module_type_dir(name), MODULES_DIR).replace(os.sep, '/')


def categories():
    '''Every category path (with its parents), sorted.'''
    out = set()
    for name in module_type_names():
        parts = category_of(name).split('/')
        for i in range(1, len(parts) + 1):
            out.add('/'.join(parts[:i]))
    return sorted(out)


def matches_selector(name, selector):
    '''--module selector: a module_type name, or a category path (prefix) such as "cell" or "cell/neurons".'''
    if selector == name:
        return True
    cat = category_of(name)
    sel = selector.strip('/')
    return cat == sel or cat.startswith(sel + '/')


def select_module_types(selectors):
    names = module_type_names()
    if not selectors:
        return names
    return [n for n in names if any(matches_selector(n, s) for s in selectors)]


# ----------------------------------------------------------------------------------------------
# config formats
# ----------------------------------------------------------------------------------------------

# Module configs come in two formats: PhLynx's (this library's) and libcuflynx's original.
# Both are read; code here uses libcuflynx's internal names. module_type means different things in
# the two, so the format is told by the other keys, never by module_type.
PHLYNX_KEYS = {'module_type': 'vessel_type', 'module_subtype': 'BC_type',
               'component_file': 'module_file', 'component_type': 'module_type'}
LIBCUFLYNX_KEYS = ('vessel_type', 'BC_type', 'module_file', 'module_type')


def is_supermodule(entry):
    return entry.get('module_format') == 'supermodule'


def config_format(entry):
    if is_supermodule(entry):
        # a supermodule has no CellML component: only its type/subtype keys tell the style
        return 'libcuflynx' if {'vessel_type', 'BC_type'} & set(entry) else 'phlynx'
    phlynx = {'module_subtype', 'component_file', 'component_type'} & set(entry)
    libcuflynx = {'vessel_type', 'BC_type', 'module_file'} & set(entry)
    if phlynx and libcuflynx:
        raise ValueError(f'config entry mixes PhLynx keys {sorted(phlynx)} with libcuflynx keys {sorted(libcuflynx)}')
    if phlynx:
        missing = set(PHLYNX_KEYS) - set(entry)
        if missing:
            raise ValueError(f'PhLynx-format config entry is missing {sorted(missing)}')
        return 'phlynx'
    missing = set(LIBCUFLYNX_KEYS) - set(entry)
    if missing:
        raise ValueError(f'config entry is missing {sorted(missing)} (neither PhLynx nor libcuflynx format)')
    return 'libcuflynx'


def normalise_config_entry(entry):
    '''A config entry in either format -> libcuflynx's names (other keys unchanged).'''
    if config_format(entry) == 'libcuflynx':
        return dict(entry)
    if is_supermodule(entry):
        out = {k: v for k, v in entry.items() if k not in ('module_type', 'module_subtype')}
        return {'vessel_type': entry['module_type'], 'BC_type': entry['module_subtype'], **out}
    out = {k: v for k, v in entry.items() if k not in PHLYNX_KEYS}
    out.update({internal: entry[k] for k, internal in PHLYNX_KEYS.items()})
    return out


def to_phlynx_entry(entry):
    '''A config entry in either format -> PhLynx's key names, in PhLynx's key order.'''
    e = normalise_config_entry(entry)
    head = {'module_type': e['vessel_type'], 'module_subtype': e['BC_type'], 'module_format': e.get('module_format'),
            'component_file': e['module_file'], 'component_type': e['module_type']}
    if head['module_format'] is None:
        del head['module_format']
    return {**head, **{k: v for k, v in e.items() if k not in LIBCUFLYNX_KEYS and k != 'module_format'}}


def read_config(path):
    with open(path) as f:
        return [normalise_config_entry(e) for e in json.load(f)]


# ----------------------------------------------------------------------------------------------
# parameters
# ----------------------------------------------------------------------------------------------

INSTANCE_COLUMNS = ['variable_name', 'units', 'value', 'data_reference', 'sourced']
# the pre-versions per-module parameters file (tools that import into the old layout still write it)
PARAMETER_COLUMNS = ['vessel_type', 'BC_type', 'variable_name', 'units', 'value', 'kind', 'data_reference',
                     'sourced']

# data_reference text that does not cite a source (assumptions, tuning, placeholders, TODOs)
UNSOURCED_PATTERN = re.compile(
    r'(^$|^none$|^test$|^known$|todo|to_do|to_be|user|defined_by|placeholder|prototype|identified|'
    r'assum|wont_be|not_used|hand.?tuned|tuned|chosen|check|calculated|calc_|from_icu|nonstiff|'
    r'initial_guess|rate_constant|cam_testing|circulatory_autogen example)', re.I)


def looks_sourced(reference):
    '''Heuristic for seeding the sourced column: does the reference cite a source?'''
    ref = (reference or '').strip()
    if ':' in ref and not ref.lower().startswith('circulatory_autogen'):
        ref = ref.split(':', 1)[1].strip()   # "<CA example model>: <reference>" from the importer
    return not UNSOURCED_PATTERN.search(ref)


@dataclass
class Parameter:
    variable_name: str
    units: str
    value: str
    data_reference: str
    sourced: str = ''
    kind: str = ''                  # from the config's variables_and_units
    proposed: bool = False          # value comes from a review proposal (not in the parameters file yet)
    model_name: str = ''            # its name in a generated model, when known (a supermodule's flattened parameters)

    @property
    def is_global(self):
        return self.kind == 'global_constant'

    @property
    def is_sourced(self):
        return self.sourced.strip().lower() in ('yes', 'true', '1')

    @property
    def is_todo(self):
        return self.value.strip().upper() == 'TODO'

    @property
    def float_value(self):
        return float(self.value)


def read_parameters(path, kinds=None):
    '''An instance's parameters (kind from the config's variables_and_units, where it has one).'''
    kinds = kinds or {}
    out = []
    if not os.path.isfile(path):
        return out
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            row = {(k or '').strip(): (v or '').strip() for k, v in row.items()}
            if not row.get('variable_name'):
                continue
            out.append(Parameter(row['variable_name'], row.get('units', ''), row.get('value', ''),
                                 row.get('data_reference', ''), row.get('sourced', ''),
                                 kinds.get(row['variable_name'], '')))
    return out


def write_parameters(path, parameters):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(INSTANCE_COLUMNS)
        for p in parameters:
            w.writerow([p.variable_name, p.units, p.value, p.data_reference, p.sourced])


# ----------------------------------------------------------------------------------------------
# the version spec: two files
# ----------------------------------------------------------------------------------------------
#
# <module_type>_<version>_verification_config.json holds everything the checks read;
# <module_type>_<version>_tests.yaml holds the review and record-keeping. They are merged into one
# spec (Version.spec); the verification config wins on the keys it owns. module_type and version
# identify the version in both. modules/README.md documents every key.

IDENTITY_KEYS = ('module_type', 'version')
# in this order in the JSON file
VERIFICATION_KEYS = ('time_label', 'sim_time', 'pre_time', 'dt', 'solver', 'solver_info', 'reference_solver_info',
                     'outputs', 'run_parameters', 'harness', 'invariants', 'bc_sweep', 'timestep', 'stability',
                     'parameter_ranges', 'validation', 'supermodule')
TESTS_KEYS = ('reviewed', 'description', 'notes', 'skip', 'known_issues', 'expected_failures',
              'calibration_in_supermodule', 'reference_proposals', 'review', 'review_scope')
# readable order of the nested keys the checks know (others follow in their own order)
NESTED_ORDER = {
    'harness': ('vessel_array', 'parameters'),
    'bc_sweep': ('sweep', 'factors', 'points', 'ranges', 'bounds', 'constraints', 'exclude', 'extra_parameters', 'plot_output', 'rationale'),
    'timestep': ('scheme', 'dts', 't_end', 'min_order', 'tol', 'cvode_tol', 'roundoff', 'wrapped'),
    'stability': ('supported', 'cvode', 'solve_ivp', 'fixed_step', 'max_step_start', 'min_step', 'time_budget', 't_end',
                  'tol'),
    'supermodule': ('globals', 'equivalent'),
    'baseline': ('status', 'kind', 'note', 'reason', 'source', 'source_short', 'data', 'time_column', 'time_offset', 'variables',
                 'parameters', 'parameter_ranges', 'targets', 'z_threshold', 'pre_time', 'sim_time', 'dt', 'metric',
                 'threshold'),
    'calibrate': ('status', 'reason', 'note', 'source', 'obs_data', 'params_for_id',
                  'prediction_window', 'initial_parameters', 'starts', 'expected_parameters', 'expected_compare',
                  'expected_rtol', 'z_threshold', 'method', 'optimiser_options', 'do_ad', 'metric', 'threshold'),
}


def _ordered(d, order):
    return {**{k: d[k] for k in order if k in d}, **{k: v for k, v in d.items() if k not in order}}


def split_spec(spec):
    '''A merged spec -> (tests.yaml dict, verification_config dict), each in its stable key order.
    A key in neither list stays in tests.yaml (the structure test reports it).'''
    ident = {k: spec[k] for k in IDENTITY_KEYS if k in spec}
    verification = dict(ident)
    for k in VERIFICATION_KEYS:
        if k not in spec:
            continue
        v = copy.deepcopy(spec[k])
        if isinstance(v, dict) and k in NESTED_ORDER:
            v = _ordered(v, NESTED_ORDER[k])
        if k == 'validation' and isinstance(v, dict):
            v = {inst: _ordered({kind: _ordered(blk, NESTED_ORDER.get(kind, ())) if isinstance(blk, dict) else blk
                                 for kind, blk in (kinds or {}).items()}, ('baseline', 'calibrate'))
                 for inst, kinds in v.items()}
        verification[k] = v
    tests = dict(ident)
    tests.update({k: spec[k] for k in TESTS_KEYS if k in spec})
    tests.update({k: v for k, v in spec.items() if k not in tests and k not in verification})
    return tests, verification


def merge_spec(tests, verification):
    return {**(tests or {}), **(verification or {})}


def read_spec_files(vdir, stem):
    '''(tests.yaml dict, verification_config dict) of a version directory ({} where a file is missing).'''
    tests, verification = {}, {}
    tp = os.path.join(vdir, f'{stem}_tests.yaml')
    vp = os.path.join(vdir, f'{stem}_verification_config.json')
    if os.path.isfile(tp):
        with open(tp) as f:
            tests = yaml.safe_load(f) or {}
    if os.path.isfile(vp):
        with open(vp) as f:
            verification = json.load(f)
    return tests, verification


def write_spec(vdir, stem, spec):
    '''Writes a merged spec as <stem>_tests.yaml and <stem>_verification_config.json.'''
    tests, verification = split_spec(spec)
    with open(os.path.join(vdir, f'{stem}_tests.yaml'), 'w') as f:
        yaml.safe_dump(tests, f, sort_keys=False, width=120, default_flow_style=False, allow_unicode=True)
    with open(os.path.join(vdir, f'{stem}_verification_config.json'), 'w') as f:
        f.write(dumps_json(verification) + '\n')


def dumps_json(data, indent=2, width=110):
    '''JSON indented by ``indent``, with lists of scalars (and short all-scalar objects) kept on one
    line when they fit in ``width``: readable, and stable for diffs.'''
    def scalar(x):
        return not isinstance(x, (dict, list))

    def enc(x, level):
        pad, inner = ' ' * (indent * level), ' ' * (indent * (level + 1))
        if isinstance(x, list) and x and all(scalar(i) for i in x):
            one = json.dumps(x, ensure_ascii=False)
            if len(one) + len(inner) <= width:
                return one
        if isinstance(x, dict) and x and all(scalar(i) for i in x.values()):
            one = json.dumps(x, ensure_ascii=False)
            if len(one) + len(inner) <= width:
                return one
        if isinstance(x, list):
            if not x:
                return '[]'
            return '[\n' + ',\n'.join(inner + enc(i, level + 1) for i in x) + '\n' + pad + ']'
        if isinstance(x, dict):
            if not x:
                return '{}'
            return '{\n' + ',\n'.join(f'{inner}{json.dumps(k, ensure_ascii=False)}: {enc(v, level + 1)}'
                                      for k, v in x.items()) + '\n' + pad + '}'
        return json.dumps(x, ensure_ascii=False)
    return enc(data, 0)


def misplaced_spec_keys(tests, verification):
    '''Keys in the wrong spec file, or in neither list.'''
    problems = []
    for k in tests:
        if k in IDENTITY_KEYS or k in TESTS_KEYS:
            continue
        where = 'verification_config.json' if k in VERIFICATION_KEYS else 'neither file (unknown key)'
        problems.append(f'tests.yaml has {k!r}, which belongs in {where}')
    for k in verification:
        if k in IDENTITY_KEYS or k in VERIFICATION_KEYS:
            continue
        where = 'tests.yaml' if k in TESTS_KEYS else 'neither file (unknown key)'
        problems.append(f'verification_config.json has {k!r}, which belongs in {where}')
    return problems


# ----------------------------------------------------------------------------------------------
# module_type / version / instance
# ----------------------------------------------------------------------------------------------

@dataclass
class ModuleType:
    name: str
    dir: str

    @property
    def category(self):
        return os.path.relpath(os.path.dirname(self.dir), MODULES_DIR).replace(os.sep, '/')

    @property
    def relpath(self):
        return os.path.relpath(self.dir, MODULES_DIR).replace(os.sep, '/')

    @property
    def versions_dir(self):
        return os.path.join(self.dir, VERSIONS)

    @property
    def html_path(self):
        return os.path.join(self.dir, f'{self.name}.html')

    def version_names(self):
        if not os.path.isdir(self.versions_dir):
            return []
        return sorted((d for d in os.listdir(self.versions_dir)
                       if os.path.isdir(os.path.join(self.versions_dir, d)) and not d.startswith('.')), key=str.lower)

    def versions(self):
        return [load_version(self.name, v) for v in self.version_names()]

    def version(self, name):
        return load_version(self.name, name)

    @property
    def reviewed(self):
        vs = self.versions()
        return bool(vs) and all(v.reviewed for v in vs)


@dataclass
class Instance:
    version: 'Version'
    name: str

    @property
    def dir(self):
        return os.path.join(self.version.dir, INSTANCES, self.name)

    @property
    def key(self):
        return f'{self.version.key}/{self.name}'

    @property
    def id(self):
        return safe_id(f'{self.version.id}__{self.name}')

    def path(self, suffix):
        return os.path.join(self.dir, f'{self.name}_{suffix}')

    @property
    def parameters_path(self):
        return self.path('parameters.csv')

    @property
    def obs_data_path(self):
        return self.path('obs_data.json')

    @property
    def params_for_id_path(self):
        return self.path('params_for_id.csv')

    @property
    def calibrated_parameters_path(self):
        return self.path('calibrated_parameters.csv')

    @property
    def calibration_path(self):
        return self.path('calibration.json')

    @property
    def has_obs_data(self):
        return os.path.isfile(self.obs_data_path)

    @property
    def obs_data_name(self):
        if not self.has_obs_data:
            return None
        with open(self.obs_data_path) as f:
            return json.load(f).get('obs_data_name')

    @property
    def is_default(self):
        return self.name == self.version.default_instance_name

    @property
    def validation(self):
        '''This instance's validation spec: {baseline: {...}, calibrate: {...}}.'''
        return ((self.version.spec.get('validation') or {}).get(self.name)) or {}

    def parameters(self):
        '''The instance's own rows, then the version's globals, TODO values filled from review proposals.'''
        return self.version._apply_proposals(read_parameters(self.parameters_path, self.version.kinds))

    @property
    def results_dir(self):
        return os.path.join(self.version.results_dir, INSTANCES, self.name)

    def data_files(self):
        '''Every file of the instance except its parameters.'''
        if not os.path.isdir(self.dir):
            return []
        return sorted(f for f in os.listdir(self.dir)
                      if f != os.path.basename(self.parameters_path) and not f.endswith('.omex'))   # .omex: generated


@dataclass
class Version:
    '''One version of a module_type: its config entry, spec and instances. It is what the tests
    run (the old Module + Component in one).'''
    mtype: ModuleType
    name: str
    config: dict               # the version's config entry, libcuflynx key names
    spec: dict
    raw_entry: dict = field(default_factory=dict)

    # --- names ------------------------------------------------------------------------------
    @property
    def dir(self):
        return os.path.join(self.mtype.versions_dir, self.name)

    @property
    def stem(self):
        return f'{self.mtype.name}_{self.name}'

    @property
    def key(self):
        '''<module_type>/<version>: the pytest id and --component value.'''
        return f'{self.mtype.name}/{self.name}'

    @property
    def id(self):
        return safe_id(f'{self.mtype.name}__{self.name}')

    @property
    def label(self):
        return f'{self.mtype.name} ({self.name})'

    @property
    def vessel_type(self):
        return self.mtype.name

    @property
    def BC_type(self):
        return self.name

    @property
    def module_type(self):
        '''The CellML component (libcuflynx's module_type; PhLynx's component_type).'''
        return self.config.get('module_type')

    @property
    def category(self):
        return self.mtype.category

    @property
    def format(self):
        return self.config.get('module_format', 'cellml')

    @property
    def is_supermodule(self):
        return self.format == 'supermodule'

    @property
    def reviewed(self):
        return bool(self.spec.get('reviewed', False))

    # --- files ------------------------------------------------------------------------------
    def path(self, suffix):
        return os.path.join(self.dir, f'{self.stem}_{suffix}')

    @property
    def cellml_path(self):
        return self.path('modules.cellml')

    @property
    def config_path(self):
        return self.path('modules_config.json')

    @property
    def units_path(self):
        return self.path('units.cellml')

    @property
    def spec_path(self):
        '''<module_type>_<version>_tests.yaml: review and record-keeping (TESTS_KEYS).'''
        return self.path('tests.yaml')

    @property
    def verification_config_path(self):
        '''<module_type>_<version>_verification_config.json: what the checks run (VERIFICATION_KEYS).'''
        return self.path('verification_config.json')

    @property
    def html_path(self):
        return os.path.join(self.dir, f'{self.stem}.html')

    @property
    def plots_dir(self):
        return os.path.join(self.dir, 'plots')

    @property
    def results_dir(self):
        return os.path.join(self.dir, 'results')

    @property
    def risk_dir(self):
        return os.path.join(self.dir, 'risk')

    @property
    def instances_dir(self):
        return os.path.join(self.dir, INSTANCES)

    # --- instances --------------------------------------------------------------------------
    @property
    def default_instance_name(self):
        return self.config.get('default_instance') or DEFAULT_INSTANCE

    def instance_names(self):
        if not os.path.isdir(self.instances_dir):
            return []
        names = sorted(d for d in os.listdir(self.instances_dir)
                       if os.path.isdir(os.path.join(self.instances_dir, d)) and not d.startswith('.'))
        # the default instance first
        return sorted(names, key=lambda n: (n != self.default_instance_name, n.lower()))

    def instances(self):
        return [Instance(self, n) for n in self.instance_names()]

    def instance(self, name=None):
        return Instance(self, name or self.default_instance_name)

    @property
    def default_instance(self):
        return self.instance()

    # --- parameters -------------------------------------------------------------------------
    @property
    def kinds(self):
        if self.is_supermodule:
            # no variables_and_units: its instance rows are {var}_{submodule} or the declared globals
            return {g: 'global_constant' for g in self.supermodule_globals}
        return {v[0]: v[3].strip() for v in self.config.get('variables_and_units') or []}

    @functools.cached_property
    def _parameters(self):
        return self.default_instance.parameters()

    def parameters(self):
        '''The default instance's parameters: the version's own constants / BCs, then globals.'''
        return self._parameters

    def _apply_proposals(self, params):
        # A TODO value with a proposed value in the spec's reference_proposals is used (marked
        # proposed) so a version prepared for review can be tested before the proposal is applied.
        for var, prop in (self.spec.get('reference_proposals') or {}).items():
            if not isinstance(prop, dict) or 'value' not in prop:
                continue
            for p in params:
                if p.variable_name == var and p.is_todo:
                    p.value, p.proposed = str(prop['value']), True
        return params

    def unsourced_parameters(self):
        return [p.variable_name for p in self.parameters() if not p.is_sourced]

    def todo_parameters(self):
        return [p.variable_name for p in self.parameters() if p.is_todo]

    def boundary_conditions(self):
        return [p for p in self.parameters() if p.kind == 'boundary_condition']

    def variables(self):
        return [v for v in self.config.get('variables_and_units') or [] if v[3].strip() == 'variable']

    def all_sourced(self):
        return all(p.is_sourced for p in self.parameters())

    # --- supermodules -----------------------------------------------------------------------
    @property
    def submodules(self):
        return self.config.get('submodules') or []

    @property
    def submodule_names(self):
        return [s['name'] for s in self.submodules]

    @property
    def supermodule_globals(self):
        '''The global constants a supermodule's instance may set (spec supermodule.globals).'''
        return list((self.spec.get('supermodule') or {}).get('globals') or SUPERMODULE_DEFAULT_GLOBALS)


def load_module_type(name):
    return ModuleType(name, module_type_dir(name))


@functools.lru_cache(maxsize=None)
def _load_version_cached(mt_name, version):
    return _load_version(mt_name, version)


def _load_version(mt_name, version):
    mtype = load_module_type(mt_name)
    vdir = os.path.join(mtype.versions_dir, version)
    stem = f'{mt_name}_{version}'
    cfg_path = os.path.join(vdir, f'{stem}_modules_config.json')
    with open(cfg_path) as f:
        raw = json.load(f)
    if len(raw) != 1:
        raise ValueError(f'{cfg_path}: a version config has exactly one entry, not {len(raw)}')
    entry = normalise_config_entry(raw[0])
    tests, verification = read_spec_files(vdir, stem)
    spec = merge_spec(tests, verification)
    if isinstance(spec.get('review'), str):
        spec['review_file'] = spec['review']
        spec['review'] = read_review(spec['review'])
    return Version(mtype, version, entry, spec, raw[0])


def read_review(rel):
    '''A shared review (tests.yaml's review: reviews/<name>_review.yaml, relative to the repo root).'''
    with open(os.path.join(REPO_ROOT, rel)) as f:
        return yaml.safe_load(f) or {}


def load_version(mt_name, version):
    '''A version, loaded fresh (specs and parameters may have been edited).'''
    return _load_version(mt_name, version)


def all_versions(selectors=None):
    return [v for n in select_module_types(selectors) for v in load_module_type(n).versions()]


def version_by_key(key):
    '''"<module_type>/<version>" (or the id "<module_type>__<version>") -> Version.'''
    if '/' in key:
        mt, v = key.split('/', 1)
        return load_version(mt, v)
    for n in module_type_names():
        for v in load_module_type(n).version_names():
            if safe_id(f'{n}__{v}') == key:
                return load_version(n, v)
    raise KeyError(key)


def version_index():
    '''(module_type, version) -> Version, for the whole library.'''
    return {(v.vessel_type, v.name): v for v in all_versions()}


@functools.lru_cache(maxsize=None)
def legacy_renames():
    '''(old module_type, old module_subtype) -> (module_type, version) for every pair that changed in
    the move to versions (ion channels, the sympathetic-neuron pieces, the supermodules), from
    tools/restructure_map.yaml. circulatory_autogen's own models still use the old pairs.'''
    path = os.path.join(REPO_ROOT, 'tools', 'restructure_map.yaml')
    if not os.path.isfile(path):
        return {}
    with open(path) as f:
        m = yaml.safe_load(f)
    out = {}
    for mts in m['categories'].values():
        for mt, node in mts.items():
            for v, info in node['versions'].items():
                r = re.match(r'^[\w/]+: \(([^,]+), ([^)]+)\)$', info['from'])
                if r:
                    old = (r[1], r[2])
                else:
                    r = re.match(r'supermodules/(\w+)', info['from'])
                    old = (r[1], 'supermodule')
                if old != (mt, v):
                    out[old] = (mt, v)
    return out


# ----------------------------------------------------------------------------------------------
# compatibility: the tests and tools that took a module name
# ----------------------------------------------------------------------------------------------

def load_module(name):
    '''A ModuleType (the old per-module object); .versions() replaces .components().'''
    return load_module_type(name)
