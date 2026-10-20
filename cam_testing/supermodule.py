"""
Supermodule versions: a version whose config entry has module_format "supermodule" -- a module made
of other modules' versions. Such a version replaces a monolithic one of the same module_type (e.g.
heart version Argus2026_v01, cell/neuron/soma version sympathetic).

  <module_type>_<version>_modules_config.json   one entry: module_format "supermodule", the
                                                submodule_versions it is built from ([{"module_type",
                                                "module_subtype"}]) and "default_parameterisation"
  <module_type>_<version>_module_array.json     its instances (submodules) and their internal
                                                connections: a module array, records naming no
                                                parameterisation (libcuflynx reads it: Version.submodules)
  <module_type>_<version>_parameters.json       its instances' parameter values: for each instance, the
                                                parameterisation of its version they started from and the
                                                values (Version.instance_parameters); calibrate apply
                                                writes every instance of a calibrated version here
  parameterisations/<p>/<p>_parameters.csv      its globals (and any <submodule>/<var> rows, read over
                                                the parameters file)
  <module_type>_<version>_tests.yaml / _verification_config.json   description; supermodule.globals and
                                                supermodule.equivalent (system models it must reproduce)

A model uses one record (an instance of the version), e.g. {"name": "heart", "module_type": "heart",
"module_subtype": "Argus2026_v01", "parameterisation": "default", "per_submodule_inputs": {"ra": ["venous_svc"]},
"per_submodule_outputs": {"aov": ["aortic_root"]}}. libcuflynx expands it: each submodule becomes
<name>__<submodule>, the host modules in per_submodule_inputs / per_submodule_outputs are coupled to
that submodule, and the parameterisation's parameters are renamed to the prefixed names. Precedence:
the model's parameters file > the supermodule's parameterisation > the supermodule's parameters file >
each submodule version's default parameterisation.

The version-level tests (tests/test_modules.py):
  supermodule_structure_test    the module array and parameters file load with libcuflynx (submodules
                                of the listed submodule_versions, internal connections naming siblings,
                                one parameters entry per instance); every submodule is a library version;
                                the parameters file names each instance's variables and a parameterisation
                                it has; the parameterisations' parameters name submodules (or are
                                declared globals)
  supermodule_equivalence_test  each supermodule.equivalent entry: a system model using the version
                                reproduces a system model with the submodules written out
"""
import csv
import os

from cam_testing import checks, harness, library


def supermodules():
    '''Every supermodule version in the library.'''
    return [v for v in library.all_versions() if v.is_supermodule]


def load_supermodule(key):
    '''A supermodule version by "<module_type>/<version>" (or its module_type, if it has one supermodule version).'''
    if '/' in key:
        return library.version_by_key(key)
    found = [v for v in supermodules() if v.vessel_type == key]
    if len(found) != 1:
        raise KeyError(f'{key}: {len(found)} supermodule versions')
    return found[0]


def parameters_file_problems(version):
    '''The supermodule's parameters file: libcuflynx loads it (one entry per instance, its schema), and each
    instance's parameters are variables of its version (not its globals, which the supermodule's
    parameterisation sets), or paths into a supermodule instance, and its parameterisation one the version has.'''
    try:
        instances = version.instance_parameters
    except ValueError as e:
        return [str(e)]
    problems = []
    for sub in version.submodules:
        name = sub['name']
        try:
            sv = library.load_version(sub['module_type'], sub['module_subtype'])
        except (OSError, ValueError, KeyError):
            continue                     # reported as not a library version
        where = f'parameters file, instance {name}'
        started = instances[name]['parameterisation']
        if started != 'default' and started not in sv.parameterisation_names():
            problems.append(f'{where}: {sv.key} has no parameterisation {started} (it has {sv.parameterisation_names()})')
        kinds = {v[0]: v[3] for v in sv.config.get('variables_and_units') or []} if not sv.is_supermodule else {}
        for var in instances[name]['parameters']:
            path, _, leaf = var.rpartition('/')
            if sv.is_supermodule:
                owner = harness.resolve_submodule_path(sv, path) if path else None
                if owner is None:
                    problems.append(f'{where}: parameter {var} is not <its submodule path>/<variable> of {sv.key}')
            elif path:
                problems.append(f'{where}: parameter {var} names a path, but {sv.key} is not a supermodule')
            elif leaf not in kinds:
                problems.append(f'{where}: parameter {var} is not a variable of {sv.key}')
            elif kinds[leaf] == 'global_constant':
                problems.append(f'{where}: parameter {var} is a global of {sv.key}; the supermodule\'s '
                                'parameterisation sets it')
    return problems


def structure_problems(version):
    try:
        # libcuflynx checks the array: names unique, connections between siblings, versions listed
        subs = version.submodules
    except ValueError as e:
        return [str(e)]
    index = library.version_index(include_libraries=True)   # submodules may be extra libraries' versions
    problems = []
    for sub in subs:
        key = (sub['module_type'], sub['module_subtype'])
        if key not in index:
            problems.append(f'{sub["name"]}: ({key[0]}, {key[1]}) is not a version in the module library')
    problems += parameters_file_problems(version)
    declared = version.supermodule_globals
    for inst in version.parameterisations():
        with open(inst.parameters_path) as f:
            for row in csv.DictReader(f):
                var = (row.get('variable_name') or '').strip()
                path = var.rpartition('/')[0]
                if path and harness.resolve_submodule_path(version, path) is None:
                    problems.append(f'parameterisation {inst.name}: parameter {var}: {path} is not a submodule path')
                elif not path and var not in declared:
                    problems.append(f'parameterisation {inst.name}: parameter {var} names no submodule and is not a declared global')
    if version.default_parameterisation_name not in version.parameterisation_names():
        problems.append(f'default_parameterisation {version.default_parameterisation_name} has no '
                        'parameterisations/ directory')
    return problems


def prefixed_output_map(names, instance, submodules):
    '''Output names of a model with the submodules flattened ('ra/u') -> the supermodule model's ('heart__ra/u').'''
    subs = set(submodules)
    out = {}
    for n in names:
        vessel, _, var = n.partition('/')
        out[n] = f'{instance}{harness.LEVEL_SEPARATOR}{vessel}/{var}' if vessel in subs else n
    return out


def equivalence(version, entry, work_dir):
    '''(rows, missing) comparing entry["model"] (uses the version) with entry["reproduces"].'''
    from cam_testing import system as systems
    model = systems.load_system(entry['model'])
    target = systems.load_system(entry['reproduces'])
    solver_info = entry.get('solver_info') or {'rtol': 1e-10, 'atol': 1e-12}
    t_ref, ref = systems.simulate(systems.generate(target, os.path.join(work_dir, 'target')), target.spec, solver_info)
    t_new, new = systems.simulate(systems.generate(model, os.path.join(work_dir, 'model')), model.spec, solver_info)
    # submodule outputs renamed <instance>__<sub> (the instance: the supermodule's record name in the
    # model); an explicit output_map in the spec (for nested
    # supermodules, or variables that moved between components) takes precedence
    output_map = {**prefixed_output_map(ref, entry.get('instance', version.vessel_type), version.submodule_names),
                  **harness.model_output_map(entry.get('output_map') or {})}
    ignore = harness.model_output_keys(entry.get('ignore') or {})
    ref = {k: v for k, v in ref.items() if k not in ignore}
    return systems.compare(ref, new, output_map, float(entry.get('tol', 1e-9)))


def equivalence_result_name(entry):
    return 'supermodule_equivalence_test__' + library.safe_id(entry['model'])


def structure_check(version):
    problems = structure_problems(version)
    if problems:
        r = checks.Result('supermodule_structure_test', checks.FAILED, f'{len(problems)} problems', details=problems)
    else:
        r = checks.Result('supermodule_structure_test', checks.PASSED,
                          f'{len(version.submodules)} submodules, all library versions; '
                          'connections internal; the parameters file names their variables; '
                          'parameters name submodules or globals')
    return checks.save(version, r)


def equivalence_check(version, entry, work_dir):
    tol = float(entry.get('tol', 1e-9))
    name = equivalence_result_name(entry)
    try:
        rows, missing = equivalence(version, entry, work_dir)
    except Exception as e:  # noqa: BLE001
        return checks.save(version, checks.Result(name, checks.FAILED, f'{type(e).__name__}: {str(e)[-1500:]}'))
    bad = [r for r in rows if not r['ok']]
    worst = max(rows, key=lambda r: r['difference']) if rows else None
    metrics = {'model': entry['model'], 'reproduces': entry['reproduces'], 'compared': len(rows),
               'missing': missing, 'tol': tol, 'worst': worst}
    if bad or missing or not rows:
        r = checks.Result(name, checks.FAILED, f'{len(bad)} of {len(rows)} outputs differ; {len(missing)} missing',
                          metrics, details=[f"{x['reference']} -> {x['model']}: {x['difference']:.3g}" for x in bad[:20]]
                          + missing[:20])
    else:
        r = checks.Result(name, checks.PASSED, f'{entry["model"]} reproduces {entry["reproduces"]}: all {len(rows)} '
                          f'outputs within {tol:g} (worst {worst["reference"]} {worst["difference"]:.1e})', metrics)
    return checks.save(version, r)
