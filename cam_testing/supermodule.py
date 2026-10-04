"""
Supermodule versions: a version whose config entry has module_format "supermodule" -- a module made
of other modules' versions. Such a version replaces a monolithic one of the same module_type (e.g.
heart version Argus2026_v01, cell/neuron/soma version sympathetic).

  <module_type>_<version>_modules_config.json   one entry: module_format "supermodule", the
                                                submodules (a vessel array of library versions, each
                                                with its "instance"; internal connections only) and
                                                "default_instance"
  instances/<instance>/<instance>_parameters.csv  its parameters, {var}_{submodule} or globals
  <module_type>_<version>_tests.yaml / _verification_config.json   description; supermodule.globals and
                                                supermodule.equivalent (system models it must reproduce)

A model uses one record, e.g. {"name": "heart", "module_type": "heart", "module_subtype":
"Argus2026_v01", "instance": "default", "per_submodule_inputs": {"ra": ["venous_svc"]},
"per_submodule_outputs": {"aov": ["aortic_root"]}}. libcuflynx expands it: each submodule becomes
<name>_<submodule>, the host modules in per_submodule_inputs / per_submodule_outputs are coupled to
that submodule, and the instance's parameters are renamed to the prefixed names. Precedence: the
model's parameters file > the supermodule's instance > the submodules' instances.

The version-level tests (tests/test_modules.py):
  supermodule_structure_test    every submodule names a library (module_type, version) and an instance
                                it has; internal connections name sibling submodules; the instance's
                                parameters name submodules (or are declared globals)
  supermodule_equivalence_test  each supermodule.equivalent entry: a system model using the version
                                reproduces a system model with the submodules written out
"""
import csv
import dataclasses
import os

from cam_testing import checks, library


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


def structure_problems(version):
    index = library.version_index()
    names = set(version.submodule_names)
    problems = []
    if len(names) != len(version.submodules):
        problems.append('duplicate submodule names')
    template = bool(version.config.get('template'))
    for sub in version.submodules:
        if template and not sub.get('module_subtype'):
            # an open slot: its module_type and each listed choice must be in the library
            fitting = [v for (mt, v) in index if mt == sub['module_type']]
            if not fitting:
                problems.append(f'{sub["name"]}: no module_type {sub["module_type"]} in the module library')
            for c in sub.get('choices') or []:
                if c not in fitting:
                    problems.append(f'{sub["name"]}: choice {c} is not a version of {sub["module_type"]}')
            continue
        key = (sub['module_type'], sub['module_subtype'])
        target = index.get(key)
        if target is None:
            problems.append(f'{sub["name"]}: ({key[0]}, {key[1]}) is not a version in the module library')
        elif sub.get('instance') and sub['instance'] not in target.instance_names():
            problems.append(f'{sub["name"]}: {key[0]}/{key[1]} has no instance {sub["instance"]} '
                            f'(it has {target.instance_names()})')
        for n in sub.get('inp_instances', []) + sub.get('out_instances', []):
            if n not in names:
                problems.append(f'{sub["name"]}: {n} is not a submodule (external connections belong in the host\'s '
                                'per_submodule_inputs / per_submodule_outputs)')
    declared = version.supermodule_globals
    shared = set(version.config.get('shared_parameters') or [])
    for key, routes in (version.config.get('routes') or {}).items():
        for port_type, subs in routes.items():
            for n in [subs] if isinstance(subs, str) else subs:
                if n not in names:
                    problems.append(f'routes.{key}.{port_type}: {n} is not a submodule')
    for inst in version.instances():
        with open(inst.parameters_path) as f:
            for row in csv.DictReader(f):
                var = (row.get('variable_name') or '').strip()
                owner = next((n for n in sorted(names, key=len, reverse=True) if var.endswith('_' + n)), None)
                if owner is None and var not in declared and var not in shared:
                    problems.append(f'instance {inst.name}: parameter {var} names no submodule and is not a declared global')
    if version.default_instance_name not in version.instance_names():
        problems.append(f'default_instance {version.default_instance_name} has no instances/ directory')
    return problems


def prefixed_output_map(names, instance, submodules):
    '''Output names of a model with the submodules flattened ('ra/u') -> the supermodule model's ('heart_ra/u').'''
    subs = set(submodules)
    out = {}
    for n in names:
        vessel, _, var = n.partition('/')
        out[n] = f'{instance}_{vessel}/{var}' if vessel in subs else n
    return out


def equivalence(version, entry, work_dir):
    '''(rows, missing) comparing entry["model"] (uses the version) with entry["reproduces"].'''
    from cam_testing import system as systems
    model = systems.load_system(entry['model'])
    target = systems.load_system(entry['reproduces'])
    solver_info = entry.get('solver_info') or {'rtol': 1e-10, 'atol': 1e-12}
    t_ref, ref = systems.simulate(systems.generate(target, os.path.join(work_dir, 'target')), target.spec, solver_info)
    t_new, new = systems.simulate(systems.generate(model, os.path.join(work_dir, 'model')), model.spec, solver_info)
    # submodule outputs renamed <instance>_<sub>; an explicit output_map in the spec (for nested
    # supermodules, or variables that moved between components) takes precedence
    output_map = {**prefixed_output_map(ref, entry.get('instance', version.vessel_type), version.submodule_names),
                  **(entry.get('output_map') or {})}
    ignore = entry.get('ignore') or {}
    # vessels replaced by supermodules (e.g. a model rebuilt with lumped vessels): only their mapped
    # outputs are compared, the rest have no single counterpart
    replaced = tuple(f'{v}/' for v in entry.get('replaced_vessels') or [])
    mapped = entry.get('output_map') or {}
    ref = {k: v for k, v in ref.items() if k not in ignore and (not k.startswith(replaced) or k in mapped)}
    return systems.compare(ref, new, output_map, float(entry.get('tol', 1e-9)), wrapped=entry.get('wrapped'),
                           reciprocal=entry.get('reciprocal'))


def equivalence_result_name(entry):
    return 'supermodule_equivalence_test__' + library.safe_id(entry['model'])


def structure_check(version):
    problems = structure_problems(version)
    if problems:
        r = checks.Result('supermodule_structure_test', checks.FAILED, f'{len(problems)} problems', details=problems)
    else:
        r = checks.Result('supermodule_structure_test', checks.PASSED,
                          f'{len(version.submodules)} submodules, all library versions with their instances; '
                          'connections internal; parameters name submodules or globals')
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


# ---- a supermodule version that reproduces a monolithic version ----------------------------------
#
# supermodule.equivalent_version in the spec:
#   version                 <module_type>/<version> of the monolithic version
#   monolithic_parameters   {variable: value}: the monolithic version's parameters for the comparison
#                           (its default instance's, with its TODOs and boundary conditions filled to
#                           match the supermodule's default instance)
#   output_map              {monolithic output: "<submodule>/<variable>"}
#   tol, solver_info
# Both are generated alone (the vessel "mod") and simulated over the supermodule's sim_time.

def version_equivalence(version, entry, work_dir):
    from cam_testing import harness, system as systems
    mt, v = entry['version'].split('/', 1)
    mono = library.load_version(mt, v)
    # the monolithic version's parameters, with the entry's values (and any the instance lacks)
    values = dict(entry.get('monolithic_parameters') or {})
    config = {r[0]: r for r in mono.config.get('variables_and_units') or []}
    params = []
    for p in mono.parameters():
        if p.variable_name in values:
            p = dataclasses.replace(p, value=str(values.pop(p.variable_name)))
        params.append(p)
    for name, value in values.items():
        units, kind = (config[name][1], config[name][3]) if name in config else ('dimensionless', 'constant')
        params.append(library.Parameter(name, units, str(value), 'supermodule.equivalent_version', 'no', kind))
    solver_info = entry.get('solver_info') or {'rtol': 1e-10, 'atol': 1e-12}
    spec = {'sim_time': version.spec.get('sim_time', 1.0), 'dt': version.spec.get('dt', 0.001), 'pre_time': 0.0}
    output_map = {f'{harness.VESSEL}/{k}': f'{harness.VESSEL}_{target}' for k, target in entry['output_map'].items()}
    mono_path = harness.generate(mono, os.path.join(work_dir, 'monolithic'), parameters=params)
    # a monolithic version tested in a network (spec harness: e.g. a flow source upstream) is compared
    # with the supermodule in the same network, in its place
    network = mono.spec.get('harness')
    if network and network.get('vessel_array'):
        rows = [[r[0], version.BC_type, version.vessel_type, *r[3:]] if r[0] == harness.VESSEL else r
                for r in network['vessel_array']]
        version = dataclasses.replace(version, spec={**version.spec, 'harness': {**network, 'vessel_array': rows}})
    super_path = harness.generate(version, os.path.join(work_dir, 'supermodule'))
    _, ref = systems.simulate(mono_path, spec, solver_info, names=list(output_map))
    _, new = systems.simulate(super_path, spec, solver_info, names=list(output_map.values()))
    missing = [k for k in output_map if k not in ref]
    rows, missing_new = systems.compare(ref, new, output_map, float(entry.get('tol', 1e-6)))
    return rows, missing + missing_new


def version_equivalence_check(version, entry, work_dir):
    tol = float(entry.get('tol', 1e-6))
    name = 'supermodule_version_equivalence_test'
    try:
        rows, missing = version_equivalence(version, entry, work_dir)
    except Exception as e:  # noqa: BLE001
        return checks.save(version, checks.Result(name, checks.FAILED, f'{type(e).__name__}: {str(e)[-1500:]}'))
    bad = [r for r in rows if not r['ok']]
    worst = max(rows, key=lambda r: r['difference']) if rows else None
    metrics = {'reproduces': entry['version'], 'compared': len(rows), 'missing': missing, 'tol': tol, 'worst': worst}
    if bad or missing or not rows:
        r = checks.Result(name, checks.FAILED, f'{len(bad)} of {len(rows)} outputs differ; {len(missing)} missing',
                          metrics, details=[f"{x['reference']} -> {x['model']}: {x['difference']:.3g}" for x in bad]
                          + missing)
    else:
        r = checks.Result(name, checks.PASSED, f'reproduces {entry["version"]}: all {len(rows)} outputs within '
                          f'{tol:g} (worst {worst["reference"]} {worst["difference"]:.1e})', metrics)
    return checks.save(version, r)
