"""
Smooths a system model module by module, then as a whole (libcuflynx.utilities.smooth_cellml).

    python -m cam_testing.smooth --module-array sys_module_array.json --parameters sys_parameters.json \\
        [--external-modules-dir DIR ...] [--allow-unchecked] [--output-dir DIR]

1. The original system is generated and run once (its values set the ad hoc harnesses below).
2. Each module of the system is checked alone:
   - a library version (directly, or through library.legacy_renames) by verification_test_smoothing,
     reusing its saved result when that is current (same fingerprint and settings);
   - a module with no library version but a config in an external modules dir by an ad hoc harness:
     the module alone, its constants at the system's values and its inputs held at their mean in the
     system run, each swept over DEFAULT_FACTORS x its value. It is reported as "checked with the
     default sweep", not as a library result;
   - a module found nowhere is "not checked".
   A failing module stops the run with its report. A module not checked stops it too, unless
   --allow-unchecked.
3. The system is probed, run, smoothed and compared with the original; smoothing_report.json merges
   the system's sites with each site's module and that module's result.
"""
import argparse
import contextlib
import glob
import io
import json
import math
import os
import shutil
import sys
import tempfile

import numpy as np

from cam_testing import checks, harness, module_array
from cam_testing.library import legacy_renames, library_dirs, version_index

DEFAULT_FACTORS = [0.5, 0.75, 1.0, 1.5, 2.0]
LIBRARY, AD_HOC, NOT_CHECKED = 'library', 'checked with the default sweep', 'not checked'
DEFAULT_SPEC = {'sim_time': 5.0, 'pre_time': 20.0, 'dt': 0.01, 'solver': 'CVODE_myokit',
                'solver_info': {'rtol': 1e-6, 'atol': 1e-8, 'MaximumStep': 1e-3}}


class SmoothingStopped(Exception):
    """A module failed its smoothing check, or was not checked: the system step is not run."""


def resolve(record, index=None, renames=None):
    """The library Version of a module-array record (directly or through legacy_renames), or None."""
    index = version_index(include_libraries=True) if index is None else index
    renames = legacy_renames() if renames is None else renames
    pair = (record.get('module_type'), record.get('module_subtype'))
    return index.get(pair) or index.get(renames.get(pair, pair))


def external_entry(record, external_dirs):
    """(config entry, its directory) of a record's module in the external modules dirs, or None.
    Both config forms are read (module_type/module_subtype, and vessel_type/BC_type)."""
    want = (record.get('module_type'), record.get('module_subtype'))
    for d in external_dirs:
        for path in sorted(glob.glob(os.path.join(d, '**', '*_config.json'), recursive=True)):
            with open(path) as f:
                try:
                    entries = json.load(f)
                except ValueError:
                    continue
            for e in entries if isinstance(entries, list) else []:
                pair = (e.get('module_type', e.get('vessel_type')), e.get('module_subtype', e.get('BC_type')))
                if (e.get('vessel_type'), e.get('BC_type')) == want or pair == want:
                    return e, os.path.dirname(path)
    return None


def parameter_values(parameters_path):
    """{model name: value} of a system's parameters file (JSON or CSV)."""
    from libcuflynx.utilities import parameters_file as pf
    from libcuflynx.utilities.parameter_names import model_name
    out = {}
    for r in pf.read_parameter_rows(parameters_path):
        try:
            out[model_name(r['variable_name'])] = float(r['value'])
        except (TypeError, ValueError):
            pass
    return out


def _generate(work_dir, prefix, records, parameters_path, external_dirs, quiet=True, builtin=False):
    """Generates a model from module-array records and a parameters file; returns its .cellml path.
    ``builtin``: also libcuflynx's bundled modules (for modules in neither library)."""
    from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture
    resources = os.path.join(work_dir, 'resources')
    os.makedirs(resources, exist_ok=True)
    module_array.write_records(os.path.join(resources, f'{prefix}_module_array.json'), records)
    ext = os.path.splitext(parameters_path)[1]
    shutil.copy(parameters_path, os.path.join(resources, f'{prefix}_parameters{ext}'))
    config = {'file_prefix': prefix, 'input_param_file': f'{prefix}_parameters{ext}', 'model_type': 'cellml',
              'solver': 'CVODE_myokit', 'resources_dir': resources,
              'generated_models_dir': os.path.join(work_dir, 'generated_models'),
              'module_library_dirs': list(library_dirs()) + list(external_dirs), 'use_builtin_modules': builtin,
              'DEBUG': False}
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(log) if quiet else contextlib.nullcontext():
            ok = generate_with_new_architecture(False, config)
    except (Exception, SystemExit) as e:
        raise harness.GenerationFailed(f'{type(e).__name__}: {e}\n{log.getvalue()[-4000:]}') from None
    if not ok:
        raise harness.GenerationFailed(f'generation returned False:\n{log.getvalue()[-4000:]}')
    return os.path.join(work_dir, 'generated_models', prefix, f'{prefix}.cellml')


def library_check(version, settings_override=None, force=False):
    """verification_test_smoothing of a library version: the saved result when current, else a new run."""
    settings = checks.smoothing_settings(version.spec, **(settings_override or {}))
    saved = checks.load(version, 'verification_test_smoothing')
    if not force and saved is not None and \
            saved.metrics.get('fingerprint') == checks.smoothing_fingerprint(version, settings):
        return saved, True
    cm = checks.ComponentModel(version)
    return checks.verification_test_smoothing(cm, **(settings_override or {})), False


def ad_hoc_check(record, entry, external_dirs, system_values, system_means, spec, settings, work_dir, out_dir,
                 exempt=()):
    """The ad hoc harness of a module with no library version: the module alone, generated as the vessel
    "mod", its constants at the system's values and its boundary conditions held at their mean in the
    system run, each swept over DEFAULT_FACTORS x its value. Returns (status, message, metrics)."""
    name = record['name']
    rows, nominal, points = [], {}, [{'label': 'nominal', 'params': {}}]
    outputs, missing = [], []
    for var, units, _access, kind, *_ in entry.get('variables_and_units') or []:
        if kind == 'variable':
            outputs.append(var)
            continue
        if kind == 'global_constant':
            pname, path, value = var, f'global/{var}', system_values.get(var)
        elif kind == 'constant':
            pname, path, value = f'{var}_{harness.VESSEL}', f'{harness.VESSEL}/{var}', system_values.get(f'{var}_{name}')
        elif kind == 'boundary_condition':
            pname, path = f'{var}_{harness.VESSEL}', f'{harness.VESSEL}/{var}'
            value = system_means.get(f'{name}/{var}', system_values.get(f'{var}_{name}'))
        else:
            continue
        if value is None:
            missing.append(var)
            continue
        rows.append((path, units, value))
        nominal['parameters/' + pname] = float(value)
        if value != 0:
            points += [{'label': f'{var}={value * f:.4g}', 'params': {'parameters/' + pname: value * f}}
                       for f in DEFAULT_FACTORS if f != 1.0]
    if missing:
        return checks.FAILED, f'no system value for {", ".join(missing)}', {'missing': missing}
    from libcuflynx.utilities import parameters_file as pf
    os.makedirs(work_dir, exist_ok=True)
    params_path = os.path.join(work_dir, f'{name}_parameters.json')
    doc = pf.new_document()
    for path, units, value in rows:
        pf.set_leaf(doc, pf.name_path(path), {'value': value, 'units': units, 'data_reference': 'system value'})
    pf.write_parameters_file(params_path, doc)
    records = [{'name': harness.VESSEL, 'module_type': record['module_type'],
                'module_subtype': record['module_subtype'], 'inp_instances': [], 'out_instances': []}]
    model_path = _generate(work_dir, f'adhoc_{name}', records, params_path, external_dirs)
    # the system's exemptions of this instance ("<name>_module/<var>"), as the module alone names them
    own = []
    for e in exempt:
        comp, _, var = e['variable'].rpartition('/')
        if comp in ('', f'{name}_module'):
            own.append({'variable': f'{harness.VESSEL}_module/{var}' if comp else var, 'reason': e.get('reason', '')})
    settings = dict(settings, exempt_entries=own)
    status, message, metrics, _figs, details = checks.compare_smoothed(
        model_path, spec, outputs or checks.harness.state_outputs(harness.simulation_helper(model_path, spec)),
        points, nominal, settings, out_dir, work_dir, f'adhoc_{name}')
    metrics['details'] = details
    return status, message, metrics


def system_means(model_path, spec, names):
    """{'instance/var': mean over the logged run} of a generated system model. An input connected to
    another module's variable is read from the variable that defines it (the solver keeps one name
    per connected set)."""
    import xml.etree.ElementTree as ET
    from cam_testing.system import simulate
    from libcuflynx.utilities import smooth_cellml as sc
    root = ET.parse(checks.flat_model_path(model_path)).getroot()
    graph = sc._ModelGraph(root)
    members = {}
    for comp in root.findall(sc._C + 'component'):
        for v in comp.findall(sc._C + 'variable'):
            members.setdefault(graph.find((comp.get('name'), v.get('name'))), []).append((comp.get('name'), v.get('name')))
    candidates = {}
    for name in names:
        instance, _, var = name.partition('/')
        for comp, v in members.get(graph.find((f'{instance}_module', var)), []):
            short = comp[:-len('_module')] if comp.endswith('_module') else comp
            candidates.setdefault(name, []).append(f'{short}/{v}')
    wanted = sorted({c for cs in candidates.values() for c in cs})
    _, out = simulate(model_path, spec, spec.get('solver_info'), names=wanted)
    means = {k: float(np.mean(v)) for k, v in out.items() if np.all(np.isfinite(v))}
    return {name: next(means[c] for c in cs if c in means) for name, cs in candidates.items()
            if any(c in means for c in cs)}


def run(module_array_path, parameters_path, external_dirs=(), allow_unchecked=False, output_dir=None,
        spec=None, settings=None, force=False, work_dir=None, exempt=None, modules_only=False, system_model=None):
    """The pipeline above. Returns the merged report; raises SmoothingStopped before the system step
    when a module fails (or is not checked, unless allow_unchecked). ``exempt``: the system's own
    exemptions ([{"variable": "[component/]glob", "reason"}]), merged with each library version's
    smoothing.exempt (as "<instance>_module/<variable>"). ``modules_only``: stop after the module checks.
    ``system_model``: a generated model of the system to use instead of generating it from this library
    (e.g. one generated with libcuflynx's built-in modules); its flat CellML must sit next to it.
    The module table is written to module_table.csv."""
    exempt = list(exempt or [])
    spec = dict(DEFAULT_SPEC, **(spec or {}))
    settings = dict(checks.DEFAULT_SMOOTHING, **(settings or {}))
    prefix = os.path.basename(module_array_path).split('_module_array')[0].split('_vessel_array')[0]
    output_dir = output_dir or os.path.join(os.path.dirname(os.path.abspath(module_array_path)), 'smoothing')
    work_dir = work_dir or tempfile.mkdtemp(prefix=f'cam_smooth_{prefix}_')
    os.makedirs(output_dir, exist_ok=True)
    records = module_array.read_records(module_array_path)
    index, renames = version_index(include_libraries=True), legacy_renames()

    # the original system, mapped onto this library's versions
    mapped = []
    for rec in records:
        version = resolve(rec, index, renames)
        mapped.append(dict(rec, module_type=version.vessel_type, module_subtype=version.name) if version else rec)
    unresolved = [r for r in records if resolve(r, index, renames) is None and not external_entry(r, external_dirs)]
    if system_model is None:
        system_model = _generate(os.path.join(work_dir, 'system'), prefix, mapped, parameters_path, external_dirs,
                                 builtin=bool(unresolved))
    values = parameter_values(parameters_path)

    modules, needs_means = [], []
    for rec in records:
        version = resolve(rec, index, renames)
        found = None if version else external_entry(rec, external_dirs)
        entry = found[0] if found else None
        modules.append((rec, version, entry))
        if entry:
            needs_means += [f'{rec["name"]}/{v[0]}' for v in entry.get('variables_and_units') or []
                            if len(v) > 3 and v[3] == 'boundary_condition']
    means = system_means(system_model, spec, needs_means) if needs_means else {}

    results, checked = {}, {}
    for rec, version, entry in modules:
        key = f'{rec.get("module_type")}/{rec.get("module_subtype")}'
        if key in checked:
            results[rec['name']] = checked[key]
            continue
        if version is not None:
            result, reused = library_check(version, force=force)
            res = {'kind': LIBRARY, 'version': version.key, 'status': result.status, 'message': result.message,
                   'reused': reused, 'worst_sites': result.metrics.get('attribution', [])[:3]}
        elif entry is not None:
            status, message, metrics = ad_hoc_check(
                rec, entry, external_dirs, values, means, spec, settings,
                os.path.join(work_dir, 'adhoc', rec['name']), os.path.join(output_dir, 'adhoc', rec['name']), exempt)
            res = {'kind': AD_HOC, 'status': status, 'message': message,
                   'worst_sites': metrics.get('attribution', [])[:3]}
        else:
            res = {'kind': NOT_CHECKED, 'status': checks.SKIPPED,
                   'message': 'no library version and no config in the external modules dirs'}
        checked[key] = results[rec['name']] = dict(res, module=key)

    # each library version's exemptions, for its instances in the system
    system_exempt = list(exempt)
    for rec, version, _entry in modules:
        if version is not None:
            system_exempt += checks.exempt_entries((version.spec.get('smoothing') or {}).get('exempt'),
                                                   component=f'{rec["name"]}_module')
    report = {'module_array': os.path.abspath(module_array_path), 'exempt': system_exempt, 'parameters': os.path.abspath(parameters_path),
              'settings': settings, 'spec': spec, 'modules': results}
    with contextlib.redirect_stdout(io.StringIO()) as ranges_out:
        try:
            from cam_testing import ranges
            report['outside_tested_ranges'] = ranges.check(module_array_path, parameters_path, out=sys.stdout)
        except Exception as e:      # an informational check; the smoothing does not depend on it
            report['outside_tested_ranges'] = f'not checked ({type(e).__name__}: {e})'
    report['ranges_check'] = ranges_out.getvalue()
    report_path = os.path.join(output_dir, 'smoothing_report.json')
    write_module_table(os.path.join(output_dir, 'module_table.csv'), results)
    if modules_only:
        report['system'] = {'status': 'not run', 'reason': 'modules only'}
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2, default=checks._json_default)
        return report

    failed = {n: r for n, r in results.items() if r['status'] == checks.FAILED}
    unchecked = {n: r for n, r in results.items() if r['kind'] == NOT_CHECKED}
    if failed or (unchecked and not allow_unchecked):
        report['system'] = {'status': 'not run', 'reason': ('modules failed: ' + ', '.join(sorted(failed))) if failed
                            else 'modules not checked (use --allow-unchecked): ' + ', '.join(sorted(unchecked))}
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2, default=checks._json_default)
        raise SmoothingStopped(report['system']['reason'] + '\n' + '\n'.join(
            f'  {n}: {r["module"]}: {r["message"]}' for n, r in {**failed, **unchecked}.items()))

    # the system
    from libcuflynx.solver_wrappers import get_simulation_helper
    with contextlib.redirect_stdout(io.StringIO()):
        helper = get_simulation_helper(model_path=system_model, solver='CVODE_myokit', model_type='cellml',
                                       dt=spec['dt'], sim_time=spec['sim_time'], solver_info=spec['solver_info'],
                                       pre_time=spec['pre_time'])
    outputs = harness.state_outputs(helper)
    status, message, metrics, _figs, details = checks.compare_smoothed(
        system_model, spec, outputs, [{'label': 'nominal', 'params': {}}], {}, dict(settings, exempt_entries=system_exempt),
        output_dir, os.path.join(work_dir, 'system'), prefix)
    with open(metrics.get('report', os.path.join(output_dir, 'smoothing_report.json'))) as f:
        system_report = json.load(f)
    by_component = {f'{n}_module': r for n, r in results.items()}
    for site in system_report.get('sites', []):
        module = by_component.get(site['component'])
        site['module'] = None if module is None else {k: module[k] for k in ('module', 'kind', 'status')}
    report['system'] = {'status': status, 'message': message, 'details': details,
                        'metrics': {k: v for k, v in metrics.items() if k != 'points'},
                        'sites': system_report.get('sites', []), 'not_smoothed': system_report.get('not_smoothed', [])}
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2, default=checks._json_default)
    return report


def write_module_table(path, results):
    """One row per module instance: its module, how it was checked, its status and message."""
    import csv
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['instance', 'module', 'checked', 'status', 'reused', 'worst_sites', 'message'])
        for name, r in results.items():
            worst = '; '.join(f'{a.get("variable")} {a.get("excess", a.get("l2", 0)):.3g}' for a in r.get('worst_sites') or [])
            w.writerow([name, r['module'], r['kind'], r['status'], r.get('reused', ''), worst, r['message']])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--module-array', required=True)
    parser.add_argument('--parameters', required=True)
    parser.add_argument('--external-modules-dir', action='append', default=[],
                        help='a directory of module configs and CellML outside the library (repeatable)')
    parser.add_argument('--exempt', help='JSON [{"variable": "[component/]glob", "reason": ...}] for the system')
    parser.add_argument('--allow-unchecked', action='store_true',
                        help='run the system step although some modules could not be checked')
    parser.add_argument('--output-dir')
    parser.add_argument('--sim-time', type=float)
    parser.add_argument('--pre-time', type=float)
    parser.add_argument('--dt', type=float)
    parser.add_argument('--rel', type=float, help='eps relative to each range for the system and ad hoc modules')
    parser.add_argument('--settings', help='JSON smoothing settings for the system and ad hoc modules, as the spec\'s '
                                           'smoothing key (e.g. {"rel": {"maxmin": 1e-3, "variables": {"v_leak": 1e-5}}})')
    parser.add_argument('--tol', type=float)
    parser.add_argument('--force', action='store_true', help='rerun the library checks even when current')
    parser.add_argument('--modules-only', action='store_true', help='stop after the module checks')
    parser.add_argument('--system-model', help='a generated model of the system to use instead of generating one')
    args = parser.parse_args(argv)
    spec = {k: v for k, v in (('sim_time', args.sim_time), ('pre_time', args.pre_time), ('dt', args.dt))
            if v is not None}
    settings = {}
    if args.settings:
        with open(args.settings) as f:
            settings.update(json.load(f))
    settings.update({k: v for k, v in (('rel', args.rel), ('tol', args.tol)) if v is not None})
    try:
        exempt = None
        if args.exempt:
            with open(args.exempt) as f:
                exempt = json.load(f)
        report = run(args.module_array, args.parameters, args.external_modules_dir, args.allow_unchecked,
                     args.output_dir, spec, settings, args.force, exempt=exempt, modules_only=args.modules_only,
                     system_model=args.system_model)
    except SmoothingStopped as e:
        print(f'stopped before the system step: {e}')
        return 2
    for name, r in report['modules'].items():
        print(f'{name:24s} {r["module"]:40s} {r["kind"]:30s} {r["status"]}: {r["message"][:120]}')
    if args.modules_only:
        return 0
    print(f'system: {report["system"]["status"]}: {report["system"]["message"]}')
    return 0 if report['system']['status'] == checks.PASSED else 1


if __name__ == '__main__':
    sys.exit(main())
