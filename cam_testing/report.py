"""
HTML reports at two levels, and a site index:

    <module_type>/<module_type>.html                        every version, linking to its page
    <module_type>/versions/<v>/<module_type>_<v>.html       the version: equations, variables,
                                                            tests, references, and its instances
                                                            with their validation / calibration

    python -m cam_testing.report                  # every module_type and its versions
    python -m cam_testing.report --module heart   # a module_type and those nested in it (or a category: --module cell)
    python -m cam_testing.report --site           # also assemble site/ as GitHub Pages serves it

Pages read the result JSON and plots the tests wrote, so a report shows the last local (or
CI) test run. Plots are referenced relatively (plots/...), so a page works opened from disk.
"""
import argparse
import datetime
import json
import os
import shutil

from jinja2 import Environment, FileSystemLoader, select_autoescape

from cam_testing import bib, checks, phlynx, ranges, risk
from cam_testing import omex as omex_mod
from cam_testing import supermodule as sm
from cam_testing.library import LICENCES, REPO_ROOT, load_module_type, module_type_names, select_module_types
from cam_testing.mathml import component_equation_targets, component_equations, component_variables

TEMPLATES = os.path.join(os.path.dirname(__file__), 'templates')
SITE_DIR = os.path.join(REPO_ROOT, 'site')

TEST_TITLES = {
    'run_test': 'Run',
    'verification_test_invariants': 'Verification: invariants',
    'verification_test_BC': 'Verification: boundary conditions',
    'verification_test_timestep': 'Verification: timestep convergence',
    'stability_test': 'Stability: solvers & settings',
    'validation_test_baseline': 'Validation: baseline data',
    'validation_test_calibrate': 'Validation: calibrate & predict',
    'version_calibration': 'Calibration (all instances)',
    'phlynx_export_test': 'PhLynx: build & export',
    'cuflynx_simulate_test': 'CUFLynx: import & simulate',
    'phlynx_equivalence_test': 'PhLynx → CUFLynx vs libcuflynx',
    'supermodule_structure_test': 'Supermodule: structure',
    'supermodule_equivalence_test': 'Supermodule: reproduces',
}
TEST_SHORT = {
    'run_test': 'Run', 'verification_test_invariants': 'Invariants', 'verification_test_BC': 'BC sweep', 'verification_test_timestep': 'Timestep',
    'stability_test': 'Stability', 'validation_test_baseline': 'Baseline', 'validation_test_calibrate': 'Calibrate', 'version_calibration': 'Calibration',
    'phlynx_export_test': 'PhLynx export', 'cuflynx_simulate_test': 'CUFLynx simulate', 'phlynx_equivalence_test': 'PhLynx equivalence',
    'supermodule_structure_test': 'Supermodule structure', 'supermodule_equivalence_test': 'Supermodule reproduces',
}
TEST_ABOUT = {
    'run_test': 'Generates the version alone with libcuflynx (every boundary condition becomes a '
                'parameter, valued from the default instance; a supermodule is expanded into its submodules, '
                'and the boundary conditions no internal connection closes become parameters), simulates it, '
                'and checks every output is finite.',
    'verification_test_invariants': 'Checks the simulation against what the version is supposed to do: '
                                    'the invariants in its spec (exact solutions, conservation laws, bounds, '
                                    'delays) evaluated at the run parameters.',
    'verification_test_BC': 'Sweeps each boundary condition and constant (bc_sweep.sweep: all, a component '
                            'version\'s default; a supermodule version sweeps its global constants and the boundary '
                            'conditions no internal connection closes, since each submodule\'s own parameters are '
                            'swept in that submodule version\'s tests) over a range and checks every run completes, '
                            'stays finite and keeps the invariants. Values outside a parameter\'s valid range (bc_sweep.bounds, '
                            'a gate initial value outside [0, 1], a negative physical quantity, or a broken '
                            'bc_sweep.constraints relation) are skipped. A failed run is rerun once with a tight '
                            'reference solver: if that fails too it is a parameter failure (the model breaks '
                            'there); if it passes, a numerical failure, and the cheapest CVODE setting that '
                            'passes every point is recommended. The test passes when there are no parameter '
                            'failures and every numerical failure passes with the recommended settings.',
    'verification_test_timestep': 'Integrates the generated right-hand side with a fixed-step scheme at '
                                  'successively halved steps. Differences between successive solutions must '
                                  'shrink at the scheme\'s order, and the finest solution must agree with '
                                  'libcuflynx\'s CVODE run at tight tolerances.',
    'stability_test': 'Runs the version with a matrix of solvers, tolerances and timesteps. A '
                      'configuration works when it finishes with finite outputs within the tolerance of a '
                      'tight-tolerance reference. Declared-supported configurations must work.',
    'validation_test_baseline': 'Compares the model, at the instance\'s parameters, with published or '
                                'experimental baseline data.',
    'validation_test_calibrate': 'Calibrates parameters to the instance\'s obs_data with libcuflynx parameter '
                                 'identification, then checks predictions against held-out data; writes the '
                                 'instance\'s calibrated parameters. Not applicable to an instance without '
                                 'obs_data ("no calibration data in this instance").',
    'version_calibration': 'The version\'s calibration over its instances: failed when no instance has '
                           'calibration data (obs_data), otherwise failed if any instance\'s calibration '
                           'failed and passed when all of them pass. Each instance\'s result is in Instances. '
                           'A version with no calibration data of its own that is a submodule of supermodule '
                           'versions is calibrated as part of them: "Pass in super" when any of those '
                           'supermodules passes calibration (plainly or, transitively, in its own supermodules), '
                           '"Fail in super" when none does (tests.yaml calibration_in_supermodule: false opts out).',
    'phlynx_export_test': 'Builds the version\'s test network in PhLynx (its own code, loaded with this '
                          'library\'s modules and parameters), checks every connection was made, and exports '
                          'the .omex PhLynx sends to CUFLynx. PhLynx has no supermodule support, so a supermodule '
                          'version fails this test (a recorded known issue).',
    'cuflynx_simulate_test': 'Imports that .omex into a released CUFLynx through its API and simulates it; '
                             'every output must be finite.',
    'phlynx_equivalence_test': 'Compares CUFLynx\'s simulation of the PhLynx-built model with libcuflynx\'s '
                               'model of the same network (outputs matched by instance and variable, '
                               'normalised difference within 1e-6).',
    'supermodule_structure_test': 'Every submodule is a library version with an existing instance, internal '
                                  'connections name sibling submodules, and the instance parameters name '
                                  'submodules or declared globals.',
    'supermodule_equivalence_test': 'A system model using this supermodule version reproduces the system model '
                                    'with its submodules written out, output for output.',
}
# The test columns of every version, component or supermodule, in this order: verification, the
# version-level calibration, the PhLynx -> CUFLynx pipeline, then the supermodule tests. A column
# that doesn't apply to a version is N/A with the reason (not_applicable_reason); one that applies
# but has no result is "Not run". The validation tests are shown per instance.
SUPERMODULE_TESTS = ['supermodule_structure_test', 'supermodule_equivalence_test']
REPORT_TESTS = checks.VERSION_TESTS + ['version_calibration'] + phlynx.PIPELINE_TESTS + SUPERMODULE_TESTS
INSTANCE_TESTS = checks.INSTANCE_TESTS
NOT_A_SUPERMODULE = 'not a supermodule'
NO_EQUIVALENT = ('no system model to reproduce: the spec lists no supermodule.equivalent entry (the version is '
                 'checked through the models that use it)')
STATUS_LABEL = {'passed': 'Passed', 'failed': 'Failed', 'skipped': 'Skipped', 'pending': 'Pending',
                'not_applicable': 'N/A', checks.PASSED_IN_SUPER: 'Pass in super', checks.FAILED_IN_SUPER: 'Fail in super',
                None: 'Not run'}
# in the status counts, a pass / fail in a supermodule counts as a pass / fail
COUNT_AS = {checks.PASSED_IN_SUPER: checks.PASSED, checks.FAILED_IN_SUPER: checks.FAILED}


def _fmt(value):
    if isinstance(value, float):
        if value == 0:
            return '0'
        if abs(value) >= 1e4 or abs(value) < 1e-3:
            return f'{value:.3e}'
        return f'{value:.4g}'
    if isinstance(value, list):
        return ', '.join(_fmt(v) for v in value)
    if isinstance(value, dict):
        return '; '.join(f'{k}: {_fmt(v)}' for k, v in value.items())
    if value is None:
        return '—'
    return str(value)


def _ports(entry):
    out = []
    for side in ('entrance_ports', 'exit_ports', 'general_ports'):
        for port in entry.get(side, []):
            out.append({'side': side.split('_')[0], 'type': port['port_type'],
                        'variables': ', '.join(port['variables']), 'multi': port.get('multi_port', '')})
    return out


def _cellml_equations(version):
    equations, unsupported, cellml_vars = [], set(), {}
    if version.format == 'cellml' and os.path.isfile(version.cellml_path):
        equations, unsupported = component_equations(version.cellml_path, version.module_type)
        cellml_vars = {v[0]: v for v in component_variables(version.cellml_path, version.module_type)}
    return equations, unsupported, cellml_vars


# a baseline whose parameters were fitted to the same data (spec kind: fit_check)
FIT_CHECK_TITLE = 'Goodness of fit: published parameters against the data they were fitted to'
FIT_CHECK_ABOUT = ('The parameters were fitted to this data, so this shows how well that fit reproduces it; it '
                   'is not an independent validation.')


def _test_entry(version, test, r, spec_block=None, title=None):
    base = test.split('__')[0]
    return {
        'key': test, 'title': title or TEST_TITLES.get(base, test), 'short': TEST_SHORT.get(base, test),
        'about': TEST_ABOUT.get(base, ''),
        'status': r.status if r else None, 'status_label': STATUS_LABEL.get(r.status if r else None, 'Not run'),
        'message': r.message if r else 'This test has not been run yet.',
        'metrics': r.metrics if r else {}, 'plots': r.plots if r else [],
        'details': r.details if r else [], 'timestamp': r.timestamp if r else '',
        'known_issue': (version.spec.get('expected_failures') or {}).get(test),
        'proposed_data': (spec_block or {}).get('status') == checks.PROPOSED,
        'links': [],          # related version pages (version_calibration in a supermodule)
    }


def instance_context(version, inst):
    '''One instance: its files, and its validation (baseline) and calibration results.'''
    tests = []
    for test in INSTANCE_TESTS:
        kind = test.rsplit('_', 1)[1]
        v = inst.validation.get(kind) or {}
        r = checks.load(version, test, inst)
        if r is None:
            # not run (e.g. slow tests excluded): what a run would record without data
            if kind == 'calibrate' and not inst.has_obs_data:
                r = checks.Result(test, checks.NOT_APPLICABLE, checks.NO_CALIBRATION_DATA)
            elif not v:
                r = checks.Result(test, checks.NOT_APPLICABLE,
                                  checks.NO_BASELINE_DATA if kind == 'baseline' else checks.NO_CALIBRATION_DATA)
            elif v.get('status', checks.PENDING) in (checks.PENDING, checks.SKIPPED, checks.NOT_APPLICABLE):
                r = checks.Result(test, v.get('status', checks.PENDING), v.get('reason', 'no validation data chosen yet'))
        entry = _test_entry(version, test, r, v, title=FIT_CHECK_TITLE if v.get('kind') == 'fit_check' else None)
        if v.get('kind') == 'fit_check':
            entry['short'] = 'Fit check'
            entry['about'] = v.get('note') or FIT_CHECK_ABOUT
        entry['anchor'] = f'{inst.name}--{test}'
        # the publication figures the data were extracted from, shown beside the test's plots
        entry['source_figures'] = inst.source_figures()
        entry['source_figures_missing'] = inst.needs_source_figures() and not entry['source_figures']
        tests.append(entry)
    calibration = None
    if os.path.isfile(inst.calibration_path):
        with open(inst.calibration_path) as f:
            calibration = json.load(f)
    from cam_testing import omex as omex_mod
    omex_file = omex_mod.omex_path(version, inst)
    omex_result = checks.load(version, f'cuflynx_instance_omex_test__{inst.name}')
    return {'name': inst.name, 'is_default': inst.is_default, 'files': inst.data_files(),
            # the generated COMBINE archive for CUFLynx (tools/build_instance_omex.py), and its CUFLynx check
            'omex': os.path.relpath(omex_file, version.dir) if os.path.isfile(omex_file) else None,
            'omex_status': omex_result.status if omex_result else None,
            'omex_message': omex_result.message if omex_result else 'not checked yet (tests/test_instance_omex.py)',
            'n_parameters': len(inst.parameters()), 'obs_data_name': inst.obs_data_name,
            'has_obs_data': inst.has_obs_data, 'calibration': calibration,
            'calibrated_file': os.path.basename(inst.calibrated_parameters_path)
            if os.path.isfile(inst.calibrated_parameters_path) else None,
            'rel_dir': os.path.relpath(inst.dir, version.dir),
            'tests': tests, 'baseline': tests[0], 'calibrate': tests[1]}


def not_applicable_reason(version, test):
    '''Why a report column doesn't apply to the version, or None when it does.'''
    if test in SUPERMODULE_TESTS and not version.is_supermodule:
        return NOT_A_SUPERMODULE
    if test == 'supermodule_equivalence_test' and not (version.spec.get('supermodule') or {}).get('equivalent'):
        return NO_EQUIVALENT
    if checks.is_cpp(version):
        if test in ('verification_test_invariants', 'verification_test_BC', 'verification_test_timestep', 'stability_test'):
            return checks.CPP_NOT_APPLICABLE
        if test in phlynx.PIPELINE_TESTS:
            return phlynx.CPP_REASON
    return None


def _equivalence_result(version):
    '''The supermodule_equivalence_test column: one result over the spec's supermodule.equivalent
    entries (each recorded on its own): failed if any failed, passed only if all ran and passed,
    otherwise not run.'''
    entries = (version.spec.get('supermodule') or {}).get('equivalent') or []
    results = [(e, checks.load(version, sm.equivalence_result_name(e))) for e in entries]
    lines = [f'{e["model"]} reproduces {e["reproduces"]}: '
             + (f'{r.status}: {r.message}' if r else 'not run') for e, r in results]
    if any(r is not None and r.status == checks.FAILED for _, r in results):
        status = checks.FAILED
    elif results and all(r is not None and r.status == checks.PASSED for _, r in results):
        status = checks.PASSED
    else:
        return None if all(r is None for _, r in results) else checks.Result(
            'supermodule_equivalence_test', checks.PENDING, 'not every system model has been compared yet', details=lines)
    plots = [p for _, r in results if r for p in r.plots]
    metrics = {'entries': [dict(r.metrics, status=r.status) for _, r in results if r]}
    if len(results) == 1:
        r = results[0][1]
        return checks.Result('supermodule_equivalence_test', status, r.message, r.metrics, plots, r.details, r.timestamp)
    return checks.Result('supermodule_equivalence_test', status, '; '.join(lines), metrics, plots, lines,
                         max((r.timestamp for _, r in results if r), default=''))


def version_tests(version):
    '''The version's report columns (REPORT_TESTS), the same for every version.'''
    tests = []
    for test in REPORT_TESTS:
        reason = not_applicable_reason(version, test)
        if reason:
            r = checks.Result(test, checks.NOT_APPLICABLE, reason)
        elif test == 'version_calibration':
            r = checks.version_calibration(version, _supermodule_index())
        elif test == 'supermodule_equivalence_test':
            r = _equivalence_result(version)
        else:
            r = checks.load(version, test)
        title = None
        entries = (version.spec.get('supermodule') or {}).get('equivalent') or []
        if test == 'supermodule_equivalence_test' and entries:
            title = 'Supermodule: reproduces ' + '; '.join(f'{e["reproduces"]} ({e["model"]})' for e in entries)
        entry = _test_entry(version, test, r, title=title)
        if r is not None and r.status in (checks.PASSED_IN_SUPER, checks.FAILED_IN_SUPER):
            entry['links'] = _supermodule_links(version, r)
        tests.append(entry)
    return tests


_index_cache = []


def _supermodule_index():
    '''The library's submodule -> supermodules index, built once per report run.'''
    if not _index_cache:
        _index_cache.append(checks.supermodule_index())
    return _index_cache[0]


def _supermodule_links(version, r):
    '''Links to the supermodule version pages a calibration in super comes from: href relative to
    this version's page, href_module relative to its module_type page.'''
    from cam_testing.library import load_version
    links = []
    for row in (r.metrics or {}).get('supermodules') or []:
        target = load_version(row['module_type'], row['version']).html_path
        links.append({'key': row['key'], 'status': row['status'],
                      'status_label': STATUS_LABEL.get(row['status'], row['status']),
                      'href': os.path.relpath(target, version.dir),
                      'href_module': os.path.relpath(target, version.mtype.dir)})
    return links


# ---- contents summary ---------------------------------------------------------------------------

CONTENTS_ABOUT = {
    'states': 'State variables: the variables with a d/dt equation in the CellML component',
    'algebraic': 'Algebraic variables: the variables an algebraic equation (x = ...) defines in the CellML component',
    'parameters': 'Parameters and constants: the config\'s variables_and_units of kind constant or global_constant',
    'boundary_conditions': 'Boundary conditions: the config\'s variables_and_units of kind boundary_condition (set by a '
                           'connected module, or a parameter when nothing connects them)',
    'ports': 'Ports: the config\'s entrance, exit and general ports, through which models connect the version',
    'equations': 'Equations: the equations of the CellML component',
    'instances': 'Instances: the parameter sets of this version (instances/)',
    'parts': 'Parts: the submodules this supermodule version is built from',
    'levels': 'Nested levels: how deep the supermodule nests (1 when every part is a component)',
}
SUPER_ABOUT_SUFFIX = ', summed over every part (all nested levels), counting each use of a version'
_contents_cache = {}


def _leaf_contents(version):
    '''Counts of one component version (cached per report run).'''
    if version.key in _contents_cache:
        return _contents_cache[version.key]
    kinds = [v[3].strip() for v in version.config.get('variables_and_units') or []]
    targets = []
    if version.format == 'cellml' and os.path.isfile(version.cellml_path):
        targets = component_equation_targets(version.cellml_path, version.module_type)
    out = {'states': len({v for v, k in targets if k == 'ode'}),
           'algebraic': len({v for v, k in targets if k == 'algebraic'}),
           'parameters': sum(k in ('constant', 'global_constant') for k in kinds),
           'boundary_conditions': kinds.count('boundary_condition'),
           'ports': len(_ports(version.config)), 'equations': len(targets)}
    _contents_cache[version.key] = out
    return out


def _sub_version(sub):
    from cam_testing.library import load_version
    try:
        return load_version(sub['module_type'], sub['module_subtype'])
    except (OSError, ValueError, KeyError):
        return None


def contents(version, _depth=0):
    '''The contents summary bar: {key: count} (states, algebraic, parameters, boundary_conditions, ports,
    equations, instances; a supermodule adds parts and levels, its other counts summed over its parts).'''
    if not version.is_supermodule:
        return dict(_leaf_contents(version), instances=len(version.instance_names()))
    if _depth > 10:
        raise ValueError(f'{version.key}: supermodules nest more than 10 levels (a cycle?)')
    total = {k: 0 for k in ('states', 'algebraic', 'parameters', 'boundary_conditions', 'ports', 'equations')}
    levels, leaves = 1, 0
    for sub in version.submodules:
        sv = _sub_version(sub)
        if sv is None:
            continue
        c = contents(sv, _depth + 1)
        for k in total:
            total[k] += c[k]
        if sv.is_supermodule:
            levels = max(levels, c['levels'] + 1)
            leaves += c['components']
        else:
            leaves += 1
    return dict(total, instances=len(version.instance_names()), parts=len(version.submodules), levels=levels,
                components=leaves)


def contents_bar(version, c, comp_id):
    '''[{key, count, label, about, href}] for the version page's summary bar (and, compact, the module_type page).'''
    sup = version.is_supermodule
    structure = f'#{comp_id}-structure'
    href = {'states': f'#{comp_id}-variables', 'algebraic': f'#{comp_id}-equations',
            'parameters': f'#{comp_id}-variables', 'boundary_conditions': f'#{comp_id}-variables',
            'ports': f'#{comp_id}-ports', 'equations': f'#{comp_id}-equations', 'instances': '#instances',
            'parts': structure, 'levels': structure}
    labels = {'states': ('state', 'states'), 'algebraic': ('algebraic', 'algebraic'),
              'parameters': ('parameter', 'parameters'), 'boundary_conditions': ('boundary condition', 'boundary conditions'),
              'ports': ('port', 'ports'), 'equations': ('equation', 'equations'), 'instances': ('instance', 'instances'),
              'parts': ('part', 'parts'), 'levels': ('nested level', 'nested levels')}
    short = {'states': 'st', 'algebraic': 'alg', 'parameters': 'par', 'boundary_conditions': 'BC', 'ports': 'ports',
             'equations': 'eq', 'instances': 'inst', 'parts': 'parts', 'levels': 'levels'}
    keys = ['states', 'algebraic', 'parameters', 'boundary_conditions', 'ports', 'equations', 'instances']
    if sup:
        keys = ['parts', 'levels'] + keys
    out = []
    for k in keys:
        n = c.get(k, 0)
        about = CONTENTS_ABOUT[k]
        if sup and k not in ('instances', 'parts', 'levels'):
            about += SUPER_ABOUT_SUFFIX + f' ({c.get("components", 0)} component versions)'
            link = structure        # a supermodule has no equations or variables of its own: its parts do
        else:
            link = href[k]
        out.append({'key': k, 'count': n, 'label': labels[k][0 if n == 1 else 1],
                    'short': 'level' if k == 'levels' and n == 1 else short[k], 'about': about,
                    'href': link})
    return out


# ---- unit consistency ---------------------------------------------------------------------------

def unit_consistency_context(version, equations):
    '''The version page's "Unit consistency" block: shown only when libcellml found failures (or could not
    check the equations); None when every equation is unit-consistent or the check doesn't apply.'''
    r = checks.load_unit_consistency(version)
    if r.get('status') not in (checks.UNITS_INCONSISTENT, checks.UNITS_NOT_CHECKED):
        return None
    groups = {}
    for f in r.get('failures') or []:
        key = f.get('equation') or f['message']
        g = groups.setdefault(key, {'equation': f.get('equation'), 'variable': f.get('variable'),
                                    'latex': equations[f['equation_index']] if f.get('equation_index') is not None
                                    and f['equation_index'] < len(equations) else None,
                                    'issues': []})
        g['issues'].append(f)
    return {'status': r['status'], 'message': r.get('message', ''), 'n_failures': r.get('n_failures', 0),
            'equations': list(groups.values()), 'definition_issues': r.get('definition_issues') or [],
            'other_issues': r.get('other_issues') or [], 'libcellml_version': r.get('libcellml_version'),
            'timestamp': r.get('timestamp')}


# ---- supermodule structure ----------------------------------------------------------------------

def _port_types(version, sides):
    if version is None or version.is_supermodule:
        return set()
    return {p['port_type'] for s in sides for p in version.config.get(s) or []}


def _inner_version(sv, inner):
    '''The version of part ``inner`` of supermodule version ``sv`` (or sv itself when inner is None).'''
    if sv is None or inner is None:
        return sv
    sub = next((s for s in sv.submodules if s['name'] == inner), None)
    return _sub_version(sub) if sub else None


def _edge_ports(src, dst):
    '''The port types libcuflynx can connect from src (exit / general) to dst (entrance / general).'''
    return sorted(_port_types(src, ('exit_ports', 'general_ports')) & _port_types(dst, ('entrance_ports', 'general_ports')))


def _global_names(version, _depth=0):
    '''The global constants a version uses: a component's global_constant variables; a supermodule's globals
    and its parts'.'''
    if not version.is_supermodule:
        return {v[0] for v in version.config.get('variables_and_units') or [] if v[3].strip() == 'global_constant'}
    out = set(version.supermodule_globals)
    if _depth < 10:
        for sub in version.submodules:
            sv = _sub_version(sub)
            if sv is not None:
                out |= _global_names(sv, _depth + 1)
    return out


def _system_model_uses(version):
    '''(model, record) for every system-model record using this version (system_models/*/*/*_vessel_array.json).'''
    import glob
    from cam_testing.library import SYSTEM_MODELS_DIR
    out = []
    for path in sorted(glob.glob(os.path.join(SYSTEM_MODELS_DIR, '*', '*', '*_vessel_array.json'))):
        try:
            with open(path) as f:
                rows = json.load(f)
        except (OSError, ValueError):
            continue
        for r in rows if isinstance(rows, list) else []:
            if isinstance(r, dict) and (r.get('module_type'), r.get('module_subtype')) == (version.vessel_type, version.name):
                out.append((os.path.relpath(os.path.dirname(path), SYSTEM_MODELS_DIR).replace(os.sep, '/'), r))
    return out


def _mermaid_id(prefix, name):
    return prefix + ''.join(ch if ch.isalnum() else '_' for ch in name)


def _mermaid_text(text):
    return text.replace('"', '#quot;')


def _structure_tree(version, from_dir, _depth=0):
    '''Nested [{name, key, href, is_supermodule, children}] of a supermodule's parts, for the collapsible tree.'''
    out = []
    for sub in version.submodules:
        sv = _sub_version(sub)
        node = {'name': sub['name'], 'key': f'{sub["module_type"]}/{sub["module_subtype"]}',
                'href': os.path.relpath(sv.html_path, from_dir) if sv else None,
                'is_supermodule': bool(sv and sv.is_supermodule), 'children': []}
        if sv is not None and sv.is_supermodule and _depth < 10:
            node['children'] = _structure_tree(sv, from_dir, _depth + 1)
        out.append(node)
    return out


def supermodule_structure(version):
    '''A supermodule version's parts and connections: {parts, edges, external, globals, mermaid, tree}.
    Edges are the internal connections (inp_instances / out_instances between parts, and a part's
    per_submodule_inputs / per_submodule_outputs naming a sibling, which couple one of that part's own
    parts), with the port types libcuflynx can connect where the configs give them. External are the
    couplings from outside: the supermodule versions that use this one and the system models whose
    records name it, through per_submodule_inputs / per_submodule_outputs.'''
    if not version.is_supermodule:
        return None
    subs = version.submodules
    versions = {s['name']: _sub_version(s) for s in subs}
    # the supermodule's own instance: rows <var>_<part> override the part's instance; others are globals
    names = sorted(versions, key=len, reverse=True)
    overrides, global_rows = {}, []
    for p in version.default_instance.parameters():
        owner = next((n for n in names if p.variable_name.endswith('_' + n)), None)
        if owner is None:
            global_rows.append({'name': p.variable_name, 'value': p.value, 'units': p.units})
        else:
            overrides.setdefault(owner, []).append({'name': p.variable_name, 'var': p.variable_name[:-len(owner) - 1],
                                                    'value': p.value, 'units': p.units})
    sup_globals = set(version.supermodule_globals)

    edges = {}

    def add_edge(src, dst, src_inner=None, dst_inner=None):
        k = (src, dst, src_inner, dst_inner)
        if k in edges or src not in versions or dst not in versions:
            return
        a = _inner_version(versions[src], src_inner)
        b = _inner_version(versions[dst], dst_inner)
        edges[k] = {'src': src, 'dst': dst, 'src_inner': src_inner, 'dst_inner': dst_inner, 'ports': _edge_ports(a, b)}

    for s in subs:
        for t in s.get('out_instances') or []:
            add_edge(s['name'], t)
        for t in s.get('inp_instances') or []:
            add_edge(t, s['name'])
        for inner, others in (s.get('per_submodule_outputs') or {}).items():
            for t in others:
                add_edge(s['name'], t, src_inner=inner)
        for inner, others in (s.get('per_submodule_inputs') or {}).items():
            for t in others:
                add_edge(t, s['name'], dst_inner=inner)
    # a plain connection to a supermodule part that a per_submodule_* entry resolves (e.g. axon's inp_instances
    # ["soma"] and soma's per_submodule_outputs {"membrane": ["axon"]}) is one coupling: keep the resolved one
    resolved = {(e['src'], e['dst']) for e in edges.values() if e['src_inner'] or e['dst_inner']}
    edges = [e for e in edges.values() if e['src_inner'] or e['dst_inner'] or (e['src'], e['dst']) not in resolved]
    for e in edges:
        e['label'] = (e['src'] + (f'/{e["src_inner"]}' if e['src_inner'] else '') + ' → '
                      + e['dst'] + (f'/{e["dst_inner"]}' if e['dst_inner'] else ''))

    # couplings from outside: supermodules using this version, and system models
    external = {}

    def add_external(part, ext, direction, where, href):
        k = (part, ext, direction)
        row = external.setdefault(k, {'part': part, 'external': ext, 'direction': direction, 'where': []})
        if all(w['name'] != where for w in row['where']):
            row['where'].append({'name': where, 'href': href})

    for parent in _supermodule_index().get((version.vessel_type, version.name)) or []:
        href = os.path.relpath(parent.html_path, version.dir)
        for rec in parent.submodules:
            if (rec['module_type'], rec['module_subtype']) != (version.vessel_type, version.name):
                continue
            for part, others in (rec.get('per_submodule_inputs') or {}).items():
                for o in others:
                    add_external(part, o, 'in', parent.key, href)
            for part, others in (rec.get('per_submodule_outputs') or {}).items():
                for o in others:
                    add_external(part, o, 'out', parent.key, href)
    for model, rec in _system_model_uses(version):
        for part, others in (rec.get('per_submodule_inputs') or {}).items():
            for o in others:
                add_external(part, o, 'in', f'system model {model}', None)
        for part, others in (rec.get('per_submodule_outputs') or {}).items():
            for o in others:
                add_external(part, o, 'out', f'system model {model}', None)
    external = sorted(external.values(), key=lambda r: (r['part'], r['direction'], r['external']))

    parts = []
    for s in subs:
        sv = versions[s['name']]
        ports = []
        if sv is not None and not sv.is_supermodule:
            ports = [f'{p["side"]} {p["type"]}' for p in _ports(sv.config)]
        uses = _global_names(sv) if sv is not None else set()
        parts.append({
            'name': s['name'], 'module_type': s['module_type'], 'version': s['module_subtype'],
            'instance': s.get('instance') or '(default)', 'is_supermodule': bool(sv and sv.is_supermodule),
            'n_parts': len(sv.submodules) if sv is not None and sv.is_supermodule else 0,
            'href': os.path.relpath(sv.html_path, version.dir) if sv is not None else None,
            'missing': sv is None, 'ports': ports,
            'inputs': [e for e in edges if e['dst'] == s['name']],
            'outputs': [e for e in edges if e['src'] == s['name']],
            'external': [r for r in external if r['part'] == s['name']],
            'overrides': overrides.get(s['name'], []),
            'globals': sorted(uses & (sup_globals | {g['name'] for g in global_rows})),
        })

    # the diagram (mermaid flowchart)
    lines = ['flowchart LR']      # many parts stack vertically: readable at page width
    for p in parts:
        nid = _mermaid_id('p_', p['name'])
        label = f'<b>{p["name"]}</b><br/>{p["module_type"]} / {p["version"]}'
        if p['is_supermodule']:
            label += f'<br/><i>supermodule, {p["n_parts"]} parts</i>'
        shape = ('[["', '"]]') if p['is_supermodule'] else ('["', '"]')
        lines.append(f'  {nid}{shape[0]}{_mermaid_text(label)}{shape[1]}')
        if p['href']:
            lines.append(f'  click {nid} "{p["href"]}" "open {p["module_type"]} {p["version"]}"')
    labelled = len(edges) <= 20
    for e in edges:
        label = []
        if e['src_inner'] or e['dst_inner']:
            label.append(e['label'])
        if e['ports'] and labelled:
            label.append(', '.join(e['ports']))
        arrow = f'-->|"{_mermaid_text(" · ".join(label))}"|' if label and (labelled or e['src_inner'] or e['dst_inner']) else '-->'
        lines.append(f'  {_mermaid_id("p_", e["src"])} {arrow} {_mermaid_id("p_", e["dst"])}')
    # one dashed node per outside module (e.g. a volume_sum coupled to all four chambers)
    xids = {}
    for r in external:
        where = sorted({w['name'] for w in r['where']})
        key = (r['external'], tuple(where))
        if key not in xids:
            xids[key] = f'x{len(xids)}'
            short = [w.rsplit('/', 1)[-1] if w.startswith('system model ') else w for w in where]   # model name only
            shown = ', '.join(short[:2]) + (f' +{len(short) - 2}' if len(short) > 2 else '')
            label = f'‹{r["external"]}›<br/><i>{shown}</i>'
            lines.append(f'  {xids[key]}(["{_mermaid_text(label)}"]):::ext')
        xid, pid = xids[key], _mermaid_id('p_', r['part'])
        lines.append(f'  {xid} -.-> {pid}' if r['direction'] == 'in' else f'  {pid} -.-> {xid}')
    lines.append('  classDef ext stroke-dasharray: 5 4')
    return {'parts': parts, 'edges': edges, 'external': external, 'labelled': labelled,
            'globals': sorted(sup_globals), 'global_rows': global_rows, 'mermaid': '\n'.join(lines),
            'tree': _structure_tree(version, version.dir)}


def component_context(version):
    '''The version page's content (what was a component's section of a module page).'''
    equations, unsupported, cellml_vars = _cellml_equations(version)
    params = {p.variable_name: p for p in version.parameters()}
    spread = ranges.validated_spread(version) if not version.is_supermodule else {}
    variables = []
    for name, units, access, kind in version.config.get('variables_and_units') or []:
        cv = cellml_vars.get(name)
        p = params.get(name)
        variables.append({
            'name': name, 'units': units, 'kind': kind.strip(), 'access': access,
            'initial': cv[2] if cv else None,
            'value': (p.value + ' (proposed)' if p.proposed else p.value) if p else None,
            'todo': bool(p and (p.is_todo or p.proposed)),
            'reference': p.data_reference if p else '',
            'sourced': p.is_sourced if p else None,
            'tested': ranges.tested_range(version, p) if p else None,
            'validated': spread.get(name),
        })
    tests = version_tests(version)
    risk_result = risk.load(version) if not version.is_supermodule else None
    if risk_result is not None and not all(os.path.isfile(os.path.join(version.dir, p)) for p in risk_result.get('plots', [])):
        risk_result['plots'] = risk.make_plots(version, risk_result)   # e.g. on a fresh clone
    proposals = version.spec.get('reference_proposals') or {}
    bib_keys = set(bib.read(bib.bib_path(version))) | set(bib.read(bib.proposed_bib_path(version)))
    references = []
    for p in version.parameters():
        key = bib.reference_key(p.data_reference)
        references.append({'name': p.variable_name, 'value': p.value + (' (proposed)' if p.proposed else ''), 'units': p.units,
                           'reference': p.data_reference, 'key': key if key in bib_keys else None,
                           'note': p.data_reference.split(';', 1)[1].strip() if ';' in p.data_reference else '',
                           'sourced': p.is_sourced, 'proposal': proposals.get(p.variable_name)})
    return {
        'risk': risk_result, 'references': references, 'has_proposals': bool(proposals),
        'id': version.id, 'key': version.key, 'label': version.label, 'vessel_type': version.vessel_type,
        'BC_type': version.BC_type, 'module_type': version.module_type, 'category': version.category,
        'format': version.format, 'config_notes': version.config.get('notes'),
        'description': version.spec.get('description'),
        'submodules': version.submodules,
        'skip': version.spec.get('skip'), 'notes': version.spec.get('notes'),
        'sweep_rationale': (version.spec.get('bc_sweep') or {}).get('rationale'),
        'equations': equations, 'unsupported': sorted(unsupported),
        'ports': _ports(version.config), 'variables': variables,
        'todo': version.todo_parameters(), 'unsourced': version.unsourced_parameters(),
        'invariants': version.spec.get('invariants') or [],
        'tests': tests, 'instances': [instance_context(version, i) for i in version.instances()],
        'phlynx_compatible': _phlynx_compatible(tests),
        'contents': contents_bar(version, contents(version), version.id),
        'unit_consistency': unit_consistency_context(version, equations),
        'structure_view': supermodule_structure(version),
    }


def _phlynx_compatible(tests):
    '''True when PhLynx builds and exports it, CUFLynx simulates it and it matches libcuflynx; None if not run.'''
    status = {t['key']: t['status'] for t in tests if t['key'] in phlynx.PIPELINE_TESTS}
    if not status or all(v is None for v in status.values()):
        return None
    return all(v == checks.PASSED for v in status.values())


def _status_counts(components):
    counts = {'passed': 0, 'failed': 0, 'skipped': 0, 'pending': 0, 'not_applicable': 0, None: 0}
    for c in components:
        for t in c['tests'] + [t for i in c['instances'] for t in i['tests']]:
            status = COUNT_AS.get(t['status'], t['status'])
            counts[status] = counts.get(status, 0) + 1
    return counts


def _structure(version):
    path = os.path.join(version.results_dir, 'structure.json')
    if os.path.isfile(path):
        with open(path) as f:
            return json.load(f)
    return None


def _libcuflynx_version():
    try:
        from importlib.metadata import version
        return version('libcuflynx')
    except Exception:
        return 'unknown'


def _emph(text):
    '''Escapes text and renders **...** (a reference's flagged difference from the literature) bold.'''
    import re
    from markupsafe import Markup, escape
    return Markup(re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', str(escape(text or ''))))


def _env():
    env = Environment(loader=FileSystemLoader(TEMPLATES), autoescape=select_autoescape(['html']))
    env.filters['fmt'] = _fmt
    env.filters['emph'] = _emph
    return env


def _now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')


def _module_phlynx(components):
    run = [c for c in components if c['phlynx_compatible'] is not None]
    ok = sum(bool(c['phlynx_compatible']) for c in run)
    total = len(components)
    return {'run': len(run), 'compatible': ok, 'total': total,
            'all': bool(run) and len(run) == total and ok == total}


def _bc_summary(comp):
    '''One line on the version's boundary-condition sweep: skipped points, failures by class and the
    recommended solver settings; None when it hasn't run.'''
    t = next((t for t in comp['tests'] if t['key'] == 'verification_test_BC'), None)
    m = (t or {}).get('metrics') or {}
    if 'n_runs' not in m or 'failures' not in m:
        return None
    parts = [f"{m['n_runs'] - m.get('n_failures', 0)}/{m['n_runs']} runs pass"]
    if m.get('skipped_out_of_range'):
        parts.append(f"{len(m['skipped_out_of_range'])} skipped (outside the valid range)")
    if m.get('n_parameter_failures'):
        parts.append(f"{m['n_parameter_failures']} parameter failure{'s' if m['n_parameter_failures'] > 1 else ''}")
    if m.get('n_numerical_failures'):
        n = m['n_numerical_failures']
        rc = m.get('recommended')
        parts.append(f"{n} numerical failure{'s' if n > 1 else ''}"
                     + (f", fixed by the recommended {rc['label']}" if rc else ', not fixed by any candidate setting'))
    return '; '.join(parts)


def version_context(version):
    comp = component_context(version)
    return {
        'name': version.label, 'key': version.key, 'module_type': version.vessel_type, 'version': version.name,
        'category': version.category, 'location': version.mtype.location, 'reviewed': version.reviewed,
        'components': [comp],
        'is_supermodule': version.is_supermodule,
        'all_sourced': version.all_sourced(),
        'bibliography': [{'key': k, 'text': bib.format_entry(f), 'link': bib.link(f), 'proposed': False}
                         for k, f in bib.read(bib.bib_path(version)).items()]
                        + [{'key': k, 'text': bib.format_entry(f), 'link': bib.link(f), 'proposed': True}
                           for k, f in bib.read(bib.proposed_bib_path(version)).items()],
        'review': version.spec.get('review'), 'review_scope': version.spec.get('review_scope'),
        'bib_keys': sorted(set(bib.read(bib.bib_path(version))) | set(bib.read(bib.proposed_bib_path(version)))),
        'bib_file': os.path.basename(bib.bib_path(version)),
        'n_sourced': sum(p.is_sourced for p in version.parameters()), 'n_parameters': len(version.parameters()),
        'counts': _status_counts([comp]), 'structure': _structure(version),
        'max_risk': comp['risk']['failure_probability'] if comp['risk'] else None,
        'known_issues': (version.spec.get('known_issues') or [])
                        + [f'{t}: {why}' for t, why in (version.spec.get('expected_failures') or {}).items()],
        'test_keys': [t['key'] for t in comp['tests']], 'test_short': {t['key']: t['short'] for t in comp['tests']},
        'phlynx': _module_phlynx([comp]), 'bc_summary': _bc_summary(comp),
        'licence': version.licence, 'licence_url': LICENCES.get(version.licence), 'creators': version.creators,
        'unit_failures': (comp['unit_consistency'] or {}).get('n_failures', 0),
        'module_type_href': f'../../{version.vessel_type}.html',
        'generated': _now(), 'libcuflynx_version': _libcuflynx_version(),
    }


def build_version(version):
    ctx = version_context(version)
    with open(version.html_path, 'w') as f:
        f.write(_env().get_template('version.html').render(**ctx))
    return version.html_path, ctx


def module_context(name, version_contexts=None):
    '''A module_type: a summary row per version (from the version pages' contexts).'''
    mtype = load_module_type(name)
    vctx = version_contexts if version_contexts is not None else [version_context(v) for v in mtype.versions()]
    rows = []
    for c in vctx:
        comp = c['components'][0]
        rows.append({'version': c['version'], 'href': f'versions/{c["version"]}/{name}_{c["version"]}.html',
                     'format': comp['format'], 'component_type': comp['module_type'], 'notes': comp['config_notes'],
                     'counts': c['counts'], 'reviewed': c['reviewed'], 'max_risk': c['max_risk'],
                     'all_sourced': c['all_sourced'], 'n_sourced': c['n_sourced'], 'n_parameters': c['n_parameters'],
                     'phlynx': c['phlynx'], 'phlynx_compatible': comp['phlynx_compatible'], 'tests': comp['tests'],
                     'instances': comp['instances'], 'submodules': comp['submodules'],
                     'known_issues': len(c['known_issues']), 'is_supermodule': c['is_supermodule'],
                     'contents': comp['contents'], 'unit_failures': c['unit_failures'],
                     'unit_status': (comp['unit_consistency'] or {}).get('status'),
                     'creators': c['creators'], 'licence': c['licence']})
    counts = {}
    for c in vctx:
        for k, n in c['counts'].items():
            counts[k] = counts.get(k, 0) + n
    return {
        'name': name, 'category': mtype.category, 'group': mtype.group, 'relpath': mtype.relpath,
        'location': mtype.location, 'parent': mtype.parent, 'nested': mtype.nested, 'versions': rows,
        'reviewed': bool(rows) and all(r['reviewed'] for r in rows),
        'counts': counts, 'n_versions': len(rows),
        'n_instances': sum(len(r['instances']) for r in rows),
        'all_sourced': all(r['all_sourced'] for r in rows),
        'n_sourced': sum(r['n_sourced'] for r in rows), 'n_parameters': sum(r['n_parameters'] for r in rows),
        'max_risk': max((r['max_risk'] for r in rows if r['max_risk'] is not None), default=None),
        'known_issues': sum(r['known_issues'] for r in rows),
        'phlynx': {'run': sum(r['phlynx']['run'] for r in rows),
                   'compatible': sum(r['phlynx']['compatible'] for r in rows),
                   'total': len(rows),
                   'all': bool(rows) and all(r['phlynx']['all'] for r in rows)},
        'ready_for_review': any(c.get('review') and not c['reviewed'] for c in vctx),
        'test_keys': REPORT_TESTS, 'test_short': TEST_SHORT,
        'generated': _now(), 'libcuflynx_version': _libcuflynx_version(),
    }


def build_module(name):
    '''Writes every version page of a module_type and the module_type page; returns (path, context).'''
    mtype = load_module_type(name)
    vctx = [build_version(v)[1] for v in mtype.versions()]
    ctx = module_context(name, vctx)
    ctx['version_contexts'] = vctx
    with open(mtype.html_path, 'w') as f:
        f.write(_env().get_template('module_type.html').render(**ctx))
    return mtype.html_path, ctx


def build_index(contexts, out_path, href):
    '''The site index: module_types grouped by category (a module_type directly under modules/ and the
    ones nested in it form their own group), each with links to its versions. A nested module_type
    is listed by its path in the group, e.g. neuron/soma under cell.'''
    groups = {}
    for c in contexts:
        groups.setdefault(c['group'], []).append(c)
    out_groups = []
    for cat, ctxs in sorted(groups.items()):
        rows = []
        for c in ctxs:
            c['display'] = c['relpath'][len(c['category']) + 1:] if c['category'] else c['relpath']
        for c in sorted(ctxs, key=lambda c: c['display'].lower()):
            rows.append({'name': c['name'], 'display': c['display'], 'href': href(c['relpath'], f'{c["name"]}.html'),
                         'n_versions': c['n_versions'], 'n_instances': c['n_instances'], 'counts': c['counts'],
                         'reviewed': c['reviewed'], 'known_issues': c['known_issues'], 'max_risk': c['max_risk'],
                         'all_sourced': c['all_sourced'], 'n_sourced': c['n_sourced'], 'n_parameters': c['n_parameters'],
                         'phlynx': c['phlynx'], 'ready_for_review': c['ready_for_review'] and not c['reviewed'],
                         'versions': [{'name': v['version'], 'href': href(c['relpath'], v['href']),
                                       'counts': v['counts'], 'reviewed': v['reviewed'],
                                       'phlynx_compatible': v['phlynx_compatible'],
                                       'is_supermodule': v['is_supermodule']} for v in c['versions']]})
        out_groups.append({'category': cat, 'anchor': cat.replace('/', '-'), 'rows': rows,
                           'n_versions': sum(r['n_versions'] for r in rows)})
    all_rows = [r for g in out_groups for r in g['rows']]
    totals = {k: sum(r['counts'].get(k, 0) for r in all_rows)
              for k in ('passed', 'failed', 'skipped', 'pending', 'not_applicable', None)}
    html_text = _env().get_template('index.html').render(
        groups=out_groups, totals=totals, n_module_types=len(all_rows), n_versions=sum(r['n_versions'] for r in all_rows),
        generated=_now(), libcuflynx_version=_libcuflynx_version())
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(html_text)
    return out_path


def build_review_queue(contexts, out_path, href, standalone=True):
    '''One page listing every version awaiting review. A review written for a whole old module was
    copied to each of its versions; such a review is shown once, with the versions it covers.'''
    queue, reviewed, by_review = [], [], {}
    for c in contexts:
        for vc in c.get('version_contexts') or []:
            if vc['reviewed']:
                reviewed.append(vc['key'])
                continue
            if not vc.get('review'):
                continue
            key = json.dumps(vc['review'], sort_keys=True)
            comp = vc['components'][0]
            item = {'key': vc['key'],
                    'href': href(c['relpath'], f'versions/{vc["version"]}/{c["name"]}_{vc["version"]}.html'),
                    'counts': vc['counts'], 'n_known': len(vc['known_issues']),
                    'n_sourced': vc['n_sourced'], 'n_parameters': vc['n_parameters'], 'max_risk': vc['max_risk'],
                    'n_proposals': sum(1 for r in comp['references'] if r['proposal'])}
            if key not in by_review:
                by_review[key] = {'review': vc['review'], 'scope': vc.get('review_scope'), 'versions': []}
                queue.append(by_review[key])
            by_review[key]['versions'].append(item)
    for q in queue:
        vs = q['versions']
        q['name'] = vs[0]['key'] if len(vs) == 1 else f'{vs[0]["key"]} and {len(vs) - 1} more'
        q['anchor'] = vs[0]['key'].replace('/', '--')
        q['href'] = vs[0]['href']
        q['n_components'] = len(vs)
        q['counts'] = {k: sum(v['counts'].get(k, 0) for v in vs) for k in ('passed', 'failed')}
        q['n_known'] = sum(v['n_known'] for v in vs)
        q['n_sourced'] = sum(v['n_sourced'] for v in vs)
        q['n_parameters'] = sum(v['n_parameters'] for v in vs)
        q['n_proposals'] = sum(v['n_proposals'] for v in vs)
        q['max_risk'] = max((v['max_risk'] for v in vs if v['max_risk'] is not None), default=None)
    html_text = _env().get_template('review_queue.html').render(
        modules=queue, reviewed=reviewed, standalone=standalone, generated=_now())
    with open(out_path, 'w') as f:
        f.write(html_text)
    return out_path


def _href(relpath, page):
    return f'modules/{relpath}/{page}'


def assemble_site(names, contexts):
    '''site/index.html + site/modules/<module_type path>/ (its page, and each version's page and
    plots; a nested module_type's directory is inside its parent's), as GitHub Pages serves it.'''
    if os.path.isdir(SITE_DIR):
        shutil.rmtree(SITE_DIR)
    for name in names:
        mtype = load_module_type(name)
        dest = os.path.join(SITE_DIR, 'modules', mtype.relpath)
        os.makedirs(dest, exist_ok=True)
        if not os.path.isfile(mtype.html_path) or any(not os.path.isfile(v.html_path) for v in mtype.versions()):
            build_module(name)
        shutil.copy2(mtype.html_path, dest)
        for v in mtype.versions():
            vdest = os.path.join(dest, 'versions', v.name)
            os.makedirs(vdest, exist_ok=True)
            shutil.copy2(v.html_path, vdest)
            if os.path.isdir(v.plots_dir):
                shutil.copytree(v.plots_dir, os.path.join(vdest, 'plots'))
            # the instances' CUFLynx archives, where they have been built (make omex)
            for inst in v.instances():
                archive = omex_mod.omex_path(v, inst)
                if os.path.isfile(archive):
                    idest = os.path.join(vdest, os.path.relpath(inst.dir, v.dir))
                    os.makedirs(idest, exist_ok=True)
                    shutil.copy2(archive, idest)
    build_review_queue(contexts, os.path.join(SITE_DIR, 'review_queue.html'), _href)
    return build_index(contexts, os.path.join(SITE_DIR, 'index.html'), _href)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--module', action='append', default=[],
                        help='module_type (with those nested in it) or path under modules/ (repeatable); default all')
    parser.add_argument('--site', action='store_true', help='also assemble site/ for GitHub Pages')
    args = parser.parse_args(argv)
    contexts = {}
    _index_cache.clear()
    _contents_cache.clear()
    for name in select_module_types(args.module):
        out, ctx = build_module(name)
        contexts[name] = ctx
        print(f'wrote {os.path.relpath(out, REPO_ROOT)} and {ctx["n_versions"]} version page(s)')
    if args.site:
        # the index always covers every module_type, even when only some pages were rebuilt
        for name in module_type_names():
            if name not in contexts:
                vctx = [version_context(v) for v in load_module_type(name).versions()]
                contexts[name] = dict(module_context(name, vctx), version_contexts=vctx)
        out = assemble_site(module_type_names(), [contexts[n] for n in module_type_names()])
        print(f'wrote {os.path.relpath(out, REPO_ROOT)}')


if __name__ == '__main__':
    main()
