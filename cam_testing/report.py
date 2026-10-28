"""
HTML reports at two levels, and a site index:

    <module_type>/<module_type>.html                        every version, linking to its page
    <module_type>/versions/<v>/<module_type>_<v>.html       the version: equations, variables,
                                                            tests, references, and its parameterisations
                                                            with their validation / calibration
    system_models/<category>/<model>/<model>.html           a system model: its figure (PhLynx screenshot, or a
                                                            generated diagram) with each instance linked to its
                                                            version page, its parameters and its system tests

    python -m cam_testing.report                  # every module_type and its versions, and every system model
    python -m cam_testing.report --module heart   # a module_type and those nested in it (or a category: --module cell)
    python -m cam_testing.report --site           # also assemble site/ as GitHub Pages serves it
    python -m cam_testing.report --zip out.zip    # site/ zipped with a README.txt, to send by email

Pages read the result JSON and plots the tests wrote, so a report shows the last local (or
CI) test run. Plots are referenced relatively (plots/...), so a page works opened from disk.
"""
import argparse
import datetime
import json
import os
import shutil
import urllib.parse

from jinja2 import Environment, FileSystemLoader, select_autoescape

from cam_testing import bib, bondgraph, checks, module_array, phlynx, ranges, risk
from cam_testing import omex as omex_mod
from cam_testing import supermodule as sm
from cam_testing import paths
from cam_testing.library import LICENCES, load_module_type, module_type_names, select_module_types
from cam_testing.mathml import component_equation_targets, component_equations, component_variables

TEMPLATES = os.path.join(os.path.dirname(__file__), 'templates')

# The apps a version page opens its files in. The CellDL editor and PhLynx take the file in the
# URL's hash, standard base64 of its bytes, which the page builds when the link is clicked
# (templates/_open_in.js); CUFLynx takes the archive's absolute URL, percent-encoded, through its
# cuflynx:// link scheme. CAM_CELLDL_EDITOR_URL, when set, replaces CELLDL_EDITOR_URL (e.g. a local
# editor build, http://127.0.0.1:8475/).
CELLDL_EDITOR_URL = 'https://physiomelinks.github.io/CellDLEditorComponent/'
CELLDL_EDITOR_ENV = 'CAM_CELLDL_EDITOR_URL'
# A note shown with a version's bond graph and its CellDL link when set (python -m cam_testing.report
# --celldl-note TEXT, or this variable): e.g. that the linked editor may not draw some element types yet.
CELLDL_NOTE_ENV = 'CAM_CELLDL_NOTE'
PHLYNX_OPEN = 'https://phlynx.com/?open=omex#'
CUFLYNX_OPEN = 'cuflynx://open?url='
PHLYNX_NOT_EXPORTED = 'not exported yet (the PhLynx pipeline, make pipeline, keeps its export with the results)'


def __getattr__(name):
    # the former constant: site/ of the repo under test (cam_testing.paths)
    if name == 'SITE_DIR':
        return paths.roots().site_dir
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

TEST_TITLES = {
    'review_test': 'Reviewed',
    'bondgraph_celldl_test': 'Bond graph: opens in CellDL',
    'bondgraph_reactions_test': 'Bond graph: reactions match Rhea',
    'run_test': 'Run',
    'verification_test_invariants': 'Verification: invariants',
    'verification_test_BC': 'Verification: boundary conditions',
    'verification_test_timestep': 'Verification: timestep convergence',
    'stability_test': 'Stability: solvers & settings',
    'validation_test_baseline': 'Validation: baseline data',
    'validation_test_calibrate': 'Validation: calibrate & predict',
    'version_calibration': 'Calibration (all parameterisations)',
    'phlynx_export_test': 'PhLynx: build & export',
    'cuflynx_simulate_test': 'CUFLynx: import & simulate',
    'phlynx_equivalence_test': 'PhLynx → CUFLynx vs libcuflynx',
    'supermodule_structure_test': 'Supermodule: structure',
    'supermodule_equivalence_test': 'Supermodule: reproduces',
    'coupled_validation_test': 'Validation: coupled, against the CellML grid',
}
TEST_SHORT = {
    'review_test': 'Reviewed', 'run_test': 'Run', 'verification_test_invariants': 'Invariants', 'verification_test_BC': 'BC sweep', 'verification_test_timestep': 'Timestep',
    'stability_test': 'Stability', 'validation_test_baseline': 'Baseline', 'validation_test_calibrate': 'Calibrate', 'version_calibration': 'Calibration',
    'phlynx_export_test': 'PhLynx export', 'cuflynx_simulate_test': 'CUFLynx simulate', 'phlynx_equivalence_test': 'PhLynx equivalence',
    'supermodule_structure_test': 'Supermodule structure', 'supermodule_equivalence_test': 'Supermodule reproduces',
    'coupled_validation_test': 'Coupled',
    'bondgraph_celldl_test': 'Opens in CellDL', 'bondgraph_reactions_test': 'Reactions vs Rhea',
}
TEST_ABOUT = {
    'review_test': 'Whether the version has been reviewed: its equations, parameter sources (each value traced '
                   'to its source), tests, data and known issues checked with its owner (tests.yaml reviewed: '
                   'true). A version not reviewed yet fails this check. It is still published and usable, but '
                   'its values and behaviour are provisional, and CI runs only its quick tests (Run, the '
                   'PhLynx -> CUFLynx pipeline) until it is reviewed, so its other columns show Not run.',
    'bondgraph_celldl_test': 'The version\'s bond graph (<version>_bondgraph.svg) is a CellDL document the editor '
                             'opens: Turtle BG-RDF with one bond-graph model, every element and connection drawn and '
                             'nothing drawn that the RDF lacks, and the same graph as <version>_bondgraph.ttl. With '
                             'BG_TOOLS set, bg2cellml also makes CellML from it. N/A without a bond graph.',
    'bondgraph_reactions_test': 'Every reaction annotated with a Rhea reaction (bqbiol:isVersionOf, or bqbiol:hasPart '
                                'for a lumped one) has Rhea\'s participants (ChEBI terms on its species) and '
                                'stoichiometry (bond counts), allowing the reverse direction and one common factor. A '
                                'difference only in H+ or H2O is a warning. Pending while some reactions are not '
                                'annotated; N/A without a bond graph.',
    'run_test': 'Generates the version alone with libcuflynx (every boundary condition becomes a '
                'parameter, valued from the default parameterisation; a supermodule is expanded into its submodules, '
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
    'validation_test_baseline': 'Compares the model, at the parameterisation\'s parameters, with published or '
                                'experimental baseline data.',
    'validation_test_calibrate': 'Calibrates parameters to the parameterisation\'s obs_data with libcuflynx '
                                 'parameter identification, then checks predictions against held-out data; '
                                 'writes the parameterisation\'s calibrated parameters. Not applicable to a '
                                 'parameterisation without obs_data ("no calibration data in this '
                                 'parameterisation").',
    'version_calibration': 'The version\'s calibration over its parameterisations: failed when no '
                           'parameterisation has calibration data (obs_data), otherwise failed if any '
                           'parameterisation\'s calibration failed and passed when all of them pass. Each '
                           'parameterisation\'s result is in Parameterisations. '
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
    'supermodule_structure_test': 'Its module array (<module_type>_<version>_module_array.json) loads with '
                                  'libcuflynx: every submodule is a library version its config lists in '
                                  'submodule_versions and internal connections name sibling submodules; each '
                                  'parameterisation (nested JSON) gives its instances their own variables, and its '
                                  'globals are declared.',
    'supermodule_equivalence_test': 'A system model using this supermodule version reproduces the system model '
                                    'with its submodules written out, output for output.',
    'coupled_validation_test': 'The version in use in the coupled system models its spec lists (coupled_systems): '
                               '0D modules generated as C++ and coupled to the FEniCS model tissue_diffusion_FEniCS, '
                               'against the same problem with the tissue as a finite-volume grid of CellML cells '
                               '(tests/test_coupled_systems.py). The heatmaps show the mid-plane of the box at a few '
                               'times, in both, and their difference on the grid.',
}
# The test columns of every version, component or supermodule, in this order: verification, the
# version-level calibration, the PhLynx -> CUFLynx pipeline, then the supermodule tests. A column
# that doesn't apply to a version is N/A with the reason (not_applicable_reason); one that applies
# but has no result is "Not run". The validation tests are shown per parameterisation.
SUPERMODULE_TESTS = ['supermodule_structure_test', 'supermodule_equivalence_test']
COUPLED_TESTS = ['coupled_validation_test']
REPORT_TESTS = (checks.VERSION_TESTS + ['version_calibration'] + phlynx.PIPELINE_TESTS + SUPERMODULE_TESTS
                + COUPLED_TESTS + bondgraph.BONDGRAPH_TESTS)
PARAMETERISATION_TESTS = checks.PARAMETERISATION_TESTS
# the test overview's parameterisation columns: each parameterisation's validation tests and its CUFLynx archive
PARAM_COLS = [TEST_SHORT[t] for t in PARAMETERISATION_TESTS] + ['CUFLynx archive']
NOT_A_SUPERMODULE = 'not a supermodule, and no supermodule version uses it'
NOT_COUPLED = 'not used in a coupled system model (the spec lists no coupled_systems)'
NO_EQUIVALENT = ('no system model to reproduce: the spec lists no supermodule.equivalent entry (the version is '
                 'checked through the models that use it)')
STATUS_LABEL = {'passed': 'Passed', 'failed': 'Failed', 'skipped': 'Skipped', 'pending': 'Pending', 'partial': 'Partial',
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


def _not_run_reason(version, test):
    '''Why a test that applies to the version has no result.'''
    from cam_testing.pytest_plugin import QUICK_TESTS
    if not getattr(version, 'reviewed', True) and test not in QUICK_TESTS:
        return ('Not run: the version is not reviewed yet, and CI runs only the quick tests (Run, the '
                'supermodule tests, the PhLynx -> CUFLynx pipeline) for such versions. It runs the full set '
                'once the version is reviewed (locally: pytest --include-unreviewed).')
    return 'This test has not been run yet.'


def _test_entry(version, test, r, spec_block=None, title=None):
    base = test.split('__')[0]
    return {
        'key': test, 'title': title or TEST_TITLES.get(base, test), 'short': TEST_SHORT.get(base, test),
        'about': TEST_ABOUT.get(base, ''),
        'status': r.status if r else None, 'status_label': STATUS_LABEL.get(r.status if r else None, 'Not run'),
        'message': r.message if r else _not_run_reason(version, base),
        'metrics': r.metrics if r else {}, 'plots': r.plots if r else [],
        'details': r.details if r else [], 'timestamp': r.timestamp if r else '',
        'known_issue': (version.spec.get('expected_failures') or {}).get(test),
        'proposed_data': (spec_block or {}).get('status') == checks.PROPOSED,
        'links': [],          # related version pages (version_calibration in a supermodule)
    }


def parameterisation_context(version, inst):
    '''One parameterisation: its files, and its validation (baseline) and calibration results.'''
    tests = []
    for test in PARAMETERISATION_TESTS:
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
    omex_result = checks.load(version, f'cuflynx_parameterisation_omex_test__{inst.name}')
    return {'name': inst.name, 'is_default': inst.is_default, 'files': inst.data_files(),
            # the generated COMBINE archive for CUFLynx (tools/build_parameterisation_omex.py), and its CUFLynx check
            'omex': os.path.relpath(omex_file, version.dir) if os.path.isfile(omex_file) else None,
            'omex_status': omex_result.status if omex_result else None,
            'omex_message': omex_result.message if omex_result else 'not checked yet (tests/test_parameterisation_omex.py)',
            'omex_known_issue': (version.spec.get('expected_failures') or {}).get('cuflynx_parameterisation_omex_test'),
            'n_parameters': len(inst.parameters()), 'obs_data_name': inst.obs_data_name,
            'has_obs_data': inst.has_obs_data, 'calibration': calibration,
            'calibrated_file': os.path.basename(inst.calibrated_parameters_path)
            if os.path.isfile(inst.calibrated_parameters_path) else None,
            'rel_dir': os.path.relpath(inst.dir, version.dir),
            'tests': tests, 'baseline': tests[0], 'calibrate': tests[1]}


def not_applicable_reason(version, test):
    '''Why a report column doesn't apply to the version, or None when it does.'''
    if test in bondgraph.BONDGRAPH_TESTS and not bondgraph.is_bondgraph(version):
        return bondgraph.NOT_A_BONDGRAPH
    if test in SUPERMODULE_TESTS and not version.is_supermodule:
        return NOT_A_SUPERMODULE
    if test == 'supermodule_equivalence_test' and not (version.spec.get('supermodule') or {}).get('equivalent'):
        return NO_EQUIVALENT
    if test in COUPLED_TESTS and not version.spec.get('coupled_systems'):
        return NOT_COUPLED
    if checks.is_non_cellml(version):
        if test in ('verification_test_invariants', 'verification_test_BC', 'verification_test_timestep', 'stability_test'):
            return checks.non_cellml_reason(version)
        if test in phlynx.PIPELINE_TESTS:
            return phlynx.EXTERNAL_REASON if checks.is_external_api(version) else phlynx.CPP_REASON
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
        if test in SUPERMODULE_TESTS and not version.is_supermodule:
            # a component: the result of the supermodule versions that use it, if any
            r = _in_supermodules(version, test) or checks.Result(test, checks.NOT_APPLICABLE, reason)
        elif reason:
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


def _in_supermodules(version, test):
    '''A component's supermodule_structure_test / supermodule_equivalence_test column: the results of
    the supermodule versions that use it as a submodule. "passed_in_super" when all of them pass,
    "failed_in_super" when any fails, pending when some haven't run, None when none has run or no
    supermodule uses it; for the equivalence only the supermodules that have a system model to
    reproduce count.'''
    users = _supermodule_index().get((version.vessel_type, version.name)) or []
    if not users:
        return None
    rows = []
    for s in users:
        if test == 'supermodule_equivalence_test':
            if not (s.spec.get('supermodule') or {}).get('equivalent'):
                continue
            r = _equivalence_result(s)
        else:
            r = checks.load(s, test)
        rows.append({'key': s.key, 'module_type': s.vessel_type, 'version': s.name,
                     'status': r.status if r else None, 'message': r.message if r else 'not run yet'})
    if not rows:
        return checks.Result(test, checks.NOT_APPLICABLE,
                             'not a supermodule; the supermodules using it (' + ', '.join(s.key for s in users)
                             + ') have no system model to reproduce')
    listed = '; '.join(f'{r["key"]}: {r["status"] or "not run"} ({r["message"]})' for r in rows)
    metrics = {'supermodules': rows}
    statuses = [r['status'] for r in rows]
    if any(st in (checks.FAILED, checks.FAILED_IN_SUPER) for st in statuses):
        return checks.Result(test, checks.FAILED_IN_SUPER, f'in its supermodules: {listed}', metrics)
    if all(st in (checks.PASSED, checks.PASSED_IN_SUPER) for st in statuses):
        return checks.Result(test, checks.PASSED_IN_SUPER, f'passes in its supermodules: {listed}', metrics)
    if all(st is None for st in statuses):
        return None
    return checks.Result(test, checks.PENDING, f'not every supermodule using it has run: {listed}', metrics)


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
    'parameterisations': 'Parameterisations: the parameter sets of this version (parameterisations/)',
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
    equations, parameterisations; a supermodule adds parts and levels, its other counts summed over its parts).'''
    if not version.is_supermodule:
        return dict(_leaf_contents(version), parameterisations=len(version.parameterisation_names()))
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
    return dict(total, parameterisations=len(version.parameterisation_names()), parts=len(version.submodules),
                levels=levels, components=leaves)


def contents_bar(version, c, comp_id):
    '''[{key, count, label, about, href}] for the version page's summary bar (and, compact, the module_type page).'''
    sup = version.is_supermodule
    structure = f'#{comp_id}-structure'
    href = {'states': f'#{comp_id}-variables', 'algebraic': f'#{comp_id}-equations',
            'parameters': f'#{comp_id}-variables', 'boundary_conditions': f'#{comp_id}-variables',
            'ports': f'#{comp_id}-ports', 'equations': f'#{comp_id}-equations', 'parameterisations': '#parameterisations',
            'parts': structure, 'levels': structure}
    labels = {'states': ('state', 'states'), 'algebraic': ('algebraic', 'algebraic'),
              'parameters': ('parameter', 'parameters'), 'boundary_conditions': ('boundary condition', 'boundary conditions'),
              'ports': ('port', 'ports'), 'equations': ('equation', 'equations'),
              'parameterisations': ('parameterisation', 'parameterisations'),
              'parts': ('part', 'parts'), 'levels': ('nested level', 'nested levels')}
    short = {'states': 'st', 'algebraic': 'alg', 'parameters': 'par', 'boundary_conditions': 'BC', 'ports': 'ports',
             'equations': 'eq', 'parameterisations': 'psets', 'parts': 'parts', 'levels': 'levels'}
    keys = ['states', 'algebraic', 'parameters', 'boundary_conditions', 'ports', 'equations', 'parameterisations']
    if sup:
        keys = ['parts', 'levels'] + keys
    out = []
    for k in keys:
        n = c.get(k, 0)
        about = CONTENTS_ABOUT[k]
        if sup and k not in ('parameterisations', 'parts', 'levels'):
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
    '''(model, record) for every system-model record using this version (system_models/*/*/*_module_array.json).'''
    import glob
    from cam_testing.library import SYSTEM_MODELS_DIR
    out = []
    for path in sorted(glob.glob(os.path.join(SYSTEM_MODELS_DIR, '*', '*', '*_module_array.json'))):
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


def _same_number(a, b):
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()


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
    # what this supermodule sets for a part: its template's values (its default parameterisation, nested by
    # part) that differ from the part's version's default parameterisation; the template's globals
    overrides, global_rows = {}, []
    try:
        instance_parameters = version.instance_parameters
    except ValueError:
        instance_parameters = {}       # the structure test reports it
    for name, inst in instance_parameters.items():
        sv = versions.get(name)
        defaults = {} if sv is None or sv.is_supermodule else \
            {q.variable_name: q.value for q in sv.default_parameterisation.parameters()}
        for var, q in inst.items():
            if var in defaults and _same_number(defaults[var], q['value']):
                continue
            overrides.setdefault(name, {})[f'{name}/{var}'] = {
                'name': f'{name}/{var}', 'var': var.rpartition('/')[2], 'value': str(q['value']),
                'units': q.get('units', '')}
    for p in version.default_parameterisation.parameters():
        if '/' not in p.variable_name:
            global_rows.append({'name': p.variable_name, 'value': p.value, 'units': p.units})
    overrides = {k: list(v.values()) for k, v in overrides.items()}
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
            'is_supermodule': bool(sv and sv.is_supermodule),
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


def required_citations_context(version):
    """The papers the version's model is built on (config "required_citations"), each linked to its
    bibliography entry, with the best-guess ones marked uncertain; for a supermodule also its
    submodules' (linked to their own version pages)."""
    entries = bib.read(bib.bib_path(version))
    unc = set(version.uncertain_citations)
    own = [{'key': k, 'uncertain': k in unc, 'href': f'#ref-{k}' if k in entries else None,
            'short': _short_citation(entries.get(k))} for k in version.required_citations]
    subs = []
    here = os.path.dirname(version.html_path)
    for sub, k, u in version.submodule_citations():
        if k in version.required_citations:
            continue
        sub_entries = bib.read(bib.bib_path(sub))
        subs.append({'key': k, 'uncertain': u, 'via': sub.key, 'short': _short_citation(sub_entries.get(k)),
                     'href': os.path.relpath(sub.html_path, here) + f'#ref-{k}'})
    # one row per key: a key several submodules cite is listed once, naming them all
    merged = {}
    for r in subs:
        m = merged.setdefault(r['key'], {**r, 'via': []})
        m['via'].append(r['via'])
        m['uncertain'] = m['uncertain'] and r['uncertain']
    return {'own': own, 'submodules': list(merged.values())}


def _short_citation(fields):
    """'Belluzzi & Sacchi 1986' style author-year label from a bib entry (None when there is none)."""
    if not fields:
        return None
    names = [a.strip() for a in (fields.get('author') or '').split(' and ') if a.strip()]
    surname = lambda a: a.split(',')[0].strip() if ',' in a else a.split()[-1]   # noqa: E731
    if not names:
        who = fields.get('title', '')[:40]
    elif len(names) == 1 or names[-1] == 'others' and len(names) == 2:
        who = surname(names[0]) + (' et al.' if len(names) == 2 else '')
    elif len(names) == 2:
        who = f'{surname(names[0])} & {surname(names[1])}'
    else:
        who = f'{surname(names[0])} et al.'
    return f"{who} {fields.get('year', 'n.d.')}".strip()


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
        'required_citations': required_citations_context(version),
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
        'tests': tests,
        'parameterisations': [parameterisation_context(version, i) for i in version.parameterisations()],
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
    counts = {'passed': 0, 'partial': 0, 'failed': 0, 'skipped': 0, 'pending': 0, 'not_applicable': 0, None: 0}
    for c in components:
        for t in c['tests'] + [t for i in c['parameterisations'] for t in i['tests']]:
            status = COUNT_AS.get(t['status'], t['status'])
            counts[status] = counts.get(status, 0) + 1
        # each parameterisation's CUFLynx archive check (shown in the overview's parameterisation table)
        for i in c['parameterisations']:
            if i.get('omex_status') is not None:
                counts[i['omex_status']] = counts.get(i['omex_status'], 0) + 1
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


def celldl_editor_url():
    '''The CellDL editor the version pages link to: $CAM_CELLDL_EDITOR_URL, else CELLDL_EDITOR_URL.'''
    return os.environ.get(CELLDL_EDITOR_ENV) or CELLDL_EDITOR_URL


def celldl_open(name):
    '''The editor URL a CellDL document's base64 follows: <editor>?open=celldl&name=<name>#.'''
    return f'{celldl_editor_url()}?open=celldl&name={urllib.parse.quote(name)}#'


def open_in_context(version):
    '''The version page's "Open in" links: {celldl, phlynx, cuflynx}. celldl (a bond graph only) and
    phlynx are {href: the file relative to the page, open: the app URL its base64 follows}; phlynx is
    {reason} where PhLynx has no form of the version and {missing} before an export is kept.
    cuflynx is the URL prefix of the per-parameterisation links. celldl also has open_by_url: the
    editor URL the document's absolute, percent-encoded URL follows, for a document too large for a
    link (the server must allow the editor's origin to fetch it: CORS).'''
    celldl = None
    if bondgraph.is_bondgraph(version):
        href = os.path.relpath(bondgraph.document_path(version), version.dir).replace(os.sep, '/')
        celldl = {'href': href, 'open': celldl_open(os.path.basename(href)), 'open_by_url': celldl_editor_url() + '?open='}
        if os.environ.get(CELLDL_NOTE_ENV):
            celldl['note'] = os.environ[CELLDL_NOTE_ENV]
    archive = phlynx.archive_path(version)
    if phlynx._cpp(version):
        phl = {'reason': phlynx._reason(version)}
    elif version.is_supermodule:
        phl = {'reason': phlynx.NO_SUPERMODULE_SUPPORT}
    elif os.path.isfile(archive):
        phl = {'href': os.path.relpath(archive, version.dir).replace(os.sep, '/'), 'open': PHLYNX_OPEN}
    else:
        phl = {'missing': PHLYNX_NOT_EXPORTED}
    return {'celldl': celldl, 'phlynx': phl, 'cuflynx': CUFLYNX_OPEN}


def version_context(version):
    comp = component_context(version)
    return {
        'open_in': open_in_context(version),
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
        'param_cols': PARAM_COLS,
        'phlynx': _module_phlynx([comp]), 'bc_summary': _bc_summary(comp),
        'licence': version.licence, 'licence_url': LICENCES.get(version.licence), 'creators': version.creators,
        'unit_failures': (comp['unit_consistency'] or {}).get('n_failures', 0),
        'module_type_href': f'../../{version.vessel_type}.html',
        'generated': _now(), 'libcuflynx_version': _libcuflynx_version(),
    }


def build_version(version, out_path=None):
    '''Writes the version's page (default its html_path; out_path: the same page elsewhere, e.g. in the
    site for an extra library's version); returns (path, context).'''
    ctx = version_context(version)
    out_path = out_path or version.html_path
    with open(out_path, 'w') as f:
        f.write(_env().get_template('version.html').render(**ctx))
    return out_path, ctx


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
                     'parameterisations': comp['parameterisations'], 'submodules': comp['submodules'],
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
        'n_parameterisations': sum(len(r['parameterisations']) for r in rows),
        'all_sourced': all(r['all_sourced'] for r in rows),
        'n_sourced': sum(r['n_sourced'] for r in rows), 'n_parameters': sum(r['n_parameters'] for r in rows),
        'max_risk': max((r['max_risk'] for r in rows if r['max_risk'] is not None), default=None),
        'known_issues': sum(r['known_issues'] for r in rows),
        'phlynx': {'run': sum(r['phlynx']['run'] for r in rows),
                   'compatible': sum(r['phlynx']['compatible'] for r in rows),
                   'total': len(rows),
                   'all': bool(rows) and all(r['phlynx']['all'] for r in rows)},
        'ready_for_review': any(c.get('review') and not c['reviewed'] for c in vctx),
        'test_keys': REPORT_TESTS, 'test_short': TEST_SHORT, 'param_cols': PARAM_COLS,
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


def build_index(contexts, out_path, href, systems=()):
    '''The site index: the system models (systems: system_rows), then module_types grouped by category (a
    module_type directly under modules/ and the ones nested in it form their own group), each with links
    to its versions. A nested module_type is listed by its path in the group, e.g. neuron/soma under cell.'''
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
                         'n_versions': c['n_versions'], 'n_parameterisations': c['n_parameterisations'],
                         'counts': c['counts'],
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
              for k in ('passed', 'partial', 'failed', 'skipped', 'pending', 'not_applicable', None)}
    html_text = _env().get_template('index.html').render(
        groups=out_groups, systems=list(systems), totals=totals, n_module_types=len(all_rows),
        n_versions=sum(r['n_versions'] for r in all_rows),
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
        q['counts'] = {k: sum(v['counts'].get(k, 0) for v in vs) for k in ('passed', 'partial', 'failed')}
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


# ---- system models ------------------------------------------------------------------------------
#
# Every system model (system_models/<category>/<model>/) gets a page, <model>.html beside it (and
# site/system_models/<category>/<model>/<model>.html), generated from its module array, parameters,
# spec and test results. Its figure is a PhLynx screenshot when the model has one, <model>_phlynx.png
# with <model>_phlynx_nodes.json giving each node's instance and pixel box in the image (the page
# overlays a link per node to that instance's version page); otherwise a diagram generated from the
# module array's connections.

SYSTEM_TESTS = ['system_run_test', 'system_invariants_test', 'system_equivalence_test']
SYSTEM_TEST_TITLES = {'system_run_test': 'Run', 'system_invariants_test': 'Verification: invariants',
                      'system_equivalence_test': 'Equivalence: reproduces its reference'}
SYSTEM_TEST_ABOUT = {
    'system_run_test': 'The model generates with libcuflynx from this library and runs; every output is finite.',
    'system_invariants_test': 'The spec\'s invariants (conservation laws, bounds) hold over the run.',
    'system_equivalence_test': 'The model reproduces its reference (circulatory_autogen\'s original, or a ready '
                               'CellML model such as a PhLynx export): every output the reference logs agrees '
                               'to the spec\'s equivalence.tol.',
}
SYSTEM_FIGURE_SUFFIX = '_phlynx.png'
SYSTEM_NODES_SUFFIX = '_phlynx_nodes.json'


def system_page_path(system):
    '''The system model's page, <model>.html beside its module array.'''
    return os.path.join(system.dir, f'{system.name}.html')


def system_figure_paths(system):
    '''(<model>_phlynx.png, <model>_phlynx_nodes.json) of a system model.'''
    return (os.path.join(system.dir, system.name + SYSTEM_FIGURE_SUFFIX),
            os.path.join(system.dir, system.name + SYSTEM_NODES_SUFFIX))


def _png_size(path):
    with open(path, 'rb') as f:
        head = f.read(24)
    if head[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError(f'{path}: not a PNG')
    return int.from_bytes(head[16:20], 'big'), int.from_bytes(head[20:24], 'big')


def system_figure(system, instance_names):
    '''The PhLynx screenshot and its node boxes, or None when the model has no <model>_phlynx.png:
    {src, width, height, nodes: [{instance, left, top, width, height} in % of the image], errors}.
    <model>_phlynx_nodes.json is {"width", "height", "source" (optional: what the screenshot shows),
    "nodes": [{"instance", "box": [x, y, w, h]}]} in the image's pixels; errors list nodes naming no instance or lying outside the image.'''
    png, nodes_path = system_figure_paths(system)
    if not os.path.isfile(png):
        return None
    width, height = _png_size(png)
    out = {'src': os.path.basename(png), 'width': width, 'height': height, 'nodes': [], 'errors': []}
    if not os.path.isfile(nodes_path):
        out['errors'].append(f'no {os.path.basename(nodes_path)}: the figure has no clickable nodes')
        return out
    with open(nodes_path) as f:
        doc = json.load(f)
    out['source'] = doc.get('source')
    if (doc.get('width'), doc.get('height')) != (width, height):
        out['errors'].append(f'{os.path.basename(nodes_path)} is for a {doc.get("width")}x{doc.get("height")} '
                             f'image; {os.path.basename(png)} is {width}x{height}')
    for node in doc.get('nodes') or []:
        name, box = node.get('instance'), node.get('box') or []
        if name not in instance_names:
            out['errors'].append(f'node {name!r}: no instance of that name in the module array')
            continue
        if len(box) != 4 or box[2] <= 0 or box[3] <= 0 or box[0] < 0 or box[1] < 0 \
                or box[0] + box[2] > width or box[1] + box[3] > height:
            out['errors'].append(f'node {name!r}: box {box} is not inside the {width}x{height} image')
            continue
        x, y, w, h = box
        out['nodes'].append({'instance': name, 'left': 100 * x / width, 'top': 100 * y / height,
                             'width': 100 * w / width, 'height': 100 * h / height})
    return out


def _system_records(system):
    with open(os.path.join(system.dir, f'{system.name}_module_array.json')) as f:
        rows = json.load(f)
    return [r for r in rows if isinstance(r, dict) and r.get('name')]


def system_connections(records):
    '''The module array's connections, [(from, to)] in order: out_instances and inp_instances, and a
    supermodule instance's per_submodule_outputs / per_submodule_inputs as <instance>/<submodule>.'''
    edges = []

    def add(a, b):
        if (a, b) not in edges:
            edges.append((a, b))
    for r in records:
        for dst in r.get('out_instances') or []:
            add(r['name'], dst)
        for src in r.get('inp_instances') or []:
            add(src, r['name'])
        for sub, dsts in (r.get('per_submodule_outputs') or {}).items():
            for dst in dsts:
                add(f'{r["name"]}/{sub}', dst)
        for sub, srcs in (r.get('per_submodule_inputs') or {}).items():
            for src in srcs:
                add(src, f'{r["name"]}/{sub}')
    return edges


def _parameter_groups(path):
    '''A nested parameters file as [{name, rows: [{name, value, units, reference}]}], one group per top key.'''
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        doc = json.load(f)
    groups = []
    for top, node in doc.items():
        if not isinstance(node, dict):
            continue
        rows = []

        def walk(n, prefix):
            for k, v in n.items():
                if isinstance(v, dict) and 'value' in v:
                    rows.append({'name': '/'.join(prefix + [k]), 'value': v.get('value'), 'units': v.get('units', ''),
                                 'reference': v.get('data_reference', '')})
                elif isinstance(v, dict):
                    walk(v, prefix + [k])
        walk(node, [])
        groups.append({'name': top, 'rows': rows})
    return groups


def _system_tests(system):
    out = []
    for test in SYSTEM_TESTS:
        path = os.path.join(system.results_dir, f'{test}.json')
        r = {}
        if os.path.isfile(path):
            with open(path) as f:
                r = json.load(f)
        known = (system.spec.get('expected_failures') or {}).get(test)
        out.append({'key': test, 'title': SYSTEM_TEST_TITLES[test], 'about': SYSTEM_TEST_ABOUT[test],
                    'status': r.get('status'), 'message': r.get('message', 'not run yet'),
                    'timestamp': r.get('timestamp'), 'details': r.get('details') or [],
                    'rows': ((r.get('metrics') or {}).get('rows') or [])[:30], 'known_issue': known})
    return out


def _system_plots(system):
    if not os.path.isdir(system.plots_dir):
        return []
    return [f'plots/{f}' for f in sorted(os.listdir(system.plots_dir)) if f.lower().endswith(('.png', '.svg'))]


def system_instance(record, href_of):
    '''One module-array record: its version (linked by href_of(version)) and, for a supermodule, its
    submodules (linked too). missing says why a record names no version of the libraries.'''
    from cam_testing.library import load_version
    out = {'name': record['name'], 'module_type': record.get('module_type'), 'version': record.get('module_subtype'),
           'parameterisation': record.get('parameterisation', ''), 'href': None, 'missing': None,
           'is_supermodule': False, 'submodules': [],
           'inputs': list(record.get('inp_instances') or []), 'outputs': list(record.get('out_instances') or [])}
    try:
        v = load_version(record['module_type'], record['module_subtype'])
    except (OSError, ValueError, KeyError) as e:
        out['missing'] = f'not a version of the module libraries ({type(e).__name__})'
        return out
    out['href'] = href_of(v)
    out['is_supermodule'] = v.is_supermodule
    out['in_primary'] = v.mtype.is_primary
    for sub in v.submodules if v.is_supermodule else []:
        sv = _sub_version(sub)
        out['submodules'].append({'name': sub['name'], 'key': f'{sub["module_type"]}/{sub["module_subtype"]}',
                                  'href': href_of(sv) if sv else None})
    return out


def _system_mermaid(instances, edges):
    '''The fallback figure: a flowchart of the instances and their connections, each node linking to its page.'''
    ids = {i['name']: _mermaid_id('n_', i['name']) for i in instances}
    lines = ['flowchart LR']
    for i in instances:
        lines.append(f'  {ids[i["name"]]}["{_mermaid_text(i["name"])}<br/><small>{_mermaid_text(i["module_type"] or "")}'
                     f' {_mermaid_text(i["version"] or "")}</small>"]')
    for a, b in edges:
        ia, ib = ids.get(a.split('/')[0]), ids.get(b.split('/')[0])
        if ia and ib:
            label = '/'.join(p for p in (a.partition('/')[2], b.partition('/')[2]) if p)
            lines.append(f'  {ia} -->|"{_mermaid_text(label)}"| {ib}' if label else f'  {ia} --> {ib}')
    for i in instances:
        if i['href']:
            lines.append(f'  click {ids[i["name"]]} href "{i["href"]}" "{_mermaid_text(i["module_type"] + "/" + i["version"])}"')
    return '\n'.join(lines)


def system_context(system, href_of=None, page_dir=None):
    '''The system page's context. href_of(version) is the version page's href from the page (default:
    relative to page_dir, the system model's directory).'''
    page_dir = page_dir or system.dir
    href_of = href_of or (lambda v: os.path.relpath(v.html_path, page_dir).replace(os.sep, '/'))
    records = _system_records(system)
    instances = [system_instance(r, href_of) for r in records]
    edges = system_connections(records)
    spec = system.spec
    notes = spec.get('notes') or []
    notes = [notes] if isinstance(notes, str) else list(notes)
    figure = system_figure(system, {i['name'] for i in instances})
    hrefs = {i['name']: i['href'] for i in instances}
    for n in (figure or {}).get('nodes') or []:
        n['href'] = hrefs.get(n['instance'])
        n['key'] = next(f'{i["module_type"]}/{i["version"]}' for i in instances if i['name'] == n['instance'])
    tests = _system_tests(system)
    counts = {}
    for t in tests:
        counts[t['status']] = counts.get(t['status'], 0) + 1
    return {
        'name': system.name, 'category': system.category, 'id': system.id,
        'description': spec.get('description') or (notes[0] if notes else ''),
        'notes': notes if spec.get('description') else notes[1:],
        'reviewed': bool(spec.get('reviewed')), 'skip': spec.get('skip'),
        'known_issues': [f'{t}: {why}' for t, why in (spec.get('expected_failures') or {}).items()],
        'sim': {k: spec.get(k) for k in ('sim_time', 'pre_time', 'dt') if spec.get(k) is not None},
        'instances': instances, 'edges': edges, 'figure': figure,
        'mermaid': None if figure else _system_mermaid(instances, edges),
        'parameters': _parameter_groups(os.path.join(system.dir, f'{system.name}_parameters.json')),
        'parameters_file': f'{system.name}_parameters.json',
        'tests': tests, 'counts': counts, 'plots': _system_plots(system),
        'generated': _now(), 'libcuflynx_version': _libcuflynx_version(),
    }


def build_system(system, out_path=None, href_of=None, index_href=None):
    '''Writes the system model's page (default <model>.html beside it; index_href: its link to the
    site index, in the site); returns (path, context).'''
    out_path = out_path or system_page_path(system)
    ctx = system_context(system, href_of, os.path.dirname(out_path))
    ctx['index_href'] = index_href
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(_env().get_template('system_model.html').render(**ctx))
    return out_path, ctx


def system_rows(contexts, href):
    '''The site index's "System models" rows, from the system pages' contexts.'''
    return [{'id': c['id'], 'name': c['name'], 'category': c['category'], 'href': href(c),
             'n_instances': len(c['instances']), 'counts': c['counts'], 'description': c['description'],
             'has_figure': bool(c['figure']), 'tests': c['tests']}
            for c in sorted(contexts, key=lambda c: c['id'].lower())]


def _href(relpath, page):
    return f'modules/{relpath}/{page}'


def _copy_version_files(v, vdest):
    '''A version's plots, bond graph, PhLynx export and parameterisation archives, beside its page in the site.'''
    if os.path.isdir(v.plots_dir):
        shutil.copytree(v.plots_dir, os.path.join(vdest, 'plots'), dirs_exist_ok=True)
    # the files the page opens in the CellDL editor and PhLynx, at the same relative paths
    for path in (bondgraph.document_path(v), phlynx.archive_path(v)):
        if os.path.isfile(path):
            fdest = os.path.join(vdest, os.path.relpath(path, v.dir))
            os.makedirs(os.path.dirname(fdest), exist_ok=True)
            shutil.copy2(path, fdest)
    # the parameterisations' CUFLynx archives, where they have been built (make omex)
    for inst in v.parameterisations():
        archive = omex_mod.omex_path(v, inst)
        if os.path.isfile(archive):
            idest = os.path.join(vdest, os.path.relpath(inst.dir, v.dir))
            os.makedirs(idest, exist_ok=True)
            shutil.copy2(archive, idest)


def _site_version_path(site_dir, v):
    return os.path.join(site_dir, 'modules', v.mtype.relpath, 'versions', v.name, f'{v.stem}.html')


def _extra_versions(systems, names):
    '''The extra libraries' versions the site's pages link to: those the system models use, and (for a
    supermodule, at any depth) their submodules, and the submodules of this library's supermodules.'''
    from cam_testing.library import load_version
    seeds = []
    for s in systems:
        for r in _system_records(s):
            try:
                seeds.append(load_version(r['module_type'], r['module_subtype']))
            except (OSError, ValueError, KeyError):
                pass
    for name in names:
        for v in load_module_type(name).versions():
            seeds += [sv for sv in (_sub_version(sub) for sub in v.submodules) if sv] if v.is_supermodule else []
    out, todo = {}, list(seeds)
    while todo:
        v = todo.pop()
        if v.key in out:
            continue
        if not v.mtype.is_primary:
            out[v.key] = v
        if v.is_supermodule:
            todo += [sv for sv in (_sub_version(sub) for sub in v.submodules) if sv]
    return sorted(out.values(), key=lambda v: v.key)


def assemble_site(names, contexts, systems=()):
    '''site/index.html + site/modules/<module_type path>/ (its page, and each version's page, plots,
    bond graph, PhLynx export and instance archives; a nested module_type's directory is inside its
    parent's), as GitHub Pages serves it; and site/system_models/<category>/<model>/ for each of
    systems (its page and PhLynx screenshot). The extra libraries' versions the pages link to get their
    pages under site/modules/ too, at their own library's paths.'''
    site_dir = paths.roots().site_dir
    if os.path.isdir(site_dir):
        shutil.rmtree(site_dir)
    for name in names:
        mtype = load_module_type(name)
        dest = os.path.join(site_dir, 'modules', mtype.relpath)
        os.makedirs(dest, exist_ok=True)
        if not os.path.isfile(mtype.html_path) or any(not os.path.isfile(v.html_path) for v in mtype.versions()):
            build_module(name)
        shutil.copy2(mtype.html_path, dest)
        for v in mtype.versions():
            vdest = os.path.join(dest, 'versions', v.name)
            os.makedirs(vdest, exist_ok=True)
            shutil.copy2(v.html_path, vdest)
            _copy_version_files(v, vdest)
    extra = _extra_versions(systems, names)
    by_type = {}
    for v in extra:
        out = _site_version_path(site_dir, v)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        by_type.setdefault(v.mtype.name, []).append(build_version(v, out)[1])
        _copy_version_files(v, os.path.dirname(out))
    for name, vctx in by_type.items():
        mtype = load_module_type(name)
        page = os.path.join(site_dir, 'modules', mtype.relpath, f'{name}.html')
        if not os.path.isfile(page):        # a module_type of this library already has its page
            os.makedirs(os.path.dirname(page), exist_ok=True)
            with open(page, 'w') as f:
                f.write(_env().get_template('module_type.html').render(**module_context(name, vctx)))
    system_ctx = []
    for s in systems:
        sdest = os.path.join(site_dir, 'system_models', s.category, s.name)
        os.makedirs(sdest, exist_ok=True)
        png = system_figure_paths(s)[0]
        if os.path.isfile(png):
            shutil.copy2(png, sdest)
        if os.path.isdir(s.plots_dir):
            shutil.copytree(s.plots_dir, os.path.join(sdest, 'plots'))
        system_ctx.append(build_system(
            s, os.path.join(sdest, f'{s.name}.html'),
            href_of=lambda v, d=sdest: os.path.relpath(_site_version_path(site_dir, v), d).replace(os.sep, '/'),
            index_href=os.path.relpath(os.path.join(site_dir, 'index.html'), sdest).replace(os.sep, '/'))[1])
    build_review_queue(contexts, os.path.join(site_dir, 'review_queue.html'), _href)
    return build_index(contexts, os.path.join(site_dir, 'index.html'), _href,
                       system_rows(system_ctx, lambda c: f'system_models/{c["category"]}/{c["name"]}/{c["name"]}.html'))


SHARE_README = '''{title}

A copy of the module library's review site, to read offline.

1. Unzip this file.
2. Open site/index.html in a web browser.
3. {start}

Every page links to the others by relative paths, so keep the folder as it is. Equations, fonts and
diagrams load from the web when you are online; the pages work without them.
'''


def zip_site(zip_path, title='Module library review site', start=None):
    '''Zips site/ (assembled first) as site/... with a README.txt at the top; returns zip_path.'''
    import zipfile
    site_dir = paths.roots().site_dir
    start = start or 'Click a system model under "System models", or a module type, to open its page.'
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('README.txt', SHARE_README.format(title=title, start=start))
        for root, _, files in os.walk(site_dir):
            for f in sorted(files):
                path = os.path.join(root, f)
                z.write(path, os.path.join('site', os.path.relpath(path, site_dir)))
    return zip_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--module', action='append', default=[],
                        help='module_type (with those nested in it) or path under modules/ (repeatable); default all')
    parser.add_argument('--site', action='store_true', help='also assemble site/ for GitHub Pages')
    parser.add_argument('--celldl-note', metavar='TEXT',
                        help=f'a note shown with every bond graph and its CellDL link (${CELLDL_NOTE_ENV})')
    parser.add_argument('--zip', metavar='PATH', help='assemble site/ and zip it, with a README.txt, to send')
    parser.add_argument('--zip-title', default='Module library review site', help="the zip README's title")
    parser.add_argument('--zip-start', help="the zip README's third step (what to open first)")
    args = parser.parse_args(argv)
    if args.celldl_note:
        os.environ[CELLDL_NOTE_ENV] = args.celldl_note
    from cam_testing import system as systems_mod
    contexts = {}
    _index_cache.clear()
    _contents_cache.clear()
    for name in select_module_types(args.module):
        out, ctx = build_module(name)
        contexts[name] = ctx
        print(f'wrote {os.path.relpath(out, paths.roots().repo_root)} and {ctx["n_versions"]} version page(s)')
    # every system model gets its page, as every version does
    systems = systems_mod.system_models()
    for s in systems:
        build_system(s)
    print(f'wrote {len(systems)} system model page(s)')
    if args.site or args.zip:
        # the index always covers every module_type, even when only some pages were rebuilt
        for name in module_type_names():
            if name not in contexts:
                vctx = [version_context(v) for v in load_module_type(name).versions()]
                contexts[name] = dict(module_context(name, vctx), version_contexts=vctx)
        out = assemble_site(module_type_names(), [contexts[n] for n in module_type_names()], systems)
        print(f'wrote {os.path.relpath(out, paths.roots().repo_root)}')
    if args.zip:
        zip_site(args.zip, args.zip_title, args.zip_start)
        print(f'wrote {args.zip} ({os.path.getsize(args.zip) / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
