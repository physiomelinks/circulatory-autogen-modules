"""
HTML reports: one page per module (modules/<name>/<name>.html) and a site index.

    python -m cam_testing.report                  # every module + site/index.html
    python -m cam_testing.report --module heart   # one module
    python -m cam_testing.report --site           # also assemble site/ as GitHub Pages serves it

Pages read the result JSON and plots the tests wrote, so a report shows the last local (or
CI) test run. Plots are referenced relatively (plots/...), so a page works opened from disk.
"""
import argparse
import datetime
import html
import json
import os
import shutil

from jinja2 import Environment, FileSystemLoader, select_autoescape

from cam_testing import bib, checks, phlynx, ranges, risk
from cam_testing.library import REPO_ROOT, load_module, module_names
from cam_testing.mathml import component_equations, component_variables

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
    'phlynx_export_test': 'PhLynx: build & export',
    'cuflynx_simulate_test': 'CUFLynx: import & simulate',
    'phlynx_equivalence_test': 'PhLynx → CUFLynx vs libcuflynx',
}
TEST_SHORT = {
    'run_test': 'Run', 'verification_test_invariants': 'Invariants', 'verification_test_BC': 'BC sweep', 'verification_test_timestep': 'Timestep',
    'stability_test': 'Stability', 'validation_test_baseline': 'Baseline', 'validation_test_calibrate': 'Calibrate',
    'phlynx_export_test': 'PhLynx', 'cuflynx_simulate_test': 'CUFLynx', 'phlynx_equivalence_test': 'PhLynx≡',
}
TEST_ABOUT = {
    'run_test': 'Generates the component alone with libcuflynx (every boundary condition becomes a '
                'parameter), simulates it, and checks every output is finite.',
    'verification_test_invariants': 'Checks the simulation against what the component is supposed to do: '
                                    'the invariants in its spec (exact solutions, conservation laws, bounds, '
                                    'delays) evaluated at the run parameters.',
    'verification_test_BC': 'Sweeps each boundary condition (or, for a self-contained component, each '
                            'constant) over a range and checks every run completes, stays finite and '
                            'keeps the invariants.',
    'verification_test_timestep': 'Integrates the generated right-hand side with a fixed-step scheme at '
                                  'successively halved steps. Differences between successive solutions must '
                                  'shrink at the scheme\'s order, and the finest solution must agree with '
                                  'libcuflynx\'s CVODE run at tight tolerances.',
    'stability_test': 'Runs the component with a matrix of solvers, tolerances and timesteps. A '
                      'configuration works when it finishes with finite outputs within the tolerance of a '
                      'tight-tolerance reference. Declared-supported configurations must work.',
    'validation_test_baseline': 'Compares the model with published or experimental baseline data.',
    'validation_test_calibrate': 'Calibrates parameters to one data set with libcuflynx parameter '
                                 'identification, then checks predictions against held-out data.',
    'phlynx_export_test': 'Builds the component\'s test network in PhLynx (its own code, loaded with this '
                          'library\'s modules and parameters), checks every connection was made, and exports '
                          'the .omex PhLynx sends to CUFLynx.',
    'cuflynx_simulate_test': 'Imports that .omex into a released CUFLynx through its API and simulates it; '
                             'every output must be finite.',
    'phlynx_equivalence_test': 'Compares CUFLynx\'s simulation of the PhLynx-built model with libcuflynx\'s '
                               'model of the same network (outputs matched by instance and variable, '
                               'normalised difference within 1e-6).',
}
# The tests shown for every component: the V&V tests, then the PhLynx -> CUFLynx pipeline
REPORT_TESTS = checks.TESTS + phlynx.PIPELINE_TESTS
STATUS_LABEL = {'passed': 'Passed', 'failed': 'Failed', 'skipped': 'Skipped', 'pending': 'Pending',
                'not_applicable': 'N/A', None: 'Not run'}


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


def component_context(component):
    module = component.module
    equations, unsupported = ([], set())
    variables = []
    if component.config.get('module_format', 'cellml') == 'cellml':
        # the component may live in another module's file (shared components)
        cellml = module.cellml_path
        equations, unsupported = component_equations(cellml, component.module_type)
        if not equations:
            for other in module_names():
                path = load_module(other).cellml_path
                eqs, uns = component_equations(path, component.module_type)
                if eqs:
                    cellml, equations, unsupported = path, eqs, uns
                    break
        cellml_vars = {v[0]: v for v in component_variables(cellml, component.module_type)}
    else:
        cellml_vars = {}
    params = {p.variable_name: p for p in component.parameters()}
    spread = ranges.validated_spread(component)
    for name, units, access, kind in component.config['variables_and_units']:
        cv = cellml_vars.get(name)
        p = params.get(name)
        variables.append({
            'name': name, 'units': units, 'kind': kind.strip(), 'access': access,
            'initial': cv[2] if cv else None,
            'value': (p.value + ' (proposed)' if p.proposed else p.value) if p else None, 'todo': bool(p and (p.is_todo or p.proposed)),
            'reference': p.data_reference if p else '',
            'sourced': p.is_sourced if p else None,
            'tested': ranges.tested_range(component, p) if p else None,
            'validated': spread.get(name),
        })

    tests = []
    validation_spec = component.spec.get('validation') or {}
    for test in REPORT_TESTS:
        r = checks.load(component, test)
        if r is None and test.startswith('validation_test_'):
            # not run (e.g. slow tests excluded): a skipped/pending spec still says why
            v = validation_spec.get(test.rsplit('_', 1)[1]) or {}
            if v.get('status', checks.PENDING) in (checks.PENDING, checks.SKIPPED, checks.NOT_APPLICABLE):
                r = checks.Result(test, v.get('status', checks.PENDING), v.get('reason', 'no validation data chosen yet'))
        tests.append({
            'key': test, 'title': TEST_TITLES[test], 'short': TEST_SHORT[test], 'about': TEST_ABOUT[test],
            'status': r.status if r else None, 'status_label': STATUS_LABEL[r.status if r else None],
            'message': r.message if r else 'This test has not been run for this component yet.',
            'metrics': r.metrics if r else {}, 'plots': r.plots if r else [],
            'details': r.details if r else [], 'timestamp': r.timestamp if r else '',
            'known_issue': (component.spec.get('expected_failures') or {}).get(test),
            'proposed_data': test.startswith('validation_test_') and
                             (validation_spec.get(test.rsplit('_', 1)[1]) or {}).get('status') == checks.PROPOSED,
        })
    risk_result = risk.load(component)
    if risk_result is not None and not all(os.path.isfile(os.path.join(module.dir, p)) for p in risk_result.get('plots', [])):
        risk_result['plots'] = risk.make_plots(component, risk_result)   # e.g. on a fresh clone
    proposals = component.spec.get('reference_proposals') or {}
    bib_keys = set(bib.read(bib.bib_path(module))) | set(bib.read(bib.proposed_bib_path(module)))
    references = []
    for p in component.parameters():
        key = bib.reference_key(p.data_reference)
        references.append({'name': p.variable_name, 'value': p.value + (' (proposed)' if p.proposed else ''), 'units': p.units,
                           'reference': p.data_reference, 'key': key if key in bib_keys else None,
                           'note': p.data_reference.split(';', 1)[1].strip() if ';' in p.data_reference else '',
                           'sourced': p.is_sourced, 'proposal': proposals.get(p.variable_name)})
    return {
        'risk': risk_result, 'references': references, 'has_proposals': bool(proposals),
        'id': component.id, 'label': component.label, 'vessel_type': component.vessel_type,
        'BC_type': component.BC_type, 'module_type': component.module_type,
        'format': component.config.get('module_format', 'cellml'),
        'skip': component.spec.get('skip'), 'notes': component.spec.get('notes'),
        'sweep_rationale': (component.spec.get('bc_sweep') or {}).get('rationale'),
        'equations': equations, 'unsupported': sorted(unsupported),
        'ports': _ports(component.config), 'variables': variables,
        'todo': component.todo_parameters(), 'unsourced': component.unsourced_parameters(), 'invariants': component.spec.get('invariants') or [],
        'validation': component.spec.get('validation') or {},
        'tests': tests,
        'phlynx_compatible': _phlynx_compatible(tests),
    }


def _phlynx_compatible(tests):
    '''True when PhLynx builds and exports it, CUFLynx simulates it and it matches libcuflynx; None if not run.'''
    status = {t['key']: t['status'] for t in tests if t['key'] in phlynx.PIPELINE_TESTS}
    if all(v is None for v in status.values()):
        return None
    return all(v == checks.PASSED for v in status.values())


def _status_counts(components):
    counts = {'passed': 0, 'failed': 0, 'skipped': 0, 'pending': 0, 'not_applicable': 0, None: 0}
    for c in components:
        for t in c['tests']:
            counts[t['status']] = counts.get(t['status'], 0) + 1
    return counts


def _structure(module):
    path = os.path.join(module.results_dir, 'structure.json')
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


def _env():
    env = Environment(loader=FileSystemLoader(TEMPLATES), autoescape=select_autoescape(['html']))
    env.filters['fmt'] = _fmt
    return env


def module_context(name):
    module = load_module(name)
    components = [component_context(c) for c in module.components()]
    return {
        'name': name, 'reviewed': module.reviewed, 'components': components,
        'all_sourced': module.all_sourced(),
        'bibliography': [{'key': k, 'text': bib.format_entry(f), 'link': bib.link(f), 'proposed': False}
                         for k, f in bib.read(bib.bib_path(module)).items()]
                        + [{'key': k, 'text': bib.format_entry(f), 'link': bib.link(f), 'proposed': True}
                           for k, f in bib.read(bib.proposed_bib_path(module)).items()],
        'review': module.spec.get('review'),
        'bib_keys': sorted(set(bib.read(bib.bib_path(module))) | set(bib.read(bib.proposed_bib_path(module)))),
        'bib_file': os.path.basename(bib.bib_path(module)),
        'n_sourced': sum(p.is_sourced for p in module.parameters), 'n_parameters': len(module.parameters),
        'counts': _status_counts(components), 'structure': _structure(module),
        'max_risk': max((c['risk']['failure_probability'] for c in components if c['risk']), default=None),
        'known_issues': (module.spec.get('known_issues') or [])
                        + [f'{c.id}: {t}: {why}' for c in module.components()
                           for t, why in (c.spec.get('expected_failures') or {}).items()],
        'test_keys': REPORT_TESTS, 'test_short': TEST_SHORT,
        'phlynx': _module_phlynx(components),
        'generated': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC'),
        'libcuflynx_version': _libcuflynx_version(),
    }


def _module_phlynx(components):
    run = [c for c in components if c['phlynx_compatible'] is not None]
    ok = sum(bool(c['phlynx_compatible']) for c in run)
    return {'run': len(run), 'compatible': ok, 'total': len(components),
            'all': bool(run) and len(run) == len(components) and ok == len(components)}


def build_module(name):
    ctx = module_context(name)
    out = os.path.join(load_module(name).dir, f'{name}.html')
    with open(out, 'w') as f:
        f.write(_env().get_template('module.html').render(**ctx))
    return out, ctx


def build_index(contexts, out_path, module_href):
    rows = [{'name': c['name'], 'href': module_href(c['name']), 'n_components': len(c['components']),
             'counts': c['counts'], 'reviewed': c['reviewed'],
             'known_issues': len(c['known_issues']), 'max_risk': c.get('max_risk'),
             'all_sourced': c['all_sourced'], 'n_sourced': c['n_sourced'], 'n_parameters': c['n_parameters'],
             'phlynx': c['phlynx'],
             'ready_for_review': bool(c.get('review')) and not c['reviewed']}
            for c in contexts]
    totals = {k: sum(r['counts'].get(k, 0) for r in rows) for k in ('passed', 'failed', 'skipped', 'pending', 'not_applicable', None)}
    html_text = _env().get_template('index.html').render(
        rows=rows, totals=totals, generated=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC'),
        libcuflynx_version=_libcuflynx_version())
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(html_text)
    return out_path


def build_review_queue(contexts, out_path, module_href, standalone=True):
    '''One page listing every module awaiting review: comment, questions, proposed fixes.'''
    queue, reviewed = [], []
    for c in contexts:
        module = load_module(c['name'])
        if module.reviewed:
            reviewed.append(c['name'])
            continue
        if not module.spec.get('review'):
            continue
        queue.append({
            'name': c['name'], 'href': module_href(c['name']), 'review': module.spec['review'],
            'n_components': len(c['components']), 'counts': c['counts'], 'n_known': len(c['known_issues']),
            'n_sourced': c['n_sourced'], 'n_parameters': c['n_parameters'], 'max_risk': c.get('max_risk'),
            'n_proposals': sum(len(comp.get('reference_proposals') or {}) for comp in module.spec.get('components', [])),
        })
    html_text = _env().get_template('review_queue.html').render(
        modules=queue, reviewed=reviewed, standalone=standalone,
        generated=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC'))
    with open(out_path, 'w') as f:
        f.write(html_text)
    return out_path


def assemble_site(names, contexts):
    '''site/index.html + site/modules/<name>/{<name>.html, plots/}, as GitHub Pages serves it.'''
    if os.path.isdir(SITE_DIR):
        shutil.rmtree(SITE_DIR)
    for name in names:
        module = load_module(name)
        dest = os.path.join(SITE_DIR, 'modules', name)
        os.makedirs(dest, exist_ok=True)
        if not os.path.isfile(os.path.join(module.dir, f'{name}.html')):
            build_module(name)
        shutil.copy2(os.path.join(module.dir, f'{name}.html'), dest)
        if os.path.isdir(module.plots_dir):
            shutil.copytree(module.plots_dir, os.path.join(dest, 'plots'))
    build_review_queue(contexts, os.path.join(SITE_DIR, 'review_queue.html'), lambda n: f'modules/{n}/{n}.html')
    return build_index(contexts, os.path.join(SITE_DIR, 'index.html'), lambda n: f'modules/{n}/{n}.html')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--module', action='append', default=[], help='module name (repeatable); default all')
    parser.add_argument('--site', action='store_true', help='also assemble site/ for GitHub Pages')
    args = parser.parse_args(argv)
    names = args.module or module_names()
    contexts = []
    for name in names:
        out, ctx = build_module(name)
        contexts.append(ctx)
        print(f'wrote {os.path.relpath(out, REPO_ROOT)}')
    if args.site:
        # the index always covers every module, even when only some pages were rebuilt
        all_contexts = {c['name']: c for c in contexts}
        for name in module_names():
            if name not in all_contexts:
                all_contexts[name] = module_context(name)
        out = assemble_site(module_names(), [all_contexts[n] for n in module_names()])
        print(f'wrote {os.path.relpath(out, REPO_ROOT)}')


if __name__ == '__main__':
    main()
