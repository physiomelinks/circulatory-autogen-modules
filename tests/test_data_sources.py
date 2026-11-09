"""
The data sources of calibration and validation data (cam_testing.sources): the SOURCES.md parser, each
kind and role of the "Data sources" table, the data_sources_test, and the table on a version page and a
system-model page. A synthetic parameterisation (tests/fixtures/sources/: the enterocyte whole-cell data's
SOURCES.md) is added to a copy of Lotka_Volterra/Lotka1925_v01 in a repo of its own. Nothing is simulated.
"""
import json
import os
import shutil
from types import SimpleNamespace

import pytest

from cam_testing import checks, paths, report, sources, system as systems
from cam_testing.library import load_version

FIXTURE = os.path.join(os.path.dirname(__file__), 'fixtures', 'sources', 'SOURCES.md')
LIBRARY = os.path.join(paths.PACKAGE_CHECKOUT, 'modules')
LV = os.path.join('benchmarks', 'Lotka_Volterra')
P = 'enterocyte_demo'


def item(name, value=1.0, std=0.1, weight=1.0, **kw):
    return dict({'data_item_name': name, 'data_type': 'constant', 'operation': 'mean', 'operands': [f'mod/{name}'],
                 'unit': 'dimensionless', 'cost_type': 'gaussian_MLE', 'weight': weight, 'value': value, 'std': std}, **kw)


def obs_data(**changes):
    '''The synthetic obs_data: the fixture's items (two terms, one without its "term of" label), an unsourced
    target, two string data_references, a held-out prediction and a prediction without data.'''
    data = [
        item('u_e_rest', -0.038, 0.003),
        item('c_Na_rest', 15, 4),
        item('u_e_with_sugar', 0, 1e15, 0.0),
        item('u_e_depolarisation_with_sugar', 0.006, 0.0015, operation='calculate_two_observable_difference',
             operation_kwargs={'x1': 'u_e_with_sugar', 'x2': 'u_e_rest'}),
        item('na_time_to_90pc_after_sugar', 120, 40),
        item('v_NKA_rest', 0, 1e15, 0.0),
        item('v_NKA_sugar', 0, 1e15, 0.0),
        item('R_te', 30, 5, data_reference='rose1971studies Table 2'),
        item('G_K_total', 2, 1),
    ]
    preds = [
        item('c_Na_peak_sugar', 20, 2),
        {'data_item_name': 'J_glucose', 'operation': 'mean', 'operands': ['mod/J_glucose'], 'unit': 'mol_per_s',
         'data_reference': 'rose1971studies; Fig. 8, digitised'},
    ]
    obs = {'obs_data_name': P, 'protocol_info': {'pre_times': [0.0], 'sim_times': [[1.0]], 'params_to_change': {}},
           'data_items': data, 'prediction_items': preds}
    obs.update(changes)
    return obs


def rows_by_name(table):
    return {r['name']: r for r in table['rows']}


# ---- the SOURCES.md parser ----------------------------------------------------------------------

def test_parser_reads_headings_labels_and_prose():
    with open(FIXTURE, encoding='utf-8') as f:
        sections, order = sources.parse_sources_md(f.read())
    # the heading is the item name: a trailing " — value" and "(proposed)" are dropped
    assert order[0] == 'enterocyte, whole cell'
    assert {'u_e_rest', 'c_Na_rest', 'c_Na_rise_with_sugar', 'u_e_with_sugar', 'v_NKA_sugar'} <= set(sections)
    s = sections['u_e_rest']
    assert s['kind'] == 'measured' and s['role'] == 'target'
    assert s['species'].startswith('rabbit ileum')
    assert s['primary_read'].startswith('apical ψ_mc −36 mV')       # unicode kept
    assert s['sigma_reason'].startswith('3 mV')
    assert s['prose'] == "The serosal composition has not been checked against the paper's methods."
    assert sources.classify_kind(sections['c_Na_rest']['kind']) == ('assumed', 'not measured in an enterocyte')
    assert sources.classify_kind(sections['u_e_with_sugar']['kind']) == ('term', 'u_e_depolarisation_with_sugar')
    assert sections['v_NKA_sugar'] == {'prose': '', '_level': 2, 'role': 'term'}


def test_heading_names_and_reference_forms():
    assert sources.heading_name('`c_Na_rest` — 15 ± 4 mM') == 'c_Na_rest'
    assert sources.heading_name('u_e_depolarisation_with_sugar – +6.0 +- 1.5 mV') == 'u_e_depolarisation_with_sugar'
    assert sources.heading_name('c_Na_peak_sugar   (held out, not fitted)') == 'c_Na_peak_sugar'
    assert sources.read_reference('smith2020 Fig. 2') == {'source': 'smith2020 Fig. 2'}
    assert sources.read_reference(None) == {} and sources.read_reference('') == {}
    # an object gives no source (the row reports it as a problem)
    assert sources.read_reference({'source': 'x'}) == {}
    # bold labels are read too
    sections, _ = sources.parse_sources_md('## a\n\n- **Kind:** assumed\n- **Source**: b2020\n')
    assert sections['a']['kind'] == 'assumed' and sections['a']['source'] == 'b2020'


def test_kinds_keep_their_detail_and_synthetic_is_a_kind():
    # parentheses around the whole detail are dropped, one inside it is kept
    assert sources.classify_kind('derived (synthetic: FitzHugh-Nagumo output at (a, b, c) = (0.2, 0.2, 3.0))') == \
        ('derived', 'synthetic: FitzHugh-Nagumo output at (a, b, c) = (0.2, 0.2, 3.0)')
    assert sources.classify_kind('assumed (from (x)) and (y)') == ('assumed', '(from (x)) and (y)')
    # model output used as data
    assert sources.classify_kind('synthetic (Liang 2009 model output, Fig. 3)') == \
        ('synthetic', 'Liang 2009 model output, Fig. 3')
    sections, _ = sources.parse_sources_md('## a\n\n- Kind: synthetic (the benchmark series)\n')
    row = sources._row('a', section=sections['a'], role='target')
    assert row['kind_label'] == 'synthetic (the benchmark series)' and not row['problem']


def test_data_reference_label_is_the_source_not_prose():
    sections, _ = sources.parse_sources_md('## a\n\n- Data reference: b2020; Fig. 2,\n  digitised\n'
                                           '\n## b\n\n- data reference: c2021; Table 1\n')
    assert sections['a']['source'] == 'b2020; Fig. 2, digitised' and sections['a']['prose'] == ''
    assert sections['a'][sources.DATA_REFERENCE] == 'b2020; Fig. 2, digitised'
    assert sections['b']['source'] == 'c2021; Table 1'


# ---- kinds and roles ----------------------------------------------------------------------------

@pytest.fixture
def table(tmp_path):
    shutil.copy(FIXTURE, tmp_path / 'SOURCES.md')
    return sources.data_table(str(tmp_path), obs=obs_data())


def test_kinds_and_roles(table):
    r = rows_by_name(table)
    # measured target from SOURCES.md
    assert (r['u_e_rest']['kind'], r['u_e_rest']['role'], r['u_e_rest']['origin']) == ('measured', 'target', 'SOURCES.md')
    # an assumed value can be a target
    assert r['c_Na_rest']['kind'] == 'assumed' and r['c_Na_rest']['role'] == 'target'
    assert r['c_Na_rest']['kind_label'] == 'assumed (not measured in an enterocyte)'
    # weight-0 terms: "term of <item>, not a measurement", from the label or from the data
    assert r['u_e_with_sugar']['kind_label'] == 'term of u_e_depolarisation_with_sugar, not a measurement'
    assert r['u_e_with_sugar']['role'] == 'term' and not r['u_e_with_sugar']['problem']
    # the fixture's Role: term matches the data, so it is not flagged
    assert r['u_e_with_sugar']['stated_role'] == '' and r['v_NKA_rest']['stated_role'] == ''
    assert sum(x['role'] == 'target' for x in table['rows']) == 6
    assert r['v_NKA_rest']['kind_label'] == 'term of v_NKA_ratio_with_sugar, not a measurement'
    assert r['v_NKA_sugar']['kind_label'] == 'term of another item, not a measurement'
    assert 'derived: term of' in r['v_NKA_sugar']['problem']
    # a measured value can be held out (a prediction_item with data)
    assert r['c_Na_peak_sugar']['kind'] == 'measured' and r['c_Na_peak_sugar']['role'] == 'held out'
    # a prediction without data, sourced by its data_reference string alone
    j = r['J_glucose']
    assert (j['role'], j['origin'], j['source'], j['kind'], j['problem']) == \
        ('prediction', 'data_reference', 'rose1971studies; Fig. 8, digitised', None, None)
    # a string data_reference, and an item with nothing
    assert r['R_te']['source'] == 'rose1971studies Table 2' and r['R_te']['kind_label'] == 'not stated'
    assert r['G_K_total']['missing'] and r['G_K_total']['problem'] == sources.MISSING
    assert (table['n_missing'], table['n_problems']) == (1, 2)
    # headings naming no item are reported (not failures)
    assert table['extra_headings'] == ['c_Na_rise_with_sugar', 'v_NKA_ratio_with_sugar', 'c_Na_sugar']


def test_only_weight_0_makes_a_term(tmp_path):
    obs = {'data_items': [item('a', 0, 5.0, 0.0), item('b', 2, 1, operation_kwargs={'x1': 'a', 'x2': 'c'}),
                          item('c', 1, 1, 0.0)], 'prediction_items': []}
    r = rows_by_name(sources.data_table(str(tmp_path), obs=obs))
    assert r['a']['kind_label'] == 'term of b, not a measurement'
    assert 'derived: term of b' in r['a']['problem']
    assert r['c']['role'] == 'term'
    # weight 0 with an ordinary std, named by no item: still a term, and its source must say so
    obs['data_items'] = [item('d', 3, 1, 0.0), item('e', 3, 1e20, 1.0)]
    r = rows_by_name(sources.data_table(str(tmp_path), obs=obs))
    assert r['d']['role'] == 'term' and r['d']['kind'] == 'term' and 'derived: term of' in r['d']['problem']
    # a std of very large order never makes a term: a weight-1 item is a target
    assert r['e']['role'] == 'target' and r['e']['kind'] is None and r['e']['problem'] == sources.MISSING


def test_an_empty_section_is_not_a_source(tmp_path):
    (tmp_path / 'SOURCES.md').write_text('## a\n\n## b\n\nRead off Fig. 2.\n')
    obs = {'data_items': [item('a'), item('b')], 'prediction_items': []}
    r = rows_by_name(sources.data_table(str(tmp_path), obs=obs))
    assert r['a']['missing'] and r['a']['origin'] is None and r['a']['problem'] == sources.MISSING
    assert r['a']['problem_kind'] == 'missing'
    # prose under the heading is a source
    assert r['b']['origin'] == 'SOURCES.md' and not r['b']['problem']


def test_a_term_names_the_item_that_uses_it(tmp_path):
    def term_problem(kind):
        (tmp_path / 'SOURCES.md').write_text(f'## a\n\n- Kind: {kind}\n')
        obs = {'data_items': [item('a', 0, 1, 0.0), item('b', 2, 1, operation_kwargs={'x1': 'a', 'x2': 'c'}),
                              item('c', 1, 1, data_reference='k2020')], 'prediction_items': []}
        return rows_by_name(sources.data_table(str(tmp_path), obs=obs))['a']
    assert term_problem('derived: term of b')['problem'] is None
    # "term of" without "derived:"
    row = term_problem('term of b')
    assert row['problem_kind'] == 'term' and 'must say "derived: term of b"' in row['problem']
    # the wrong item, naming the one that uses it
    row = term_problem('derived: term of zzz')
    assert row['problem_kind'] == 'term'
    assert 'names zzz, but b uses it' in row['problem'] and 'must say "derived: term of b"' in row['problem']


def test_data_reference_and_section_combine(tmp_path):
    (tmp_path / 'SOURCES.md').write_text('## a\n\n- Kind: measured\n- Data reference: k2020; Fig. 1\n'
                                         '- Figure: Fig. 1B\n\n## b\n\n- Data reference: k2020; Fig. 2\n'
                                         '\n## c\n\n- Kind: assumed\n- Source: an older note\n- Data reference: k2021\n- Figure: Table 2\n'
                                         '\n## e\n\n- Kind: measured\n')
    obs = {'data_items': [item('a', data_reference='k2020; Fig. 1'), item('b', data_reference='k2020; Fig. 3'),
                          item('c', data_reference='k2021'),
                          item('d', data_reference='k2022'), item('e', data_reference='k2023; Fig. 4')],
           'prediction_items': []}
    obs_path = str(tmp_path / 'x_obs_data.json')
    t = sources.data_table(str(tmp_path), obs=obs, obs_path=obs_path)
    r = rows_by_name(t)
    # section fields plus data_reference, the data_reference winning over a Source line; the origin says both
    assert (r['a']['origin'], r['a']['kind'], r['a']['figure'], r['a']['problem']) == \
        ('data_reference + SOURCES.md', 'measured', 'Fig. 1B', None)
    assert (r['c']['source'], r['c']['figure'], r['c']['kind'], r['c']['problem']) == ('k2021', 'Table 2', 'assumed', None)
    assert r['c']['origin'] == 'data_reference + SOURCES.md' and not r['c']['mismatch']
    assert r['d']['origin'] == 'data_reference'
    # the two strings differ: a problem naming both strings and both files
    b = r['b']
    assert b['mismatch'] and not b['missing']
    assert all(s in b['problem'] for s in ('"k2020; Fig. 3"', '"k2020; Fig. 2"', 'x_obs_data.json', 'SOURCES.md'))
    # a data_reference and a section without its Data reference line: a problem naming both files
    e = r['e']
    assert e['problem_kind'] == 'reference' and not e['mismatch'] and not e['missing']
    assert all(s in e['problem'] for s in ('"k2023; Fig. 4"', 'x_obs_data.json', 'SOURCES.md',
                                           '"- Data reference: k2023; Fig. 4"'))
    assert t['n_problems'] == 2 and t['problems'] == {'missing': 0, 'reference': 2, 'not_string': 0, 'term': 0}


def test_an_object_data_reference_is_a_problem(tmp_path):
    '''data_reference is a string, as in libcuflynx: an object (even with a section that would source the item)
    is a problem naming the file and saying where its parts go, on obs_data items and baseline targets and blocks.'''
    (tmp_path / 'SOURCES.md').write_text('## a\n\n- Kind: measured\n- Data reference: k2020; Fig. 1\n')
    obs = {'data_items': [item('a', data_reference={'source': 'k2020; Fig. 1', 'figure': 'Fig. 1B'}),
                          item('b', data_reference='k2021')],
           'prediction_items': [{'data_item_name': 'p', 'operands': ['mod/p'], 'data_reference': {'source': 'k'}}]}
    t = sources.data_table(str(tmp_path), obs=obs, obs_path=str(tmp_path / 'x_obs_data.json'))
    r = rows_by_name(t)
    for name in ('a', 'p'):
        assert r[name]['problem_kind'] == 'not_string' and r[name]['source'] != 'k'
        assert r[name]['problem'] == ('data_reference in x_obs_data.json is an object: it must be a string; put the '
                                      'citation as "<BibTeX key>; <short note>" in data_reference and the detail '
                                      "in the item's section of SOURCES.md")
    assert not r['b']['problem']
    assert t['problems'] == {'missing': 0, 'reference': 0, 'not_string': 2, 'term': 0}
    assert sources.summary(t['problems']) == '2 data_reference not a string'
    # the data_sources_test's details name the item
    assert sources.table_problems(t, 'demo: ', '')[0].startswith('demo: a: data_reference in x_obs_data.json is an object')
    # baseline targets and a baseline block's own data_reference: the same rule
    block = {'status': 'active', 'targets': [{'name': 'EDV', 'expr': 'np.max(q)', 'value': 142, 'std': 21,
                                              'data_reference': {'source': 'maceira2006normalized'}}],
             'data': 'p/x.csv', 'variables': {'mod/x': 'x_col'}, 'data_reference': {'source': 'ref2000'}}
    rows = sources.data_table(str(tmp_path), baseline=block)['rows']
    assert [x['problem_kind'] for x in rows] == ['not_string', 'not_string']


def test_baseline_targets_and_columns(tmp_path):
    (tmp_path / 'SOURCES.md').write_text('## EDV\n\n- Kind: measured\n- Source: maceira2006normalized Table 8\n')
    block = {'status': 'active', 'targets': [{'name': 'EDV', 'expr': 'np.max(q)', 'value': 142, 'std': 21},
                                             {'name': 'ESV', 'expr': 'np.min(q)', 'value': 47, 'std': 10}]}
    r = rows_by_name(sources.data_table(str(tmp_path), baseline=block))
    assert r['EDV']['role'] == 'validation' and not r['EDV']['missing'] and r['ESV']['missing']
    # a pending baseline is not read by the validation test: no items
    assert sources.data_table(str(tmp_path), baseline=dict(block, status='pending')) is None
    series = {'status': 'active', 'data': 'p/x.csv', 'variables': {'mod/x': 'x_col'}, 'data_reference': 'ref2000'}
    row = sources.data_table(str(tmp_path), baseline=series)['rows'][0]
    assert (row['name'], row['variable'], row['source']) == ('x_col', 'mod/x', 'ref2000')


# ---- the version, its test and its page ---------------------------------------------------------

@pytest.fixture
def repo(tmp_path):
    root = tmp_path / 'repo'
    shutil.copytree(os.path.join(LIBRARY, LV), root / 'modules' / LV,
                    ignore=shutil.ignore_patterns('results', 'plots', 'risk', '*.html', '*.omex'))
    shutil.copy(os.path.join(LIBRARY, 'README.md'), root / 'modules' / 'README.md')
    vdir = root / 'modules' / LV / 'versions' / 'Lotka1925_v01'
    pdir = vdir / 'parameterisations' / P
    pdir.mkdir()
    shutil.copy(vdir / 'parameterisations' / 'default' / 'default_parameters.json', pdir / f'{P}_parameters.json')
    (pdir / f'{P}_obs_data.json').write_text(json.dumps(obs_data()))
    shutil.copy(FIXTURE, pdir / 'SOURCES.md')
    report._index_cache.clear()
    report._contents_cache.clear()
    paths.configure(repo_root=str(root), library_dirs=[], export=False)
    try:
        yield SimpleNamespace(root=root, version=load_version('Lotka_Volterra', 'Lotka1925_v01'), pdir=pdir)
    finally:
        paths.reset()
        report._index_cache.clear()
        report._contents_cache.clear()


def _strip_lv_data(repo):
    '''Gives the copied version's own hudson_bay_lynx_hare data a data_reference each, so only the
    synthetic parameterisation decides the test.'''
    path = repo.version.parameterisation('hudson_bay_lynx_hare').obs_data_path
    with open(path) as f:
        obs = json.load(f)
    for it in obs['data_items'] + obs['prediction_items']:
        it['data_reference'] = 'hewitt1921conservation (via carpenter2018predator)'
    with open(path, 'w') as f:
        json.dump(obs, f)
    # without the SOURCES.md, whose Data reference lines would differ from that string
    md = os.path.join(os.path.dirname(path), 'SOURCES.md')
    if os.path.exists(md):
        os.remove(md)
    # carpenter2018's baseline data columns, each section repeating the block's data_reference
    carpenter = repo.version.parameterisation('carpenter2018')
    ref = carpenter.validation['baseline']['data_reference']
    with open(os.path.join(carpenter.dir, 'SOURCES.md'), 'w') as f:
        f.write(f'## Hare\n\n- Kind: measured\n- Data reference: {ref}\n\n## Lynx\n\n- Kind: measured\n'
                f'- Data reference: {ref}\n')


def test_data_sources_test_fails_on_unsourced_items_and_lists_where_to_add_them(repo):
    _strip_lv_data(repo)
    r = sources.check(repo.version)
    assert r.status == checks.FAILED
    assert r.message == '2 of 22 data items have a problem: 1 no source, 1 term not labelled'
    assert (r.metrics['n_missing'], r.metrics['n_reference'], r.metrics['n_term']) == (1, 0, 1)
    assert any(d.startswith(f'{P}: G_K_total: source missing') and f'parameterisations/{P}/SOURCES.md' in d
               for d in r.details)
    assert any(d.startswith(f'{P}: v_NKA_sugar: weight-0 term') for d in r.details)
    # the extra headings are listed, not counted
    assert any('"c_Na_sugar" names no data item' in d for d in r.details)


def test_data_sources_test_passes_once_every_item_is_sourced(repo):
    _strip_lv_data(repo)
    with open(repo.pdir / 'SOURCES.md', 'a', encoding='utf-8') as f:
        f.write('- Kind: derived: term of v_NKA_ratio_with_sugar\n\n## G_K_total\n\n- Kind: assumed\n')
    assert sources.check(repo.version).status == checks.PASSED


def test_data_sources_test_reports_a_data_reference_mismatch(repo):
    _strip_lv_data(repo)
    with open(repo.pdir / 'SOURCES.md', 'a', encoding='utf-8') as f:
        f.write('- Kind: derived: term of v_NKA_ratio_with_sugar\n\n## G_K_total\n\n- Kind: assumed\n'
                '\n## R_te\n\n- Data reference: rose1971studies Table 3\n')
    r = sources.check(repo.version)
    assert r.status == checks.FAILED
    assert r.message == '1 of 22 data items have a problem: 1 Data reference mismatch/missing'
    d = next(d for d in r.details if d.startswith(f'{P}: R_te: '))
    assert '"rose1971studies Table 2" in ' + f'parameterisations/{P}/{P}_obs_data.json' in d
    assert f'"rose1971studies Table 3" in parameterisations/{P}/SOURCES.md' in d


def test_data_sources_test_unsourced_and_not_applicable(repo, tmp_path):
    # the copied version's own data, without their SOURCES.md or data_reference keys, are unsourced
    hb = repo.version.parameterisation('hudson_bay_lynx_hare')
    os.remove(os.path.join(hb.dir, 'SOURCES.md'))
    with open(hb.obs_data_path) as f:
        obs = json.load(f)
    for it in obs['data_items'] + obs['prediction_items']:
        it.pop('data_reference', None)
    with open(hb.obs_data_path, 'w') as f:
        json.dump(obs, f)
    r = sources.check(repo.version)
    assert r.status == checks.FAILED and 'hudson_bay_lynx_hare: mod/x: source missing' in '\n'.join(r.details)
    shutil.rmtree(repo.pdir)
    for p in ('hudson_bay_lynx_hare', 'carpenter2018'):
        shutil.rmtree(repo.version.parameterisation(p).dir)
    assert sources.check(load_version('Lotka_Volterra', 'Lotka1925_v01')).status == checks.NOT_APPLICABLE


def test_a_reviewed_version_with_unsourced_data_fails_ci(repo, tmp_path):
    '''The suite test fails, and tools/ci_reviewed_failures.py makes that blocking for a reviewed version
    (Lotka1925_v01 is reviewed) and a warning for one not reviewed yet.'''
    import subprocess
    import sys
    from cam_testing.suite.test_modules import data_sources_test
    assert repo.version.reviewed
    with pytest.raises(pytest.fail.Exception, match='1 no source'):
        data_sources_test(('Lotka_Volterra', 'Lotka1925_v01'))
    junit = tmp_path / 'junit.xml'
    junit.write_text('<testsuites><testsuite><testcase name="data_sources_test[Lotka_Volterra/Lotka1925_v01]">'
                     '<failure message="x"/></testcase></testsuite></testsuites>')
    env = dict(os.environ, CAM_REPO_ROOT=str(repo.root), CAM_MODULE_LIBRARY_DIRS='')
    tool = os.path.join(paths.PACKAGE_CHECKOUT, 'tools', 'ci_reviewed_failures.py')
    out = subprocess.run([sys.executable, tool, str(junit)], env=env, capture_output=True, text=True)
    assert out.returncode == 1 and '::error::data_sources_test[Lotka_Volterra/Lotka1925_v01] failed' in out.stdout
    spec = repo.version.spec_path
    with open(spec) as f:
        text = f.read()
    with open(spec, 'w') as f:
        f.write(text.replace('reviewed: true', 'reviewed: false'))
    out = subprocess.run([sys.executable, tool, str(junit)], env=env, capture_output=True, text=True)
    assert out.returncode == 0 and '::warning::' in out.stdout


def test_version_page_shows_the_table_and_the_count(repo):
    _strip_lv_data(repo)
    ctx = report.version_context(repo.version)
    assert ctx['n_data_items'] == 22 and ctx['n_data_problems'] == 2
    html = report._env().get_template('version.html').render(**ctx)
    assert 'id="data-sources"' in html and 'id="data-sources-enterocyte_demo"' in html
    assert 'source missing' in html and 'term of u_e_depolarisation_with_sugar, not a measurement' in html
    assert '2 of 22 data items with a problem: 1 no source, 1 term not labelled' in html
    assert '2 with a problem: 1 no source, 1 term not labelled' in html      # the table's own count
    index = report._env().get_template('index.html').render(
        groups=[], systems=[], totals={}, n_module_types=0, n_versions=0, generated='', libcuflynx_version='',
        data_sources={'items': 22, 'problems': 2, 'counts': ctx['data_problems'], 'systems': 1,
                      'system_counts': {'missing': 0, 'reference': 1, 'not_string': 0, 'term': 0}})
    assert ('2 of 22 data items with a problem (1 no source, 1 term not labelled), and 1 in system models '
            '(1 Data reference mismatch/missing)') in index
    # bib keys of the version link to its References
    assert '<a href="#ref-hewitt1921conservation">hewitt1921conservation</a>' in html
    page = report.parameterisation_context(repo.version, repo.version.parameterisation(P))['sources']
    assert page['sources_md'] == f'parameterisations/{P}/SOURCES.md'


def test_cite_links_known_keys_only():
    html = str(report._cite('rose1971studies p. 648; <b>x</b> rose1971studiesX', ['rose1971studies']))
    assert html.startswith('<a href="#ref-rose1971studies">rose1971studies</a> p. 648; &lt;b&gt;')
    assert html.endswith('rose1971studiesX')


def test_system_page_shows_its_obs_data_sources(tmp_path):
    root = tmp_path / 'repo'
    (root / 'modules').mkdir(parents=True)
    sdir = root / 'system_models' / 'demo' / 'chain'
    sdir.mkdir(parents=True)
    (sdir / 'chain_module_array.json').write_text(json.dumps([]))
    (sdir / 'chain_system.yaml').write_text('model: chain\ncategory: demo\n')
    (sdir / 'chain_obs_data.json').write_text(json.dumps(obs_data()))
    shutil.copy(FIXTURE, sdir / 'SOURCES.md')
    paths.configure(repo_root=str(root), library_dirs=[], export=False)
    try:
        system = systems.load_system('demo/chain')
        ctx = report.system_context(system)
        assert ctx['sources']['n_problems'] == 2 and ctx['sources']['sources_md'] == 'SOURCES.md'
        html = report._env().get_template('system_model.html').render(**ctx)
        assert 'id="data-sources-h"' in html and 'source missing' in html
        assert report.system_rows([ctx], lambda c: 'x.html')[0]['n_data_problems'] == 2
    finally:
        paths.reset()


def test_a_reviewed_system_model_with_unsourced_data_fails_ci(tmp_path):
    '''system_data_sources_test fails for a system model with an unsourced item, and
    tools/ci_reviewed_failures.py makes that blocking when its spec has reviewed: true, a warning otherwise.'''
    import subprocess
    import sys
    from cam_testing.suite.test_systems import system_data_sources_test
    root = tmp_path / 'repo'
    (root / 'modules').mkdir(parents=True)
    sdir = root / 'system_models' / 'demo' / 'chain'
    sdir.mkdir(parents=True)
    (sdir / 'chain_module_array.json').write_text(json.dumps([]))
    spec = sdir / 'chain_system.yaml'
    spec.write_text('model: chain\ncategory: demo\nreviewed: true\n')
    (sdir / 'chain_obs_data.json').write_text(json.dumps({'data_items': [item('a', data_reference='k2020'), item('b')],
                                                          'prediction_items': []}))
    paths.configure(repo_root=str(root), library_dirs=[], export=False)
    try:
        with pytest.raises(pytest.fail.Exception, match=r'1 of 2 data items have a problem: 1 no source'):
            system_data_sources_test(systems.load_system('demo/chain'))
    finally:
        paths.reset()
    junit = tmp_path / 'junit.xml'
    junit.write_text('<testsuites><testsuite><testcase name="system_data_sources_test[demo/chain]">'
                     '<failure message="x"/></testcase></testsuite></testsuites>')
    env = dict(os.environ, CAM_REPO_ROOT=str(root), CAM_MODULE_LIBRARY_DIRS='')
    tool = os.path.join(paths.PACKAGE_CHECKOUT, 'tools', 'ci_reviewed_failures.py')
    out = subprocess.run([sys.executable, tool, str(junit)], env=env, capture_output=True, text=True)
    assert out.returncode == 1 and '::error::system_data_sources_test[demo/chain] failed' in out.stdout
    spec.write_text('model: chain\ncategory: demo\nreviewed: false\n')
    out = subprocess.run([sys.executable, tool, str(junit)], env=env, capture_output=True, text=True)
    assert out.returncode == 0 and '::warning::' in out.stdout and 'system model not reviewed yet' in out.stdout
