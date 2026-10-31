"""
Parameter derivations (cam_testing.derivations): the data_sources_test's derivation checks on a version, the
system_derivations_test on a system model, and the "Derived parameters" table on both pages.

The toy (tests/fixtures/derivations/derivations/): alpha = births / N by a CellML model, births an obs_data
input_item and N a data_item, then gamma = d0 + alpha by a Python function, chained on alpha. It is added as a
parameterisation of a copy of Lotka_Volterra/Lotka1925_v01 (reviewed) in a repo of its own, and as a system
model. Nothing is simulated.
"""
import json
import os
import shutil
from types import SimpleNamespace

import pytest

from cam_testing import checks, derivations, paths, report, sources, system as systems
from cam_testing.library import load_version

pytest.importorskip('libcuflynx.derivations', reason='libcuflynx has no derivations '
                    '(physiomelinks/circulatory_autogen#621)')

FIXTURE = os.path.join(os.path.dirname(__file__), 'fixtures', 'derivations', 'derivations')
LIBRARY = os.path.join(paths.PACKAGE_CHECKOUT, 'modules')
LV = os.path.join('benchmarks', 'Lotka_Volterra')
P = 'derived_demo'
REF = 'carpenter2018predator'
# a data_item's data_reference needs physiomelinks/circulatory_autogen#618: N_hares is sourced here
SOURCES_MD = '## N_hares\n\n- Kind: measured\n- Source: hewitt1921conservation\n'


def obs_data(name):
    return {'obs_data_name': name,
            'protocol_info': {'pre_times': [0.0], 'sim_times': [[1.0]], 'params_to_change': {}},
            'data_items': [{'data_item_name': 'N_hares', 'data_type': 'constant', 'operation': 'mean',
                            'operands': ['Lotka_Volterra/x'], 'unit': 'dimensionless', 'value': 30.0, 'std': 5.0,
                            'weight': 1.0}],
            'input_items': [{'data_item_name': 'births', 'value': 16.35, 'std': 2.0, 'unit': 'dimensionless',
                             'source': 'pelt counts', 'data_reference': 'hewitt1921conservation'},
                            {'data_item_name': 'd0', 'value': 0.258, 'unit': 'dimensionless',
                             'data_reference': REF}]}


def derivations_doc(prefix=''):
    return {'derivations': [
        {'name': 'gamma', 'inputs': {'d0': {'input_item': 'd0'}, 'alpha': {'parameter': f'{prefix}alpha'}},
         'outputs': {'gamma': f'{prefix}gamma'}, 'python': 'derivations/derivation_funcs.py:death_rate'},
        {'name': 'alpha', 'inputs': {'births': {'input_item': 'births'}, 'N': {'data_item': 'N_hares'}},
         'outputs': {'alpha': f'{prefix}alpha'}, 'cellml': 'derivations/birth_rate.cellml',
         'comment': 'births per hare'}]}


def leaves():
    leaf = {'units': 'dimensionless', 'data_reference': REF}
    return {'x_init': dict(leaf, value=33.956), 'y_init': dict(leaf, value=5.933),
            'alpha': {'value': 0.545, 'units': 'dimensionless', 'data_reference': 'derived; derivation alpha',
                      'derivation': 'alpha'},
            'beta': dict(leaf, value=0.028),
            'delta': dict(leaf, value=0.024),
            'gamma': {'value': 0.803, 'units': 'dimensionless', 'data_reference': 'derived; derivation gamma',
                      'derivation': 'gamma'}}


def _write(path, doc):
    with open(path, 'w') as f:
        json.dump(doc, f, indent=1)


def _read(path):
    with open(path) as f:
        return json.load(f)


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / 'repo'
    shutil.copytree(os.path.join(LIBRARY, LV), root / 'modules' / LV,
                    ignore=shutil.ignore_patterns('results', 'plots', 'risk', '*.html', '*.omex'))
    shutil.copy(os.path.join(LIBRARY, 'README.md'), root / 'modules' / 'README.md')
    vdir = root / 'modules' / LV / 'versions' / 'Lotka1925_v01'
    pdir = vdir / 'parameterisations' / P
    pdir.mkdir()
    _write(pdir / f'{P}_parameters.json', dict({'schema_version': 1}, **leaves()))
    _write(pdir / f'{P}_obs_data.json', obs_data(P))
    _write(pdir / f'{P}_derivations.json', derivations_doc())
    shutil.copytree(FIXTURE, pdir / 'derivations', ignore=shutil.ignore_patterns('__pycache__'))
    (pdir / 'SOURCES.md').write_text(SOURCES_MD)
    # the copied version's own data, sourced, so only the toy decides the tests
    obs_path = vdir / 'parameterisations' / 'hudson_bay_lynx_hare' / 'hudson_bay_lynx_hare_obs_data.json'
    obs = _read(obs_path)
    for it in obs['data_items'] + obs['prediction_items']:
        it['data_reference'] = 'hewitt1921conservation'
    _write(obs_path, obs)
    (vdir / 'parameterisations' / 'carpenter2018' / 'SOURCES.md').write_text(
        '## Hare\n\n- Kind: measured\n- Source: x\n\n## Lynx\n\n- Kind: measured\n- Source: x\n')
    # a system model with the same toy, its parameters under an instance
    sdir = root / 'system_models' / 'demo' / 'chain'
    sdir.mkdir(parents=True)
    _write(sdir / 'chain_module_array.json', [])
    (sdir / 'chain_system.yaml').write_text('model: chain\ncategory: demo\nreviewed: true\n')
    _write(sdir / 'chain_parameters.json', {'schema_version': 1, 'lv': leaves()})
    _write(sdir / 'chain_obs_data.json', obs_data('chain'))
    _write(sdir / 'chain_derivations.json', derivations_doc('lv/'))
    shutil.copytree(FIXTURE, sdir / 'derivations', ignore=shutil.ignore_patterns('__pycache__'))
    (sdir / 'SOURCES.md').write_text(SOURCES_MD)
    report._index_cache.clear()
    report._contents_cache.clear()
    paths.configure(repo_root=str(root), library_dirs=[], export=False)
    try:
        yield SimpleNamespace(root=root, pdir=pdir, sdir=sdir,
                              version=lambda: load_version('Lotka_Volterra', 'Lotka1925_v01'),
                              system=lambda: systems.load_system('demo/chain'))
    finally:
        paths.reset()
        report._index_cache.clear()
        report._contents_cache.clear()


def _edit(path, change):
    doc = _read(path)
    change(doc)
    _write(path, doc)


# ---- the checks ---------------------------------------------------------------------------------

def test_a_version_and_a_system_whose_derivations_check_pass(repo):
    r = sources.check(repo.version())
    assert r.status == checks.PASSED, r.details
    assert r.message.endswith('; 1 derivations files check') and r.metrics['warnings'] == []
    # the input_items are rows of the Data sources table, role input
    rows = {x['name']: x for x in sources.parameterisation_table(repo.version().parameterisation(P))['rows']}
    assert rows['births']['role'] == 'input' and rows['births']['comment'] == 'pelt counts'
    assert derivations.system_check(repo.system()) == (
        'passed', 'every derivation input resolves and is sourced; every derived value matches its cache', [])


@pytest.mark.parametrize('where', ['version', 'system'])
def test_an_input_item_without_a_source_fails(repo, where):
    path = repo.pdir / f'{P}_obs_data.json' if where == 'version' else repo.sdir / 'chain_obs_data.json'
    _edit(path, lambda d: d['input_items'][1].pop('data_reference'))
    if where == 'version':
        r = sources.check(repo.version())
        assert r.status == checks.FAILED and r.message.startswith('1 of ') and r.message.endswith(' data items have no source')
        assert any(d.startswith(f'{P}: d0: source missing') for d in r.details)
    else:
        status, msg, details = derivations.system_check(repo.system())
        assert status == 'failed' and details[0].startswith('d0: input_item source missing')
    # its section in SOURCES.md sources it
    with open(path.parent / 'SOURCES.md', 'a') as f:
        f.write('\n## d0\n\n- Kind: assumed\n')
    assert (sources.check(repo.version()).status if where == 'version'
            else derivations.system_check(repo.system())[0]) == 'passed'


@pytest.mark.parametrize('where', ['version', 'system'])
def test_a_stale_cache_fails(repo, where):
    path, key = ((repo.pdir / f'{P}_parameters.json', None) if where == 'version'
                 else (repo.sdir / 'chain_parameters.json', 'lv'))
    _edit(path, lambda d: (d[key] if key else d)['gamma'].update(value=0.81))
    if where == 'version':
        r = sources.check(repo.version())
        assert r.status == checks.FAILED and r.message == '1 derivation problems'
        detail = r.details[0]
    else:
        status, msg, details = derivations.system_check(repo.system())
        assert status == 'failed'
        detail = details[0]
    assert 'gamma: cached 0.81' in detail and "derivation 'gamma' gives 0.803" in detail


def test_inputs_that_do_not_resolve_and_an_invalid_cellml_model_fail(repo):
    # two bad inputs of one derivation: each is reported
    _edit(repo.pdir / f'{P}_derivations.json',
          lambda d: d['derivations'][1]['inputs'].update(N={'data_item': 'no_such_item'},
                                                         births={'input_item': 'no_such_input'}))
    with open(repo.pdir / 'derivations' / 'birth_rate.cellml') as f:
        text = f.read()
    with open(repo.pdir / 'derivations' / 'birth_rate.cellml', 'w') as f:
        f.write(text.replace('units="dimensionless" initial_value="1"', 'units="no_such_units" initial_value="1"'))
    r = sources.check(repo.version())
    assert r.status == checks.FAILED
    text = '\n'.join(r.details)
    assert "input 'N' (data_item 'no_such_item')" in text and 'no_such_units' in text
    assert "input 'births' (input_item 'no_such_input')" in text
    assert r.message == '3 derivation problems'


def test_a_derived_parameter_without_a_derivation_warns_in_a_reviewed_version_only(repo):
    _edit(repo.pdir / f'{P}_parameters.json', lambda d: d['beta'].update(data_reference='derived from pelt ratios'))
    r = sources.check(repo.version())
    assert r.status == checks.PASSED
    assert r.metrics['warnings'] == [f'{P}: beta: data_reference says "derived" but it names no derivation (add it '
                                     f'to {P}_derivations.json)']
    assert f'warning: {r.metrics["warnings"][0]}' in r.details
    # the suite test passes and warns
    from cam_testing.suite.test_modules import data_sources_test
    with pytest.warns(UserWarning, match='beta: data_reference says "derived"'):
        data_sources_test(('Lotka_Volterra', 'Lotka1925_v01'))
    spec = repo.version().spec_path
    with open(spec) as f:
        text = f.read()
    with open(spec, 'w') as f:
        f.write(text.replace('reviewed: true', 'reviewed: false'))
    assert sources.check(repo.version()).metrics['warnings'] == []


def test_without_libcuflynx_derivations_the_check_is_skipped(repo, monkeypatch):
    monkeypatch.setattr(derivations, 'unavailable', lambda: 'no libcuflynx.derivations')
    r = sources.check(repo.version())
    assert r.status == checks.SKIPPED and r.message.endswith('derivations not checked: no libcuflynx.derivations')
    assert derivations.system_check(repo.system())[:2] == ('skipped', 'derivations not checked: no libcuflynx.derivations')
    page = report.parameterisation_context(repo.version(), repo.version().parameterisation(P))['derived']
    assert page['skipped'] == 'no libcuflynx.derivations' and page['rows'] == []


def test_a_system_without_derivations_is_not_applicable(repo):
    os.remove(repo.sdir / 'chain_derivations.json')
    assert derivations.system_check(repo.system())[0] == 'not_applicable'


# ---- the pages ----------------------------------------------------------------------------------

def _check_table(tbl):
    rows = {r['name']: r for r in tbl['rows']}
    assert tbl['problems'] == [] and tbl['skipped'] is None
    alpha, gamma = rows['alpha'], rows['gamma']
    assert alpha['equations'] == [r'\alpha = \frac{\mathrm{births}}{N}']
    assert (gamma['function'], gamma['doc']) == ('death_rate', 'The predator death rate: a baseline rate plus the '
                                                               'prey birth rate.')
    births = next(i for i in alpha['inputs'] if i['symbol'] == 'births')
    assert (births['value'], births['std'], births['units'], births['source'], births['missing']) == \
        (16.35, 2.0, 'dimensionless', 'hewitt1921conservation', False)
    a = next(i for i in gamma['inputs'] if i['symbol'] == 'alpha')
    assert (a['kind'], a['source']) == ('parameter', 'derived; derivation alpha')
    assert gamma['outputs'][0]['value'] == pytest.approx(0.803) and not gamma['outputs'][0]['stale']


def test_version_page_shows_the_derived_parameters_table(repo):
    ctx = report.version_context(repo.version())
    inst = next(i for i in ctx['components'][0]['parameterisations'] if i['name'] == P)
    _check_table(inst['derived'])
    assert inst['derived']['file'] == f'parameterisations/{P}/{P}_derivations.json'
    html = report._env().get_template('version.html').render(**ctx)
    assert 'id="derived-parameters"' in html and 'id="derivation-alpha"' in html
    assert 'not declarative' in html and 'death_rate()' in html and r'\[ \alpha' in html
    _edit(repo.pdir / f'{P}_parameters.json', lambda d: d['gamma'].update(value=0.81))
    html = report._env().get_template('version.html').render(**report.version_context(repo.version()))
    assert 'cache stale: 0.81' in html


def test_system_page_shows_the_derived_parameters_table(repo):
    ctx = report.system_context(repo.system())
    _check_table(ctx['derived'])
    assert ctx['derived']['file'] == 'chain_derivations.json'
    assert [r['role'] for r in ctx['sources']['rows']] == ['target', 'input', 'input']
    html = report._env().get_template('system_model.html').render(**ctx)
    assert 'id="derived-h"' in html and 'tex-svg.js' in html and 'id="derivation-gamma"' in html
    assert 'Derived parameters' in [t['title'] for t in ctx['tests']]
