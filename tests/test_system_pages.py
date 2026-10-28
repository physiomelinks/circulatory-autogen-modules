"""
The system model pages (cam_testing.report: build_system, the site's system_models/ and its index, the
shareable zip) on a small synthetic system model in a repo of its own, built from this library's versions:
an inlet, an artery, the heart supermodule and a record naming no version. No model is generated or simulated.
"""
import json
import os
import re
import struct
import zipfile
import zlib
from types import SimpleNamespace

import pytest

from cam_testing import paths, report, system as systems
from cam_testing.library import load_version

pytest.importorskip('libcuflynx')

LIBRARY = os.path.join(paths.PACKAGE_CHECKOUT, 'modules')
RECORDS = [
    {'name': 'flow_in', 'module_type': 'inlet_flow', 'module_subtype': 'constant', 'parameterisation': 'default',
     'inp_instances': [], 'out_instances': ['aorta']},
    {'name': 'aorta', 'module_type': 'arterial_simple', 'module_subtype': 'vp', 'parameterisation': 'default',
     'inp_instances': ['flow_in'], 'out_instances': []},
    {'name': 'heart', 'module_type': 'heart', 'module_subtype': 'Argus2026_v01', 'parameterisation': 'default',
     'per_submodule_inputs': {'ra': ['aorta']}, 'per_submodule_outputs': {'aov': ['aorta']}},
    {'name': 'ghost', 'module_type': 'no_such_type', 'module_subtype': 'v1', 'parameterisation': 'default'},
]
PARAMETERS = {'schema_version': 1,
              'global': {'rho': {'value': 1050, 'units': 'Js2_per_m5', 'data_reference': 'known'}},
              'aorta': {'R': {'value': 1.5e6, 'units': 'Js_per_m6', 'data_reference': 'Blanco_2013'}},
              'heart': {'lv': {'E_es': {'value': 2.5e8, 'units': 'J_per_m6', 'data_reference': 'fit'}}}}


def _png(path, width, height):
    '''A valid width x height white PNG.'''
    def chunk(kind, data):
        return struct.pack('>I', len(data)) + kind + data + struct.pack('>I', zlib.crc32(kind + data) & 0xffffffff)
    raw = b''.join(b'\x00' + b'\xff' * (3 * width) for _ in range(height))
    with open(path, 'wb') as f:
        f.write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', width, height, 8, 2, 0, 0, 0))
                + chunk(b'IDAT', zlib.compress(raw)) + chunk(b'IEND', b''))


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / 'repo'
    (root / 'modules').mkdir(parents=True)
    sdir = root / 'system_models' / 'demo' / 'chain'
    (sdir / 'results').mkdir(parents=True)
    (sdir / 'chain_module_array.json').write_text(json.dumps(RECORDS))
    (sdir / 'chain_parameters.json').write_text(json.dumps(PARAMETERS))
    (sdir / 'chain_system.yaml').write_text('model: chain\ncategory: demo\nsim_time: 1.0\n'
                                            'description: An inlet feeding an artery, with a heart.\n'
                                            'notes: [a note]\n')
    (sdir / 'results' / 'system_run_test.json').write_text(json.dumps(
        {'test': 'system_run_test', 'status': 'passed', 'message': 'generated and ran 1.0 s', 'metrics': {},
         'details': [], 'timestamp': '2026-10-28 00:00 UTC'}))
    # the report's per-run caches (the supermodule index) belong to the repo they were built for
    report._index_cache.clear()
    report._contents_cache.clear()
    paths.configure(repo_root=str(root), library_dirs=[LIBRARY], export=False)
    try:
        yield SimpleNamespace(root=root, system=systems.load_system('demo/chain'), dir=sdir)
    finally:
        paths.reset()
        report._index_cache.clear()
        report._contents_cache.clear()


def _hrefs(html):
    return re.findall(r'(?:href|src)="([^"#:]+)"', html)


def test_page_without_a_screenshot_draws_a_diagram(repo):
    out, ctx = report.build_system(repo.system)
    assert out == str(repo.dir / 'chain.html')
    html = open(out).read()
    assert ctx['figure'] is None and 'class="mermaid"' in html
    assert 'An inlet feeding an artery, with a heart.' in html and 'a note' in html
    by_name = {i['name']: i for i in ctx['instances']}
    # each record links to its version's page, relative to the system page
    aorta = load_version('arterial_simple', 'vp')
    assert by_name['aorta']['href'] == os.path.relpath(aorta.html_path, repo.dir)
    assert f'click n_aorta href "{by_name["aorta"]["href"]}"' in ctx['mermaid']
    # a supermodule instance lists its submodules, each linked
    heart = by_name['heart']
    assert heart['is_supermodule'] and {'ra', 'lv', 'aov'} <= {s['name'] for s in heart['submodules']}
    assert all(s['href'] for s in heart['submodules'])
    assert by_name['ghost']['href'] is None and 'not a version of the module libraries' in html
    # connections, parameters (collapsible, nested paths) and the tests
    assert ('flow_in', 'aorta') in ctx['edges'] and ('heart/aov', 'aorta') in ctx['edges'] \
        and ('aorta', 'heart/ra') in ctx['edges']
    assert '<details><summary><span class="mono">heart</span> · 1 value</summary>' in html
    assert 'lv/E_es' in html
    statuses = {t['key']: t['status'] for t in ctx['tests']}
    assert statuses == {'system_run_test': 'passed', 'system_invariants_test': None, 'system_equivalence_test': None}


def test_screenshot_nodes_link_to_the_version_pages(repo):
    _png(repo.dir / 'chain_phlynx.png', 200, 100)
    (repo.dir / 'chain_phlynx_nodes.json').write_text(json.dumps({
        'width': 200, 'height': 100, 'source': 'chain.omex',
        'nodes': [{'instance': 'aorta', 'box': [20, 10, 50, 40]}, {'instance': 'heart', 'box': [100, 50, 100, 50]}]}))
    out, ctx = report.build_system(repo.system)
    html = open(out).read()
    assert ctx['mermaid'] is None and 'src="chain_phlynx.png"' in html and not ctx['figure']['errors']
    links = re.findall(r'<a class="node" href="([^"]+)" data-instance="([^"]+)"[^>]*style="([^"]+)"', html)
    heart = load_version('heart', 'Argus2026_v01')
    assert [(h, i) for h, i, _ in links] == [
        (os.path.relpath(load_version('arterial_simple', 'vp').html_path, repo.dir), 'aorta'),
        (os.path.relpath(heart.html_path, repo.dir), 'heart')]
    # positioned in % of the image, so the boxes scale with it
    assert links[0][2] == 'left:10.000%;top:10.000%;width:25.000%;height:40.000%'


def test_screenshot_node_errors(repo):
    _png(repo.dir / 'chain_phlynx.png', 200, 100)
    (repo.dir / 'chain_phlynx_nodes.json').write_text(json.dumps({
        'width': 200, 'height': 100,
        'nodes': [{'instance': 'nobody', 'box': [0, 0, 10, 10]}, {'instance': 'aorta', 'box': [190, 90, 20, 20]}]}))
    ctx = report.system_context(repo.system)
    assert ctx['figure']['nodes'] == []
    assert any("'nobody': no instance" in e for e in ctx['figure']['errors'])
    assert any("'aorta': box" in e for e in ctx['figure']['errors'])
    (repo.dir / 'chain_phlynx_nodes.json').write_text(json.dumps({'width': 300, 'height': 100, 'nodes': []}))
    assert any('300x100' in e for e in report.system_context(repo.system)['figure']['errors'])


def test_site_and_zip_link_every_page(repo, monkeypatch):
    '''The site holds the system page, its screenshot, and the pages of the extra library's versions it links
    to (the heart's submodules too); every relative link of the system page resolves, and the zip has it all.'''
    _png(repo.dir / 'chain_phlynx.png', 200, 100)
    (repo.dir / 'chain_phlynx_nodes.json').write_text(json.dumps({
        'width': 200, 'height': 100, 'nodes': [{'instance': 'heart', 'box': [0, 0, 100, 50]}]}))
    monkeypatch.setattr(report, 'build_review_queue', lambda *a, **k: None)
    index = report.assemble_site([], [], [repo.system])
    site = repo.root / 'site'
    page = site / 'system_models' / 'demo' / 'chain' / 'chain.html'
    assert (page.parent / 'chain_phlynx.png').is_file()
    assert 'href="system_models/demo/chain/chain.html"' in open(index).read()
    html = page.read_text()
    targets = _hrefs(html)
    assert len(targets) >= 5
    for href in targets:
        assert (page.parent / href).resolve().is_file(), href
    chamber = load_version('chamber', next(s['module_subtype'] for s in load_version('heart', 'Argus2026_v01').submodules
                                           if s['name'] == 'lv'))
    assert (site / 'modules' / chamber.mtype.relpath / 'versions' / chamber.name / f'{chamber.stem}.html').is_file()
    assert (site / 'modules' / 'heart' / 'heart.html').is_file()
    z = report.zip_site(str(repo.root / 'share.zip'), title='Chain', start='Click System models, then chain.')
    with zipfile.ZipFile(z) as f:
        names = f.namelist()
        readme = f.read('README.txt').decode()
    assert 'site/index.html' in names and 'site/system_models/demo/chain/chain.html' in names
    assert readme.startswith('Chain') and 'Click System models, then chain.' in readme


def test_celldl_note_is_a_report_option(tmp_path, monkeypatch):
    '''--celldl-note (CAM_CELLDL_NOTE) puts a note under every bond graph and beside its CellDL link.'''
    vdir = tmp_path / 'v'
    vdir.mkdir()
    v = SimpleNamespace(name='v', stem='mt_v', dir=str(vdir), is_supermodule=False, config={'module_format': 'cellml'},
                        spec={}, results_dir=str(vdir / 'results'))
    v.path = lambda suffix: str(vdir / f'mt_v_{suffix}')
    (vdir / 'mt_v_bondgraph.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
    monkeypatch.setattr(report.phlynx, '_cpp', lambda version: False)
    ctx = report.version_context(load_version('Lotka_Volterra', 'Lotka1925_v01'))
    ctx['open_in'] = report.open_in_context(v)
    assert 'note' not in ctx['open_in']['celldl']
    assert 'id="celldl-note"' not in report._env().get_template('version.html').render(**ctx)
    note = 'The public editor may not draw signal blocks yet.'
    monkeypatch.setenv(report.CELLDL_NOTE_ENV, note)
    ctx['open_in'] = report.open_in_context(v)
    html = report._env().get_template('version.html').render(**ctx)
    assert re.search(r'<p class="celldl-note" id="celldl-note"[^>]*>' + re.escape(note) + '</p>', html)
    assert '<a href="#celldl-note"' in html
