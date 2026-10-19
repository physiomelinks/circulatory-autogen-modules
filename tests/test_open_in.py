"""
The version page's "Open in" links (CellDL for a bond graph, PhLynx for the version's PhLynx export,
CUFLynx per instance archive), the files the site copies for them, the PhLynx export kept with the
results, the calibration result in an instance's CUFLynx archive, and the page's base64 encoder
(cam_testing/templates/_open_in.js, run under node). No model is generated or simulated.
"""
import base64
import io
import json
import os
import re
import shutil
import subprocess
import zipfile
from types import SimpleNamespace

import pytest

from cam_testing import bondgraph, omex, phlynx, report
from cam_testing.library import load_version

OPEN_IN_JS = os.path.join(report.TEMPLATES, '_open_in.js')


def _fake_version(tmp_path, name='v', supermodule=False, module_format='cellml'):
    vdir = tmp_path / name
    vdir.mkdir(parents=True, exist_ok=True)
    stem = f'mt_{name}'
    v = SimpleNamespace(name=name, vessel_type='mt', key=f'mt/{name}', id=f'mt__{name}', stem=stem, dir=str(vdir),
                        results_dir=str(vdir / 'results'), plots_dir=str(vdir / 'plots'),
                        html_path=str(vdir / f'{stem}.html'), is_supermodule=supermodule,
                        config={'module_format': module_format}, spec={})
    v.path = lambda suffix: str(vdir / f'{stem}_{suffix}')
    v.instances = lambda: []
    return v


def _render(ctx):
    return report._env().get_template('version.html').render(**ctx)


# ---- which links a version gets -----------------------------------------------------------------

def test_celldl_link_only_for_a_bond_graph(tmp_path):
    v = _fake_version(tmp_path)
    assert report.open_in_context(v)['celldl'] is None
    with open(bondgraph.document_path(v), 'w') as f:
        f.write('<svg xmlns="http://www.w3.org/2000/svg"/>')
    assert report.open_in_context(v)['celldl'] == {'href': 'mt_v_bondgraph.svg', 'open': report.CELLDL_OPEN}
    assert report.CELLDL_OPEN == 'https://finbarargus.github.io/CellDLEditorComponent/?open=celldl#'


def test_phlynx_link_to_the_kept_export(tmp_path):
    v = _fake_version(tmp_path)
    assert report.open_in_context(v)['phlynx'] == {'missing': report.PHLYNX_NOT_EXPORTED}
    os.makedirs(v.results_dir)
    with open(phlynx.archive_path(v), 'wb') as f:
        f.write(b'PK')
    assert report.open_in_context(v)['phlynx'] == {'href': 'results/mt_v_phlynx.omex', 'open': 'https://phlynx.com/?open=omex#'}


def test_phlynx_not_applicable_to_a_supermodule_or_cpp(tmp_path):
    sup = _fake_version(tmp_path, 'sup', supermodule=True)
    assert report.open_in_context(sup)['phlynx'] == {'reason': phlynx.NO_SUPERMODULE_SUPPORT}
    cpp = _fake_version(tmp_path, 'cpp', module_format='cpp')
    assert report.open_in_context(cpp)['phlynx'] == {'reason': phlynx.CPP_REASON}


# ---- the version page -------------------------------------------------------------------------

def test_page_of_a_version_without_a_bond_graph_has_no_celldl_link():
    ctx = report.version_context(load_version('Lotka_Volterra', 'Lotka1925_v01'))
    assert ctx['open_in']['celldl'] is None
    html = _render(ctx)
    assert 'id="open-celldl"' not in html and report.CELLDL_EDITOR_URL not in html
    assert 'class="open-in"' in html


def test_page_of_a_bond_graph_has_the_celldl_and_phlynx_links(tmp_path):
    ctx = report.version_context(load_version('Lotka_Volterra', 'Lotka1925_v01'))
    v = _fake_version(tmp_path)
    with open(bondgraph.document_path(v), 'w') as f:
        f.write('<svg xmlns="http://www.w3.org/2000/svg"/>')
    os.makedirs(v.results_dir)
    with open(phlynx.archive_path(v), 'wb') as f:
        f.write(b'PK')
    ctx['open_in'] = report.open_in_context(v)
    html = _render(ctx)
    celldl = re.search(r'<a [^>]*id="open-celldl"[^>]*>', html).group(0)
    assert 'data-open-file="mt_v_bondgraph.svg"' in celldl and 'href="mt_v_bondgraph.svg"' in celldl
    assert f'data-open-url="{report.CELLDL_OPEN}"' in celldl
    ph = re.search(r'<a [^>]*id="open-phlynx"[^>]*>', html).group(0)
    assert 'data-open-file="results/mt_v_phlynx.omex"' in ph and 'data-open-url="https://phlynx.com/?open=omex#"' in ph
    # the encoder is inline, so the page needs nothing else to open the file
    assert 'camOpenIn' in html and 'bytesToBase64' in html


def test_supermodule_page_shows_phlynx_not_applicable():
    ctx = report.version_context(load_version('neuron', 'sympathetic'))
    assert ctx['open_in']['phlynx'] == {'reason': phlynx.NO_SUPERMODULE_SUPPORT}
    html = _render(ctx)
    assert 'id="open-phlynx"' not in html and 'id="phlynx-na"' in html
    assert 'PhLynx: N/A' in html and phlynx.NO_SUPERMODULE_SUPPORT.replace('"', '&#34;') in html


def test_cuflynx_link_per_instance_archive():
    ctx = report.version_context(load_version('Lotka_Volterra', 'Lotka1925_v01'))
    instances = ctx['components'][0]['instances']
    assert instances
    for i in instances:
        i['omex'] = f'instances/{i["name"]}/Lotka_Volterra_Lotka1925_v01_{i["name"]}.omex'
    # an instance without a built archive gets no link
    instances[0]['omex'] = None
    html = _render(ctx)
    links = re.findall(r'<a [^>]*data-cuflynx-omex="([^"]*)"[^>]*>', html)
    assert links == [i['omex'] for i in instances[1:]]
    tags = re.findall(r'<a [^>]*data-cuflynx-omex="[^"]*"[^>]*>', html)
    assert all('data-cuflynx-open="cuflynx://open?url="' in t and 'cuflynx://' in t for t in tags)
    assert 'needs a CUFLynx that handles cuflynx:// links' in html


# ---- the site ---------------------------------------------------------------------------------

def test_assemble_site_copies_the_bond_graph_and_phlynx_export(tmp_path, monkeypatch):
    v = _fake_version(tmp_path / 'lib')
    for path, data in ((v.html_path, b'<html></html>'), (bondgraph.document_path(v), b'<svg/>'),
                       (phlynx.archive_path(v), b'PK\x03\x04')):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'wb') as f:
            f.write(data)
    plain = _fake_version(tmp_path / 'lib', 'plain')
    with open(plain.html_path, 'w') as f:
        f.write('<html></html>')
    mtype = SimpleNamespace(relpath='cat/mt', html_path=str(tmp_path / 'lib' / 'mt.html'), versions=lambda: [v, plain])
    with open(mtype.html_path, 'w') as f:
        f.write('<html></html>')
    site = tmp_path / 'site'
    monkeypatch.setattr(report, 'load_module_type', lambda name: mtype)
    monkeypatch.setattr(report.paths, 'roots', lambda: SimpleNamespace(site_dir=str(site)))
    monkeypatch.setattr(report, 'build_review_queue', lambda *a, **k: None)
    monkeypatch.setattr(report, 'build_index', lambda *a, **k: None)
    report.assemble_site(['mt'], [])
    vdest = site / 'modules' / 'cat' / 'mt' / 'versions' / 'v'
    # at the paths the page links them by
    assert (vdest / 'mt_v_bondgraph.svg').read_bytes() == b'<svg/>'
    assert (vdest / 'results' / 'mt_v_phlynx.omex').read_bytes() == b'PK\x03\x04'
    assert sorted(os.listdir(site / 'modules' / 'cat' / 'mt' / 'versions' / 'plain')) == ['mt_plain.html']


# ---- the PhLynx export kept with the results ----------------------------------------------------

def test_keep_archives_copies_exports_and_drops_stale_ones(tmp_path):
    ok, failed = _fake_version(tmp_path, 'ok'), _fake_version(tmp_path, 'failed')
    exported = tmp_path / 'work' / 'mt__ok.omex'
    exported.parent.mkdir()
    exported.write_bytes(b'PK new')
    os.makedirs(failed.results_dir)
    with open(phlynx.archive_path(failed), 'wb') as f:
        f.write(b'PK old')
    phlynx.keep_archives([ok, failed], {'mt__ok': {'ok': True, 'omex': str(exported)},
                                        'mt__failed': {'ok': False, 'omex': None}})
    assert phlynx.archive_path(ok) == os.path.join(ok.results_dir, 'mt_ok_phlynx.omex')
    with open(phlynx.archive_path(ok), 'rb') as f:
        assert f.read() == b'PK new'
    assert not os.path.exists(phlynx.archive_path(failed))


# ---- the calibration result in an instance's CUFLynx archive -------------------------------------

def _omex_inputs(tmp_path, monkeypatch, calibrated=True):
    vdir = tmp_path / 'v'
    idir = vdir / 'instances' / 'base'
    idir.mkdir(parents=True)
    stem = 'mt_v'
    version = SimpleNamespace(vessel_type='mt', name='v', id='mt__v', key='mt/v', is_supermodule=False,
                              cellml_path=str(vdir / f'{stem}_modules.cellml'), config_path=str(vdir / f'{stem}_modules_config.json'),
                              units_path=str(vdir / f'{stem}_units.cellml'),
                              verification_config_path=str(vdir / f'{stem}_verification_config.json'))
    inst = SimpleNamespace(name='base', dir=str(idir), parameters=lambda: [])
    for attr, suffix in (('obs_data_path', 'obs_data.json'), ('params_for_id_path', 'params_for_id.csv'),
                         ('parameters_path', 'parameters.csv'), ('calibration_path', 'calibration.json'),
                         ('calibrated_parameters_path', 'calibrated_parameters.csv')):
        setattr(inst, attr, str(idir / f'base_{suffix}'))
    files = {inst.obs_data_path: '{"obs_data_name": "base"}', inst.params_for_id_path: 'vessel_name,param_name\n',
             inst.parameters_path: 'variable_name,value\n', str(idir / 'SOURCES.md'): '# sources\n',
             str(idir / 'a_raw_data.csv'): 't,x\n', version.cellml_path: '<model/>', version.config_path: '[]'}
    if calibrated:
        files[inst.calibration_path] = json.dumps({'calibrated_parameters': {'k': 1.5}, 'cost': 0.1})
        files[inst.calibrated_parameters_path] = 'variable_name,value\nk,1.5\n'
    for path, text in files.items():
        with open(path, 'w') as f:
            f.write(text)

    def generate(v, work_dir, **kw):
        model = os.path.join(work_dir, 'model.cellml')
        with open(model[:-len('.cellml')] + '_flat.cellml', 'w') as f:
            f.write('<model name="flat"/>')
        return model
    from cam_testing import harness
    monkeypatch.setattr(harness, 'generate', generate)
    return version, inst


def test_omex_carries_the_calibration_result_after_the_study(tmp_path, monkeypatch):
    version, inst = _omex_inputs(tmp_path, monkeypatch)
    work = tmp_path / 'work'
    work.mkdir()
    names = [n for n, _, _ in omex.members(version, inst, str(work))]
    assert names[:6] == ['mt_v_base.cellml', 'base_obs_data.json', 'base_params_for_id.csv', 'base_parameters.csv',
                         'base_calibration.json', 'base_calibrated_parameters.csv']
    assert names.count('base_calibration.json') == 1 and names.count('base_calibrated_parameters.csv') == 1
    # CUFLynx takes the first JSON with "obs" and the first CSV with "param" in its name
    assert next(n for n in names if n.endswith('.json') and 'obs' in n) == 'base_obs_data.json'
    assert next(n for n in names if n.endswith('.csv') and 'param' in n) == 'base_params_for_id.csv'

    path = omex.build(version, inst, str(work), str(tmp_path / 'out.omex'))
    with zipfile.ZipFile(path) as z:
        manifest = z.read('manifest.xml').decode()
        assert json.loads(z.read('base_calibration.json'))['calibrated_parameters'] == {'k': 1.5}
    assert ('location="./base_calibration.json" format="http://purl.org/NET/mediatypes/application/json"' in manifest)
    assert ('location="./base_calibrated_parameters.csv" format="http://purl.org/NET/mediatypes/text/csv"' in manifest)


def test_omex_without_a_calibration_result(tmp_path, monkeypatch):
    version, inst = _omex_inputs(tmp_path, monkeypatch, calibrated=False)
    work = tmp_path / 'work'
    work.mkdir()
    names = [n for n, _, _ in omex.members(version, inst, str(work))]
    assert not any('calibrat' in n for n in names)
    assert names[:4] == ['mt_v_base.cellml', 'base_obs_data.json', 'base_params_for_id.csv', 'base_parameters.csv']


# ---- the page's encoder, under node --------------------------------------------------------------

NODE_SCRIPT = r'''
const fs = require('fs'), vm = require('vm');
const [js, data, out] = process.argv.slice(2);
const ctx = vm.createContext({ btoa, URL, encodeURIComponent });
vm.runInContext(fs.readFileSync(js, 'utf8'), ctx);
const bytes = new Uint8Array(fs.readFileSync(data));
const b64 = ctx.camOpenIn.bytesToBase64(bytes);
const back = Buffer.from(b64, 'base64');
fs.writeFileSync(out, JSON.stringify({
  b64, roundtrip: back.equals(Buffer.from(bytes)),
  page: ctx.camOpenIn.cuflynxUrl('cuflynx://open?url=', 'instances/a b/x.omex',
                                 'https://physiomelinks.github.io/cam/modules/v/x.html'),
  disk: ctx.camOpenIn.cuflynxUrl('cuflynx://open?url=', 'instances/i/x.omex', 'file:///home/u/modules/v/x.html'),
}));
'''


def test_page_base64_round_trips_a_binary_file(tmp_path):
    node = shutil.which('node')
    if not node:
        pytest.skip('node is not installed')
    # every byte value, and more than several of the encoder's chunks (an odd length, so padding is used)
    data = bytes(range(256)) * 8 + os.urandom(3 * 0x8000 + 7)
    zbuf = io.BytesIO()
    with zipfile.ZipFile(zbuf, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('a.bin', data)
    blob = tmp_path / 'blob.omex'
    blob.write_bytes(zbuf.getvalue() + data)
    script, out = tmp_path / 'enc.js', tmp_path / 'out.json'
    script.write_text(NODE_SCRIPT)
    subprocess.run([node, str(script), OPEN_IN_JS, str(blob), str(out)], check=True, timeout=60)
    r = json.loads(out.read_text())
    assert r['roundtrip']
    assert r['b64'] == base64.b64encode(blob.read_bytes()).decode()
    assert r['page'] == 'cuflynx://open?url=https%3A%2F%2Fphysiomelinks.github.io%2Fcam%2Fmodules%2Fv%2Finstances%2Fa%2520b%2Fx.omex'
    assert r['disk'] == 'cuflynx://open?url=file%3A%2F%2F%2Fhome%2Fu%2Fmodules%2Fv%2Finstances%2Fi%2Fx.omex'
