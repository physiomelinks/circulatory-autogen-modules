"""The CellDL bridge (cam_testing/bridges/celldl/convert.mjs) with BGRDF_DIR: a BG-RDF checkout's ontology and
templates replace the bundled ones of the same name, and bundled templates it doesn't have are kept. Converts the RC
circuit in tests/fixtures/celldl; skipped without node or the bridge's packages (make celldl-setup)."""
import json
import os
import shutil
import subprocess

import pytest

from cam_testing import bondgraph

DOCUMENT = os.path.join(os.path.dirname(__file__), 'fixtures', 'celldl', 'rc_circuit_bondgraph.svg')
BRIDGE_DIR = os.environ.get('CELLDL_BRIDGE_DIR') or os.path.dirname(bondgraph.CELLDL_BRIDGE)
# the bundled framework, written out as a BG-RDF checkout's files
DUMP = '''import fs from 'node:fs'
const [, , entry, out] = process.argv
const { bgRdfFramework } = await import(entry)
const f = bgRdfFramework()
fs.mkdirSync(`${out}/schema`, { recursive: true })
fs.mkdirSync(`${out}/templates`, { recursive: true })
fs.writeFileSync(`${out}/schema/ontology.ttl`, f.ontology)
for (const [name, text] of Object.entries(f.templates)) fs.writeFileSync(`${out}/templates/${name}`, text)
'''


def convert(tmp_path, name):
    out = tmp_path / f'{name}.cellml'
    try:
        report = bondgraph.convert(DOCUMENT, str(out))
    except bondgraph.ConverterUnavailable as e:
        pytest.skip(f'the CellDL bridge cannot run: {e}')
    assert report['ok'], report
    return out.read_text()


@pytest.fixture
def bundled(tmp_path, monkeypatch):
    '''The bundled framework as files in a BG-RDF checkout's layout, and the CellML made with the bundle.'''
    monkeypatch.delenv('BGRDF_DIR', raising=False)
    cellml = convert(tmp_path, 'bundled')
    package = os.path.join(BRIDGE_DIR, 'node_modules', '@celldl', 'editor')
    with open(os.path.join(package, 'package.json')) as f:
        entry = os.path.join(package, json.load(f)['exports']['./bgrdf']['import'])
    script = tmp_path / 'dump.mjs'
    script.write_text(DUMP)
    out = tmp_path / 'bundled-bgrdf'
    subprocess.run([shutil.which('node'), str(script), entry, str(out)], check=True, timeout=120)
    return out, cellml


def test_a_bgrdf_dir_template_replaces_the_bundled_one(bundled, tmp_path, monkeypatch):
    files, bundled_cellml = bundled
    electrical = files / 'templates' / 'electrical.ttl'
    text = electrical.read_text()
    start = text.index('bgf:ElectricalResistor\n')
    law = text.index('<divide/>', start)
    electrical.write_text(text[:law] + '<times/>' + text[law + len('<divide/>'):])     # v = u R, not v = u / R
    monkeypatch.setenv('BGRDF_DIR', str(files))
    cellml = convert(tmp_path, 'bgrdf_dir')
    assert cellml != bundled_cellml
    assert cellml.count('<divide/>') == bundled_cellml.count('<divide/>') - 1
    assert cellml.count('<times/>') == bundled_cellml.count('<times/>') + 1


def test_bundled_templates_a_bgrdf_dir_lacks_are_kept(bundled, tmp_path, monkeypatch):
    files, bundled_cellml = bundled
    for template in (files / 'templates').iterdir():
        if template.name != 'chemical.ttl':
            template.unlink()
    monkeypatch.setenv('BGRDF_DIR', str(files))
    assert convert(tmp_path, 'bgrdf_dir') == bundled_cellml
