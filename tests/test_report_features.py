"""
Unit tests of the report's informational features: the licence / creator keys, the contents summary
bar, the libcellml unit-consistency check (mapped to equations) and the supermodule structure view.
No model is generated or simulated.
"""
import json
import os
from types import SimpleNamespace

from cam_testing import checks, report
from cam_testing.library import (DEFAULT_LICENCE, load_version, normalise_config_entry, to_phlynx_entry,
                                 with_record_keys)


# ---- licence and creator ------------------------------------------------------------------------

def test_record_keys_inserted_after_default_instance_and_carried_through():
    entry = {'module_type': 'm', 'module_subtype': 'v', 'component_file': 'f.cellml', 'component_type': 'c',
             'default_instance': 'default', 'entrance_ports': []}
    out = with_record_keys(entry)
    assert list(out) == ['module_type', 'module_subtype', 'component_file', 'component_type', 'default_instance',
                         'licence', 'creator', 'entrance_ports']
    assert out['licence'] == DEFAULT_LICENCE and out['creator'] == []
    # existing values are kept
    assert with_record_keys(dict(entry, licence='MIT', creator=['A. Author']))['creator'] == ['A. Author']
    # both config formats carry them
    assert normalise_config_entry(out)['licence'] == DEFAULT_LICENCE
    assert to_phlynx_entry(normalise_config_entry(out))['creator'] == []


def test_version_page_shows_licence_and_unset_creator():
    v = load_version('Lotka_Volterra', 'nn')
    assert v.licence == 'CC0-1.0' and v.creators == []
    ctx = report.version_context(v)
    assert ctx['licence_url'] == 'https://spdx.org/licenses/CC0-1.0.html'
    html = report._env().get_template('version.html').render(**ctx)
    assert 'creator: not set (to be given in review)'.lower() in html.lower()
    assert 'https://spdx.org/licenses/CC0-1.0.html' in html


# ---- contents summary ---------------------------------------------------------------------------

def test_contents_of_a_small_version():
    v = load_version('Lotka_Volterra', 'nn')
    c = report.contents(v)
    assert {k: c[k] for k in ('states', 'algebraic', 'parameters', 'boundary_conditions', 'ports', 'equations')} == \
        {'states': 2, 'algebraic': 0, 'parameters': 6, 'boundary_conditions': 0, 'ports': 0, 'equations': 2}
    assert c['instances'] == len(v.instance_names())
    bar = report.contents_bar(v, c, v.id)
    assert [b['key'] for b in bar] == ['states', 'algebraic', 'parameters', 'boundary_conditions', 'ports', 'equations',
                                       'instances']
    assert all(b['about'] and b['href'].startswith('#') for b in bar)


def test_contents_of_a_supermodule_sum_its_parts():
    v = load_version('neuron', 'sympathetic')
    c = report.contents(v)
    assert c['parts'] == 3 and c['levels'] == 2
    parts = [report.contents(load_version(s['module_type'], s['module_subtype'])) for s in v.submodules]
    assert c['states'] == sum(p['states'] for p in parts) > 0
    assert [b['key'] for b in report.contents_bar(v, c, v.id)][:2] == ['parts', 'levels']


# ---- unit consistency ---------------------------------------------------------------------------

TINY_CELLML = '''<?xml version='1.0' encoding='UTF-8'?>
<model name="tiny" xmlns="http://www.cellml.org/cellml/1.1#" xmlns:cellml="http://www.cellml.org/cellml/1.1#">
    <component name="tiny">
        <variable name="t" public_interface="in" units="second"/>
        <variable name="k" public_interface="in" units="per_second"/>
        <variable name="V" public_interface="in" units="volt"/>
        <variable name="x" initial_value="1" public_interface="out" units="metre"/>
        <variable name="y" public_interface="out" units="metre"/>
        <variable name="z" public_interface="out" units="metre"/>
        <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>x</ci></apply>
                <apply><times/><ci>k</ci><ci>x</ci></apply></apply>
            <apply><eq/><ci>y</ci><apply><plus/><ci>x</ci><ci>V</ci></apply></apply>
            <apply><eq/><ci>z</ci><ci>x</ci></apply>
        </math>
    </component>
</model>
'''
TINY_UNITS = '''<?xml version='1.0' encoding='UTF-8'?>
<model name="Units" xmlns="http://www.cellml.org/cellml/1.1#" xmlns:cellml="http://www.cellml.org/cellml/1.1#">
    <units name="per_second"><unit exponent="-1" units="second"/></units>
</model>
'''


def _tiny_version(tmp_path, cellml=TINY_CELLML):
    (tmp_path / 'tiny.cellml').write_text(cellml)
    (tmp_path / 'tiny_units.cellml').write_text(TINY_UNITS)
    (tmp_path / 'tiny_config.json').write_text('[]')
    return SimpleNamespace(cellml_path=str(tmp_path / 'tiny.cellml'), units_path=str(tmp_path / 'tiny_units.cellml'),
                           config_path=str(tmp_path / 'tiny_config.json'), results_dir=str(tmp_path / 'results'),
                           is_supermodule=False, format='cellml', module_type='tiny')


def test_unit_consistency_maps_the_failing_equation(tmp_path):
    v = _tiny_version(tmp_path)
    r = checks.unit_consistency(v)
    assert r['status'] == checks.UNITS_INCONSISTENT
    assert r['n_failures'] == 1 and r['n_equations'] == 1
    (f,) = r['failures']
    assert f['variable'] == 'y' and f['equation_index'] == 1      # y = x + V: metre + volt
    assert f['kind'] == 'not equivalent'
    assert {s['units'] for s in f['sides']} == {'metre', 'volt'}
    with open(os.path.join(v.results_dir, 'unit_consistency.json')) as fh:
        assert json.load(fh)['n_failures'] == 1
    # the page block lists only the failing equation, rendered as the page renders equations
    ctx = report.unit_consistency_context(v, ['eq0', 'eq1', 'eq2'])
    assert [g['latex'] for g in ctx['equations']] == ['eq1']


def test_unit_consistency_nothing_shown_when_consistent(tmp_path):
    consistent = TINY_CELLML.replace('<apply><plus/><ci>x</ci><ci>V</ci></apply>', '<ci>x</ci>')
    v = _tiny_version(tmp_path, consistent)
    assert checks.unit_consistency(v)['status'] == checks.UNITS_CONSISTENT
    assert report.unit_consistency_context(v, []) is None


def test_unit_consistency_not_applicable_to_a_supermodule():
    assert checks.unit_consistency(load_version('neuron', 'sympathetic'), save=False)['status'] == checks.NOT_APPLICABLE


# ---- supermodule structure ----------------------------------------------------------------------

def test_supermodule_structure_of_the_sympathetic_neuron():
    v = load_version('neuron', 'sympathetic')
    sv = report.supermodule_structure(v)
    parts = {p['name']: p for p in sv['parts']}
    assert list(parts) == ['soma', 'axon', 'varicosity']
    assert parts['soma']['is_supermodule'] and parts['varicosity']['is_supermodule']
    assert parts['axon']['href'].endswith('axon_sympathetic_monolithic_v01.html')
    # soma's membrane drives the axon, the axon drives the varicosity's membrane (resolved per_submodule_*)
    edges = {(e['src'], e['src_inner'], e['dst'], e['dst_inner']): e for e in sv['edges']}
    assert set(edges) == {('soma', 'membrane', 'axon', None), ('axon', None, 'varicosity', 'membrane')}
    assert edges[('soma', 'membrane', 'axon', None)]['ports'] == ['VI_port']
    # the supermodule's instance overrides the axon's parameters (rows <var>_axon); R, F, T are globals
    assert {o['name'] for o in parts['axon']['overrides']} >= {'C_axon', 'R_axon'}
    assert {'R', 'F', 'T'} <= set(sv['globals'])
    assert 'click p_soma' in sv['mermaid'] and 'p_soma -->' in sv['mermaid']
    assert [n['name'] for n in sv['tree']] == ['soma', 'axon', 'varicosity']
    assert sv['tree'][0]['children']


def test_nested_supermodule_shows_its_outside_coupling():
    sv = report.supermodule_structure(load_version('soma', 'sympathetic'))
    assert any(x['part'] == 'membrane' and x['external'] == 'axon' and x['direction'] == 'out'
               and x['where'][0]['name'] == 'neuron/sympathetic' for x in sv['external'])
    assert report.supermodule_structure(load_version('Lotka_Volterra', 'nn')) is None
