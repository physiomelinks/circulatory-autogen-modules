"""
Detailed metabolism modules (tools/build_metabolism_detailed.py) against the lumped modules of the generic adaptive cell.

    pytest tests/test_metabolism_detailed.py --include-unreviewed -s

cell/glycolysis detailed_ArgusUNPUBLISHED_v01 runs glucose to lactate in eleven steps; energy_metabolism
adaptive_ArgusUNPUBLISHED_v01 lumps them into two (hexokinase, then the lower glycolysis with Hill kinetics in G6P).
Checks: the bond graph (every step balanced in elements and charge, thermodynamically consistent; its BG-RDF and
diagram); the detailed module rests at the fermenter's resting flux and conserves carbon; and, with the cytosol's
metabolites clamped alike, how the two respond to glucose and to the ATP/ADP ratio (printed with -s). They agree at
rest (calibrated so) and as glucose rises to about five times the resting flux; beyond, the detailed pathway saturates
lower (hexose phosphates build up behind aldolase's equilibrium) and it is more sensitive to ADP/ATP (the PGK and
pyruvate kinase equilibria), which the lumped kinetics do not have.
"""
import contextlib
import io
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from scipy.optimize import brentq

from cam_testing import checks, harness
from cam_testing.library import load_version

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'tools'))
import build_metabolism_detailed as detailed  # noqa: E402
from bondgraph import BondGraph, R_GAS  # noqa: E402

pytestmark = pytest.mark.module_vv
VERSION = detailed.VERSION
VDIR = os.path.join(REPO, 'modules', 'cell', 'glycolysis', 'versions', VERSION)
INTERNAL = [s.name for s in detailed.GLY_SPECIES if not s.boundary]


@pytest.fixture(scope='module')
def graph():
    bg = BondGraph('glycolysis', 'Detailed glycolysis', detailed.GLY_SPECIES, detailed.GLY_REACTIONS,
                   temperature=detailed.T)
    return bg


@pytest.fixture(scope='module')
def det(tmp_path_factory):
    return checks.ComponentModel(load_version('glycolysis', VERSION), work_dir=str(tmp_path_factory.mktemp('det')))


@pytest.fixture(scope='module')
def lump(tmp_path_factory):
    return checks.ComponentModel(load_version('energy_metabolism', 'adaptive_ArgusUNPUBLISHED_v01'),
                                 work_dir=str(tmp_path_factory.mktemp('lump')))


def _run(cm, params, t_end, dt, outputs):
    with contextlib.redirect_stdout(io.StringIO()):
        helper = harness.simulation_helper(cm.model_path, dict(cm.spec, sim_time=t_end, dt=dt))
        helper.reset_states()
        if params:
            helper.set_param_vals(['parameters/' + k + '_mod' for k in params], list(params.values()), change_states=True)
        assert helper.run() is not False
    t = np.asarray(helper.get_time(), dtype=float)
    return t, {o: np.asarray(helper.get_results([[harness.output_name('mod/' + o)]], flatten=True)[0], dtype=float).ravel()
               for o in outputs}


def _detailed_flux(det, glc, atp, adp):
    _, r = _run(det, {'C_glc': glc, 'C_ATP': atp, 'C_ADP': adp}, 600.0, 6.0, ['v_HK', 'v_LDH'])
    return r['v_HK'][-1], 0.5 * r['v_LDH'][-1]


def _lumped_flux(lump, glc, atp, adp):
    '''The lumped cytosol's glycolytic flux with its glucose, lactate, ATP, ADP and Pi clamped: hexokinase and the lower
    glycolysis at the G6P where they balance, from the module's own rate expressions.'''
    m = lump.helper.model.clone()
    st = {s.qname().split('.')[-1]: s for s in m.states()}
    var = {v.qname().split('.')[-1]: v for v in m.variables(deep=True)}
    for k, x in (('C_glc', glc), ('C_ATP', atp), ('C_ADP', adp), ('C_Pi', detailed.GLY_REST['Pi']),
                 ('C_lac', detailed.GLY_REST['lac'])):
        st[k].set_initial_value(x)

    def gap(g6p):
        st['C_G6P'].set_initial_value(g6p)
        return var['v_gly'].eval() - var['v_lg'].eval()
    st['C_G6P'].set_initial_value(brentq(gap, 1e-9, 1e3))
    return var['v_gly'].eval()


def test_every_step_balances_and_the_thermodynamics_are_consistent(graph):
    '''Each reaction balances C, H, O, N, P and charge; the species' standard potentials reproduce every textbook
    standard Gibbs energy (no cycle to break them: residual round-off); the internal species are steady at the resting
    fluxes, and at rest every step runs down its Gibbs energy.'''
    assert graph.element_imbalance() == {} and graph.charge_imbalance() == {}
    assert graph.fit_residual < 1e-6
    for r in graph.reactions:
        assert graph.dG0_fitted(r) == pytest.approx(r.dG0, abs=1e-6)
    steady = graph.steady(detailed.GLY_FLUX)
    assert all(abs(steady[s]) < 1e-15 for s in INTERNAL)
    rest = detailed.near_equilibrium_rest(graph)
    for r in graph.reactions:
        assert graph.reaction_gibbs(r, rest) < 0.0, r.name
        if r.name in detailed.NEAR_EQ:
            assert graph.reaction_gibbs(r, rest) == pytest.approx(R_GAS * detailed.T * np.log(1 - 1 / detailed.EXCHANGE),
                                                                  rel=1e-9)


def test_its_bond_graph_is_written_as_bg_rdf_and_drawn(graph):
    '''The BG-RDF (reduced form, as the CellDL editor writes it) has the model, a storage node per species and a reaction
    node per reaction, and a bond per unit of stoichiometry; the diagram is well-formed SVG carrying that RDF where
    CellDL keeps it, with an element per species and reaction named as in the bond graph.'''
    ttl = open(os.path.join(VDIR, f'glycolysis_{VERSION}_bondgraph.ttl')).read()
    assert 'a bgf:BondgraphModel' in ttl and 'bgf:usesTemplate tpl:chemical.ttl' in ttl
    assert ttl.count('a bgf:ChemicalStorageNode') == len(graph.species)
    assert ttl.count('a bgf:ChemicalReactionNode') == len(graph.reactions)
    assert ttl.count('bgf:hasSource') == sum(abs(n) for r in graph.reactions for n in r.stoich.values())
    svg = ET.parse(os.path.join(VDIR, f'glycolysis_{VERSION}_bondgraph.svg')).getroot()
    ns = '{http://www.w3.org/2000/svg}'
    meta = svg.find(ns + 'metadata')
    assert meta is not None and meta.get('id') == 'celldl-rdf-metadata' and meta.get('data-content-type') == 'text/turtle'
    assert meta.text.strip() == ttl.strip()
    ids = {g.get('id') for g in svg.iter(ns + 'g')}
    assert {f'q_{s.name}' for s in graph.species} | {f'Re_{r.name}' for r in graph.reactions} <= ids


def test_it_rests_at_the_fermenters_glycolytic_flux(det):
    '''At its default (the fermenter's rest) every intermediate stays put for a minute, glucose is used at half the
    resting ATP turnover and ATP made at twice that, and the reactions dissipate free energy.'''
    outs = [f'C_{s}' for s in INTERNAL] + ['v_HK', 'j_ATP', 'j_glc', 'P_diss']
    _, r = _run(det, {}, 60.0, 1.0, outs)
    for s in INTERNAL:
        np.testing.assert_allclose(r[f'C_{s}'], r[f'C_{s}'][0], rtol=1e-7, err_msg=s)
    assert r['v_HK'][-1] == pytest.approx(detailed.J_REST, rel=1e-6)
    assert r['j_glc'][-1] == pytest.approx(-detailed.J_REST, rel=1e-6)
    assert r['j_ATP'][-1] == pytest.approx(2.0 * detailed.J_REST, rel=1e-6)
    assert np.all(r['P_diss'] > 0.0)


def test_carbon_is_conserved(det):
    '''Ten times the resting glucose for a second, on a fine grid: the carbon the intermediates gain is what the
    chemostats give (relative 1e-4 of the turnover).'''
    carbon = {s.name: s.formula.get('C', 0) for s in detailed.GLY_SPECIES}
    bound = [s.name for s in detailed.GLY_SPECIES if s.boundary and carbon[s.name]]
    outs = [f'C_{s}' for s in INTERNAL] + [f'j_{s}' for s in bound]
    t, r = _run(det, {'C_glc': 10.0 * detailed.GLY_REST['glc']}, 1.0, 1e-4, outs)
    inner = sum(carbon[s] * r[f'C_{s}'] for s in INTERNAL)
    given = sum(carbon[s] * r[f'j_{s}'] for s in bound)
    turnover = np.trapz(np.abs(given), t)
    assert abs(inner[-1] - inner[0] + np.trapz(given, t)) <= 1e-4 * turnover


def test_compared_with_the_lumped_glycolysis(det, lump, graph):
    '''With the cytosol's glucose, lactate, ATP, ADP and Pi clamped alike: the two carry the same flux at rest; as glucose
    rises to five times rest they stay within 20%; at saturating glucose the detailed pathway carries less (its lower
    half saturates, hexose phosphates build up), and as ATP falls and ADP rises it slows sooner. Their overall standard
    Gibbs energies of glucose -> 2 lactate + 2 ATP differ by the lumped dG0_ferm (textbook -196 kJ/mol; BioNumbers
    -183.6) against the steps' sum.'''
    rest = detailed.GLY_REST
    rows = []
    for glc in (0.25 * rest['glc'], rest['glc'], 5.0 * rest['glc'], 2500.0 * rest['glc']):
        hk, low = _detailed_flux(det, glc, rest['ATP'], rest['ADP'])
        rows.append(('glucose', glc, _lumped_flux(lump, glc, rest['ATP'], rest['ADP']), hk, low))
    for atp, adp in ((6.0, 2.0), (3.0, 4.0)):
        hk, low = _detailed_flux(det, rest['glc'], atp, adp)
        rows.append((f'ATP {atp} ADP {adp}', rest['glc'], _lumped_flux(lump, rest['glc'], atp, adp), hk, low))
    print('\ncondition            glucose [mM]   lumped [mM/s]   detailed: hexokinase, lower half [mM/s]')
    for name, glc, l, hk, low in rows:
        print(f'{name:20s} {glc:12.4g} {l:15.5f} {hk:12.5f} {low:12.5f}')
    (_, _, l0, h0, d0), (_, _, l1, h1, d1), (_, _, l5, h5, d5), (_, _, ls, hs, ds) = rows[:4]
    assert h1 == pytest.approx(l1, rel=1e-6) and d1 == pytest.approx(l1, rel=1e-6), 'calibrated alike at rest'
    assert d0 == pytest.approx(l0, rel=0.4) and h0 == pytest.approx(l0, rel=0.2)
    assert d5 == pytest.approx(l5, rel=0.2) and h5 == pytest.approx(l5, rel=0.1)
    assert ds < 0.7 * ls, 'the detailed pathway saturates lower'
    (_, _, la, _, da), (_, _, lb, _, db) = rows[4:]
    assert da < 0.8 * la and db < la and db < da, 'and slows sooner as ATP falls'

    lumped = detailed.CYTO['dG0_ferm'] - 2.0 * detailed.CYTO['dG0_hyd']
    weights = {r.name: (1 if r.name in ('HK', 'PGI', 'PFK', 'ALD', 'TPI') else 2) for r in graph.reactions}
    steps = sum(weights[r.name] * r.dG0 for r in graph.reactions)
    print(f'glucose -> 2 lactate + 2 ATP: lumped {lumped / 1e3:.1f} kJ/mol, the steps {steps / 1e3:.1f} kJ/mol')
    assert steps < -100e3 and lumped < -100e3 and abs(steps - lumped) < 15e3
