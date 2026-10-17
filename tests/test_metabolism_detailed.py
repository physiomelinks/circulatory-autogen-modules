"""
Detailed metabolism modules (tools/build_metabolism_detailed.py) against the lumped modules of the generic adaptive cell,
and the bond-graph descriptions of both (tools/bondgraph.py, tools/build_bondgraph_lumped.py).

    pytest tests/test_metabolism_detailed.py --include-unreviewed -s

Each detailed module is a bond graph: every step balanced in elements and charge and thermodynamically consistent, its
BG-RDF giving the module's own rates, its diagram carrying that RDF as CellDL does. Comparisons (printed with -s), with
the cytosol clamped alike:
  glycolysis (eleven steps) against energy_metabolism adaptive's two: alike at rest and up to about five times the
    resting flux; beyond, the detailed pathway saturates lower and is more sensitive to ADP/ATP;
  the mitochondrion against energy_metabolism mitochondrial's lumped respiration: alike at rest, respiratory control
    through the matrix potential; complex IV respires down to a few uM O2, the lumped K_O2 (a tissue value) does not,
    and at high ADP the detailed one saturates sooner;
  the glycogen store (ADP-glucose route) against the lumped store: alike at the comparison point, less responsive to
    G6P (pyrophosphate builds up), no ceiling;
  gluconeogenesis: its bypasses close futile cycles with glycolysis' kinases worth one ATP each, and glycerol enters only
    where NADH/NAD+ is low.
energy_metabolism adaptive's own bond graph has its stoichiometry; its phosphorylase's equilibrium leaves out glycogen's
activity (a strict xfail).
"""
import contextlib
import io
import json
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
import build_bondgraph_lumped as lumped  # noqa: E402
import build_celldl as celldl_build  # noqa: E402
import celldl  # noqa: E402
import compare_celldl as celldl_compare  # noqa: E402
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
def mito():
    return detailed.mitochondrion_graph()


@pytest.fixture(scope='module')
def mdet(tmp_path_factory):
    return checks.ComponentModel(load_version('mitochondrion', VERSION), work_dir=str(tmp_path_factory.mktemp('mdet')))


@pytest.fixture(scope='module')
def mlump(tmp_path_factory):
    return checks.ComponentModel(load_version('energy_metabolism', 'mitochondrial_ArgusUNPUBLISHED_v01'),
                                 work_dir=str(tmp_path_factory.mktemp('mlump')))


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


MODULES = {'glycolysis': (detailed.GLY_FLUX, detailed.CYTO['V_c']), 'mitochondrion': (detailed.MITO_FLUX, detailed.MITO['V_c']),
           'glycogen': (detailed.GLYC_FLUX, detailed.CYTO['V_c']),
           'gluconeogenesis': (detailed.GNG_FLUX, detailed.MITO['V_c'])}


@pytest.fixture(scope='module')
def built():
    '''Each module's bond graph as built (saturations set) and its rate constants.'''
    with contextlib.redirect_stdout(io.StringIO()):
        return {m: detailed.BUILDERS[m]() for m in MODULES}


def _document(module):
    """The version's CellDL document (_bondgraph.svg) and its Turtle (_bondgraph.ttl)."""
    vdir, stem = celldl_build.version_dir(module)
    with open(os.path.join(vdir, f'{stem}_bondgraph.ttl')) as f:
        ttl = f.read()
    return ET.parse(os.path.join(vdir, f'{stem}_bondgraph.svg')).getroot(), ttl


@pytest.mark.parametrize('module', list(MODULES))
def test_its_bond_graph_is_a_celldl_document(module, built):
    """The version's bond graph is a CellDL document: its Turtle (the _bondgraph.ttl) in CellDL's metadata element, a
    celldl:Document and bgf:BondgraphModel; in its diagram layer a component per species, reaction and chemostat copy
    (one beside each reaction that uses it) and a bond per species of each reaction; each element the template of its
    constitutive relation (saturating where the module saturates)."""
    graph = built[module][0]
    svg, ttl = _document(module)
    ns = '{http://www.w3.org/2000/svg}'
    meta = svg.find(ns + 'metadata')
    assert meta is not None and meta.get('id') == 'celldl-rdf-metadata' and meta.get('data-content-type') == 'text/turtle'
    assert meta.text == ttl
    assert '<> a celldl:Document, bgf:BondgraphModel ;' in ttl
    layer = [g for g in svg.iter(ns + 'g') if g.get('id') == 'celldl-diagram-layer'][0]
    components = {g.get('id') for g in layer.iter(ns + 'g') if g.get('class') == 'celldl-Component'}
    internal = [sp for sp in graph.species if not sp.boundary]
    boundary = {sp.name for sp in graph.species if sp.boundary}
    copies = {f'Se_{sp}_{r.name}' for r in graph.reactions for sp in r.stoich if sp in boundary}
    assert components == {f'q_{sp.name}' for sp in internal} | {f'Re_{r.name}' for r in graph.reactions} | copies
    bonds = [p for p in layer.iter(ns + 'path') if 'celldl-Connection' in (p.get('class') or '')]
    assert len(bonds) == sum(len(r.stoich) for r in graph.reactions)
    for r in graph.reactions:
        law = 'SaturatingChemicalReactionNode' if r.sigma else 'ChemicalReactionNode'
        assert f':Re_{r.name} a celldl:Component, bgf:{law} ;' in ttl
    for sp in internal:
        law = 'SaturatingChemicalStorageNode' if sp.k_sat is not None else 'ChemicalStorageNode'
        assert f':q_{sp.name} a celldl:Component, bgf:{law} ;' in ttl


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


MITO_INTERNAL = [s.name for s in detailed.MITO_SPECIES if not s.boundary]


def test_the_mitochondrion_balances_and_its_thermodynamics_are_consistent(mito):
    """Each reaction balances C, H, O, N, P, S and charge; the standard potentials reproduce the textbook Gibbs energies
    and midpoint potentials; at the resting fluxes every matrix species is steady and no net charge crosses the inner
    membrane; at rest every step runs down its electrochemical Gibbs energy, the near-equilibrium ones at the set
    exchange, and the proton-motive force is the matrix potential plus the pH difference."""
    assert mito.element_imbalance() == {} and mito.charge_imbalance() == {}
    assert mito.fit_residual < 1e-6
    steady = mito.steady(detailed.MITO_FLUX)
    assert all(abs(steady[s]) < 1e-15 for s in MITO_INTERNAL)
    assert mito.current('m', detailed.MITO_FLUX) == pytest.approx(0.0, abs=1e-15)
    rest = detailed.mitochondrion_rest(mito)
    for r in mito.reactions:
        assert mito.reaction_gibbs(r, rest) < 0.0, r.name
        if r.name in detailed.MITO_NEAR_EQ:
            assert mito.reaction_gibbs(r, rest) == pytest.approx(R_GAS * detailed.T_M * np.log(1 - 1 / detailed.EXCHANGE),
                                                                 rel=1e-9)
    leak = {r.name: r for r in mito.reactions}['LEAK']
    dp = -detailed.PSI_REST + R_GAS * detailed.T_M / detailed.FARADAY * np.log(10.0) * (detailed.PH_M - detailed.PH_C)
    assert -mito.reaction_gibbs(leak, rest) == pytest.approx(detailed.FARADAY * dp, rel=1e-9)


def test_the_mitochondrion_rests_at_the_lumped_respirations_o2_use(mdet):
    """At its default every matrix species and the membrane potential stay put for a minute, O2 is used at the lumped
    respiration's resting rate and 4.5 ATP are exported per O2 (P/O 2.25)."""
    outs = [f'C_{s}' for s in MITO_INTERNAL] + ['psi_m', 'j_O2', 'j_ATP_c', 'P_diss']
    _, r = _run(mdet, {}, 60.0, 1.0, outs)
    for s in MITO_INTERNAL + ['psi_m']:
        key = s if s == 'psi_m' else f'C_{s}'
        np.testing.assert_allclose(r[key], r[key][0], rtol=1e-6, err_msg=s)
    assert -r['j_O2'][-1] == pytest.approx(detailed.O2_REST, rel=1e-6)
    assert r['j_ATP_c'][-1] / -r['j_O2'][-1] == pytest.approx((detailed.N_ATPS + 1.0) / 2.5, rel=1e-6)
    assert np.all(r['P_diss'] > 0.0)


def _mito_o2(mdet, params):
    _, r = _run(mdet, params, 120.0, 12.0, ['j_O2', 'j_ATP_c', 'psi_m'])
    return -r['j_O2'][-1], r['j_ATP_c'][-1], r['psi_m'][-1]


def _lumped_o2(mlump, o2, adp):
    """The lumped respiration's O2 use with O2 and ADP (ATP = A_tot - ADP) clamped: 6 r_ox from its rate expression."""
    m = mlump.helper.model.clone()
    var = {v.qname().split('.')[-1]: v for v in m.variables(deep=True)}
    st = {s.qname().split('.')[-1]: s for s in m.states()}
    var['C_O2_mod'].set_rhs(o2)
    st['ATP'].set_initial_value(detailed.MITO['A_tot'] - adp)
    return 6.0 * var['r_ox'].eval()


def test_compared_with_the_lumped_respiration(mdet, mlump):
    """With the cytosol clamped alike (ATP + ADP as the lumped module's A_tot): the two use O2 at the same rate at rest;
    both slow as O2 falls and speed up with ADP (respiratory control), the detailed one through the membrane potential
    (the synthase and the pumps), the lumped one through Michaelis-Menten factors. They differ in O2: complex IV keeps
    the detailed one respiring at nearly its full rate down to a few uM (as isolated mitochondria do), while the lumped
    K_O2 is a tissue value, which includes diffusion to the mitochondria; and at high ADP the detailed one saturates
    sooner (its dehydrogenases and pumps can carry 4 times their resting flux)."""
    a_tot, o2_rest, adp_rest = detailed.MITO['A_tot'], detailed.MITO['C_O2'], detailed._ADP
    rows = []
    for o2 in (0.25 * detailed.MITO['K_O2'], detailed.MITO['K_O2'], o2_rest, 3.0 * o2_rest):
        d, atp, psi = _mito_o2(mdet, {'C_O2': o2})
        rows.append(('O2', o2, adp_rest, _lumped_o2(mlump, o2, adp_rest), d, atp, psi))
    for adp in (0.1, adp_rest, 1.0, 2.0):
        d, atp, psi = _mito_o2(mdet, {'C_ADP_c': adp, 'C_ATP_c': a_tot - adp})
        rows.append(('ADP', o2_rest, adp, _lumped_o2(mlump, o2_rest, adp), d, atp, psi))
    print('\nvaried   O2 [mM]   ADP [mM]   lumped O2 use   detailed O2 use [mM/s]   ATP exported   psi_m [mV]')
    for name, o2, adp, l, d, atp, psi in rows:
        print(f'{name:6s} {o2:9.4f} {adp:9.3f} {l:15.5f} {d:15.5f} {atp:18.5f} {psi * 1e3:10.1f}')
    rest_row = rows[2]
    assert rest_row[4] == pytest.approx(rest_row[3], rel=1e-5), 'calibrated alike at rest'
    o2_rows, adp_rows = rows[:4], rows[4:]
    for group in (o2_rows, adp_rows):
        for k in (3, 4):
            assert all(a[k] < b[k] for a, b in zip(group, group[1:])), 'both rise with O2 and with ADP'
    assert all(a[6] < b[6] for a, b in zip(adp_rows, adp_rows[1:])), 'ADP depolarises the matrix'
    (_, _, _, l_low, d_low, _, _), (_, _, _, l_k, d_k, _, _) = o2_rows[:2]
    assert d_low > 0.8 * rest_row[4] and l_low < 0.3 * rest_row[3], \
        'cytochrome oxidase keeps respiring down to a few uM O2; the lumped K_O2 (a tissue P50) does not'
    (_, _, _, l_hi, d_hi, _, _) = adp_rows[-1]
    assert d_hi < 0.8 * l_hi, 'the detailed one saturates sooner with ADP (its dehydrogenases and pumps at 4x rest)'


def _celldl_rates(ttl):
    """Each store's rate of change [mol/s] in a CellDL document's bond graph, as the BG-RDF templates define it: a store's
    exp(u/RT) is K q (times k_sat/(k_sat + q) for a saturating one), a chemical potential source's exp(u/RT) its value's;
    a reaction's Pf and Pr are the products over its substrates and products of exp(u/RT) to its bonds' bondCount, and
    it flows kappa (Pf - Pr), over 1 + (Pf + Pr)/sigma when saturating."""
    import rdflib
    bgf = rdflib.Namespace('https://bg-rdf.org/ontologies/bondgraph-framework#')
    g = rdflib.Graph().parse(data=ttl, format='turtle', publicID='https://example.org/document')

    def number(literal):
        return float(str(literal).split()[0])

    def parameter(e, name):
        for b in g.objects(e, bgf.parameterValue):
            if str(g.value(b, bgf.varName)) == name:
                return number(g.value(b, bgf.hasValue))

    kind = {e: str(t).rsplit('#', 1)[1] for e, t in g.subject_objects(rdflib.RDF.type) if str(t).startswith(str(bgf))}
    activity = {}
    for e, k in kind.items():
        if k in ('ChemicalStorageNode', 'SaturatingChemicalStorageNode'):
            q = number(g.value(e, bgf.hasValue))
            activity[e] = parameter(e, 'K') * q
            if k == 'SaturatingChemicalStorageNode':
                activity[e] *= parameter(e, 'k_sat') / (parameter(e, 'k_sat') + q)
        elif k == 'ChemicalPotential':
            activity[e] = np.exp(number(g.value(e, bgf.hasValue)) / celldl.TEMPLATE_RT)
    sides = {e: ([], []) for e, k in kind.items() if 'Reaction' in k}
    for b in g.subjects(bgf.hasSource, None):
        a, t = g.value(b, bgf.hasSource), g.value(b, bgf.hasTarget)
        n = int(g.value(b, bgf.bondCount) or 1)
        if t in sides:
            sides[t][0].append((a, n))
        else:
            sides[a][1].append((t, n))
    rate = {e: 0.0 for e, k in kind.items() if 'Storage' in k}
    for r, (sub, prod) in sides.items():
        pf = np.prod([activity[e] ** n for e, n in sub])
        pr = np.prod([activity[e] ** n for e, n in prod])
        v = parameter(r, 'kappa') * (pf - pr)
        if kind[r] == 'SaturatingChemicalReactionNode':
            v /= 1.0 + (pf + pr) / parameter(r, 'sigma')
        for e, n in sub:
            if e in rate:
                rate[e] -= n * v
        for e, n in prod:
            if e in rate:
                rate[e] += n * v
    return {str(e).rsplit('#', 1)[1][2:]: x for e, x in rate.items()}


@pytest.mark.parametrize('module', list(MODULES))
def test_the_celldl_document_has_the_modules_resting_rates(module, built):
    """Evaluated as the BG-RDF templates define it (what CellDL's bg2cellml generates), the CellDL document gives every
    store the module's own resting rate of change, S v times the volume (to 1e-9 of the largest flow): its stoichiometry,
    the stores' K, the chemostats' potentials and the rate constants (with saturation, regulation and psi_m folded in)
    say the same as the module at rest."""
    bg = built[module][0]
    flux, volume = MODULES[module]
    want = {s.name: x * volume for s, x in zip(bg.species, bg.S() @ [flux[r.name] for r in bg.reactions])
            if not s.boundary}
    got = _celldl_rates(_document(module)[1])
    assert set(got) == set(want)
    scale = max(abs(x) for x in flux.values()) * volume   # at rest the species are steady: scale by the turnover
    assert max(abs(got[s] - want[s]) for s in want) <= 1e-9 * scale


BG_TOOLS = os.environ.get('BG_TOOLS')   # a BG-Tools checkout (with constitutive-law templates) and its .venv


@pytest.mark.skipif(not BG_TOOLS, reason='set BG_TOOLS to a BG-Tools checkout with its .venv to run bg2cellml')
@pytest.mark.parametrize('module', ['glycogen', 'gluconeogenesis'])
def test_bg2cellml_makes_the_module_from_its_celldl_document(module, tmp_path):
    """CellDL's bg2cellml turns the module's bond graph into CellML that runs as the module does, at rest and from the
    verification config's run_parameters, to 1e-8 of each species' scale (tools/compare_celldl.py). Glycolysis and the
    mitochondrion match at rest only: their allosteric activation and dynamic matrix potential are not yet expressible
    in the templates."""
    import subprocess
    for name, spec in celldl_build.specs(module).items():
        with open(tmp_path / f'{name}.json', 'w') as f:
            json.dump(spec, f)
        (tmp_path / f'{name}.ttl').write_text(spec['turtle'])
        script = ('import asyncio, sys; from pathlib import Path; sys.path.insert(0, sys.argv[1]); '
                  'from bg2cellml.bondgraph.framework import get_framework\n'
                  'async def main():\n    fw = await get_framework(); src = Path(sys.argv[2])\n'
                  '    m = fw.make_bondgraph_model(src.as_uri(), src.read_text()); assert not m.has_issues\n'
                  '    Path(sys.argv[3]).write_text(m.make_cellml_model().to_xml())\nasyncio.run(main())')
        subprocess.run([os.path.join(BG_TOOLS, '.venv', 'bin', 'python'), '-c', script, BG_TOOLS,
                        str(tmp_path / f'{name}.ttl'), str(tmp_path / f'{name}.cellml')], check=True, cwd=BG_TOOLS)
    with contextlib.redirect_stdout(io.StringIO()):
        result = celldl_compare.compare(str(tmp_path), module)
    for name, r in result.items():
        assert r['max_relative_difference'] < 1e-8, (name, r['species'])


@pytest.fixture(scope='module')
def gdet(tmp_path_factory):
    return checks.ComponentModel(load_version('glycogen', VERSION), work_dir=str(tmp_path_factory.mktemp('gdet')))


@pytest.fixture(scope='module')
def ndet(tmp_path_factory):
    return checks.ComponentModel(load_version('gluconeogenesis', VERSION), work_dir=str(tmp_path_factory.mktemp('ndet')))


def test_the_glycogen_steps_balance_and_sum_to_the_lumped_store():
    """Each step balances in elements and charge; their standard Gibbs energies sum to the lumped store's (G6P + ATP ->
    a glucose unit + ADP + 2 Pi, dG0_store - dG0_hk); at the comparison point the intermediates are steady and every
    step runs down its Gibbs energy (phosphoglucomutase backwards: the store is being spent there)."""
    bg = detailed.glycogen_graph()
    assert bg.element_imbalance() == {} and bg.charge_imbalance() == {}
    store = sum(r.dG0 for r in bg.reactions if r.name in ('PGM', 'GLGC', 'GLGA', 'PPA'))
    assert store == pytest.approx(detailed.CYTO['dG0_store'] - detailed.CYTO['dG0_hk'], abs=1.0)
    steady = bg.steady(detailed.GLYC_FLUX)
    assert all(abs(steady[s]) < 1e-15 for s in ('G1P', 'ADPglc', 'PPi'))
    rest = detailed.glycogen_rest(bg)
    for r in bg.reactions:
        assert bg.reaction_gibbs(r, rest) * detailed.GLYC_FLUX[r.name] < 0.0, r.name


def test_the_glycogen_store_rests_at_the_lumped_rates(gdet):
    """At the comparison point the intermediates stay put for a second and glycogen changes at the lumped store's net
    rate (synthesis less phosphorylase)."""
    outs = ['C_G1P', 'C_ADPglc', 'C_PPi', 'C_Gly', 'v_GLGA', 'v_GLGP']
    t, r = _run(gdet, {}, 1.0, 0.01, outs)
    for k in ('C_G1P', 'C_ADPglc', 'C_PPi'):
        np.testing.assert_allclose(r[k], r[k][0], rtol=1e-4, err_msg=k)
    assert r['v_GLGA'][0] == pytest.approx(detailed.V_STORE, rel=1e-6)
    assert r['v_GLGP'][0] == pytest.approx(detailed.V_MOB, rel=1e-6)
    assert (r['C_Gly'][-1] - r['C_Gly'][0]) / (t[-1] - t[0]) == pytest.approx(detailed.V_STORE - detailed.V_MOB, rel=1e-3)


def _glycogen_net(gdet, glyc, g6p, atp, adp):
    '''Net storage after ten minutes (pyrophosphate and ADP-glucose settle in that time; glycogen moves a few mM).'''
    _, r = _run(gdet, {'C_Gly_init': glyc, 'C_G6P': g6p, 'C_ATP': atp, 'C_ADP': adp}, 600.0, 60.0, ['v_GLGA', 'v_GLGP'])
    return r['v_GLGA'][-1] - r['v_GLGP'][-1]


def test_compared_with_the_lumped_glycogen_store(gdet):
    """With the cytosol clamped alike (the lumped store fully allocated), after ten minutes: the same net rate at the
    comparison point; both store more with more G6P and spend the store as ATP falls (the phosphorylase's ATP
    inhibition, in both), and slow as glycogen runs out. They differ in how they get there: the detailed store
    responds to G6P less (pyrophosphate builds up behind its pyrophosphorylase, as its pyrophosphatase is mass action,
    and the phosphorylase runs closer to equilibrium); the lumped store slows near its ceiling ((C/200 mM)^8) while
    the detailed steps, which have no ceiling, slow only as glycogen's activity saturates."""
    c = detailed.CYTO
    rows = []
    for name, glyc, g6p, atp, adp in (('rest', 50.0, c['C_G6P_init'], c['C_ATP_init'], c['C_ADP_init']),
                                      ('G6P x5', 50.0, 5.0, c['C_ATP_init'], c['C_ADP_init']),
                                      ('G6P x20', 50.0, 20.0, c['C_ATP_init'], c['C_ADP_init']),
                                      ('ATP 2, ADP 3', 50.0, c['C_G6P_init'], 2.0, 3.0),
                                      ('little glycogen', 2.0, c['C_G6P_init'], c['C_ATP_init'], c['C_ADP_init']),
                                      ('near the ceiling, G6P x20', 190.0, 20.0, c['C_ATP_init'], c['C_ADP_init'])):
        store, mob = detailed.lumped_store_rates(g6p, atp, adp, c['C_Pi_init'], glyc)
        rows.append((name, store - mob, _glycogen_net(gdet, glyc, g6p, atp, adp)))
    print('\ncondition                     lumped net storage   detailed [mM/s]')
    for name, l, d in rows:
        print(f'{name:28s} {l:18.5f} {d:12.5f}')
    net = {name: (l, d) for name, l, d in rows}
    assert net['rest'][1] == pytest.approx(net['rest'][0], rel=0.02), 'calibrated alike (glycogen drifts in ten minutes)'
    for k in (0, 1):
        assert net['rest'][k] < net['G6P x5'][k] < net['G6P x20'][k], 'both store more with G6P'
        assert net['ATP 2, ADP 3'][k] < net['rest'][k], 'and spend it as ATP falls'
        assert abs(net['little glycogen'][k]) < abs(net['rest'][k]), 'the phosphorylase slows as glycogen runs out'
    assert net['G6P x20'][1] < net['G6P x20'][0], 'the detailed store responds to G6P less'
    assert net['near the ceiling, G6P x20'][0] < 0.7 * net['G6P x20'][0], 'the lumped store slows near its ceiling'
    assert net['near the ceiling, G6P x20'][1] > 0.75 * net['G6P x20'][1], \
        'the detailed one has no ceiling (it slows only as glycogen\'s activity saturates)'


def test_the_gluconeogenic_bypasses_close_futile_cycles_with_glycolysis(built):
    """Each step balances; fructose 1,6-bisphosphatase with phosphofructokinase, and glucose 6-phosphatase with
    hexokinase (the glycolysis module's), sum to ATP hydrolysis (the cycles they would form if both ran); at rest every
    step runs down its Gibbs energy; and the dehydrogenase takes glycerol 3-phosphate to DHAP only where NADH/NAD+ is
    low: at the fermenter's NADH/NAD+ and pH its Gibbs energy that way is positive."""
    bg = built['gluconeogenesis'][0]
    assert bg.element_imbalance() == {} and bg.charge_imbalance() == {}
    gng = {r.name: r.dG0 for r in bg.reactions}
    gly = {r.name: r.dG0 for r in detailed.GLY_REACTIONS}
    hydrolysis = detailed.CYTO['dG0_hyd']
    assert gly['PFK'] + gng['FBPASE'] == pytest.approx(hydrolysis, abs=200.0)
    assert gly['HK'] + gng['G6PASE'] == pytest.approx(hydrolysis, abs=200.0)
    rest = detailed.gluconeogenesis_rest(bg)
    for r in bg.reactions:
        assert bg.reaction_gibbs(r, rest) < 0.0, r.name
    g3pdh = {r.name: r for r in bg.reactions}['G3PDH']
    fermenter = dict(rest, NAD=detailed.GLY_REST['NAD'], NADH=detailed.GLY_REST['NADH'], H=detailed.GLY_REST['H'],
                     DHAP=detailed.GLY_REST['DHAP'])
    assert bg.reaction_gibbs(g3pdh, fermenter) > 0.0


def test_gluconeogenesis_rests_and_takes_up_more_glycerol_when_given_it(ndet):
    """At its default glycerol 3-phosphate stays put for a minute and glucose is made at the set rate; with ten times the
    glycerol, more glycerol is taken up (and passed on as DHAP)."""
    _, r = _run(ndet, {}, 60.0, 1.0, ['C_Gly3P', 'j_glc', 'j_glycerol', 'j_DHAP'])
    np.testing.assert_allclose(r['C_Gly3P'], r['C_Gly3P'][0], rtol=1e-7)
    assert r['j_glc'][-1] == pytest.approx(detailed.J_GNG, rel=1e-6)
    _, r10 = _run(ndet, {'C_glycerol': 10.0 * detailed.GNG_REST['glycerol']}, 60.0, 1.0, ['j_glycerol', 'j_DHAP'])
    assert -r10['j_glycerol'][-1] > -r['j_glycerol'][-1] and r10['j_DHAP'][-1] == pytest.approx(-r10['j_glycerol'][-1],
                                                                                               rel=1e-6)



@pytest.fixture(scope='module')
def cytosol_point():
    return lumped.evaluate(*lumped.CYTOSOL, lumped.CYTOSOL_EXPORT, lumped.CYTOSOL_RATES, lumped.CYTOSOL_STATES)


def test_the_lumped_cytosols_bond_graph_is_its_stoichiometry(cytosol_point):
    """The bond graph drawn for energy_metabolism adaptive has the module's own stoichiometry: at the export point (every
    reaction running) each state's rate of change in the module, less its transport (the harness's membrane face), is
    S v from the module's own reaction rates; and its BG-RDF, evaluated as the chemical template defines it, gives those rates of change."""
    v, deriv, conc = cytosol_point
    bg = lumped.cytosol_graph()
    net = bg.steady(v)
    names = {f'{s}:{n}': n for s, ns in lumped.CYTOSOL_TRANSPORT.items() for n in ns}
    j, _, _ = lumped.evaluate(*lumped.CYTOSOL, lumped.CYTOSOL_EXPORT, names, {})
    for s, d in deriv.items():
        transport = sum(j[f'{s}:{n}'] for n in lumped.CYTOSOL_TRANSPORT[s]) / detailed.CYTO['V_c']
        assert d - transport == pytest.approx(net[s], rel=1e-9, abs=1e-15), s
    got = _celldl_rates(_document('energy_metabolism_adaptive')[1])
    stores = {s.name for s in bg.species if not s.boundary}
    assert set(got) == stores
    scale = max(abs(net[s]) for s in stores)
    assert max(abs(got[s] - net[s] * detailed.CYTO['V_c']) for s in stores) <= 1e-9 * scale * detailed.CYTO['V_c']


def test_the_lumped_cytosols_affinities_are_its_bond_graphs(cytosol_point):
    """Each reaction's affinity in the module (A_ = RT (ln K_eq - ln Gamma)) is the bond graph's at the export point
    (the hexokinase step, the lower glycolysis, ATP hydrolysis and the store; the matrix export is irreversible)."""
    _, _, conc = cytosol_point
    bg = lumped.cytosol_graph()
    affinity = _cytosol_affinities()
    for r in bg.reactions:
        if r.name in ('HK', 'LG', 'HYD', 'STORE'):
            assert affinity[r.name] == pytest.approx(-bg.reaction_gibbs(r, dict(conc, mx=1.0)), rel=1e-6), r.name


@pytest.mark.xfail(strict=True, reason='energy_metabolism adaptive: the phosphorylase\'s equilibrium (Gamma_mob = G6P/Pi) '
                                       'leaves out glycogen\'s activity a_glyc, which the store\'s includes, so the store '
                                       'and phosphorylase cycle does not have ATP hydrolysis\' equilibrium while a_glyc < 1')
def test_the_lumped_phosphorylases_affinity_is_its_bond_graphs(cytosol_point):
    """The phosphorylase's affinity in the module is the bond graph's (glycogen with the activity the store gives it),
    so the store and phosphorylase cycle sums to ATP hydrolysis (Wegscheider)."""
    _, _, conc = cytosol_point
    mob = {r.name: r for r in lumped.cytosol_graph().reactions}['MOB']
    assert _cytosol_affinities()['MOB'] == pytest.approx(-lumped.cytosol_graph().reaction_gibbs(mob, dict(conc, mx=1.0)),
                                                         rel=1e-6)


def _cytosol_affinities():
    names = {'HK': 'A_hk', 'LG': 'A_lg', 'HYD': 'A_hyd', 'STORE': 'A_store', 'MOB': 'A_mob'}
    a, _, _ = lumped.evaluate(*lumped.CYTOSOL, lumped.CYTOSOL_EXPORT, names, {})
    return a


def test_cyclic_pathways_are_drawn_as_rings(built):
    """The CellDL layout (tools/celldl.py) draws the longest cycle of reactions linked by their own intermediates as a
    ring, from its earliest reaction: the TCA cycle from citrate synthase, the glycogen store's pyrophosphorylase, synthase and phosphorylase;
    glycolysis' NAD+/NADH loop is a cofactor's and not a ring, and the lumped cytosol has none."""
    assert built['mitochondrion'][0].cycle() == ['CS', 'ACO', 'IDH', 'OGDH', 'SCS', 'CII', 'FUM', 'MDH']
    assert built['glycogen'][0].cycle() == ['GLGC', 'GLGA', 'GLGP']
    assert built['glycolysis'][0].cycle() == [] and built['gluconeogenesis'][0].cycle() == []
    assert lumped.cytosol_graph().cycle() == []
