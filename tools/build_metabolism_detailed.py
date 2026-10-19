"""Builds detailed metabolism modules: each step of a pathway as its own reaction, to compare against the lumped
modules of the generic adaptive cell (tools/build_adaptive_cell.py), which make simplifying assumptions (energy_metabolism
adaptive: glycolysis as two steps). The pathways follow the carbohydrate-metabolism model the designer supplied
(CarbohydrateMetabolism.cellml: glycolysis, glycogen, gluconeogenesis, the mitochondrial carriers, the TCA cycle and the
electron transport chain, element-balanced), here split into modules with physical units, thermodynamically consistent
kinetics and a bond-graph description each (tools/bondgraph.py: the module's equations, its BG-RDF and its diagram).

    python tools/build_metabolism_detailed.py

Module versions (all detailed_ArgusUNPUBLISHED_v01):

  cell/glycolysis        glucose to lactate in eleven steps (hexokinase to lactate dehydrogenase), NAD+/NADH inside;
                         the cytosol's glucose, lactate, ATP, ADP, Pi, H+ and water as chemostats; calibrated to the
                         fermenter's resting flux and concentrations (compared with energy_metabolism adaptive)
  cell/mitochondrion     cytosolic pyruvate to CO2: the pyruvate and phosphate carriers, pyruvate dehydrogenase, the TCA
                         cycle, complexes I-IV with the Q and cytochrome c pools, ATP synthase, the translocase and a
                         proton leak, the matrix potential a state on the inner membrane's capacitance; calibrated to
                         the lumped respiration's resting O2 use (compared with energy_metabolism mitochondrial)
  cell/glycogen          the glycogen store by the ADP-glucose route (phosphoglucomutase, GlgC, GlgA, pyrophosphatase,
                         GlgP), glycogen with a saturating activity; calibrated to the lumped store's rates (compared
                         with energy_metabolism adaptive's store)
  cell/gluconeogenesis   the gluconeogenic bypasses (fructose 1,6-bisphosphatase, glucose 6-phosphatase) and glycerol's
                         way in (glycerol kinase, glycerol 3-phosphate dehydrogenase), in a respiring cytosol

Each version also gets <module>_<version>_bondgraph.ttl (BG-RDF) and _bondgraph.svg (a CellDL-style diagram carrying
that RDF). tools/build_bondgraph_lumped.py writes the same for lumped modules.

These modules change nothing in the game's cell: they are references to test the lumped ones against
(tests/test_metabolism_detailed.py).
"""
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import module_writer  # noqa: E402
from module_writer import REPO, register_bib, write_version  # noqa: E402
from bondgraph import BondGraph, Reaction, Species, FARADAY, R_GAS  # noqa: E402

VERSION = 'detailed_ArgusUNPUBLISHED_v01'
PH = 'unreferenced_placeholder'
TEXTBOOK = f'{PH}; textbook transformed standard Gibbs energy at pH 7 (Lehninger, from memory, unchecked)'
module_writer.UNITS.update({'W_per_m3': [('watt', None, 1, 1), ('metre', None, -3, 1)],
                            'F_per_m3': [('farad', None, 1, 1), ('metre', None, -3, 1)]})

register_bib({
    'argus2026unpublished': dict(
        type='unpublished', author='Argus, Finbar', title='Unpublished modules of the circulatory-autogen module library',
        year='2026', note='Placeholder citation (required_citations, 2026-10-06): replace with the real reference when it '
                          'exists: the owner\'s unpublished or modified models, no publication'),
    'bennett2009absolute': dict(type='article',
                                author='Bennett, Bryson D. and Kimball, Elizabeth H. and Gao, Melissa and Osterhout, '
                                       'Robin and Van Dien, Stephen J. and Rabinowitz, Joshua D.',
                                title='Absolute metabolite concentrations and implied enzyme active site occupancy in '
                                      'Escherichia coli',
                                journal='Nature Chemical Biology', year='2009', volume='5', number='8', pages='593--599',
                                doi='10.1038/nchembio.186'),
})


def csv_values(path):
    with open(os.path.join(REPO, path)) as f:
        return {r['variable_name']: float(r['value']) for r in csv.DictReader(f) if r['value'] not in ('', None)}


# ------------------------------------------------------------------ the fermenter's resting cytosol (the lumped modules)
CYTO = csv_values('modules/cell/energy_metabolism/versions/adaptive_ArgusUNPUBLISHED_v01/parameterisations/default/default_parameters.csv')
BIO = csv_values('modules/cell/biosynthesis/versions/adaptive_ArgusUNPUBLISHED_v01/parameterisations/default/default_parameters.csv')
T = CYTO['T']
PH_I = 7.6   # the fermenter's cytosol (generic_adaptive_cell parameters pH_i)
J_REST = BIO['J_ref'] / CYTO['V_c'] / 2.0   # mM/s of glucose: the resting ATP turnover, two ATP per glucose

F = lambda **k: k   # noqa: E731
GLY_SPECIES = [
    Species('glc', 'glucose', F(C=6, H=12, O=6), boundary=True),
    Species('G6P', 'glucose 6-phosphate', F(C=6, H=11, O=9, P=1), charge=-2),
    Species('F6P', 'fructose 6-phosphate', F(C=6, H=11, O=9, P=1), charge=-2),
    Species('F16BP', 'fructose 1,6-bisphosphate', F(C=6, H=10, O=12, P=2), charge=-4),
    Species('DHAP', 'dihydroxyacetone phosphate', F(C=3, H=5, O=6, P=1), charge=-2),
    Species('GAP', 'glyceraldehyde 3-phosphate', F(C=3, H=5, O=6, P=1), charge=-2),
    Species('BPG13', '1,3-bisphosphoglycerate', F(C=3, H=4, O=10, P=2), charge=-4),
    Species('PG3', '3-phosphoglycerate', F(C=3, H=4, O=7, P=1), charge=-3),
    Species('PG2', '2-phosphoglycerate', F(C=3, H=4, O=7, P=1), charge=-3),
    Species('PEP', 'phosphoenolpyruvate', F(C=3, H=2, O=6, P=1), charge=-3),
    Species('Pyr', 'pyruvate', F(C=3, H=3, O=3), charge=-1),
    Species('lac', 'lactate', F(C=3, H=5, O=3), boundary=True, charge=-1),
    Species('NAD', 'NAD+', F(C=21, H=26, N=7, O=14, P=2), charge=-1),
    Species('NADH', 'NADH', F(C=21, H=27, N=7, O=14, P=2), charge=-2),
    Species('ATP', 'ATP', F(C=10, H=12, N=5, O=13, P=3), boundary=True, charge=-4),
    Species('ADP', 'ADP', F(C=10, H=12, N=5, O=10, P=2), boundary=True, charge=-3),
    Species('Pi', 'phosphate (HPO4 2-)', F(H=1, O=4, P=1), boundary=True, charge=-2),
    Species('H', 'H+', F(H=1), boundary=True, charge=1, c_std=1e-4),
    Species('H2O', 'water', F(H=2, O=1), boundary=True, c_std=55500.0),
]
GLY_REACTIONS = [
    Reaction('HK', 'hexokinase', {'glc': -1, 'ATP': -1, 'G6P': 1, 'ADP': 1, 'H': 1}, -16.7e3, TEXTBOOK),
    Reaction('PGI', 'glucose-6-phosphate isomerase', {'G6P': -1, 'F6P': 1}, 1.7e3, TEXTBOOK),
    Reaction('PFK', 'phosphofructokinase', {'F6P': -1, 'ATP': -1, 'F16BP': 1, 'ADP': 1, 'H': 1}, -14.2e3, TEXTBOOK),
    Reaction('ALD', 'aldolase', {'F16BP': -1, 'DHAP': 1, 'GAP': 1}, 23.8e3, TEXTBOOK),
    Reaction('TPI', 'triose-phosphate isomerase', {'DHAP': -1, 'GAP': 1}, 7.5e3, TEXTBOOK),
    Reaction('GAPDH', 'glyceraldehyde-3-P dehydrogenase', {'GAP': -1, 'Pi': -1, 'NAD': -1, 'BPG13': 1, 'NADH': 1, 'H': 1},
             6.3e3, TEXTBOOK),
    Reaction('PGK', 'phosphoglycerate kinase', {'BPG13': -1, 'ADP': -1, 'PG3': 1, 'ATP': 1}, -18.5e3, TEXTBOOK),
    Reaction('PGM', 'phosphoglycerate mutase', {'PG3': -1, 'PG2': 1}, 4.4e3, TEXTBOOK),
    Reaction('ENO', 'enolase', {'PG2': -1, 'PEP': 1, 'H2O': 1}, 1.8e3,
             f'{PH}; textbook +1.8 kJ/mol (Berg, from memory; Lehninger\'s +7.5 cannot carry the flux at the cell\'s rest)'),
    # activated feed-forward by fructose 1,6-bisphosphate (E. coli PykF), so the lower pathway keeps pace as flux rises
    Reaction('PK', 'pyruvate kinase', {'PEP': -1, 'ADP': -1, 'H': -1, 'Pyr': 1, 'ATP': 1}, -31.4e3, TEXTBOOK,
             activators={'F16BP': 0.2}),
    Reaction('LDH', 'lactate dehydrogenase', {'Pyr': -1, 'NADH': -1, 'H': -1, 'lac': 1, 'NAD': 1}, -25.1e3, TEXTBOOK),
]
# resting concentrations [mM]: the cytosol's from the lumped module; F1,6BP, pyruvate and NAD+/NADH E. coli-like
# (Bennett et al. 2009, from memory); the other intermediates from the near-equilibrium steps (near_equilibrium_rest)
NEAR_EQ = ('PGI', 'ALD', 'TPI', 'GAPDH', 'PGK', 'PGM', 'ENO')
EXCHANGE = 1000.0   # a near-equilibrium step's forward rate over its net rate at rest: Gibbs energy RT ln(1 - 1/1000)
GLY_REST = dict(glc=CYTO['C_glc_init'], G6P=CYTO['C_G6P_init'], F16BP=15.2, Pyr=0.37, lac=CYTO['C_lac_init'], NAD=2.55,
                NADH=0.0834, ATP=CYTO['C_ATP_init'], ADP=CYTO['C_ADP_init'], Pi=CYTO['C_Pi_init'],
                H=1e3 * 10 ** -PH_I, H2O=55500.0)
GLY_REST_REF = {
    'glc': 'the fermenter\'s resting cytosol (energy_metabolism adaptive)', 'G6P': 'the fermenter\'s resting cytosol',
    'lac': 'the fermenter\'s resting cytosol', 'ATP': 'the fermenter\'s resting cytosol', 'ADP': 'the fermenter\'s '
    'resting cytosol', 'Pi': 'the fermenter\'s resting cytosol', 'H': 'the fermenter\'s cytosol pH 7.6',
    'H2O': 'definitional; pure water, 55.5 M', 'F16BP': f'{PH}; Bennett 2009 (from memory): 15.2 mM',
    'GAP': 'calibrated: aldolase and triose-phosphate isomerase near equilibrium', 'BPG13': 'calibrated: GAPDH near '
    'equilibrium', 'PG3': 'calibrated: PGK near equilibrium', 'PG2': 'calibrated: PGM near equilibrium',
    'PEP': 'calibrated: enolase near equilibrium', 'F6P': 'calibrated: glucose-6-phosphate isomerase near equilibrium',
    'DHAP': 'calibrated: aldolase and triose-phosphate isomerase near equilibrium',
    'Pyr': f'{PH}; E. coli-like (from memory)', 'NAD': f'{PH}; Bennett 2009 (from memory): 2.55 mM',
    'NADH': f'{PH}; Bennett 2009 (from memory): 0.083 mM'}
GLY_FLUX = {r.name: (J_REST if r.name in ('HK', 'PGI', 'PFK', 'ALD', 'TPI') else 2.0 * J_REST) for r in GLY_REACTIONS}
# far from equilibrium, these saturate: each can carry HEADROOM times its resting flux at saturating substrates
# (sigma = H (Pf - Pr) - (Pf + Pr) at rest gives kappa sigma = H v_rest), as the lumped glycolysis can (its hexokinase
# step about 53 times the resting flux at saturating glucose)
SATURATED = ('HK', 'PFK', 'PK', 'LDH')
HEADROOM = 50.0


def near_equilibrium_rest(bg):
    """The intermediates the near-equilibrium steps set at rest, down the chain from G6P and F1,6BP: each step's
    forward rate is EXCHANGE times its net rate."""
    c = dict(GLY_REST)
    rx = {r.name: r for r in bg.reactions}
    c['F6P'] = bg.near_equilibrium(rx['PGI'], c, 'F6P', EXCHANGE)
    # aldolase and the isomerase together: DHAP GAP and GAP / DHAP from their two ratios
    c['DHAP'] = 1.0
    gap_per_dhap = bg.near_equilibrium(rx['TPI'], c, 'GAP', EXCHANGE)
    c['GAP'] = 1.0
    product = bg.near_equilibrium(rx['ALD'], dict(c, DHAP=1.0), 'GAP', EXCHANGE)   # DHAP x GAP with DHAP = 1
    c['DHAP'] = math.sqrt(product / gap_per_dhap)
    c['GAP'] = gap_per_dhap * c['DHAP']
    for name, x in (('GAPDH', 'BPG13'), ('PGK', 'PG3'), ('PGM', 'PG2'), ('ENO', 'PEP')):
        c[x] = bg.near_equilibrium(rx[name], c, x, EXCHANGE)
    return c


def unswept_potentials(bg):
    """The standard potentials, left out of the boundary-condition sweep: only their combinations (the reactions'
    Gibbs energies) are physical, and a relative change of one depends on the gauge they are fitted in."""
    return [f'mu0_{s.name}' for s in bg.species]


def module_parts(bg, rest, rest_ref, kappa, gibbs, where, flux_ref, sigma_ref, ka_ref, t_ref):
    """The variables, roles and parameters of a bond-graph module: R, T (and F, the compartments' potentials and
    capacitances when it has any), per species mu0, c_std and either its concentration from a port (a chemostat, with its
    net production j) or its state, per reaction kappa, Pf, Pr, v, A (sigma and K_a when it saturates or is activated).
    rest_ref: species -> source of its rest value (and k_sat_<species>, psi_<c>, C_m_<c> where used); where(species):
    where a chemostat's concentration comes from; flux_ref: what each kappa is calibrated to; sigma_ref(r) and
    ka_ref(r, species) (activation and inhibition constants): sources."""
    variables = [('t', 'second', 'in', None, None), ('R', 'J_per_mol_K', 'in', None, None),
                 ('T', 'kelvin', 'in', None, None)]
    roles = {'R': 'global_constant', 'T': 'global_constant'}
    params = [('R', R_GAS, 'definitional; gas constant', True), ('T', bg.temperature, t_ref, False)]
    if bg.potentials:
        variables.append(('F', 'C_per_mol', 'in', None, None))
        roles['F'] = 'global_constant'
        params.append(('F', FARADAY, 'definitional; Faraday constant', True))
    for comp, psi in bg.potentials.items():
        if comp in bg.capacitance:
            variables += [(f'psi_{comp}_init', 'volt', 'in', None, None),
                          (f'psi_{comp}', 'volt', 'out', f'psi_{comp}_init', f'electrical potential of compartment {comp}'),
                          (f'C_m_{comp}', 'F_per_m3', 'in', None, 'membrane capacitance per volume'),
                          (f'I_{comp}', 'millimolar_per_s', 'out', None, f'charge moved into {comp} (elementary charges)')]
            roles.update({f'psi_{comp}_init': 'constant', f'psi_{comp}': 'variable', f'C_m_{comp}': 'constant',
                          f'I_{comp}': 'variable'})
            params += [(f'psi_{comp}_init', psi, rest_ref[f'psi_{comp}'], False),
                       (f'C_m_{comp}', bg.capacitance[comp], rest_ref[f'C_m_{comp}'], False)]
        else:
            variables.append((f'psi_{comp}', 'volt', 'in', None, f'electrical potential of compartment {comp}'))
            roles[f'psi_{comp}'] = 'boundary_condition'
            params.append((f'psi_{comp}', psi, rest_ref[f'psi_{comp}'], False))
    for s in bg.species:
        variables.append((f'mu0_{s.name}', 'J_per_mol', 'in', None, f'standard chemical potential of {s.label} (from the '
                          'reactions\' Gibbs energies, tools/bondgraph.py)'))
        variables.append((f'c_std_{s.name}', 'millimolar', 'in', None, None))
        roles[f'mu0_{s.name}'] = roles[f'c_std_{s.name}'] = 'constant'
        params += [(f'mu0_{s.name}', bg.mu0[s.name], 'calibrated: fitted to the reactions\' standard Gibbs energies (along '
                    'the conserved moieties, centred on the rest state)', False),
                   (f'c_std_{s.name}', s.c_std, 'definitional; standard state (1 M; H+ at pH 7; water itself)', True)]
        if s.boundary:
            variables.append((f'C_{s.name}', 'millimolar', 'in', None, f'{s.label} {where(s)}'))
            variables.append((f'j_{s.name}', 'millimolar_per_s', 'out', None, f'net production of {s.label} here'))
            roles[f'C_{s.name}'] = 'boundary_condition'
            roles[f'j_{s.name}'] = 'variable'
            params.append((f'C_{s.name}', rest[s.name], rest_ref[s.name], False))
        else:
            variables.append((f'C_{s.name}_init', 'millimolar', 'in', None, None))
            variables.append((f'C_{s.name}', 'millimolar', 'out', f'C_{s.name}_init', s.label))
            roles[f'C_{s.name}_init'] = 'constant'
            roles[f'C_{s.name}'] = 'variable'
            params.append((f'C_{s.name}_init', rest[s.name], rest_ref[s.name], False))
    for r in bg.reactions:
        variables += [(f'kappa_{r.name}', 'millimolar_per_s', 'in', None, f'{r.enzyme}: rate constant'),
                      (f'Pf_{r.name}', 'dimensionless', 'out', None, 'exp(forward affinity / RT)'),
                      (f'Pr_{r.name}', 'dimensionless', 'out', None, 'exp(reverse affinity / RT)'),
                      (f'v_{r.name}', 'millimolar_per_s', 'out', None, f'{r.enzyme}'),
                      (f'A_{r.name}', 'J_per_mol', 'out', None, 'affinity: minus the reaction\'s Gibbs energy')]
        roles.update({f'kappa_{r.name}': 'constant', f'Pf_{r.name}': 'variable', f'Pr_{r.name}': 'variable',
                      f'v_{r.name}': 'variable', f'A_{r.name}': 'variable'})
        params.append((f'kappa_{r.name}', kappa[r.name], f'calibrated: carries {flux_ref} '
                       f'(Gibbs energy at rest {gibbs[r.name] / 1e3:+.2f} kJ/mol; dG0 {r.dG0 / 1e3:+.1f} kJ/mol, {r.ref})',
                       False))
        for s, k in r.activators.items():
            variables.append((f'Ka_{s}_{r.name}', 'millimolar', 'in', None, f'activation of {r.enzyme} by {s}'))
            roles[f'Ka_{s}_{r.name}'] = 'constant'
            params.append((f'Ka_{s}_{r.name}', k, ka_ref(r, s), False))
        for s, k in r.inhibitors.items():
            variables.append((f'Ki_{s}_{r.name}', 'millimolar', 'in', None, f'inhibition of {r.enzyme} by {s}'))
            roles[f'Ki_{s}_{r.name}'] = 'constant'
            params.append((f'Ki_{s}_{r.name}', k, ka_ref(r, s), False))
        if r.sigma:
            variables.append((f'sigma_{r.name}', 'dimensionless', 'in', None, 'saturation (half-saturated at Pf = sigma)'))
            roles[f'sigma_{r.name}'] = 'constant'
            params.append((f'sigma_{r.name}', r.sigma, sigma_ref(r), False))
    for s in bg.species:
        if s.k_sat is not None:
            variables.append((f'k_sat_{s.name}', 'millimolar', 'in', None, f'saturation of {s.label}\'s activity'))
            roles[f'k_sat_{s.name}'] = 'constant'
            params.append((f'k_sat_{s.name}', s.k_sat, rest_ref[f'k_sat_{s.name}'], False))
    return variables, roles, params


def build_glycolysis():
    global GLY_REST
    bg = BondGraph('glycolysis', 'Detailed glycolysis', GLY_SPECIES, GLY_REACTIONS, temperature=T,
                   cofactors=('NAD', 'NADH'))
    GLY_REST = near_equilibrium_rest(bg)
    bg.centre_potentials(GLY_REST)
    assert not bg.element_imbalance() and not bg.charge_imbalance(), (bg.element_imbalance(), bg.charge_imbalance())
    steady = bg.steady(GLY_FLUX)
    internal = [s.name for s in bg.species if not s.boundary]
    assert all(abs(steady[s]) < 1e-15 for s in internal), steady
    for r in bg.reactions:
        if r.name in SATURATED:
            pf, pr = bg.affinity_terms(r, GLY_REST)
            r.sigma = HEADROOM * (pf - pr) - (pf + pr)
    kappa = bg.calibrate(GLY_REST, GLY_FLUX)
    gibbs = {r.name: bg.reaction_gibbs(r, GLY_REST) for r in bg.reactions}

    variables, roles, params = module_parts(
        bg, GLY_REST, GLY_REST_REF, kappa, gibbs, lambda s: 'in the cytosol (metabolite port)',
        'the fermenter\'s resting glycolytic flux',
        lambda r: f'{PH}; saturates at {HEADROOM:g} times its resting flux (the lumped glycolysis\' headroom)',
        lambda r, s: f'{PH}; feed-forward activation constant (E. coli PykF by F1,6BP, of order 0.1-1 mM, from memory)',
        'operating point: 20 C (the fermenter\'s)')
    variables.append(('P_diss', 'W_per_m3', 'out', None, 'free energy dissipated, sum v A (>= 0, second law)'))
    roles['P_diss'] = 'variable'
    eqs = bg.equations()

    invariants = [
        {'description': 'second law: every reaction runs down its affinity, v A >= 0, so P_diss >= 0',
         'expr': ' and '.join(f'np.all(v_{r.name}*A_{r.name} >= -1e-12*np.abs(v_{r.name}).max()*np.abs(A_{r.name}).max())'
                              for r in bg.reactions) + ' and np.all(P_diss >= -1e-9*np.abs(P_diss).max())'},
        {'description': 'NAD+ and NADH are made and used only here: their sum is conserved (relative 1e-9)',
         'expr': 'np.allclose(C_NAD + C_NADH, C_NAD[0] + C_NADH[0], rtol=1e-9, atol=0)'},
        {'description': 'the ports report each chemostat\'s net production: j_glc = -v_HK, j_lac = v_LDH, j_ATP = '
                        '-v_HK - v_PFK + v_PGK + v_PK (exact)',
         'expr': 'np.allclose(j_glc, -v_HK, rtol=1e-12, atol=0) and np.allclose(j_lac, v_LDH, rtol=1e-12, atol=0) and '
                 'np.allclose(j_ATP, -v_HK - v_PFK + v_PGK + v_PK, rtol=1e-12, atol=1e-30)'},
    ]
    outputs = [f'C_{s}' for s in internal] + [f'v_{r.name}' for r in bg.reactions] + \
        [f'j_{s.name}' for s in bg.species if s.boundary] + [f'A_{r.name}' for r in bg.reactions] + ['P_diss']
    ports = [{'port_type': 'cytosol_metabolites',
              'variables': [f'C_{s.name}' for s in bg.species if s.boundary] + [f'j_{s.name}' for s in bg.species
                                                                                 if s.boundary]}]
    notes = (
        'Glycolysis step by step, from glucose to lactate (homolactic fermentation), as a bond graph: each metabolite a '
        'store (C) on a 0-junction, each enzyme a reaction (Re) between 1-junctions, with every reaction balanced in C, '
        'H, O, N and P and in charge (species as at pH 7: G6P 2-, ATP 4-, ...; H+ and water explicit). The reactions\' '
        'standard Gibbs energies are textbook values (transformed, pH 7, 1 M); the species\' standard potentials mu0 are '
        'fitted to them (tools/bondgraph.py), so the thermodynamics are consistent, and H+ has its standard state at pH '
        '7 so the cell\'s pH 7.6 enters through its concentration. Each reaction flows at v = kappa (Pf - Pr) / (1 + '
        '(Pf + Pr)/sigma), Pf and Pr the exponentials of its forward and reverse affinities over RT: mass action, as the '
        'bond-graph chemical reaction, for the near-equilibrium steps, and saturating for the four far from equilibrium '
        '(hexokinase, phosphofructokinase, pyruvate kinase, lactate dehydrogenase; each saturates at 50 times its '
        'resting flux, the lumped glycolysis\' headroom); pyruvate kinase is also activated by fructose 1,6-bisphosphate, '
        'rate x (1 + F16BP/K_a), as E. coli\'s PykF is, so the lower pathway keeps pace as the flux rises (it changes the '
        'rate, not the equilibrium). The seven near-equilibrium steps run forwards at 1000 times their net rate at '
        'rest, which sets their intermediates (F6P, DHAP, GAP, 1,3-BPG, 3PG, 2PG, PEP) from G6P and F1,6BP: they come '
        'out close to E. coli\'s. Calibrated to '
        'the fermenter\'s rest: the cytosol\'s glucose, lactate, ATP, ADP and Pi (energy_metabolism adaptive), E. coli-like '
        'intermediates where they let every step run forwards and chosen ones where they do not (GAP, 1,3-BPG, 3PG, 2PG, '
        'PEP), and its glycolytic flux (half the resting ATP turnover): each kappa carries that flux at rest. NAD+/NADH '
        'cycle inside (GAPDH makes NADH, LDH uses it). Glucose, lactate, ATP, ADP, Pi, H+ and water are chemostats on '
        'the cytosol_metabolites port, which reports each one\'s net production. Compared against the lumped glycolysis '
        'of energy_metabolism adaptive (two steps) in tests/test_metabolism_detailed.py. Its bond graph is in '
        f'glycolysis_{VERSION}_bondgraph.ttl (BG-RDF) and drawn in glycolysis_{VERSION}_bondgraph.svg.')
    write_version('cell', 'glycolysis', 'glycolysis_detailed', variables, eqs,
                  {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
                   'general_ports': ports, 'roles': roles},
                  params,
                  {'sim_time': 600.0, 'pre_time': 0.0, 'dt': 0.1, 'solver': 'CVODE_myokit',
                   'solver_info': {'rtol': 1e-10, 'atol': 1e-12},
                   'outputs': outputs,
                   # off rest: ten times the glucose, half the ATP
                   'run_parameters': {'C_glc': 10.0 * GLY_REST['glc'], 'C_ATP': 0.5 * GLY_REST['ATP']},
                   'invariants': invariants,
                   'bc_sweep': {'ranges': {'C_glc': [1e-4, 10.0], 'C_ATP': [1.0, 12.0], 'C_ADP': [0.1, 3.0],
                                           'C_Pi': [1.0, 20.0], 'C_lac': [0.01, 20.0], 'C_H': [1e-5, 1.6e-4]},
                                'exclude': ['R', 'T', 'C_H2O'] + unswept_potentials(bg),
                                'rationale': 'glucose from starving to a meal; ATP and ADP from a de-energised to a '
                                             'resting cell; Pi and lactate over an order of magnitude; pH 6.8 to 8'},
                   'timestep': {'scheme': 'rk4', 'dts': [4e-5, 2e-5, 1e-5], 't_end': 0.01, 'min_order': 0.8,
                                'tol': 1e-3, 'cvode_tol': 1e-2},
                   # the near-equilibrium steps settle in milliseconds: ten seconds show any instability (a longer
                   # horizon made this the slowest test, 8 min)
                   'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-9}], 't_end': 10.0,
                                 'time_budget': 10}},
                  {'reviewed': False, 'notes': notes},
                  ['argus2026unpublished', 'bennett2009absolute'], VERSION)
    d = os.path.join(REPO, 'modules', 'cell', 'glycolysis', 'versions', VERSION)
    rdf = bg.bg_rdf(GLY_REST, kappa, CYTO['V_c'])
    with open(os.path.join(d, f'glycolysis_{VERSION}_bondgraph.ttl'), 'w') as f:
        f.write(rdf)
    with open(os.path.join(d, f'glycolysis_{VERSION}_bondgraph.svg'), 'w') as f:
        f.write(bg.svg(rdf))
    print(f'glycolysis: J_rest {J_REST:.5f} mM/s; fit residual {bg.fit_residual:.2e} J/mol; rest ' +
          ', '.join(f'{k} {GLY_REST[k]:.3g}' for k in ('F6P', 'DHAP', 'GAP', 'BPG13', 'PG3', 'PG2', 'PEP')) + ' mM')
    for r in bg.reactions:
        print(f'  {r.name:6s} dG0 {r.dG0 / 1e3:+6.1f}  dG at rest {gibbs[r.name] / 1e3:+7.2f} kJ/mol  kappa {kappa[r.name]:.4g} '
              f'mM/s' + (f'  sigma {r.sigma:.4g}' if r.sigma else ''))
    return bg, kappa



# ------------------------------------------------------------------ the mitochondrion (pyruvate to CO2, oxidative phosphorylation)
MITO = csv_values('modules/cell/energy_metabolism/versions/mitochondrial_ArgusUNPUBLISHED_v01/parameterisations/default/'
                  'default_parameters.csv')
T_M = MITO['T']
# the lumped respiration's resting O2 use (6 r_ox at its default point): this module is calibrated to it
_ADP = MITO['A_tot'] - MITO['ATP_init']
O2_REST = 6.0 * MITO['V_ox'] * MITO['C_O2'] / (MITO['K_O2'] + MITO['C_O2']) * MITO['C_glc'] / (MITO['K_glc'] + MITO['C_glc']) \
    * _ADP / (MITO['K_ADP'] + _ADP)
J_PYR = O2_REST / 2.5   # mM/s of pyruvate oxidised: 2.5 O2 per pyruvate
MEMORY = f'{PH}; textbook midpoint potentials (Lehninger, from memory, unchecked)'


def _fe(e_acceptor, e_donor, n):
    """dG0 of moving n electrons from a couple of midpoint potential e_donor to one of e_acceptor [V]."""
    return -n * FARADAY * (e_acceptor - e_donor)


MITO_SPECIES = [
    Species('Pyr_c', 'pyruvate (cytosol)', F(C=3, H=3, O=3), True, -1, compartment='c'),
    Species('Pi_c', 'phosphate (cytosol)', F(H=1, O=4, P=1), True, -2, compartment='c'),
    Species('H_c', 'H+ (cytosol and intermembrane space)', F(H=1), True, 1, c_std=1e-4, compartment='c'),
    Species('ATP_c', 'ATP (cytosol)', F(C=10, H=12, N=5, O=13, P=3), True, -4, compartment='c'),
    Species('ADP_c', 'ADP (cytosol)', F(C=10, H=12, N=5, O=10, P=2), True, -3, compartment='c'),
    Species('O2', 'O2', F(O=2), True, 0, compartment='m'),
    Species('CO2', 'CO2', F(C=1, O=2), True, 0, compartment='m'),
    Species('H2O', 'water', F(H=2, O=1), True, 0, c_std=55500.0, compartment='m'),
    Species('H_m', 'H+ (matrix)', F(H=1), True, 1, c_std=1e-4, compartment='m'),
    Species('Pyr_m', 'pyruvate', F(C=3, H=3, O=3), False, -1, compartment='m'),
    Species('AcCoA', 'acetyl-CoA', F(C=23, H=34, N=7, O=17, P=3, S=1), False, -4, compartment='m'),
    Species('CoA', 'coenzyme A', F(C=21, H=32, N=7, O=16, P=3, S=1), False, -4, compartment='m'),
    Species('Cit', 'citrate', F(C=6, H=5, O=7), False, -3, compartment='m'),
    Species('Icit', 'isocitrate', F(C=6, H=5, O=7), False, -3, compartment='m'),
    Species('OG', '2-oxoglutarate', F(C=5, H=4, O=5), False, -2, compartment='m'),
    Species('SucCoA', 'succinyl-CoA', F(C=25, H=35, N=7, O=19, P=3, S=1), False, -5, compartment='m'),
    Species('Suc', 'succinate', F(C=4, H=4, O=4), False, -2, compartment='m'),
    Species('Fum', 'fumarate', F(C=4, H=2, O=4), False, -2, compartment='m'),
    Species('Mal', 'malate', F(C=4, H=4, O=5), False, -2, compartment='m'),
    Species('OAA', 'oxaloacetate', F(C=4, H=2, O=5), False, -2, compartment='m'),
    Species('NAD_m', 'NAD+', F(C=21, H=26, N=7, O=14, P=2), False, -1, compartment='m'),
    Species('NADH_m', 'NADH', F(C=21, H=27, N=7, O=14, P=2), False, -2, compartment='m'),
    Species('ATP_m', 'ATP (matrix)', F(C=10, H=12, N=5, O=13, P=3), False, -4, compartment='m'),
    Species('ADP_m', 'ADP (matrix)', F(C=10, H=12, N=5, O=10, P=2), False, -3, compartment='m'),
    Species('Pi_m', 'phosphate (matrix)', F(H=1, O=4, P=1), False, -2, compartment='m'),
    Species('Q', 'ubiquinone', F(C=59, H=90, O=4), False, 0, compartment='m'),
    Species('QH2', 'ubiquinol', F(C=59, H=92, O=4), False, 0, compartment='m'),
    # cytochrome c's haem iron: Fe3+/Fe2+ (the protein, unchanged, not counted), in the intermembrane space
    Species('cytc_ox', 'cytochrome c (oxidised)', F(), False, 3, compartment='p'),
    Species('cytc_red', 'cytochrome c (reduced)', F(), False, 2, compartment='p'),
]
MITO_REACTIONS = [
    Reaction('PYRC', 'mitochondrial pyruvate carrier (H+ symport)', {'Pyr_c': -1, 'H_c': -1, 'Pyr_m': 1, 'H_m': 1}, 0.0,
             'definitional; transport'),
    Reaction('PIC', 'phosphate carrier (HPO4 2- with 2 H+, electroneutral)', {'Pi_c': -1, 'H_c': -2, 'Pi_m': 1, 'H_m': 2},
             0.0, 'definitional; transport'),
    Reaction('ANT', 'adenine nucleotide translocase (ATP4- out, ADP3- in)', {'ATP_m': -1, 'ADP_c': -1, 'ATP_c': 1,
                                                                             'ADP_m': 1}, 0.0, 'definitional; transport'),
    Reaction('PDH', 'pyruvate dehydrogenase', {'Pyr_m': -1, 'CoA': -1, 'NAD_m': -1, 'AcCoA': 1, 'NADH_m': 1, 'CO2': 1},
             -33.4e3, TEXTBOOK),
    Reaction('CS', 'citrate synthase', {'OAA': -1, 'AcCoA': -1, 'H2O': -1, 'Cit': 1, 'CoA': 1, 'H_m': 1}, -32.2e3, TEXTBOOK),
    Reaction('ACO', 'aconitase', {'Cit': -1, 'Icit': 1}, 13.3e3, TEXTBOOK),
    Reaction('IDH', 'isocitrate dehydrogenase', {'Icit': -1, 'NAD_m': -1, 'OG': 1, 'NADH_m': 1, 'CO2': 1}, -20.9e3, TEXTBOOK),
    Reaction('OGDH', '2-oxoglutarate dehydrogenase', {'OG': -1, 'CoA': -1, 'NAD_m': -1, 'SucCoA': 1, 'NADH_m': 1, 'CO2': 1},
             -33.5e3, TEXTBOOK),
    Reaction('SCS', 'succinyl-CoA synthetase (ATP-forming, as E. coli SucCD)',
             {'SucCoA': -1, 'ADP_m': -1, 'Pi_m': -1, 'Suc': 1, 'CoA': 1, 'ATP_m': 1}, -2.9e3, TEXTBOOK),
    Reaction('CII', 'complex II (succinate:ubiquinone oxidoreductase)', {'Suc': -1, 'Q': -1, 'Fum': 1, 'QH2': 1},
             _fe(0.045, 0.031, 2), MEMORY + ': fumarate/succinate +0.031 V, Q/QH2 +0.045 V'),
    Reaction('FUM', 'fumarase', {'Fum': -1, 'H2O': -1, 'Mal': 1}, -3.8e3, TEXTBOOK),
    Reaction('MDH', 'malate dehydrogenase', {'Mal': -1, 'NAD_m': -1, 'OAA': 1, 'NADH_m': 1, 'H_m': 1}, 29.7e3, TEXTBOOK),
    Reaction('CI', 'complex I (NADH:ubiquinone oxidoreductase, 4 H+ pumped)',
             {'NADH_m': -1, 'Q': -1, 'H_m': -5, 'NAD_m': 1, 'QH2': 1, 'H_c': 4}, _fe(0.045, -0.320, 2),
             MEMORY + ': NAD+/NADH -0.320 V, Q/QH2 +0.045 V'),
    Reaction('CIII', 'complex III (Q cycle: 4 H+ out, 2 taken up per QH2)',
             {'QH2': -1, 'cytc_ox': -2, 'H_m': -2, 'Q': 1, 'cytc_red': 2, 'H_c': 4}, _fe(0.254, 0.045, 2),
             MEMORY + ': cytochrome c +0.254 V'),
    Reaction('CIV', 'complex IV (cytochrome c oxidase, 4 H+ pumped per O2)',
             {'cytc_red': -4, 'O2': -1, 'H_m': -8, 'cytc_ox': 4, 'H2O': 2, 'H_c': 4}, _fe(0.816, 0.254, 4),
             MEMORY + ': O2/H2O +0.816 V'),
    Reaction('ATPS', 'ATP synthase (3 H+ per ATP)', {'ADP_m': -1, 'Pi_m': -1, 'H_c': -3, 'ATP_m': 1, 'H2O': 1, 'H_m': 2},
             30.5e3, TEXTBOOK + ': ATP hydrolysis -30.5 kJ/mol'),
    Reaction('LEAK', 'proton leak across the inner membrane', {'H_c': -1, 'H_m': 1}, 0.0, 'definitional; transport'),
]
PSI_REST = -0.150   # V, matrix against cytosol
PH_C, PH_M = 7.2, 7.7
C_M = 652.0         # F/m3
MITO_REST = dict(Pyr_c=0.1, Pi_c=2.0, H_c=1e3 * 10 ** -PH_C, ATP_c=MITO['ATP_init'], ADP_c=_ADP, O2=MITO['C_O2'],
                 CO2=MITO['C_CO2'], H2O=55500.0, H_m=1e3 * 10 ** -PH_M, AcCoA=0.07, CoA=0.1, Cit=0.4, OG=0.2, Mal=0.5,
                 NAD_m=2.5, NADH_m=1.0, ATP_m=1.0, ADP_m=1.0, Q=0.5, QH2=0.05, cytc_ox=0.2, cytc_red=0.05, SucCoA=0.05,
                 Suc=0.5)
MITO_NEAR_EQ = ('PYRC', 'PIC', 'ACO', 'FUM', 'MDH')
MITO_SOLVED = ['Pyr_m', 'Pi_m', 'Icit', 'Fum', 'OAA']
MITO_SATURATED = ('PDH', 'CS', 'IDH', 'OGDH', 'SCS', 'CII', 'CIII', 'CIV', 'ATPS', 'ANT')
MITO_HEADROOM = 4.0
# per pyruvate: 4 NADH and 1 QH2 (from succinate) pump 4 x 10 + 6 = 46 charges out; ATP synthase brings 3 back per ATP,
# the translocase 1 per ATP exported, the leak the rest: 4 s + 1 + L = 46, with L = 4 (about 9% of the pumped
# protons) s = 10.25 ATP from the synthase, 11.25 exported with succinyl-CoA's
N_ATPS, N_LEAK = 10.25, 4.0
MITO_FLUX = {n: J_PYR for n in ('PYRC', 'PDH', 'CS', 'ACO', 'IDH', 'OGDH', 'SCS', 'CII', 'FUM', 'MDH')}
MITO_FLUX.update(CI=4.0 * J_PYR, CIII=5.0 * J_PYR, CIV=2.5 * J_PYR, ATPS=N_ATPS * J_PYR, ANT=(N_ATPS + 1.0) * J_PYR,
                 PIC=(N_ATPS + 1.0) * J_PYR, LEAK=N_LEAK * J_PYR)
MITO_REST_REF = {
    'Pyr_c': f'{PH}; cytosolic pyruvate of order 0.1 mM', 'Pi_c': f'{PH}; cytosolic phosphate of order 1-5 mM',
    'H_c': f'{PH}; cytosol pH 7.2', 'ATP_c': 'energy_metabolism mitochondrial ATP_init (traut1994physiological)',
    'ADP_c': 'energy_metabolism mitochondrial A_tot - ATP_init', 'O2': 'energy_metabolism mitochondrial C_O2',
    'CO2': 'energy_metabolism mitochondrial C_CO2', 'H2O': 'definitional; pure water 55.5 M',
    'H_m': f'{PH}; matrix pH 7.7 (0.5 above the cytosol), held by its buffers',
    'AcCoA': f'{PH}; matrix acetyl-CoA of order 0.1 mM', 'SucCoA': f'{PH}; succinyl-CoA of order 0.05 mM',
    'Suc': f'{PH}; succinate of order 0.5 mM', 'CoA': f'{PH}; free CoA of order 0.1 mM',
    'Cit': f'{PH}; citrate of order 0.5 mM', 'OG': f'{PH}; 2-oxoglutarate of order 0.2 mM',
    'Mal': f'{PH}; malate of order 0.5 mM', 'NAD_m': f'{PH}; matrix NAD+ of order 2.5 mM',
    'NADH_m': f'{PH}; matrix NADH of order 1 mM (NADH/NAD+ 0.4)', 'ATP_m': f'{PH}; matrix ATP/ADP about 1',
    'ADP_m': f'{PH}; matrix ADP of order 1 mM', 'Q': f'{PH}; ubiquinone pool of order 0.5 mM',
    'QH2': f'{PH}; ubiquinol: the pool 10% reduced', 'cytc_ox': f'{PH}; cytochrome c of order 0.25 mM',
    'cytc_red': f'{PH}; cytochrome c 20% reduced',
    'psi_m': f'{PH}; mitochondrial membrane potential -150 mV',
    'C_m_m': f'{PH}; inner-membrane capacitance (Beard 2005 6.756e-6 mol/V per litre of mitochondria, from memory)',
}
for _n in MITO_SOLVED:
    MITO_REST_REF[_n] = 'calibrated: set by the near-equilibrium steps at rest (each forwards at 1000 times its net rate)'


def mitochondrion_rest(bg):
    return bg.solve_near_equilibrium({n: EXCHANGE for n in MITO_NEAR_EQ}, MITO_REST, MITO_SOLVED, EXCHANGE)


def mitochondrion_graph():
    return BondGraph('mitochondrion', 'Detailed mitochondrion', MITO_SPECIES, MITO_REACTIONS, temperature=T_M,
                     potentials={'m': PSI_REST}, capacitance={'m': C_M},
                     cofactors=('NAD_m', 'NADH_m', 'ATP_m', 'ADP_m', 'Q', 'QH2', 'cytc_ox', 'cytc_red', 'CoA'))


def build_mitochondrion():
    bg = mitochondrion_graph()
    rest = mitochondrion_rest(bg)
    bg.centre_potentials(rest)
    assert not bg.element_imbalance() and not bg.charge_imbalance(), (bg.element_imbalance(), bg.charge_imbalance())
    steady = bg.steady(MITO_FLUX)
    internal = [s.name for s in bg.species if not s.boundary]
    assert all(abs(steady[s]) < 1e-15 for s in internal), {s: steady[s] for s in internal if abs(steady[s]) >= 1e-15}
    assert abs(bg.current('m', MITO_FLUX)) < 1e-15, bg.current('m', MITO_FLUX)
    for r in bg.reactions:
        if r.name in MITO_SATURATED:
            pf, pr = bg.affinity_terms(r, rest)
            r.sigma = MITO_HEADROOM * (pf - pr) - (pf + pr)
            assert r.sigma > 0, r.name
    kappa = bg.calibrate(rest, MITO_FLUX)
    gibbs = {r.name: bg.reaction_gibbs(r, rest) for r in bg.reactions}
    variables, roles, params = module_parts(
        bg, rest, MITO_REST_REF, kappa, gibbs,
        lambda s: 'in the matrix, held by its buffers' if s.name == 'H_m' else 'in the cytosol (metabolite port)',
        'its share of the lumped respiration\'s resting O2 use (energy_metabolism mitochondrial)',
        lambda r: f'{PH}; saturates at {MITO_HEADROOM:g} times its resting flux (a respiratory reserve of a few times rest)',
        lambda r, s: PH, 'operating point: 37 C (energy_metabolism mitochondrial)')
    variables.append(('P_diss', 'W_per_m3', 'out', None, 'free energy dissipated, sum v A (>= 0, second law)'))
    roles['P_diss'] = 'variable'
    eqs = bg.equations()

    invariants = [
        {'description': 'second law: every reaction runs down its affinity, v A >= 0, so P_diss >= 0',
         'expr': ' and '.join(f'np.all(v_{r.name}*A_{r.name} >= -1e-12*np.abs(v_{r.name}).max()*np.abs(A_{r.name}).max())'
                              for r in bg.reactions) + ' and np.all(P_diss >= -1e-9*np.abs(P_diss).max())'},
        {'description': 'conserved moieties (relative 1e-9): NAD+ + NADH, CoA + acetyl-CoA + succinyl-CoA, matrix ATP + '
                        'ADP, Q + QH2, cytochrome c, and the TCA cycle\'s intermediates (citrate to oxaloacetate)',
         'expr': ' and '.join(f'np.allclose({e}, ({e})[0], rtol=1e-9, atol=0)' for e in (
             'C_NAD_m + C_NADH_m', 'C_CoA + C_AcCoA + C_SucCoA', 'C_ATP_m + C_ADP_m', 'C_Q + C_QH2',
             'C_cytc_ox + C_cytc_red',
             'C_Cit + C_Icit + C_OG + C_SucCoA + C_Suc + C_Fum + C_Mal + C_OAA'))},
        {'description': 'the ports report each chemostat\'s net production: j_O2 = -v_CIV, j_Pyr_c = -v_PYRC, j_ATP_c = '
                        'v_ANT, j_CO2 = v_PDH + v_IDH + v_OGDH (exact)',
         'expr': 'np.allclose(j_O2, -v_CIV, rtol=1e-12, atol=0) and np.allclose(j_Pyr_c, -v_PYRC, rtol=1e-12, atol=0) '
                 'and np.allclose(j_ATP_c, v_ANT, rtol=1e-12, atol=0) and '
                 'np.allclose(j_CO2, v_PDH + v_IDH + v_OGDH, rtol=1e-12, atol=1e-30)'},
    ]
    outputs = [f'C_{s}' for s in internal] + ['psi_m', 'I_m'] + [f'v_{r.name}' for r in bg.reactions] + \
        [f'j_{s.name}' for s in bg.species if s.boundary] + [f'A_{r.name}' for r in bg.reactions] + ['P_diss']
    port = [s.name for s in bg.species if s.boundary and s.name != 'H_m']
    ports = [{'port_type': 'cytosol_metabolites', 'variables': [f'C_{n}' for n in port] + [f'j_{n}' for n in port]}]
    notes = (
        'A mitochondrion step by step, from cytosolic pyruvate to CO2 with oxidative phosphorylation, as a bond graph: '
        'the pyruvate and phosphate carriers, pyruvate dehydrogenase, the eight steps of the TCA cycle, complexes I to IV '
        'and the ubiquinone and cytochrome c pools, ATP synthase, the adenine nucleotide translocase and a proton leak. '
        'Every reaction is balanced in C, H, O, N, P, S and in charge (species as at pH 7, H+ and water explicit); '
        'standard Gibbs energies are textbook values (transformed, pH 7, 1 M) and midpoint potentials, and the species\' '
        'standard potentials are fitted to them (tools/bondgraph.py), so the thermodynamics are consistent. The matrix '
        'has an electrical potential psi_m, a state on the inner membrane\'s capacitance (C_m dpsi/dt = F I_m): the '
        'complexes pump charge out (4, 2 and 8 per turnover), ATP synthase (3 H+ per ATP), the translocase (ATP4- out '
        'for ADP3- in) and the leak let it back, and each species in the matrix has the electrochemical potential mu0 + '
        'RT ln(C/c_std) + z F psi_m. So respiration is controlled as in a mitochondrion: with little ADP the synthase '
        'slows, psi_m rises and the pumps slow. Each reaction flows at v = kappa (Pf - Pr) / (1 + (Pf + Pr)/sigma): mass '
        'action (the bond-graph chemical reaction) for the near-equilibrium steps, complex I and the leak; saturating at '
        '4 times their resting flux for the dehydrogenases, citrate synthase, succinyl-CoA synthetase, complexes II, III '
        'and IV, ATP synthase and the translocase. The intermembrane space shares the cytosol\'s H+; the matrix pH is held (buffered). Calibrated to '
        'the lumped respiration\'s resting O2 use (energy_metabolism mitochondrial, 37 C) with typical concentrations '
        '(placeholders) where they let every step run forwards: the near-equilibrium steps run forwards at 1000 times '
        'their net rate, which sets pyruvate and phosphate in the matrix, isocitrate, fumarate and oxaloacetate. Per '
        'pyruvate '
        '11.25 ATP are exported (P/O 2.25; the leak takes 9% of the pumped protons). Pyruvate, phosphate, H+, ATP, ADP, '
        'O2, CO2 and water are chemostats on the cytosol_metabolites port. Compared against energy_metabolism '
        'mitochondrial in tests/test_metabolism_detailed.py. Its bond graph is in '
        f'mitochondrion_{VERSION}_bondgraph.ttl (BG-RDF; psi_m folded into the stores\' K) and drawn in '
        f'mitochondrion_{VERSION}_bondgraph.svg.')
    write_version('cell', 'mitochondrion', 'mitochondrion_detailed', variables, eqs,
                  {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
                   'general_ports': ports, 'roles': roles},
                  params,
                  {'sim_time': 600.0, 'pre_time': 0.0, 'dt': 0.1, 'solver': 'CVODE_myokit',
                   'solver_info': {'rtol': 1e-10, 'atol': 1e-12},
                   'outputs': outputs,
                   # off rest: a working cell, its ATP down and ADP up
                   'run_parameters': {'C_ATP_c': 2.5, 'C_ADP_c': 1.0},
                   'invariants': invariants,
                   'bc_sweep': {'ranges': {'C_O2': [1e-3, 0.2], 'C_ATP_c': [1.0, 5.0], 'C_ADP_c': [0.05, 2.0],
                                           'C_Pi_c': [0.5, 10.0], 'C_Pyr_c': [0.01, 1.0], 'C_CO2': [0.5, 5.0],
                                           'C_H_c': [2.5e-5, 1.6e-4], 'C_H_m': [1e-5, 4e-5]},
                                'exclude': ['R', 'T', 'F', 'C_H2O'] + unswept_potentials(bg),
                                'rationale': 'O2 from near anoxia to air-saturated; ATP and ADP from a de-energised to '
                                             'a resting cell; phosphate, pyruvate and CO2 over an order of magnitude; '
                                             'cytosol pH 6.8 to 7.6, matrix 7.4 to 8'},
                   # oxaloacetate (6e-5 mM, malate dehydrogenase forwards at 18 mM/s) turns over in microseconds
                   'timestep': {'scheme': 'rk4', 'dts': [4e-7, 2e-7, 1e-7], 't_end': 1e-4, 'min_order': 0.8,
                                'tol': 1e-3, 'cvode_tol': 1e-2},
                   # stiff (microsecond pools, a millisecond membrane potential): implicit solvers only, and a second
                   # is past every transient
                   'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-9}], 't_end': 1.0,
                                 'solve_ivp': ['BDF', 'LSODA'], 'time_budget': 10}},
                  {'reviewed': False, 'notes': notes},
                  ['argus2026unpublished'], VERSION)
    d = os.path.join(REPO, 'modules', 'cell', 'mitochondrion', 'versions', VERSION)
    rdf = bg.bg_rdf(rest, kappa, MITO['V_c'])
    with open(os.path.join(d, f'mitochondrion_{VERSION}_bondgraph.ttl'), 'w') as f:
        f.write(rdf)
    with open(os.path.join(d, f'mitochondrion_{VERSION}_bondgraph.svg'), 'w') as f:
        f.write(bg.svg(rdf))
    print(f'mitochondrion: J_pyr {J_PYR:.5f} mM/s (O2 {O2_REST:.4f}); fit residual {bg.fit_residual:.2e} J/mol; rest ' +
          ', '.join(f'{k} {rest[k]:.3g}' for k in MITO_SOLVED) + ' mM')
    for r in bg.reactions:
        print(f'  {r.name:6s} dG0 {r.dG0 / 1e3:+7.1f}  dG at rest {gibbs[r.name] / 1e3:+8.2f} kJ/mol  kappa '
              f'{kappa[r.name]:.4g} mM/s' + (f'  sigma {r.sigma:.4g}' if r.sigma else ''))
    return bg, kappa



# ------------------------------------------------------------------ glycogen (the fermenter's store, ADP-glucose route)
GLYC_REST_C = 50.0   # mM glucose units: the comparison point, a quarter of the lumped store's ceiling
_RT = R_GAS * T
_K_EQ_STORE = math.exp(-(CYTO['dG0_store'] - CYTO['dG0_hk']) / _RT)
_K_EQ_MOB = math.exp(-(CYTO['dG0_hk'] + CYTO['dG0_hyd'] - CYTO['dG0_store']) / _RT)


def lumped_store_rates(g6p, atp, adp, pi, glyc):
    """energy_metabolism adaptive's v_store and v_mob [mM/s] with the store fully allocated (store_alloc 1)."""
    a = glyc / (CYTO['K_glyc_mob'] + glyc)
    store = CYTO['V_store_max'] * (g6p * atp - a * adp * pi ** 2 / (CYTO['c_std'] * _K_EQ_STORE)) \
        / ((CYTO['K_G6P_store'] + g6p) * (CYTO['K_ATP_store'] + atp)) / (1.0 + (glyc / CYTO['C_glyc_max']) ** 8)
    mob = CYTO['V_mob'] * a * (pi - g6p / _K_EQ_MOB) / (CYTO['K_Pi_mob'] + pi) * CYTO['K_i_ATP_mob'] \
        / (CYTO['K_i_ATP_mob'] + atp)
    return store, mob


V_STORE, V_MOB = lumped_store_rates(CYTO['C_G6P_init'], CYTO['C_ATP_init'], CYTO['C_ADP_init'], CYTO['C_Pi_init'],
                                    GLYC_REST_C)
STORE_HEADROOM = CYTO['V_store_max'] / V_STORE
MOB_HEADROOM = CYTO['V_mob'] * GLYC_REST_C / (CYTO['K_glyc_mob'] + GLYC_REST_C) * CYTO['K_i_ATP_mob'] \
    / (CYTO['K_i_ATP_mob'] + CYTO['C_ATP_init']) / V_MOB
GLYC_SPECIES = [
    Species('G6P', 'glucose 6-phosphate', F(C=6, H=11, O=9, P=1), True, -2),
    Species('ATP', 'ATP', F(C=10, H=12, N=5, O=13, P=3), True, -4),
    Species('ADP', 'ADP', F(C=10, H=12, N=5, O=10, P=2), True, -3),
    Species('Pi', 'phosphate', F(H=1, O=4, P=1), True, -2),
    Species('H', 'H+', F(H=1), True, 1, c_std=1e-4),
    Species('H2O', 'water', F(H=2, O=1), True, 0, c_std=55500.0),
    Species('G1P', 'glucose 1-phosphate', F(C=6, H=11, O=9, P=1), False, -2),
    Species('ADPglc', 'ADP-glucose', F(C=16, H=23, N=5, O=15, P=2), False, -2),
    Species('PPi', 'pyrophosphate', F(H=1, O=7, P=2), False, -3),
    # glucose units in glycogen: activity C/(K + C) (its chain ends), standard state saturated (activity 1)
    Species('Gly', 'glycogen (glucose units)', F(C=6, H=10, O=5), False, 0, c_std=CYTO['K_glyc_mob'],
            k_sat=CYTO['K_glyc_mob']),
]
GLYC_REACTIONS = [
    Reaction('PGM', 'phosphoglucomutase', {'G6P': -1, 'G1P': 1}, 7.3e3, TEXTBOOK),
    Reaction('GLGC', 'ADP-glucose pyrophosphorylase (GlgC)', {'G1P': -1, 'ATP': -1, 'H': -1, 'ADPglc': 1, 'PPi': 1}, 0.0,
             f'{PH}; nucleotide-sugar pyrophosphorylase, near 0 (UDP-glucose pyrophosphorylase, from memory)'),
    Reaction('GLGA', 'glycogen synthase (GlgA)', {'ADPglc': -1, 'ADP': 1, 'Gly': 1, 'H': 1}, -13.4e3,
             TEXTBOOK + ' (glycogen synthase, UDP-glucose)'),
    Reaction('PPA', 'inorganic pyrophosphatase', {'PPi': -1, 'H2O': -1, 'Pi': 2, 'H': 1}, -19.2e3, TEXTBOOK),
    Reaction('GLGP', 'glycogen phosphorylase (GlgP)', {'Gly': -1, 'Pi': -1, 'G1P': 1}, 3.1e3, TEXTBOOK,
             inhibitors={'ATP': CYTO['K_i_ATP_mob']}),
]
GLYC_FLUX = dict(PGM=V_STORE - V_MOB, GLGC=V_STORE, GLGA=V_STORE, PPA=V_STORE, GLGP=V_MOB)
GLYC_REST = dict(G6P=CYTO['C_G6P_init'], ATP=CYTO['C_ATP_init'], ADP=CYTO['C_ADP_init'], Pi=CYTO['C_Pi_init'],
                 H=1e3 * 10 ** -PH_I, H2O=55500.0, ADPglc=0.05, PPi=0.5, Gly=GLYC_REST_C)
GLYC_REST_REF = {
    'G6P': 'energy_metabolism adaptive C_G6P_init', 'ATP': 'energy_metabolism adaptive C_ATP_init',
    'ADP': 'energy_metabolism adaptive C_ADP_init', 'Pi': 'energy_metabolism adaptive C_Pi_init',
    'H': 'the fermenter\'s cytosol pH 7.6 (generic_adaptive_cell pH_i)', 'H2O': 'definitional; pure water 55.5 M',
    'ADPglc': f'{PH}; ADP-glucose of order 0.05 mM', 'PPi': f'{PH}; pyrophosphate of order 0.5 mM',
    'Gly': f'{PH}; the comparison point: glycogen at a quarter of the lumped store\'s ceiling',
    'G1P': 'calibrated: phosphoglucomutase near equilibrium at rest (forwards at 1000 times the synthesis flux)',
    'k_sat_Gly': 'energy_metabolism adaptive K_glyc_mob (the phosphorylase\'s glycogen affinity)',
}


def glycogen_graph():
    return BondGraph('glycogen', 'Detailed glycogen store', GLYC_SPECIES, GLYC_REACTIONS, temperature=T)


def glycogen_rest(bg):
    pgm = {r.name: r for r in bg.reactions}['PGM']
    ratio = EXCHANGE * V_STORE / GLYC_FLUX['PGM']
    return dict(GLYC_REST, G1P=bg.near_equilibrium(pgm, GLYC_REST, 'G1P', ratio))


def build_glycogen():
    bg = glycogen_graph()
    rest = glycogen_rest(bg)
    bg.centre_potentials(rest)
    assert not bg.element_imbalance() and not bg.charge_imbalance(), (bg.element_imbalance(), bg.charge_imbalance())
    steady = bg.steady(GLYC_FLUX)
    assert all(abs(steady[s]) < 1e-15 for s in ('G1P', 'ADPglc', 'PPi')), steady
    for r in bg.reactions:
        head = {'GLGC': STORE_HEADROOM, 'GLGA': STORE_HEADROOM, 'GLGP': MOB_HEADROOM}.get(r.name)
        if head:
            pf, pr = bg.affinity_terms(r, rest)
            r.sigma = head * (pf - pr) - (pf + pr)
            assert r.sigma > 0, r.name
    kappa = bg.calibrate(rest, GLYC_FLUX)
    gibbs = {r.name: bg.reaction_gibbs(r, rest) for r in bg.reactions}
    variables, roles, params = module_parts(
        bg, rest, GLYC_REST_REF, kappa, gibbs, lambda s: 'in the cytosol (metabolite port)',
        'the lumped store\'s synthesis and phosphorylase rates at the comparison point (energy_metabolism adaptive, '
        'store_alloc 1, glycogen 50 mM)',
        lambda r: 'calibrated: saturates at the lumped ' + ('V_store_max' if r.name != 'GLGP' else 'V_mob') +
                  ' (energy_metabolism adaptive)',
        lambda r, s: 'energy_metabolism adaptive K_i_ATP_mob (unreferenced_placeholder)',
        'operating point: 20 C (the fermenter\'s)')
    variables.append(('P_diss', 'W_per_m3', 'out', None, 'free energy dissipated, sum v A (>= 0, second law)'))
    roles['P_diss'] = 'variable'
    invariants = [
        {'description': 'second law: every reaction runs down its affinity, v A >= 0, so P_diss >= 0',
         'expr': ' and '.join(f'np.all(v_{r.name}*A_{r.name} >= -1e-12*np.abs(v_{r.name}).max()*np.abs(A_{r.name}).max())'
                              for r in bg.reactions) + ' and np.all(P_diss >= -1e-9*np.abs(P_diss).max())'},
        {'description': 'glycogen is never negative (its activity vanishes as it runs out)',
         'expr': 'np.all(C_Gly >= -1e-9*np.abs(C_Gly).max())'},
        {'description': 'the ports report each chemostat\'s net production: j_G6P = -v_PGM, j_ATP = -v_GLGC, j_ADP = '
                        'v_GLGA, j_Pi = 2 v_PPA - v_GLGP (exact)',
         'expr': 'np.allclose(j_G6P, -v_PGM, rtol=1e-12, atol=0) and np.allclose(j_ATP, -v_GLGC, rtol=1e-12, atol=0) and '
                 'np.allclose(j_ADP, v_GLGA, rtol=1e-12, atol=0) and '
                 'np.allclose(j_Pi, 2.0*v_PPA - v_GLGP, rtol=1e-12, atol=1e-30)'},
    ]
    internal = [s.name for s in bg.species if not s.boundary]
    outputs = [f'C_{s}' for s in internal] + [f'v_{r.name}' for r in bg.reactions] + \
        [f'j_{s.name}' for s in bg.species if s.boundary] + [f'A_{r.name}' for r in bg.reactions] + ['P_diss']
    port = [s.name for s in bg.species if s.boundary]
    notes = (
        'The glycogen store step by step, by the bacterial ADP-glucose route (E. coli GlgC, GlgA, GlgP): '
        'phosphoglucomutase (G6P to G1P), ADP-glucose pyrophosphorylase (G1P + ATP to ADP-glucose + PPi), glycogen '
        'synthase (a glucose unit onto glycogen, ADP out), pyrophosphatase (pulls the pyrophosphorylase forwards) and '
        'glycogen phosphorylase (glycogen + Pi to G1P; inhibited by ATP, as the lumped store). The designer\'s '
        'CarbohydrateMetabolism.cellml has the animal route (UTP, UDP-glucose): the same bond graph with UTP and UDP in '
        'place of ATP and ADP. Every reaction is balanced in elements and charge; the textbook standard Gibbs energies '
        'sum to the lumped store\'s (G6P + ATP + H2O to a glucose unit + ADP + 2 Pi + H+: -25.3 kJ/mol). Glycogen is a '
        'store whose activity saturates, C/(K + C) with K the lumped phosphorylase\'s glycogen affinity (its chain ends; '
        'the lumped a_glyc), so the phosphorylase stops as it runs out. Calibrated at the comparison point (the '
        'fermenter\'s rest with the store fully allocated and glycogen at 50 mM) to the lumped store\'s synthesis and '
        'phosphorylase rates, saturating at the lumped maxima. Unlike the lumped store it has no ceiling (glycogen\'s '
        'size is limited by branching and particle structure, not modelled). Compared against energy_metabolism adaptive '
        'in tests/test_metabolism_detailed.py; its bond graph is in '
        f'glycogen_{VERSION}_bondgraph.ttl (BG-RDF) and drawn in glycogen_{VERSION}_bondgraph.svg.')
    write_version('cell', 'glycogen', 'glycogen_detailed', variables, bg.equations(),
                  {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
                   'general_ports': [{'port_type': 'cytosol_metabolites',
                                      'variables': [f'C_{n}' for n in port] + [f'j_{n}' for n in port]}],
                   'roles': roles},
                  params,
                  {'sim_time': 600.0, 'pre_time': 0.0, 'dt': 0.1, 'solver': 'CVODE_myokit',
                   'solver_info': {'rtol': 1e-10, 'atol': 1e-12}, 'outputs': outputs,
                   # off rest: ATP low, so the store is spent
                   'run_parameters': {'C_ATP': 2.0, 'C_ADP': 3.0},
                   'invariants': invariants,
                   'bc_sweep': {'ranges': {'C_G6P': [0.05, 10.0], 'C_ATP': [1.0, 12.0], 'C_ADP': [0.1, 3.0],
                                           'C_Pi': [1.0, 20.0], 'C_H': [1e-5, 1.6e-4]},
                                'exclude': ['R', 'T', 'C_H2O'] + unswept_potentials(bg),
                                'rationale': 'G6P from starving to a meal; ATP and ADP from a de-energised to a resting '
                                             'cell; Pi over an order of magnitude; pH 6.8 to 8'},
                   'timestep': {'scheme': 'rk4', 'dts': [4e-5, 2e-5, 1e-5], 't_end': 0.01, 'min_order': 0.8,
                                'tol': 1e-3, 'cvode_tol': 1e-2},
                   'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-9}], 't_end': 10.0,
                                 'time_budget': 10}},
                  {'reviewed': False, 'notes': notes}, ['argus2026unpublished'], VERSION)
    _write_bond_graph(bg, 'glycogen', rest, kappa, CYTO['V_c'])
    _report('glycogen', bg, kappa, gibbs, f'store {V_STORE:.5f}, phosphorylase {V_MOB:.5f} mM/s; G1P {rest["G1P"]:.3g} mM')
    return bg, kappa


# ------------------------------------------------------------------ gluconeogenesis (GLY17-20 of the designer's model)
GNG_SPECIES = [
    Species('F16BP', 'fructose 1,6-bisphosphate', F(C=6, H=10, O=12, P=2), True, -4),
    Species('F6P', 'fructose 6-phosphate', F(C=6, H=11, O=9, P=1), True, -2),
    Species('G6P', 'glucose 6-phosphate', F(C=6, H=11, O=9, P=1), True, -2),
    Species('glc', 'glucose', F(C=6, H=12, O=6), True, 0),
    Species('DHAP', 'dihydroxyacetone phosphate', F(C=3, H=5, O=6, P=1), True, -2),
    Species('glycerol', 'glycerol', F(C=3, H=8, O=3), True, 0),
    Species('NAD', 'NAD+', F(C=21, H=26, N=7, O=14, P=2), True, -1),
    Species('NADH', 'NADH', F(C=21, H=27, N=7, O=14, P=2), True, -2),
    Species('ATP', 'ATP', F(C=10, H=12, N=5, O=13, P=3), True, -4),
    Species('ADP', 'ADP', F(C=10, H=12, N=5, O=10, P=2), True, -3),
    Species('Pi', 'phosphate', F(H=1, O=4, P=1), True, -2),
    Species('H', 'H+', F(H=1), True, 1, c_std=1e-4),
    Species('H2O', 'water', F(H=2, O=1), True, 0, c_std=55500.0),
    Species('Gly3P', 'glycerol 3-phosphate', F(C=3, H=7, O=6, P=1), False, -2),
]
GNG_REACTIONS = [
    Reaction('FBPASE', 'fructose 1,6-bisphosphatase', {'F16BP': -1, 'H2O': -1, 'F6P': 1, 'Pi': 1}, -16.3e3, TEXTBOOK),
    Reaction('G6PASE', 'glucose 6-phosphatase', {'G6P': -1, 'H2O': -1, 'glc': 1, 'Pi': 1}, -13.8e3, TEXTBOOK),
    Reaction('GK', 'glycerol kinase', {'glycerol': -1, 'ATP': -1, 'Gly3P': 1, 'ADP': 1, 'H': 1}, -21.3e3,
             f'{PH}; ATP hydrolysis -30.5 less glycerol 3-phosphate hydrolysis -9.2 kJ/mol (textbook, from memory)'),
    Reaction('G3PDH', 'glycerol 3-phosphate dehydrogenase (towards DHAP)',
             {'Gly3P': -1, 'NAD': -1, 'DHAP': 1, 'NADH': 1, 'H': 1}, _fe(-0.320, -0.190, 2),
             MEMORY + ': DHAP/glycerol 3-phosphate -0.190 V, NAD+/NADH -0.320 V'),
]
J_GNG = 0.01   # mM/s of glucose made
GNG_FLUX = dict(FBPASE=J_GNG, G6PASE=J_GNG, GK=0.4 * J_GNG, G3PDH=0.4 * J_GNG)
GNG_HEADROOM = 4.0
GNG_REST = dict(F16BP=0.02, F6P=0.1, G6P=0.2, glc=8.0, DHAP=0.03, glycerol=0.1, NAD=0.5, NADH=0.0005,
                ATP=MITO['ATP_init'], ADP=_ADP, Pi=5.0, H=1e3 * 10 ** -PH_C, H2O=55500.0)
GNG_REST_REF = {
    'F16BP': f'{PH}; liver fructose 1,6-bisphosphate of order 0.02 mM', 'F6P': f'{PH}; fructose 6-phosphate of order 0.1 mM',
    'G6P': f'{PH}; liver glucose 6-phosphate of order 0.2 mM', 'glc': f'{PH}; liver glucose near plasma, 8 mM',
    'DHAP': f'{PH}; DHAP of order 0.03 mM', 'glycerol': f'{PH}; glycerol of order 0.1 mM',
    'NAD': f'{PH}; free cytosolic NAD+ of order 0.5 mM',
    'NADH': f'{PH}; free cytosolic NADH/NAD+ about 1/1000 (a respiring cell)',
    'ATP': 'energy_metabolism mitochondrial ATP_init (traut1994physiological)',
    'ADP': 'energy_metabolism mitochondrial A_tot - ATP_init', 'Pi': f'{PH}; cytosolic phosphate of order 5 mM',
    'H': f'{PH}; cytosol pH 7.2', 'H2O': 'definitional; pure water 55.5 M',
    'Gly3P': 'calibrated: glycerol 3-phosphate dehydrogenase near equilibrium at rest (forwards at 1000 times its net rate)',
}


def gluconeogenesis_graph():
    return BondGraph('gluconeogenesis', 'Detailed gluconeogenic bypasses and glycerol entry', GNG_SPECIES, GNG_REACTIONS,
                     temperature=T_M)


def gluconeogenesis_rest(bg):
    g3pdh = {r.name: r for r in bg.reactions}['G3PDH']
    return dict(GNG_REST, Gly3P=bg.near_equilibrium(g3pdh, GNG_REST, 'Gly3P', EXCHANGE))


def build_gluconeogenesis():
    bg = gluconeogenesis_graph()
    rest = gluconeogenesis_rest(bg)
    bg.centre_potentials(rest)
    assert not bg.element_imbalance() and not bg.charge_imbalance(), (bg.element_imbalance(), bg.charge_imbalance())
    assert abs(bg.steady(GNG_FLUX)['Gly3P']) < 1e-15
    for r in bg.reactions:
        if r.name != 'G3PDH':
            pf, pr = bg.affinity_terms(r, rest)
            r.sigma = GNG_HEADROOM * (pf - pr) - (pf + pr)
            assert r.sigma > 0, r.name
    kappa = bg.calibrate(rest, GNG_FLUX)
    gibbs = {r.name: bg.reaction_gibbs(r, rest) for r in bg.reactions}
    variables, roles, params = module_parts(
        bg, rest, GNG_REST_REF, kappa, gibbs, lambda s: 'in the cytosol (metabolite port)',
        f'a placeholder gluconeogenic flux ({J_GNG:g} mM/s of glucose, glycerol 20% of its carbon)',
        lambda r: f'{PH}; saturates at {GNG_HEADROOM:g} times its resting flux', lambda r, s: PH,
        'operating point: 37 C (a liver-like respiring cell)')
    variables.append(('P_diss', 'W_per_m3', 'out', None, 'free energy dissipated, sum v A (>= 0, second law)'))
    roles['P_diss'] = 'variable'
    invariants = [
        {'description': 'second law: every reaction runs down its affinity, v A >= 0, so P_diss >= 0',
         'expr': ' and '.join(f'np.all(v_{r.name}*A_{r.name} >= -1e-12*np.abs(v_{r.name}).max()*np.abs(A_{r.name}).max())'
                              for r in bg.reactions) + ' and np.all(P_diss >= -1e-9*np.abs(P_diss).max())'},
        {'description': 'the ports report each chemostat\'s net production: j_glc = v_G6PASE, j_F6P = v_FBPASE, '
                        'j_glycerol = -v_GK, j_DHAP = v_G3PDH (exact)',
         'expr': 'np.allclose(j_glc, v_G6PASE, rtol=1e-12, atol=0) and np.allclose(j_F6P, v_FBPASE, rtol=1e-12, atol=0) '
                 'and np.allclose(j_glycerol, -v_GK, rtol=1e-12, atol=0) and np.allclose(j_DHAP, v_G3PDH, rtol=1e-12, '
                 'atol=0)'},
    ]
    outputs = ['C_Gly3P'] + [f'v_{r.name}' for r in bg.reactions] + \
        [f'j_{s.name}' for s in bg.species if s.boundary] + [f'A_{r.name}' for r in bg.reactions] + ['P_diss']
    port = [s.name for s in bg.species if s.boundary]
    notes = (
        'The steps of gluconeogenesis that glycolysis cannot run backwards, and glycerol\'s way in, from the designer\'s '
        'CarbohydrateMetabolism.cellml (GLY17-20): fructose 1,6-bisphosphatase and glucose 6-phosphatase (the bypasses '
        'of phosphofructokinase and hexokinase: with them each forms a cycle that hydrolyses one ATP, their standard '
        'Gibbs energies summing to ATP hydrolysis\'), glycerol kinase and glycerol 3-phosphate dehydrogenase (glycerol to '
        'DHAP, which the glycolysis module\'s aldolase and isomerase take on). Every reaction is balanced in elements and '
        'charge, with textbook standard Gibbs energies and midpoint potentials. The metabolites it shares with '
        'glycolysis are chemostats on the cytosol_metabolites port; glycerol 3-phosphate is its own. Calibrated, with '
        'placeholder concentrations, to a respiring (liver-like) cytosol at 37 C making glucose at a placeholder rate: '
        'the dehydrogenase can carry glycerol to DHAP only where NADH/NAD+ is low (about 1/1000 here); at the '
        'fermenter\'s NADH/NAD+ (1/30) it runs the other way (tests/test_metabolism_detailed.py). Its bond graph is in '
        f'gluconeogenesis_{VERSION}_bondgraph.ttl (BG-RDF) and drawn in gluconeogenesis_{VERSION}_bondgraph.svg.')
    write_version('cell', 'gluconeogenesis', 'gluconeogenesis_detailed', variables, bg.equations(),
                  {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
                   'general_ports': [{'port_type': 'cytosol_metabolites',
                                      'variables': [f'C_{n}' for n in port] + [f'j_{n}' for n in port]}],
                   'roles': roles},
                  params,
                  {'sim_time': 600.0, 'pre_time': 0.0, 'dt': 0.1, 'solver': 'CVODE_myokit',
                   'solver_info': {'rtol': 1e-10, 'atol': 1e-12}, 'outputs': outputs,
                   # off rest: ten times the glycerol
                   'run_parameters': {'C_glycerol': 1.0},
                   'invariants': invariants,
                   'bc_sweep': {'ranges': {'C_glycerol': [0.01, 5.0], 'C_NADH': [1e-4, 0.05], 'C_F16BP': [0.002, 0.2],
                                           'C_G6P': [0.02, 2.0], 'C_glc': [1.0, 20.0], 'C_ATP': [1.0, 5.0],
                                           'C_Pi': [0.5, 10.0]},
                                'exclude': ['R', 'T', 'C_H2O'] + unswept_potentials(bg),
                                'rationale': 'glycerol from trace to a fat-mobilising load; NADH/NAD+ from respiring to '
                                             'fermenting; the sugar phosphates and glucose over an order of magnitude'},
                   'timestep': {'scheme': 'rk4', 'dts': [4e-5, 2e-5, 1e-5], 't_end': 0.01, 'min_order': 0.8,
                                'tol': 1e-3, 'cvode_tol': 1e-2},
                   'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-9}], 't_end': 10.0,
                                 'time_budget': 10}},
                  {'reviewed': False, 'notes': notes}, ['argus2026unpublished'], VERSION)
    _write_bond_graph(bg, 'gluconeogenesis', rest, kappa, MITO['V_c'])
    _report('gluconeogenesis', bg, kappa, gibbs, f'Gly3P {rest["Gly3P"]:.3g} mM')
    return bg, kappa


def _write_bond_graph(bg, module, rest, kappa, volume):
    d = os.path.join(REPO, 'modules', 'cell', module, 'versions', VERSION)
    rdf = bg.bg_rdf(rest, kappa, volume)
    with open(os.path.join(d, f'{module}_{VERSION}_bondgraph.ttl'), 'w') as f:
        f.write(rdf)
    with open(os.path.join(d, f'{module}_{VERSION}_bondgraph.svg'), 'w') as f:
        f.write(bg.svg(rdf))


def _report(module, bg, kappa, gibbs, line):
    print(f'{module}: {line}; fit residual {bg.fit_residual:.2e} J/mol')
    for r in bg.reactions:
        print(f'  {r.name:6s} dG0 {r.dG0 / 1e3:+7.1f}  dG at rest {gibbs[r.name] / 1e3:+8.2f} kJ/mol  kappa '
              f'{kappa[r.name]:.4g} mM/s' + (f'  sigma {r.sigma:.4g}' if r.sigma else ''))


BUILDERS = {'glycolysis': (lambda: build_glycolysis()), 'mitochondrion': (lambda: build_mitochondrion()),
            'glycogen': (lambda: build_glycogen()), 'gluconeogenesis': (lambda: build_gluconeogenesis())}


if __name__ == '__main__':
    for build in BUILDERS.values():
        build()
    print('written')
