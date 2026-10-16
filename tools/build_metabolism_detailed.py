"""Builds detailed metabolism modules: each step of a pathway as its own reaction, to compare against the lumped
modules of the generic adaptive cell (tools/build_adaptive_cell.py), which make simplifying assumptions (energy_metabolism
adaptive: glycolysis as two steps). The pathways follow the carbohydrate-metabolism model the designer supplied
(CarbohydrateMetabolism.cellml: glycolysis, glycogen, gluconeogenesis, the mitochondrial carriers, the TCA cycle and the
electron transport chain, element-balanced), here split into modules with physical units, thermodynamically consistent
kinetics and a bond-graph description each (tools/bondgraph.py: the module's equations, its BG-RDF and its diagram).

    python tools/build_metabolism_detailed.py

Module versions:

  cell/glycolysis/detailed_ArgusUNPUBLISHED_v01   glucose to lactate in eleven steps (hexokinase to lactate
                                                  dehydrogenase), NAD+/NADH inside; glucose, lactate, ATP, ADP, Pi, H+
                                                  and water from a cytosol (chemostats), calibrated to the fermenter's
                                                  resting flux and concentrations

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
from bondgraph import BondGraph, Reaction, Species, R_GAS  # noqa: E402

VERSION = 'detailed_ArgusUNPUBLISHED_v01'
PH = 'unreferenced_placeholder'
TEXTBOOK = f'{PH}; textbook transformed standard Gibbs energy at pH 7 (Lehninger, from memory, unchecked)'
module_writer.UNITS.update({'W_per_m3': [('watt', None, 1, 1), ('metre', None, -3, 1)]})

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
CYTO = csv_values('modules/cell/energy_metabolism/versions/adaptive_ArgusUNPUBLISHED_v01/instances/default/default_parameters.csv')
BIO = csv_values('modules/cell/biosynthesis/versions/adaptive_ArgusUNPUBLISHED_v01/instances/default/default_parameters.csv')
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


def build_glycolysis():
    global GLY_REST
    bg = BondGraph('glycolysis', 'Detailed glycolysis', GLY_SPECIES, GLY_REACTIONS, temperature=T)
    GLY_REST = near_equilibrium_rest(bg)
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

    variables = [('t', 'second', 'in', None, None), ('R', 'J_per_mol_K', 'in', None, None),
                 ('T', 'kelvin', 'in', None, None)]
    roles = {'R': 'global_constant', 'T': 'global_constant'}
    params = [('R', R_GAS, 'definitional; gas constant', True), ('T', T, 'operating point: 20 C (the fermenter\'s)', False)]
    for s in bg.species:
        variables.append((f'mu0_{s.name}', 'J_per_mol', 'in', None, f'standard chemical potential of {s.label} (from the '
                          'reactions\' Gibbs energies, tools/bondgraph.py)'))
        variables.append((f'c_std_{s.name}', 'millimolar', 'in', None, None))
        roles[f'mu0_{s.name}'] = roles[f'c_std_{s.name}'] = 'constant'
        params += [(f'mu0_{s.name}', bg.mu0[s.name], 'calibrated: min-norm fit to the reactions\' standard Gibbs energies',
                    False),
                   (f'c_std_{s.name}', s.c_std, 'definitional; standard state (1 M; H+ at pH 7; water itself)', True)]
        if s.boundary:
            variables.append((f'C_{s.name}', 'millimolar', 'in', None, f'{s.label} in the cytosol (metabolite port)'))
            variables.append((f'j_{s.name}', 'millimolar_per_s', 'out', None, f'net production of {s.label} here'))
            roles[f'C_{s.name}'] = 'boundary_condition'
            roles[f'j_{s.name}'] = 'variable'
            params.append((f'C_{s.name}', GLY_REST[s.name], GLY_REST_REF[s.name], False))
        else:
            variables.append((f'C_{s.name}_init', 'millimolar', 'in', None, None))
            variables.append((f'C_{s.name}', 'millimolar', 'out', f'C_{s.name}_init', s.label))
            roles[f'C_{s.name}_init'] = 'constant'
            roles[f'C_{s.name}'] = 'variable'
            params.append((f'C_{s.name}_init', GLY_REST[s.name], GLY_REST_REF[s.name], False))
    for r in bg.reactions:
        variables += [(f'kappa_{r.name}', 'millimolar_per_s', 'in', None, f'{r.enzyme}: rate constant'),
                      (f'Pf_{r.name}', 'dimensionless', 'out', None, 'exp(forward affinity / RT)'),
                      (f'Pr_{r.name}', 'dimensionless', 'out', None, 'exp(reverse affinity / RT)'),
                      (f'v_{r.name}', 'millimolar_per_s', 'out', None, f'{r.enzyme}'),
                      (f'A_{r.name}', 'J_per_mol', 'out', None, 'affinity: minus the reaction\'s Gibbs energy')]
        roles.update({f'kappa_{r.name}': 'constant', f'Pf_{r.name}': 'variable', f'Pr_{r.name}': 'variable',
                      f'v_{r.name}': 'variable', f'A_{r.name}': 'variable'})
        params.append((f'kappa_{r.name}', kappa[r.name], f'calibrated: carries the fermenter\'s resting glycolytic flux '
                       f'(Gibbs energy at rest {gibbs[r.name] / 1e3:+.2f} kJ/mol; dG0 {r.dG0 / 1e3:+.1f} kJ/mol, {r.ref})',
                       False))
        for s, k in r.activators.items():
            variables.append((f'Ka_{s}_{r.name}', 'millimolar', 'in', None, f'activation of {r.enzyme} by {s}'))
            roles[f'Ka_{s}_{r.name}'] = 'constant'
            params.append((f'Ka_{s}_{r.name}', k, f'{PH}; feed-forward activation constant (E. coli PykF by F1,6BP, of order '
                           '0.1-1 mM, from memory)', False))
        if r.sigma:
            variables.append((f'sigma_{r.name}', 'dimensionless', 'in', None, 'saturation (half-saturated at Pf = sigma)'))
            roles[f'sigma_{r.name}'] = 'constant'
            params.append((f'sigma_{r.name}', r.sigma, f'{PH}; saturates at {HEADROOM:g} times its resting flux (the lumped '
                           'glycolysis\' headroom)', False))
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
                                'exclude': ['R', 'T', 'C_H2O'],
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


if __name__ == '__main__':
    build_glycolysis()
    print('written')
