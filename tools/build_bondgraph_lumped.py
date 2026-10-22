"""Bond-graph descriptions of existing (lumped) modules: their stoichiometry and standard Gibbs energies as a bond graph
(tools/bondgraph.py), whose CellDL document tools/build_celldl.py writes next to the version. The module's own
kinetics are not mass action, so the bond graph's rate constants are the module's rates at an export point (each
reaction's flux over its Pf - Pr there): the same flows there, mass action away from it. The module's equations are
unchanged.

    python tools/build_bondgraph_lumped.py      (prints the export point's rates)

  cell/energy_metabolism/adaptive_ArgusUNPUBLISHED_v01   the fermenter's cytosol: hexokinase, the lower glycolysis,
                                                         ATP hydrolysis, the glycogen store and phosphorylase, matrix
                                                         export; exported with the store and matrix allocated and
                                                         glycogen at 50 mM, so every reaction runs

Checked in tests/test_metabolism_detailed.py: the diagram's stoichiometry gives the module's own rates of change.
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bondgraph import BondGraph, Reaction, Species  # noqa: E402
from build_metabolism_detailed import CYTO, GLYC_REST_C, F  # noqa: E402

CYTOSOL = ('energy_metabolism', 'adaptive_ArgusUNPUBLISHED_v01')
# the export point: the module's default, with the store and the matrix allocated and some glycogen
CYTOSOL_EXPORT = {'store_alloc': 1.0, 'mx_alloc': 1.0, 'C_glyc_init': GLYC_REST_C}
# the module's rate variable for each reaction, and its state variable for each species
CYTOSOL_RATES = {'HK': 'v_gly', 'LG': 'v_lg', 'HYD': 'v_hyd', 'STORE': 'v_store', 'MOB': 'v_mob', 'MX': 'v_mx'}
# the chemostats' concentrations [mM]: H+ at pH 7 and water, each its standard concentration (as in the detailed
# modules), and the exported matrix, whose activity does not enter (the module's export is irreversible)
CYTOSOL_CHEMOSTATS = {'H': 1e-4, 'H2O': 55500.0, 'mx': 1.0}
# the module's transport into each state [mol/s] (faces and carriers; not reactions of this bond graph)
CYTOSOL_TRANSPORT = {'glc': ('J_glc_in', 'J_glc_out'), 'G6P': ('j_G6P_ext',), 'lac': ('J_lac_in', 'J_lac_out'),
                     'ATP': ('J_ATP_in', 'J_ATP_out', 'j_ATP_ext'), 'ADP': ('J_ADP_in', 'J_ADP_out', 'j_ADP_ext'),
                     'Pi': ('J_Pi_in', 'J_Pi_out', 'j_Pi_ext'), 'glyc': ()}
CYTOSOL_STATES = {'glc': 'C_glc', 'G6P': 'C_G6P', 'lac': 'C_lac', 'ATP': 'C_ATP', 'ADP': 'C_ADP', 'Pi': 'C_Pi',
                  'glyc': 'C_glyc'}


def cytosol_graph():
    """energy_metabolism adaptive with H+ and water as chemostats (the module leaves them implicit): H+ at pH 7 and
    water at its own concentration, each at its standard concentration (C/c_std = 1), so the module's dG0s at pH 7
    apply and the rates are the module's. The species have the formulae and charges of their major microspecies at
    pH 7, so every reaction balances in elements and charge (the matrix as its glucan unit).
    Glycogen has the module's activity C/(K_glyc_mob + C)."""
    n = int(CYTO['n_ATP_mx'])
    assert n == CYTO['n_ATP_mx'], 'n_ATP_mx must be a whole number for a stoichiometry'
    species = [
        Species('glc', 'glucose', F(C=6, H=12, O=6)),
        Species('G6P', 'glucose 6-phosphate', F(C=6, H=11, O=9, P=1), charge=-2),
        Species('lac', 'lactate', F(C=3, H=5, O=3), charge=-1),
        Species('ATP', 'ATP', F(C=10, H=12, N=5, O=13, P=3), charge=-4),
        Species('ADP', 'ADP', F(C=10, H=12, N=5, O=10, P=2), charge=-3),
        Species('Pi', 'phosphate', F(H=1, O=4, P=1), charge=-2),
        Species('glyc', 'glycogen (glucose units)', F(C=6, H=10, O=5), c_std=CYTO['K_glyc_mob'],
                k_sat=CYTO['K_glyc_mob']),
        Species('mx', 'matrix exported (glucose units)', F(C=6, H=10, O=5), boundary=True),
        Species('H', 'H+', F(H=1), boundary=True, charge=1, c_std=CYTOSOL_CHEMOSTATS['H']),
        Species('H2O', 'water', F(H=2, O=1), boundary=True, c_std=CYTOSOL_CHEMOSTATS['H2O']),
    ]
    hk, hyd, ferm, store = CYTO['dG0_hk'], CYTO['dG0_hyd'], CYTO['dG0_ferm'], CYTO['dG0_store']
    reactions = [
        Reaction('HK', 'hexokinase step', {'glc': -1, 'ATP': -1, 'G6P': 1, 'ADP': 1, 'H': 1}, hk, 'dG0_hk'),
        Reaction('LG', 'lower glycolysis', {'G6P': -1, 'ADP': -3, 'Pi': -2, 'H': -1, 'lac': 2, 'ATP': 3, 'H2O': 2},
                 ferm - 2 * hyd - hk, 'dG0_ferm - 2 dG0_hyd - dG0_hk'),
        Reaction('HYD', 'ATP hydrolysis (maintenance)', {'ATP': -1, 'H2O': -1, 'ADP': 1, 'Pi': 1, 'H': 1}, hyd,
                 'dG0_hyd'),
        Reaction('STORE', 'glycogen synthesis', {'G6P': -1, 'ATP': -1, 'H2O': -1, 'glyc': 1, 'ADP': 1, 'Pi': 2, 'H': 1},
                 store - hk, 'dG0_store - dG0_hk'),
        Reaction('MOB', 'glycogen phosphorylase', {'glyc': -1, 'Pi': -1, 'G6P': 1}, hk + hyd - store,
                 'dG0_hk + dG0_hyd - dG0_store'),
        Reaction('MX', 'matrix synthesis and export', {'G6P': -1, 'ATP': -n, 'H2O': -n, 'mx': 1, 'ADP': n, 'Pi': n + 1,
                                                       'H': n},
                 CYTO['dG0_mx'], 'dG0_mx'),
    ]
    return BondGraph('energy_metabolism_adaptive', 'Lumped cytosol (energy_metabolism adaptive)', species, reactions,
                     temperature=CYTO['T'], gas_constant=CYTO['R'])


def evaluate(component, version, values, rates, states):
    """The module's rates and state derivatives at its default with `values` set: ({reaction: v}, {species: d/dt},
    {species: C}) [mM/s, mM]."""
    from cam_testing import checks
    from cam_testing.library import load_version
    with contextlib.redirect_stdout(io.StringIO()):
        cm = checks.ComponentModel(load_version(component, version), work_dir=os.path.join(
            os.environ.get('TMPDIR', '/tmp'), f'bondgraph_{component}'))
        m = cm.helper.model.clone()
    var = {v.qname().split('.')[-1]: v for v in m.variables(deep=True)}
    st = {s.qname().split('.')[-1]: s for s in m.states()}
    for k, x in values.items():
        if k.endswith('_init') and k[:-5] in st:
            st[k[:-5]].set_initial_value(x)
        else:
            var[k + '_mod'].set_rhs(x)
    deriv = dict(zip([s.qname().split('.')[-1] for s in m.states()], m.evaluate_derivatives()))
    return ({r: var[v].eval() for r, v in rates.items()}, {s: deriv[c] for s, c in states.items()},
            {s: st[c].initial_value(as_float=True) for s, c in states.items()})


def build_cytosol():
    bg = cytosol_graph()
    v, _, conc = evaluate(*CYTOSOL, CYTOSOL_EXPORT, CYTOSOL_RATES, CYTOSOL_STATES)
    conc = {**{k: float(x) for k, x in conc.items()}, **CYTOSOL_CHEMOSTATS}
    kappa = bg.calibrate(conc, v)
    print(f'{CYTOSOL[0]} {CYTOSOL[1]}: ' + ', '.join(f'{r} {x:.4g}' for r, x in v.items()) + ' mM/s at the export point')
    return bg, kappa, v, conc


if __name__ == '__main__':
    build_cytosol()
    print('done')
