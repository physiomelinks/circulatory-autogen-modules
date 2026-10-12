"""
The glucose carriers AdaptOrDeity's fermenter can evolve after level 1 (transport/glucose_carrier: PTS, the
group-translocating version that delivers glucose 6-phosphate; ABC, the ATP version; GalP, the H_symport version), alone
against their analytic rate laws, and in the fermenter cell (cell/generic_adaptive_cell fermenter): the ATP each nets
per glucose, G6P kept in the cell, the membrane shares their pools take, and the level-1 behaviour, which they and the
split glycolysis (hexokinase step and lower glycolysis) must keep while their pools are at 0.

    pytest tests/test_sugar_transporters.py --include-unreviewed

tests/data/generic_adaptive_cell_fermenter_level1.json holds the fermenter's traces from before the carriers were added
(the record() function below wrote it): every output of the level-1 cell at 0, 5, ..., 60 s in three scenarios. The
cell then had one lumped glycolysis and expression targets m_target for every pool; the scenarios are translated to the
membrane proteins' shares (s_target = m_target x the share at rest).
"""
import json
import math
import os

import numpy as np
import pytest

from cam_testing import checks, harness
from cam_testing.library import load_version

pytestmark = pytest.mark.module_vv

HERE = os.path.dirname(os.path.abspath(__file__))
LEVEL1 = os.path.join(HERE, 'data', 'generic_adaptive_cell_fermenter_level1.json')
GC = 'glucose_carrier'
CARRIERS = {   # instance key: (version, instance, ATP equivalents per glucose, H+ per glucose)
    'PTS': ('PTS_ArgusUNPUBLISHED_v01', 'default', 1.0, 0.0),
    'ABC': ('ATP_ArgusUNPUBLISHED_v01', 'default', 2.0, 0.0),
    'GalP': ('H_symport_ArgusUNPUBLISHED_v01', 'default', 0.0, 1.0),
}
INSIDE = {'PTS': 'C_G6P_mod', 'ABC': 'C_glc_mod', 'GalP': 'C_glc_mod'}   # what the carrier delivers into the cell
# relative tolerance of the level-1 traces on each output's scale: the split glycolysis holds ~0.3 uM G6P at rest
LEVEL1_RTOL = 1e-3
SCENARIOS = {   # the level-1 runs recorded in LEVEL1 (fermenter parameters by their generated-model names)
    'rest': {},
    'famine': {'C_glc_e_mod_membrane_face': 0.01},
    'feast': {'C_glc_e_mod_membrane_face': 10.0, 'm_target_mod_pool_H_ATPase': 2.0, 'm_target_mod_pool_Glf': 2.0,
              'm_target_mod_pool_Flagellum': 3.0},
}


# ---- the carriers alone ------------------------------------------------------------------------------------------

def _carrier_model(key, tmp_root):
    version, instance, _, _ = CARRIERS[key]
    v = load_version(GC, version)
    inst = next(i for i in v.instances() if i.name == instance)
    return checks.ComponentModel(v, work_dir=str(tmp_root / key), instance=inst)


@pytest.fixture(scope='module')
def carrier(tmp_path_factory):
    models = {}

    def get(key):
        if key not in models:
            models[key] = _carrier_model(key, tmp_path_factory.mktemp('carriers'))
        return models[key]
    return get


def _flux(cm, params):
    t, out = cm.run(params)
    return {k.split('/')[-1]: float(v[-1]) for k, v in out.items()}


def _p(cm):
    '''Nominal parameters by variable name.'''
    return cm.param_env()


def _coupling(key, p):
    '''Free energy (J/mol) the coupled reaction releases per glucose: n_P A_ATP, or n_H F (E_H - V).'''
    if key == 'GalP':
        return p['n_H'] * p['F'] * (p['E_H'] - p['V']) * 1e-3
    if key == 'PTS':   # glc_e + ATP -> G6P_i + ADP: the hexokinase reaction across the membrane
        return -p['dG0_hk'] + p['R'] * p['T'] * math.log(p['ATP'] / p['ADP'])
    a_atp = -p['dG0_hyd'] - p['R'] * p['T'] * math.log(p['ADP'] * p['Pi'] / (p['ATP'] * p['c_std']))
    return p['n_P'] * a_atp


@pytest.mark.parametrize('key', sorted(CARRIERS))
def test_michaelis_menten_uptake(key, carrier):
    '''With no glucose inside, uptake is jmax m A_mem C/(K + C) (times (a/(1 + a))^n_P, a = ATP/K_m_ATP, for the ATP
    carriers): checked from far below to far above Km.'''
    cm = carrier(key)
    p = _p(cm)
    sat = (p['ATP'] / p['K_m_ATP'] / (1 + p['ATP'] / p['K_m_ATP'])) ** p.get('n_P', 1.0) if 'K_m_ATP' in p else 1.0
    if key == 'PTS':   # the PEP supply: the lower glycolysis' saturation in ADP and Pi
        sat *= p['ADP'] / (p['K_ADP_PEP'] + p['ADP']) * p['Pi'] / (p['K_Pi_PEP'] + p['Pi'])
    jmax = p['jmax_glc'] * p['m_expr'] * p['A_mem'] * sat
    K = p['K_glc_t']
    empty = {INSIDE[key]: 0.0}
    for c in (1e-3 * K, 0.1 * K, K, 10 * K, 1e4 * K):
        J = _flux(cm, {'C_glc_e_mod': c, **empty})['J_glc']
        assert J == pytest.approx(jmax * c / (K + c), rel=1e-9), f'{key} at {c} mM'
    # the limits: linear far below Km, saturated far above
    assert _flux(cm, {'C_glc_e_mod': 1e-3 * K, **empty})['J_glc'] == pytest.approx(jmax * 1e-3, rel=2e-3)
    assert _flux(cm, {'C_glc_e_mod': 1e4 * K, **empty})['J_glc'] == pytest.approx(jmax, rel=2e-4)


# Conditions where the equilibrium gradient is moderate: little ATP free energy for the ATP carriers (ATP 1 uM, ADP and
# Pi 10 mM: A_ATP ~ 19 kJ/mol), the resting proton motive force for GalP
EQ_CONDITIONS = {'PTS': {'ATP_mod': 1e-3, 'ADP_mod': 10.0, 'Pi_mod': 10.0},
                 'ABC': {'ATP_mod': 1e-3, 'ADP_mod': 10.0, 'Pi_mod': 10.0},
                 'GalP': {}}


@pytest.mark.parametrize('key', sorted(CARRIERS))
def test_zero_net_flux_at_equilibrium(key, carrier):
    '''At C_i = C_e exp(coupling free energy / RT) (G6P inside for the PTS), computed here from the inputs, the net flux
    is zero (relative to the one-way flux); a little more inside runs the carrier backwards, a little less forwards.'''
    cm = carrier(key)
    cond = EQ_CONDITIONS[key]
    p = {**_p(cm), **{k[:-4]: v for k, v in cond.items()}}
    c_e, inside = 1e-3, INSIDE[key]
    c_eq = c_e * math.exp(_coupling(key, p) / (p['R'] * p['T']))
    one_way = abs(_flux(cm, {**cond, 'C_glc_e_mod': c_e, inside: 0.0})['J_glc'])
    at_eq = _flux(cm, {**cond, 'C_glc_e_mod': c_e, inside: c_eq})
    assert abs(at_eq['J_glc']) <= 1e-9 * one_way, f'{key}: J = {at_eq["J_glc"]:.3g} at equilibrium ({c_eq:.4g} mM)'
    assert abs(at_eq['A_glc']) <= 1e-6 * p['R'] * p['T']
    assert _flux(cm, {**cond, 'C_glc_e_mod': c_e, inside: 1.01 * c_eq})['J_glc'] < 0
    assert _flux(cm, {**cond, 'C_glc_e_mod': c_e, inside: 0.99 * c_eq})['J_glc'] > 0


@pytest.mark.parametrize('key', sorted(CARRIERS))
def test_cost_per_glucose(key, carrier):
    '''PTS: 1 ATP equivalent (its PEP) per glucose, delivered as G6P (ATP -> ADP, no Pi: the phosphoryl goes onto the
    glucose); ABC: 2 ATP hydrolysed; GalP: 1 H+ (and its charge) per glucose.'''
    cm = carrier(key)
    _, _, n_atp, n_h = CARRIERS[key]
    out = _flux(cm, {'C_glc_e_mod': 0.05})
    assert out['J_glc'] > 0
    if key == 'PTS':
        assert out['J_ATP_cost'] == out['J_G6P'] == out['J_glc'] == -out['J_glc_e']
        assert out['j_ATP'] == -out['J_glc'] and out['j_ADP'] == out['J_glc'] and out['j_Pi'] == 0
        return
    if n_atp:
        assert out['J_ATP_cost'] / out['J_glc'] == pytest.approx(n_atp, rel=1e-12)
        assert out['j_ATP'] == -out['J_ATP_cost'] and out['j_ADP'] == out['j_Pi'] == out['J_ATP_cost']
    else:
        assert out['j_H'] / out['J_glc'] == pytest.approx(n_h, rel=1e-12)
        assert out['i_glc_H'] == pytest.approx(-_p(cm)['F'] * out['j_H'] * 1e9, rel=1e-12)
    assert out['J_glc_e'] == -out['J_glc_c'] == -out['J_glc']


# ---- in the fermenter ------------------------------------------------------------------------------------------------

@pytest.fixture(scope='module')
def fermenter(tmp_path_factory):
    v = load_version('generic_adaptive_cell', 'fermenter')
    return checks.ComponentModel(v, work_dir=str(tmp_path_factory.mktemp('fermenter')))


SHARE_POOLS = ['H_ATPase', 'NhaA', 'KefC', 'Trk', 'NCX', 'Glf', 'LacT', 'Hydrolase', 'PTS', 'ABC',
               'GalP']   # protein_pool membrane_share


def _glucose(c):
    '''The same medium glucose for the Glf face and every carrier.'''
    return {f'C_glc_e_mod_{s}': c for s in ('membrane_face', 'glc_PTS', 'glc_ABC', 'glc_GalP')}


def _ceiling(cm):
    '''The membrane's protein ceiling (um2): protein_expression's A_protein_max, from the parameters.'''
    p = cm.nominal
    a_mem = 2 * math.pi * p['R_cell_mod_membrane'] * p['h_cell_mod_membrane']
    return p['area_max_fraction_mod_protein'] * p['surface_factor_mod_protein'] * a_mem


def _share_of(cm, key, m):
    '''The share of the ceiling m units of a membrane protein take: m n_0 footprint / A_ceiling.'''
    p = cm.nominal
    return m * p[f'n_0_mod_pool_{key}'] * p[f'footprint_mod_pool_{key}'] * 1e-6 / _ceiling(cm)


def _express(cm, key, m=1.0):
    '''Start the pool at m and keep it there (a membrane protein's target is the share m takes).'''
    if key in SHARE_POOLS:
        return {f'm_init_mod_pool_{key}': m, f's_target_mod_pool_{key}': _share_of(cm, key, m)}
    return {f'm_init_mod_pool_{key}': m, f'm_target_mod_pool_{key}': m}


def _translated(cm, params):
    '''Level-1 parameters with the membrane proteins' m_target as shares.'''
    out = {}
    for k, v in params.items():
        key = k.replace('m_target_mod_pool_', '')
        if k.startswith('m_target_mod_pool_') and key in SHARE_POOLS:
            out[f's_target_mod_pool_{key}'] = v * cm.nominal[f's_target_mod_pool_{key}']
        else:
            out[k] = v
    return out


def _run(cm, params, outputs=None):
    outputs = outputs or cm.outputs()
    try:
        t, out = harness.run(cm.helper, outputs, params={'parameters/' + k: v for k, v in params.items()})
    finally:   # back to the nominal parameters, as ComponentModel.run does
        cm.helper.set_param_vals(['parameters/' + k for k in params], [cm.nominal[k] for k in params])
    return t, {**out, **{k.replace('__', '/'): v for k, v in out.items()}}   # model names and module/variable


def record(path=LEVEL1):
    '''Writes LEVEL1 from the current fermenter (run at the commit before the carriers were added).'''
    import tempfile
    cm = checks.ComponentModel(load_version('generic_adaptive_cell', 'fermenter'), work_dir=tempfile.mkdtemp())
    ref = {'times': list(range(0, 61, 5)), 'scenarios': {}}
    for name, p in SCENARIOS.items():
        t, out = cm.run(p)
        idx = [int(np.argmin(np.abs(t - x))) for x in ref['times']]
        ref['scenarios'][name] = {'parameters': p, 'outputs': {k: [float(out[k][i]) for i in idx] for k in sorted(out)}}
    json.dump(ref, open(path, 'w'), indent=1)


def level1_differences(cm, scenario):
    '''max |new - level 1| / max |level 1| of every recorded output in a scenario, largest first.'''
    ref = json.load(open(LEVEL1))
    case = ref['scenarios'][scenario]
    t, out = cm.run(_translated(cm, case['parameters']))
    idx = [int(np.argmin(np.abs(t - x))) for x in ref['times']]
    diffs = []
    for name, values in case['outputs'].items():
        assert name in out, f'{name} is no longer an output of the fermenter'
        r = np.asarray(values)
        scale = max(np.max(np.abs(r)), 1e-300)
        diffs.append((float(np.max(np.abs(np.asarray(out[name])[idx] - r))) / scale, name, scale))
    return sorted(diffs, reverse=True)


@pytest.mark.parametrize('scenario', sorted(SCENARIOS))
def test_level1_unchanged(scenario, fermenter):
    '''With the carriers' pools at 0 the fermenter reproduces its level-1 traces within LEVEL1_RTOL of each output's
    scale (outputs that are zero in level 1, round-off, are not compared). The surface hydrolase, added after the
    reference, takes its share of the membrane ceiling: f_area is higher by exactly that share.'''
    diffs = level1_differences(fermenter, scenario)
    bad = [f'{name}: differs by {d:.3g} of its scale ({scale:.3g})' for d, name, scale in diffs
           if d > LEVEL1_RTOL and scale > 1e-30 and name != 'mod_protein__f_area']
    assert not bad, '\n'.join(bad)
    ref = json.load(open(LEVEL1))
    case = ref['scenarios'][scenario]
    t, out = fermenter.run(_translated(fermenter, case['parameters']))
    idx = [int(np.argmin(np.abs(t - x))) for x in ref['times']]
    hydrolase_share = out['mod_pool_Hydrolase__share'][idx]
    assert np.allclose(np.asarray(out['mod_protein__f_area'])[idx] - case['outputs']['mod_protein__f_area'],
                       hydrolase_share, rtol=LEVEL1_RTOL, atol=0)


def test_PTS_uptake_far_above_Glf_at_10uM(fermenter):
    '''At 10 uM glucose (below Glf's Km of 4.1 mM, above the PTS's 5 uM) one unit of PTS takes up glucose tens of
    times faster than the cell's Glf, and the cell with it takes up several times more glucose over a minute: as much
    as its ATP use needs (its PEP supply slows as ATP builds up), where Glf alone falls short and ATP runs down.'''
    outs = ['mod_glc_PTS/J_glc', 'mod_membrane_face/J_glc', 'mod_cytosol/C_ATP']
    t, with_pts = _run(fermenter, {**_glucose(0.01), **_express(fermenter, 'PTS')}, outs)
    _, glf_only = _run(fermenter, _glucose(0.01), outs)
    assert with_pts['mod_glc_PTS/J_glc'][0] > 50 * with_pts['mod_membrane_face/J_glc'][0]
    taken_up = np.trapz(with_pts['mod_glc_PTS/J_glc'] + with_pts['mod_membrane_face/J_glc'], t)
    assert taken_up > 5 * np.trapz(glf_only['mod_membrane_face/J_glc'], t)
    assert with_pts['mod_cytosol/C_ATP'][-1] >= with_pts['mod_cytosol/C_ATP'][0]
    assert glf_only['mod_cytosol/C_ATP'][-1] < 0.9 * glf_only['mod_cytosol/C_ATP'][0]


def _cell_invariant_failures(cm, params):
    t, out = _run(cm, params)
    return checks._check_invariants(cm.spec, t, out, params=cm.param_env(params), where='sweep')


@pytest.mark.parametrize('key,c_e', [('PTS', 0.01), ('ABC', 0.01), ('GalP', 1.0)])
def test_fermenter_invariants_with_carrier(key, c_e, fermenter):
    '''The fermenter's invariants (charge, adenine and phosphoryl conservation, the ATP budget, carbon through every
    carrier, the carriers' stoichiometry, the second law) hold with each carrier expressed.'''
    failures = _cell_invariant_failures(fermenter, {**_glucose(c_e), **_express(fermenter, key)})
    assert not failures, '\n'.join(failures)


def test_ABC_costs_two_ATP_per_glucose(fermenter):
    '''The ABC importer spends 2 ATP per glucose and its ATP is in the cell's budget: the cytosol's ATP use beyond its
    maintenance ATPase is exactly what the H+-ATPase, protein synthesis, biosynthesis and the carriers take.'''
    outs = ['mod_glc_ABC/J_glc', 'mod_glc_ABC/J_ATP_cost', 'mod_cytosol/J_ATP_turnover', 'mod_cytosol/J_ATP_maint',
            'mod_i_H_ATPase/v_ATP_H', 'mod_protein/J_ATP_total', 'mod_bio/J_ATP_bio']
    _, o = _run(fermenter, {**_glucose(0.01), **_express(fermenter, 'ABC')}, outs)
    assert np.all(o['mod_glc_ABC/J_glc'] > 0)
    assert np.allclose(o['mod_glc_ABC/J_ATP_cost'] / o['mod_glc_ABC/J_glc'], 2.0, rtol=1e-12, atol=0)
    assert np.allclose(o['mod_cytosol/J_ATP_turnover'] - o['mod_cytosol/J_ATP_maint'],
                       o['mod_i_H_ATPase/v_ATP_H'] + o['mod_protein/J_ATP_total'] + o['mod_bio/J_ATP_bio']
                       + o['mod_glc_ABC/J_ATP_cost'], rtol=1e-9, atol=1e-30)


H_FLUXES = ['mod_i_leak_H/j_leak_H', 'mod_i_H_ATPase/j_H_ATPase', 'mod_i_NaH/j_NaH_H', 'mod_i_KH/j_KH_H',
            'mod_membrane_face/J_H_lac', 'mod_i_flagellum/j_H_flag', 'mod_cytosol/j_H_met']


def test_H_symport_adds_one_H_per_glucose(fermenter):
    '''GalP brings one H+ in per glucose, and the cell's proton balance closes only with it: Vol d(H_i - beta pH_i) =
    int (every H+ flux, GalP's included) dt.'''
    outs = H_FLUXES + ['mod_glc_GalP/j_H', 'mod_glc_GalP/J_glc', 'mod_H/H_total_rel', 'mod_membrane/Vol']
    t, o = _run(fermenter, {**_glucose(1.0), **_express(fermenter, 'GalP')}, outs)
    assert np.all(o['mod_glc_GalP/J_glc'] > 0)
    assert np.allclose(o['mod_glc_GalP/j_H'], o['mod_glc_GalP/J_glc'], rtol=1e-12, atol=0)
    stored = o['mod_membrane/Vol'][0] * (o['mod_H/H_total_rel'][-1] - o['mod_H/H_total_rel'][0])
    others = np.trapz(sum(o[f] for f in H_FLUXES), t)
    galp = np.trapz(o['mod_glc_GalP/j_H'], t)
    scale = np.trapz(sum(np.abs(o[f]) for f in H_FLUXES + ['mod_glc_GalP/j_H']), t)
    assert abs(stored - (others + galp)) <= 1e-5 * scale
    assert abs(galp) > 1e3 * abs(stored - (others + galp))   # GalP's load is resolved, not lost in the balance


# ---- energy per glucose, G6P, membrane shares ------------------------------------------------------------------------

def _alone(cm, key, m=1.0):
    '''The cell with one glucose carrier only: Glf removed (its pool at 0) and the carrier at m.'''
    if key == 'Glf':
        return _express(cm, 'Glf', m)
    return {**_express(cm, 'Glf', 0.0), **_express(cm, key, m)}


@pytest.mark.parametrize('key,net', [('Glf', 2.0), ('PTS', 2.0), ('ABC', 0.0), ('GalP', 2.0)])
def test_net_ATP_per_glucose(key, net, fermenter):
    '''With one carrier and glucose at 1 mM, per glucose taken up and used (less what accumulates as glucose or G6P
    inside, from the logged pools): glycolysis makes 2 ATP net per free glucose (3 from the lower glycolysis less the
    hexokinase ATP) and 3 per G6P the PTS delivers, of which the PTS spent 1 (its PEP); the ABC importer spends 2. Net:
    Glf 2, PTS 2, ABC 0, GalP 2 in ATP, plus 1 H+ that the H+-ATPase pumps out for 3/10 ATP (10/3 H+ per ATP), 1.7 in
    all.'''
    carrier = {'Glf': 'membrane_face', 'PTS': 'glc_PTS', 'ABC': 'glc_ABC', 'GalP': 'glc_GalP'}[key]
    outs = [f'mod_{carrier}/J_glc', 'mod_cytosol/J_ATP_gly', 'mod_cytosol/J_ATP_hk', 'mod_cytosol/J_ATP_lg',
            'mod_cytosol/C_G6P', 'mod_cytosol/C_glc', 'mod_membrane/Vol']
    outs += [f'mod_{carrier}/J_ATP_cost'] if key in ('PTS', 'ABC') else []
    outs += ['mod_glc_GalP/j_H'] if key == 'GalP' else []
    t, o = _run(fermenter, {**_glucose(1.0), **_alone(fermenter, key, 4.0)}, outs)
    assert o[f'mod_{carrier}/J_glc'][-1] > 0
    # glucose used at the end: taken up less what accumulates inside (as glucose or G6P)
    vol = o['mod_membrane/Vol'][0]
    J = o[f'mod_{carrier}/J_glc'][-1] - vol * (np.gradient(o['mod_cytosol/C_glc'], t)[-1] +
                                               np.gradient(o['mod_cytosol/C_G6P'], t)[-1])
    cost = o[f'mod_{carrier}/J_ATP_cost'][-1] if key in ('PTS', 'ABC') else 0.0
    assert (o['mod_cytosol/J_ATP_gly'][-1] - cost) / J == pytest.approx(net, rel=1e-3, abs=1e-3)
    if key == 'PTS':   # no hexokinase: the G6P comes from the carrier
        assert o['mod_cytosol/J_ATP_hk'][-1] <= 1e-6 * o['mod_cytosol/J_ATP_lg'][-1]
        assert o['mod_cytosol/J_ATP_lg'][-1] / J == pytest.approx(3.0, rel=1e-4)
    if key == 'GalP':
        n_H_pump = fermenter.nominal['n_H_mod_i_H_ATPase']
        h_per_glc = o['mod_glc_GalP/j_H'][-1] / o['mod_glc_GalP/J_glc'][-1]
        assert h_per_glc == pytest.approx(1.0, rel=1e-12)
        assert net - h_per_glc / n_H_pump == pytest.approx(1.7, rel=1e-12)


def test_PTS_delivers_G6P_and_skips_hexokinase(fermenter):
    '''In the PTS-only cell at 10 uM glucose the glucose arrives as G6P: free glucose inside runs down, the hexokinase
    step stops, and G6P is held above its resting level by the PTS's flux.'''
    outs = ['mod_glc_PTS/J_glc', 'mod_glc_PTS/J_G6P', 'mod_cytosol/C_glc', 'mod_cytosol/C_G6P', 'mod_cytosol/v_gly',
            'mod_cytosol/v_lg']
    t, o = _run(fermenter, {**_glucose(0.01), **_alone(fermenter, 'PTS', 4.0)}, outs)
    assert np.array_equal(o['mod_glc_PTS/J_G6P'], o['mod_glc_PTS/J_glc'])
    assert o['mod_cytosol/C_glc'][-1] < 1e-3 * o['mod_cytosol/C_glc'][0]
    assert o['mod_cytosol/v_gly'][-1] < 1e-6 * o['mod_cytosol/v_lg'][-1]
    assert o['mod_cytosol/C_G6P'][-1] > fermenter.nominal['C_G6P_init_mod_cytosol']


def test_G6P_is_trapped(fermenter):
    '''A cell loaded with 1 mM G6P in a medium without glucose spends it as lactate: G6P has no way out, and the
    glucose Glf lets out is only what reverse hexokinase makes at equilibrium (G6P ADP/(K_eq_hk ATP), ~5e-5 mM).'''
    outs = ['mod_membrane_face/J_glc', 'mod_membrane_face/J_lac', 'mod_cytosol/C_G6P', 'mod_cytosol/C_glc',
            'mod_cytosol/C_lac', 'mod_cytosol/C_glyc', 'mod_membrane/Vol']
    t, o = _run(fermenter, {**_glucose(0.0), 'C_G6P_init_mod_cytosol': 1.0, 'C_glc_init_mod_cytosol': 0.0}, outs)
    vol = o['mod_membrane/Vol'][0]
    g6p_spent = vol * (o['mod_cytosol/C_G6P'][0] - o['mod_cytosol/C_G6P'][-1])
    glc_out = -np.trapz(np.minimum(o['mod_membrane_face/J_glc'], 0.0), t)
    assert g6p_spent > 0.9 * vol * 1.0
    assert glc_out <= 1e-3 * g6p_spent
    # the carbon of the G6P left as lactate (2 per G6P), not as G6P
    lac_out = -np.trapz(o['mod_membrane_face/J_lac'], t) + vol * (o['mod_cytosol/C_lac'][-1] - o['mod_cytosol/C_lac'][0])
    assert lac_out >= 2 * g6p_spent * (1 - 2e-3)


def test_membrane_shares(fermenter):
    '''A membrane protein's target is its share of the membrane's protein ceiling: copies = s A_ceiling / footprint.
    At rest every transporter holds its share at rest (m = 1); with the three glucose carriers given the rest of the
    ceiling the shares sum to 1 and the copies fill the ceiling exactly (the flagellum, counted by number, comes on top).'''
    nominal = fermenter.nominal
    shares = {k: nominal[f's_target_mod_pool_{k}'] for k in SHARE_POOLS}
    for k in SHARE_POOLS:
        assert shares[k] == pytest.approx(_share_of(fermenter, k, 0.0 if k in CARRIERS else 1.0), rel=1e-12, abs=0)
    free = 1.0 - sum(shares.values())
    assert 0 < free < 1
    fill = {k: free / 3 for k in CARRIERS}
    params = {**{f's_target_mod_pool_{k}': s for k, s in fill.items()},
              **{f'm_init_mod_pool_{k}': s * _ceiling(fermenter) / (nominal[f'n_0_mod_pool_{k}'] *
                                                                   nominal[f'footprint_mod_pool_{k}'] * 1e-6)
                 for k, s in fill.items()}}
    outs = [f'mod_pool_{k}/{v}' for k in SHARE_POOLS for v in ('share', 'm', 'm_target')] + [
        'mod_protein/A_protein_max', 'mod_protein/A_protein_total', 'mod_pool_Flagellum/A_protein']
    t, o = _run(fermenter, params, outs)
    assert np.allclose(o['mod_protein/A_protein_max'], _ceiling(fermenter), rtol=1e-12, atol=0)
    total = sum(o[f'mod_pool_{k}/share'] for k in SHARE_POOLS)
    assert np.allclose(total, 1.0, rtol=1e-9, atol=0)
    assert np.allclose(o['mod_protein/A_protein_total'] - o['mod_pool_Flagellum/A_protein'],
                       o['mod_protein/A_protein_max'], rtol=1e-9, atol=0)
    for k in SHARE_POOLS:
        assert np.allclose(o[f'mod_pool_{k}/share'], {**shares, **fill}[k], rtol=1e-9, atol=1e-15), k


@pytest.mark.parametrize('c_e,share', [(1.0, 0.25), (50.0, 1.0)])
def test_PTS_cell_keeps_its_phosphate(c_e, share, fermenter):
    '''A PTS-only cell on rich glucose takes up what its ATP use needs: G6P settles where the PTS and lower glycolysis
    capacities balance, far below the cell's phosphate, and ATP stays up (no sugar-phosphate trap), up to a PTS
    filling the whole membrane ceiling in 50 mM glucose.'''
    outs = ['mod_cytosol/C_G6P', 'mod_cytosol/C_ATP', 'mod_cytosol/C_Pi', 'mod_glc_PTS/J_glc']
    m = share / _share_of(fermenter, 'PTS', 1.0)
    t, o = _run(fermenter, {**_glucose(c_e), **_alone(fermenter, 'PTS', m)}, outs)
    assert np.max(o['mod_cytosol/C_G6P']) < 0.1
    assert np.min(o['mod_cytosol/C_ATP']) > 0.9 * o['mod_cytosol/C_ATP'][0]
    assert np.min(o['mod_cytosol/C_Pi']) > 0.9 * o['mod_cytosol/C_Pi'][0]


# ---- the cell-surface polysaccharide hydrolase -------------------------------------------------------------------------

@pytest.mark.parametrize('m', [0.0, 1.0, 3.0])
def test_hydrolase_capacity_scales_with_expression(m, fermenter):
    '''The hydrolase's capacity V_hydrolysis_max = kcat n_0 m / N_A follows its expression; at the level-1 baseline it
    is hydrolysis_capacity (10) times the cell's resting glucose uptake.'''
    outs = ['mod_hydrolase/V_hydrolysis_max', 'mod_hydrolase/n_copies', 'mod_pool_Hydrolase/m',
            'mod_membrane_face/J_glc']
    t, o = _run(fermenter, _express(fermenter, 'Hydrolase', m), outs)
    p = fermenter.nominal
    assert np.allclose(o['mod_pool_Hydrolase/m'], m, rtol=1e-9, atol=1e-12)
    assert np.allclose(o['mod_hydrolase/V_hydrolysis_max'], p['kcat_mod_hydrolase'] * p['n_0_mod_hydrolase'] *
                       o['mod_pool_Hydrolase/m'] / p['N_A_mod_hydrolase'], rtol=1e-12, atol=0)
    if m == 1.0:   # the calibrated cell at rest
        assert o['mod_hydrolase/V_hydrolysis_max'][0] == pytest.approx(10.0 * o['mod_membrane_face/J_glc'][0], rel=1e-9)


def test_hydrolase_pool_costs(fermenter):
    '''Keeping the hydrolase costs its turnover, n L (c_poly + c_deg + c_recycle) k_deg / N_A, and quadrupling it costs
    building copies, n_0 L (c_syn + c_poly) / N_A per unit of m: raised to 1.5, what is built is paid for.'''
    outs = ['mod_pool_Hydrolase/J_keep', 'mod_pool_Hydrolase/J_build', 'mod_pool_Hydrolase/m', 'mod_pool_Hydrolase/f_ATP',
            'mod_pool_Hydrolase/M_build']
    p = {k.replace('_mod_pool_Hydrolase', ''): v for k, v in fermenter.nominal.items() if k.endswith('_mod_pool_Hydrolase')}
    t, o = _run(fermenter, {'s_target_mod_pool_Hydrolase': 1.5 * _share_of(fermenter, 'Hydrolase', 1.0)}, outs)
    keep = o['mod_pool_Hydrolase/m'] * p['n_0'] * p['L'] * (p['c_poly'] + p['c_deg'] + p['c_recycle']) * p['k_deg'] / \
        fermenter.nominal['N_A_mod_pool_Hydrolase'] * o['mod_pool_Hydrolase/f_ATP']
    assert np.allclose(o['mod_pool_Hydrolase/J_keep'], keep, rtol=1e-12, atol=0)
    unit = p['n_0'] * p['L'] * (p['c_syn'] + p['c_poly']) / fermenter.nominal['N_A_mod_pool_Hydrolase']
    built = o['mod_pool_Hydrolase/m'][-1] - o['mod_pool_Hydrolase/m'][0]
    assert built > 0.4
    assert np.trapz(o['mod_pool_Hydrolase/J_build'], t) == pytest.approx(unit * built, rel=1e-4)
    assert o['mod_pool_Hydrolase/M_build'][-1] == pytest.approx(unit * o['mod_pool_Hydrolase/m'][-1], rel=1e-12)


def test_fermenter_invariants_with_more_hydrolase(fermenter):
    '''The hydrolase imports nothing: with four times the enzyme the cell's charge, mass, carbon, ATP-budget and
    second-law invariants hold as before (it costs only protein).'''
    failures = _cell_invariant_failures(fermenter, {**_glucose(1.0), **_express(fermenter, 'Hydrolase', 4.0)})
    assert not failures, '\n'.join(failures)
