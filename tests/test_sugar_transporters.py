"""
The glucose carriers AdaptOrDeity's fermenter can evolve after level 1 (transport/glucose_carrier: PTS and ABC, instances
of the ATP version; GalP, the H_symport version), alone against their analytic rate laws, and in the fermenter cell
(cell/generic_adaptive_cell fermenter), whose level-1 behaviour they must not change while their pools are at 0.

    pytest tests/test_sugar_transporters.py --include-unreviewed

tests/data/generic_adaptive_cell_fermenter_level1.json holds the fermenter's traces from before the carriers were added
(the record() function below wrote it): every output of the level-1 cell at 0, 5, ..., 60 s in three scenarios.
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
    'PTS': ('ATP_ArgusUNPUBLISHED_v01', 'PTS', 1.0, 0.0),
    'ABC': ('ATP_ArgusUNPUBLISHED_v01', 'ABC', 2.0, 0.0),
    'GalP': ('H_symport_ArgusUNPUBLISHED_v01', 'default', 0.0, 1.0),
}
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
    a_atp = -p['dG0_hyd'] - p['R'] * p['T'] * math.log(p['ADP'] * p['Pi'] / (p['ATP'] * p['c_std']))
    return p['n_P'] * a_atp


@pytest.mark.parametrize('key', sorted(CARRIERS))
def test_michaelis_menten_uptake(key, carrier):
    '''With no glucose inside, uptake is jmax m A_mem C/(K + C) (times (a/(1 + a))^n_P, a = ATP/K_m_ATP, for the ATP
    carriers): checked from far below to far above Km.'''
    cm = carrier(key)
    p = _p(cm)
    sat = (p['ATP'] / p['K_m_ATP'] / (1 + p['ATP'] / p['K_m_ATP'])) ** p['n_P'] if 'n_P' in p else 1.0
    jmax = p['jmax_glc'] * p['m_expr'] * p['A_mem'] * sat
    K = p['K_glc_t']
    for c in (1e-3 * K, 0.1 * K, K, 10 * K, 1e4 * K):
        J = _flux(cm, {'C_glc_e_mod': c, 'C_glc_mod': 0.0})['J_glc']
        assert J == pytest.approx(jmax * c / (K + c), rel=1e-9), f'{key} at {c} mM'
    # the limits: linear far below Km, saturated far above
    assert _flux(cm, {'C_glc_e_mod': 1e-3 * K, 'C_glc_mod': 0.0})['J_glc'] == pytest.approx(jmax * 1e-3, rel=2e-3)
    assert _flux(cm, {'C_glc_e_mod': 1e4 * K, 'C_glc_mod': 0.0})['J_glc'] == pytest.approx(jmax, rel=2e-4)


# Conditions where the equilibrium gradient is moderate: little ATP free energy for the ATP carriers (ATP 1 uM, ADP and
# Pi 10 mM: A_ATP ~ 19 kJ/mol), the resting proton motive force for GalP
EQ_CONDITIONS = {'PTS': {'ATP_mod': 1e-3, 'ADP_mod': 10.0, 'Pi_mod': 10.0},
                 'ABC': {'ATP_mod': 1e-3, 'ADP_mod': 10.0, 'Pi_mod': 10.0},
                 'GalP': {}}


@pytest.mark.parametrize('key', sorted(CARRIERS))
def test_zero_net_flux_at_equilibrium(key, carrier):
    '''At C_i = C_e exp(coupling free energy / RT), computed here from the inputs, the net flux is zero (relative to the
    one-way flux); a little more glucose inside runs the carrier backwards, a little less forwards.'''
    cm = carrier(key)
    cond = EQ_CONDITIONS[key]
    p = {**_p(cm), **{k[:-4]: v for k, v in cond.items()}}
    c_e = 1e-3
    c_eq = c_e * math.exp(_coupling(key, p) / (p['R'] * p['T']))
    one_way = abs(_flux(cm, {**cond, 'C_glc_e_mod': c_e, 'C_glc_mod': 0.0})['J_glc'])
    at_eq = _flux(cm, {**cond, 'C_glc_e_mod': c_e, 'C_glc_mod': c_eq})
    assert abs(at_eq['J_glc']) <= 1e-9 * one_way, f'{key}: J = {at_eq["J_glc"]:.3g} at equilibrium ({c_eq:.4g} mM)'
    assert abs(at_eq['A_glc']) <= 1e-6 * p['R'] * p['T']
    assert _flux(cm, {**cond, 'C_glc_e_mod': c_e, 'C_glc_mod': 1.01 * c_eq})['J_glc'] < 0
    assert _flux(cm, {**cond, 'C_glc_e_mod': c_e, 'C_glc_mod': 0.99 * c_eq})['J_glc'] > 0


@pytest.mark.parametrize('key', sorted(CARRIERS))
def test_cost_per_glucose(key, carrier):
    '''PTS: 1 ATP equivalent (its PEP), ABC: 2 ATP, GalP: 1 H+ (and its charge) per glucose.'''
    cm = carrier(key)
    _, _, n_atp, n_h = CARRIERS[key]
    out = _flux(cm, {'C_glc_e_mod': 0.05})
    assert out['J_glc'] > 0
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


def _glucose(c):
    '''The same medium glucose for the Glf face and every carrier.'''
    return {f'C_glc_e_mod_{s}': c for s in ('membrane_face', 'glc_PTS', 'glc_ABC', 'glc_GalP')}


def _express(key, m=1.0):
    return {f'm_init_mod_pool_{key}': m, f'm_target_mod_pool_{key}': m}


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


@pytest.mark.parametrize('scenario', sorted(SCENARIOS))
def test_level1_unchanged(scenario, fermenter):
    '''With the carriers' pools at 0 the fermenter reproduces its level-1 traces (relative 1e-7 of each output's
    range; solver steps differ slightly with the added states).'''
    ref = json.load(open(LEVEL1))
    case = ref['scenarios'][scenario]
    t, out = fermenter.run(case['parameters'])
    idx = [int(np.argmin(np.abs(t - x))) for x in ref['times']]
    bad = []
    for name, values in case['outputs'].items():
        assert name in out, f'{name} is no longer an output of the fermenter'
        r = np.asarray(values)
        new = np.asarray(out[name])[idx]
        scale = max(np.max(np.abs(r)), 1e-300)
        if np.max(np.abs(new - r)) > 1e-7 * scale + 1e-30:
            bad.append(f'{name}: differs by {np.max(np.abs(new - r)) / scale:.3g} of its scale')
    assert not bad, '\n'.join(bad)


def test_PTS_uptake_far_above_Glf_at_10uM(fermenter):
    '''At 10 uM glucose (below Glf's Km of 4.1 mM, above the PTS's 5 uM) one unit of PTS takes up glucose tens of
    times faster than the cell's Glf, and the cell with it takes up far more glucose over a minute.'''
    outs = ['mod_glc_PTS/J_glc', 'mod_membrane_face/J_glc']
    t, with_pts = _run(fermenter, {**_glucose(0.01), **_express('PTS')}, outs)
    _, glf_only = _run(fermenter, _glucose(0.01), outs)
    assert with_pts['mod_glc_PTS/J_glc'][0] > 50 * with_pts['mod_membrane_face/J_glc'][0]
    taken_up = np.trapz(with_pts['mod_glc_PTS/J_glc'] + with_pts['mod_membrane_face/J_glc'], t)
    assert taken_up > 20 * np.trapz(glf_only['mod_membrane_face/J_glc'], t)


def _cell_invariant_failures(cm, params):
    t, out = _run(cm, params)
    return checks._check_invariants(cm.spec, t, out, params=cm.param_env(params), where='sweep')


@pytest.mark.parametrize('key,c_e', [('PTS', 0.01), ('ABC', 0.01), ('GalP', 1.0)])
def test_fermenter_invariants_with_carrier(key, c_e, fermenter):
    '''The fermenter's invariants (charge, adenine and phosphoryl conservation, the ATP budget, carbon through every
    carrier, the carriers' stoichiometry, the second law) hold with each carrier expressed.'''
    failures = _cell_invariant_failures(fermenter, {**_glucose(c_e), **_express(key)})
    assert not failures, '\n'.join(failures)


def test_ABC_costs_two_ATP_per_glucose(fermenter):
    '''The ABC importer spends 2 ATP per glucose and its ATP is in the cell's budget: the cytosol's ATP use beyond its
    maintenance ATPase is exactly what the H+-ATPase, protein synthesis, biosynthesis and the carriers take.'''
    outs = ['mod_glc_ABC/J_glc', 'mod_glc_ABC/J_ATP_cost', 'mod_cytosol/J_ATP_turnover', 'mod_cytosol/J_ATP_maint',
            'mod_i_H_ATPase/v_ATP_H', 'mod_protein/J_ATP_total', 'mod_bio/J_ATP_bio']
    _, o = _run(fermenter, {**_glucose(0.01), **_express('ABC')}, outs)
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
    t, o = _run(fermenter, {**_glucose(1.0), **_express('GalP')}, outs)
    assert np.all(o['mod_glc_GalP/J_glc'] > 0)
    assert np.allclose(o['mod_glc_GalP/j_H'], o['mod_glc_GalP/J_glc'], rtol=1e-12, atol=0)
    stored = o['mod_membrane/Vol'][0] * (o['mod_H/H_total_rel'][-1] - o['mod_H/H_total_rel'][0])
    others = np.trapz(sum(o[f] for f in H_FLUXES), t)
    galp = np.trapz(o['mod_glc_GalP/j_H'], t)
    scale = np.trapz(sum(np.abs(o[f]) for f in H_FLUXES + ['mod_glc_GalP/j_H']), t)
    assert abs(stored - (others + galp)) <= 1e-5 * scale
    assert abs(galp) > 1e3 * abs(stored - (others + galp))   # GalP's load is resolved, not lost in the balance
