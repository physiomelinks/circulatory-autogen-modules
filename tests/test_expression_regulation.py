"""
Regulation by falling metabolites (cell/expression_regulation ArgusUNPUBLISHED_v01) and protein turnover in
AdaptOrDeity's fermenter (cell/generic_adaptive_cell fermenter): every protein_pool's target, growth, and the maintenance
(protein synthesis, RNA and lipid upkeep, repair, the maintenance ATPase) are each multiplied by R = prod_s (1 - reg_on_s
+ reg_on_s f_s), f_s = 1/(1 + u^n), u = max(0, 1 - x)/(1 - x_half), x = C/C_rest for ATP, G6P and glucose inside and
x = H_rest/H_i for H+; the switches reg_on_s are 0 by default, and f = 1 at and above rest.

    pytest tests/test_expression_regulation.py --include-unreviewed

Checks: the switches off leave the cell exactly as it was, and on they change nothing at rest; a switch makes a target
follow its species; proteins are degraded slowly (k_deg) and built no faster than their mRNA and the ribosomes allow; and
a cell on half food whose maintenance is regulated by ATP holds its ATP and recovers when food returns, where the same
cell without regulation, or with only its proteins' targets regulated, runs its ATP out and stays there (hexokinase
needs ATP to start glycolysis again: the ATP trap).
"""
import contextlib
import io

import numpy as np
import pytest

from cam_testing import checks, harness
from cam_testing.library import load_version

pytestmark = pytest.mark.module_vv

SPECIES = ['ATP', 'G6P', 'glc', 'H']
# the pools the game offers as sliders whose regulation by ATP spares the cell's energy: the pumps, exchangers, flagellum
# and hydrolase (glycolysis, glucose uptake and lactate export stay unregulated, so food can be used when it returns)
REGULATED = ['H_ATPase', 'NhaA', 'KefC', 'Trk', 'NCX', 'Flagellum', 'Hydrolase']
CARRIER_INPUTS = ('membrane_face', 'glc_PTS', 'glc_ABC', 'glc_GalP')
# the maintenance the game can regulate: (module, switch suffix, regulation factor)
MAINTENANCE = [('protein', '', 'R_reg'), ('bio', '_upkeep', 'R_reg_upkeep'), ('bio', '_repair', 'R_reg_repair'),
               ('cytosol', '', 'R_reg')]


def _maintenance_on(species='ATP'):
    return {f'reg_on_{species}{suffix}_mod_{mod}': 1.0 for mod, suffix, _ in MAINTENANCE}


@pytest.fixture(scope='module')
def fermenter(tmp_path_factory):
    v = load_version('generic_adaptive_cell', 'fermenter')
    return checks.ComponentModel(v, work_dir=str(tmp_path_factory.mktemp('fermenter')))


def _glucose(c):
    return {f'C_glc_e_mod_{s}': c for s in CARRIER_INPUTS}


def _results(helper, outputs):
    return {o: np.asarray(helper.get_results([[harness.output_name(o)]], flatten=True)[0], dtype=float).ravel()
            for o in outputs}


def _phases(cm, phases, outputs, dt=10.0):
    '''Runs the cell through phases [(duration, parameters)], each continuing from the last; returns each phase's
    (t, outputs).'''
    out = []
    helper = harness.simulation_helper(cm.model_path, dict(cm.spec, sim_time=phases[0][0], dt=dt))
    helper.reset_states()
    start = 0.0
    for i, (duration, params) in enumerate(phases):
        if i:   # continue from the last phase's end state over this phase's duration
            helper.update_times(dt, start, duration, 0.0)
        helper.set_param_vals(['parameters/' + k for k in params], list(params.values()), change_states=i == 0)
        start += duration
        with contextlib.redirect_stdout(io.StringIO()):
            assert helper.run() is not False
        out.append((np.asarray(helper.get_time(), dtype=float), _results(helper, outputs)))
    return out


def test_switches_off_leave_every_target_exactly(fermenter):
    '''By default every switch is 0: each pool's R_reg, growth's and the maintenance's are exactly 1 and every pool goes
    to its target (bitwise); the signals themselves are exactly 1 at rest.'''
    pools = [k[len('reg_on_ATP_mod_'):] for k in fermenter.nominal if k.startswith('reg_on_ATP_mod_pool_')]
    assert len(pools) == 14
    assert all(fermenter.nominal[f'reg_on_{s}_mod_{p}'] == 0.0 for s in SPECIES for p in pools + ['bio'])
    assert all(fermenter.nominal[f'reg_on_{s}{suffix}_mod_{mod}'] == 0.0 for s in SPECIES for mod, suffix, _ in MAINTENANCE)
    outs = ([f'mod_{p}/R_reg' for p in pools] + [f'mod_{p}/m_goal' for p in pools] + [f'mod_{p}/m_target' for p in pools]
            + ['mod_bio/R_reg'] + [f'mod_{mod}/{r}' for mod, _, r in MAINTENANCE] + [f'mod_reg/f_{s}_reg' for s in SPECIES])
    _, r = _phases(fermenter, [(60.0, {})], outs, dt=1.0)[0]
    for p in pools + ['bio']:
        assert np.all(r[f'mod_{p}/R_reg'] == 1.0), p
    for mod, _, name in MAINTENANCE:
        assert np.all(r[f'mod_{mod}/{name}'] == 1.0), (mod, name)
    for p in pools:
        assert np.all(r[f'mod_{p}/m_goal'] == r[f'mod_{p}/m_target']), p   # m_target may be logged once
    for s in SPECIES:
        assert r[f'mod_reg/f_{s}_reg'][0] == 1.0, s


def test_switches_on_change_nothing_at_rest(fermenter):
    '''Every switch on (every species, every pool, growth and the maintenance) leaves the resting cell where it was for
    ten minutes: f = 1 at the resting levels, so no target moves (relative 1e-9 on each output's scale).'''
    outs = ['mod_cytosol/C_ATP', 'mod_membrane/V', 'mod_H/pH_i', 'mod_pool_H_ATPase/m', 'mod_pool_proteome/m',
            'mod_protein/J_keep_total', 'mod_cytosol/J_ATP_maint', 'mod_bio/J_RNA', 'mod_bio/R_repair']
    on = {k: 1.0 for k in fermenter.nominal if k.startswith('reg_on_')}
    assert len(on) == 4 * (14 + 1 + len(MAINTENANCE))
    _, a = _phases(fermenter, [(600.0, on)], outs)[0]
    _, b = _phases(fermenter, [(600.0, {})], outs)[0]
    for o in outs:
        np.testing.assert_allclose(a[o], b[o], rtol=0, atol=1e-9 * np.max(np.abs(b[o])), err_msg=o)


@pytest.mark.parametrize('species', SPECIES)
def test_a_switch_makes_the_target_follow_its_species(species, fermenter):
    '''Switching one species on for the H+-ATPase makes its goal its target times that species' f (exact), and growth's
    J_growth = growth_level J_ref throttle f.'''
    params = {f'reg_on_{species}_mod_pool_H_ATPase': 1.0, f'reg_on_{species}_mod_bio': 1.0,
              'growth_level_mod_bio': 1.0, **_glucose(0.01)}
    outs = ['mod_pool_H_ATPase/m_goal', 'mod_pool_H_ATPase/m_target', f'mod_reg/f_{species}_reg', 'mod_bio/J_growth',
            'mod_bio/throttle', 'mod_cytosol/C_ATP']
    _, r = _phases(fermenter, [(60.0, params)], outs, dt=1.0)[0]
    f = r[f'mod_reg/f_{species}_reg']
    np.testing.assert_allclose(r['mod_pool_H_ATPase/m_goal'], r['mod_pool_H_ATPase/m_target'] * f, rtol=1e-12)
    j_ref = fermenter.nominal['J_ref_mod_bio']
    np.testing.assert_allclose(r['mod_bio/J_growth'], j_ref * r['mod_bio/throttle'] * f, rtol=1e-12, atol=1e-35)
    assert np.all((f > 0) & (f <= 1))


def test_targets_fall_as_ATP_falls(fermenter):
    '''Starving (no glucose), with ATP regulation on the pumps, growth and the maintenance: as ATP falls the H+-ATPase's
    goal, growth and the maintenance ATPase fall with it (f_ATP is monotonic in ATP), and ATP takes longer to run out than
    without regulation (below 0.1 mM at ~170 s instead of ~100 s).'''
    outs = ['mod_cytosol/C_ATP', 'mod_reg/f_ATP_reg', 'mod_pool_H_ATPase/m_goal', 'mod_bio/R_reg', 'mod_cytosol/R_reg']
    base = {**_glucose(0.0), 'growth_level_mod_bio': 1.0}
    reg = {**base, 'reg_on_ATP_mod_bio': 1.0, **{f'reg_on_ATP_mod_pool_{p}': 1.0 for p in REGULATED}, **_maintenance_on()}
    t, on = _phases(fermenter, [(300.0, reg)], outs, dt=1.0)[0]
    _, off = _phases(fermenter, [(300.0, base)], outs, dt=1.0)[0]
    atp = on['mod_cytosol/C_ATP']
    assert np.all(np.diff(on['mod_reg/f_ATP_reg']) * np.diff(atp) >= -1e-15)
    for o in ('mod_pool_H_ATPase/m_goal', 'mod_bio/R_reg', 'mod_cytosol/R_reg'):
        np.testing.assert_allclose(on[o], on['mod_reg/f_ATP_reg'], rtol=1e-12, err_msg=o)
    assert on['mod_reg/f_ATP_reg'][-1] < 0.5 * on['mod_reg/f_ATP_reg'][0]

    def gone(a):
        return t[np.argmax(a < 0.1)]
    assert gone(atp) > 1.5 * gone(off['mod_cytosol/C_ATP'])


def test_maintenance_regulation_avoids_the_ATP_trap_on_half_food(fermenter):
    '''On half the resting medium glucose (0.05 mM) for an hour, the unregulated cell runs its ATP out and cannot use full
    food when it returns: hexokinase needs ATP, so ATP stays near 0. Regulating the pumps, exchangers, flagellum and
    hydrolase by ATP no longer helps: their copies fall only as they are degraded (about 0.1% in the hour). With the
    maintenance regulated by ATP (protein synthesis, RNA and lipid upkeep, repair, the maintenance ATPase) ATP holds above
    4 mM and the membrane near -100 mV; when food returns the cell first rebuilds the protein it lost while its synthesis
    was cut (about 0.5% of the proteome, throttled as ATP is low), so ATP rises only slowly (to ~4.7 mM in the hour) while V
    recovers to about -101 mV.'''
    outs = ['mod_cytosol/C_ATP', 'mod_membrane/V', 'mod_pool_H_ATPase/m', 'mod_pool_proteome/m']
    expression = {f'reg_on_ATP_mod_pool_{p}': 1.0 for p in REGULATED}
    fed = _glucose(fermenter.nominal['C_glc_e_mod_membrane_face'])
    (_, off_starve), (_, off_fed) = _phases(fermenter, [(3600.0, _glucose(0.05)), (3600.0, fed)], outs)
    (_, ex_starve), (_, ex_fed) = _phases(fermenter, [(3600.0, {**_glucose(0.05), **expression}), (3600.0, fed)], outs)
    (_, on_starve), (_, on_fed) = _phases(fermenter, [(3600.0, {**_glucose(0.05), **_maintenance_on()}), (3600.0, fed)],
                                          outs)
    assert off_starve['mod_cytosol/C_ATP'][-1] < 1e-3 and off_fed['mod_cytosol/C_ATP'][-1] < 1e-3
    assert ex_fed['mod_cytosol/C_ATP'][-1] < 1e-3
    assert ex_starve['mod_pool_H_ATPase/m'][-1] > 0.995
    assert on_starve['mod_cytosol/C_ATP'].min() > 4.0 and on_starve['mod_membrane/V'][-1] < -99.0
    assert on_fed['mod_cytosol/C_ATP'][-1] > on_starve['mod_cytosol/C_ATP'][-1]
    assert on_fed['mod_membrane/V'][-1] < -100.0
    lost = 1.0 - on_starve['mod_pool_proteome/m'][-1]
    assert 0.003 < lost < 0.007
    assert on_fed['mod_pool_proteome/m'][-1] > on_starve['mod_pool_proteome/m'][-1]


def test_a_protein_above_its_target_falls_only_as_it_is_degraded(fermenter):
    '''Halving the glycolytic enzymes' target: the copies are degraded at k_deg (about 1%/h) and no longer replaced above
    the target, so m = 1/2 + 1/2 exp(-k_deg t) (relative 1e-6 after an hour), and the turnover cost falls with the copies.'''
    k = fermenter.nominal['k_deg_mod_pool_Glyc']
    assert k == pytest.approx(0.01 / 3600, rel=1e-9)
    (t, r), = _phases(fermenter, [(3600.0, {'m_target_mod_pool_Glyc': 0.5})], ['mod_pool_Glyc/m'], dt=60.0)
    np.testing.assert_allclose(r['mod_pool_Glyc/m'], 0.5 + 0.5 * np.exp(-k * t), rtol=1e-6)


def test_a_protein_rises_with_its_mRNA_lifetime_and_within_the_ribosomes_capacity(fermenter):
    '''On food (10 mM glucose), doubling the H+-ATPase (its share of the membrane): a small protein rises with time constant tau = t_half_mRNA / ln 2
    (about 7 min) since its synthesis is well within the ribosomes' capacity: m(tau) = 2 - exp(-1) within 2%. Raising the
    proteome by a fifth asks for more than the ribosomes can make: synthesis spends at most J_cap = build_capacity J_ref
    (about 9 times the resting turnover) and the rise takes hours.'''
    tau = fermenter.nominal['tau_mod_pool_H_ATPase']
    assert tau == pytest.approx(300.0 / np.log(2), rel=1e-9)
    outs = ['mod_pool_H_ATPase/m', 'mod_protein/J_build_total', 'mod_protein/J_cap', 'mod_protein/s_build']
    double = {**_glucose(10.0), 's_target_mod_pool_H_ATPase': 2.0 * fermenter.nominal['s_target_mod_pool_H_ATPase']}
    (t, r), = _phases(fermenter, [(tau, double)], outs, dt=tau / 50)
    assert r['mod_pool_H_ATPase/m'][-1] == pytest.approx(2.0 - np.exp(-1.0), rel=0.02)
    assert np.all(r['mod_protein/s_build'] > 0.999)
    outs = ['mod_pool_proteome/m', 'mod_protein/J_build_total', 'mod_protein/J_cap', 'mod_protein/s_build']
    (t, r), = _phases(fermenter, [(3600.0, {**_glucose(10.0), 'm_target_mod_pool_proteome': 1.2})], outs, dt=60.0)
    assert r['mod_protein/s_build'][1] < 0.5
    assert np.all(r['mod_protein/J_build_total'] <= r['mod_protein/J_cap'] * (1 + 1e-9) + 1e-30)
    assert 1.01 < r['mod_pool_proteome/m'][-1] < 1.15
