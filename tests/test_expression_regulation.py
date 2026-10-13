"""
Regulation by falling metabolites (cell/expression_regulation ArgusUNPUBLISHED_v01) and protein turnover in
AdaptOrDeity's fermenter (cell/generic_adaptive_cell fermenter): every protein_pool's target, growth, and the maintenance
(protein synthesis, RNA and lipid upkeep, repair, the maintenance ATPase) are each multiplied by R = prod_s (1 - reg_on_s
+ reg_on_s f_s), f_s = 1/(1 + u^n), u = max(0, 1 - x)/(1 - x_half), x = C/C_rest for ATP, G6P and glucose inside and
x = H_rest/H_i for H+; the switches reg_on_s are 0 by default, and f = 1 at and above rest.

    pytest tests/test_expression_regulation.py --include-unreviewed

Checks: the switches off leave the cell exactly as it was, and on they change nothing at rest; a switch makes a target
follow its species; proteins turn over slowly (k_deg) and are built no faster than their mRNA and the ribosomes allow;
copies above a lowered target (by the player or by regulation) are removed by proteolysis within minutes, within the
proteases' capacity, and their amino acids are recycled into later building, which then costs less; and a cell on half
food whose maintenance is regulated by ATP holds its ATP and recovers when food returns, where the same cell without
regulation, or with only its proteins' targets regulated, runs its ATP out and stays there (hexokinase needs ATP to start
glycolysis again: the ATP trap).
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
    recycling = ['mod_protein/N_aa', 'mod_protein/J_prot_total', 'mod_protein/J_saved_total', 'mod_protein/J_cat',
                 'mod_protein/J_build_total']
    on = {k: 1.0 for k in fermenter.nominal if k.startswith('reg_on_')}
    assert len(on) == 4 * (14 + 1 + len(MAINTENANCE))
    _, a = _phases(fermenter, [(600.0, on)], outs + recycling)[0]
    _, b = _phases(fermenter, [(600.0, {})], outs + recycling)[0]
    for o in outs:
        np.testing.assert_allclose(a[o], b[o], rtol=0, atol=1e-9 * np.max(np.abs(b[o])), err_msg=o)
    for o in recycling:   # nothing is removed, recycled or built at rest
        assert np.all(a[o] == 0.0) and np.all(b[o] == 0.0), o


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


def test_regulation_by_ATP_removes_copies_within_minutes(fermenter):
    '''Starving (no glucose), with ATP regulation on the pumps: as ATP falls the H+-ATPase's goal falls with it and its
    copies follow within minutes (below 0.95 in a minute, about 0.89 when ATP is spent at ~100 s), removed by proteolysis
    at about 1 ATP per amino acid; it all but stops when ATP is gone (the proteases need ATP: below 2% of its peak rate
    under 1 uM ATP), and the amino acids freed go to
    the recycled pool.'''
    outs = ['mod_cytosol/C_ATP', 'mod_pool_H_ATPase/m', 'mod_pool_H_ATPase/m_goal', 'mod_protein/N_aa',
            'mod_pool_H_ATPase/r_prot']
    params = {**_glucose(0.0), 'growth_level_mod_bio': 1.0, **{f'reg_on_ATP_mod_pool_{p}': 1.0 for p in REGULATED}}
    t, r = _phases(fermenter, [(300.0, params)], outs, dt=1.0)[0]
    m = r['mod_pool_H_ATPase/m']
    assert m[60] < 0.95 and 0.85 < m[-1] < 0.92
    assert np.all(np.diff(m) <= 1e-12)
    spent = r['mod_cytosol/C_ATP'] < 1e-3
    assert spent.any() and np.all(r['mod_pool_H_ATPase/r_prot'][spent] < 2e-2 * r['mod_pool_H_ATPase/r_prot'].max())
    assert r['mod_protein/N_aa'][-1] > 1e6


def test_a_lowered_target_is_removed_within_minutes_and_its_amino_acids_recovered(fermenter):
    '''Halving the surface hydrolase: the copies above the target are removed by proteolysis with tau_prot =
    t_half_proteolysis / ln 2 (about 7 min) on top of turnover (k_deg, about 1%/h), so m = 1/2 + 1/2 exp(-(1/tau_prot +
    k_deg) t) while ATP stays near rest (2e-3). Every amino acid of the copies removed goes to the recycled pool, and each
    costs c_deg ATP to remove; the pool's amino acids are conserved (released = in the pool + reused + catabolised), and
    the ATP returned is J_saved + y_aa_cat per amino acid catabolised (relative 2e-3, trapezoid integration).'''
    nom = fermenter.nominal
    k, tau_p = nom['k_deg_mod_pool_Hydrolase'], nom['tau_prot_mod_pool_Hydrolase']
    assert k == pytest.approx(0.01 / 3600, rel=1e-9)
    assert tau_p == pytest.approx(300.0 / np.log(2), rel=1e-9)
    outs = ['mod_pool_Hydrolase/m', 'mod_pool_Hydrolase/r_prot', 'mod_pool_Hydrolase/J_prot', 'mod_protein/N_aa',
            'mod_protein/Q_aa_rel_total', 'mod_protein/Q_aa_use_total', 'mod_protein/r_aa_cat', 'mod_protein/J_cat',
            'mod_protein/J_saved_total', 'mod_protein/J_recycled', 'mod_protein/J_prot_total']
    half = {'s_target_mod_pool_Hydrolase': 0.5 * nom['s_target_mod_pool_Hydrolase']}
    (t, r), = _phases(fermenter, [(1800.0, half)], outs, dt=10.0)
    m = r['mod_pool_Hydrolase/m']
    np.testing.assert_allclose(m, 0.5 + 0.5 * np.exp(-(1 / tau_p + k) * t), rtol=2e-3)
    assert m[90] < 0.6, 'within 15 minutes'
    n_0, L, N_A = nom['n_0_mod_pool_Hydrolase'], nom['L_mod_pool_Hydrolase'], nom['N_A_mod_protein']
    released = np.trapz(r['mod_protein/Q_aa_rel_total'], t)
    assert released == pytest.approx(n_0 * L * (m[0] - m[-1]), rel=2e-3)
    np.testing.assert_allclose(r['mod_pool_Hydrolase/J_prot'],
                               nom['c_deg_mod_pool_Hydrolase'] * n_0 * L * r['mod_pool_Hydrolase/r_prot'] / N_A, rtol=1e-12)
    np.testing.assert_array_equal(r['mod_protein/J_prot_total'], r['mod_pool_Hydrolase/J_prot'])
    reused, catabolised = np.trapz(r['mod_protein/Q_aa_use_total'], t), np.trapz(r['mod_protein/r_aa_cat'], t)
    assert released == pytest.approx(r['mod_protein/N_aa'][-1] + reused + catabolised, rel=2e-3)
    np.testing.assert_allclose(r['mod_protein/J_cat'], nom['y_aa_cat_mod_protein'] * r['mod_protein/r_aa_cat'] / N_A,
                               rtol=1e-12)
    np.testing.assert_allclose(r['mod_protein/J_recycled'], r['mod_protein/J_saved_total'] + r['mod_protein/J_cat'],
                               rtol=1e-12, atol=1e-35)


def test_a_raised_target_costs_ATP_and_less_with_recycled_amino_acids(fermenter):
    '''On food (10 mM glucose), doubling the H+-ATPase costs the ATP to build its copies, (c_syn + c_poly) per amino acid
    from glucose. Doing it while the glycolytic enzymes are cut to 0.8 (their copies removed by proteolysis, at the
    proteases' capacity) builds the same copies with most of their amino acids recycled: the ATP spent building falls by
    more than half, and the build plus what recycling saved is still the full cost (relative 1e-3).'''
    nom = fermenter.nominal
    outs = ['mod_pool_H_ATPase/m', 'mod_protein/J_build_total', 'mod_protein/J_saved_total', 'mod_protein/phi_aa',
            'mod_protein/s_prot']
    double = {**_glucose(10.0), 's_target_mod_pool_H_ATPase': 2.0 * nom['s_target_mod_pool_H_ATPase']}
    (t, plain), = _phases(fermenter, [(1800.0, double)], outs, dt=10.0)
    (_, swap), = _phases(fermenter, [(1800.0, {**double, 'm_target_mod_pool_Glyc': 0.8})], outs, dt=10.0)
    n_0, L, N_A = nom['n_0_mod_pool_H_ATPase'], nom['L_mod_pool_H_ATPase'], nom['N_A_mod_protein']
    full = (nom['c_syn_mod_pool_H_ATPase'] + nom['c_poly_mod_pool_H_ATPase']) * n_0 * L / N_A
    for r in (plain, swap):
        m = r['mod_pool_H_ATPase/m']
        assert m[-1] > 1.98
        assert np.trapz(r['mod_protein/J_build_total'] + r['mod_protein/J_saved_total'], t) == pytest.approx(
            full * (m[-1] - m[0]), rel=1e-3)
    assert np.all(plain['mod_protein/J_saved_total'] == 0.0)
    assert np.trapz(swap['mod_protein/J_build_total'], t) < 0.5 * np.trapz(plain['mod_protein/J_build_total'], t)
    assert swap['mod_protein/phi_aa'].max() > 0.9 and swap['mod_protein/s_prot'][1] < 0.5


def test_holding_an_amount_costs_its_upkeep_every_second(fermenter):
    '''Keeping copies costs their turnover each second: J_keep = n L (c_poly + c_deg + c_recycle) k_deg min(f_ATP, 1) / N_A
    at or below the goal (exact), so a doubled H+-ATPase costs twice the upkeep once built.'''
    nom = fermenter.nominal
    outs = ['mod_pool_H_ATPase/m', 'mod_pool_H_ATPase/J_keep', 'mod_pool_H_ATPase/f_ATP']
    double = {**_glucose(10.0), 's_target_mod_pool_H_ATPase': 2.0 * nom['s_target_mod_pool_H_ATPase']}
    (t, r), = _phases(fermenter, [(3600.0, double)], outs, dt=60.0)
    c = sum(nom[f'{k}_mod_pool_H_ATPase'] for k in ('c_poly', 'c_deg', 'c_recycle'))
    per_unit = nom['n_0_mod_pool_H_ATPase'] * nom['L_mod_pool_H_ATPase'] * c * nom['k_deg_mod_pool_H_ATPase'] / nom[
        'N_A_mod_protein']
    m, f = r['mod_pool_H_ATPase/m'], np.minimum(r['mod_pool_H_ATPase/f_ATP'], 1.0)
    np.testing.assert_allclose(r['mod_pool_H_ATPase/J_keep'], per_unit * m * f, rtol=1e-9)
    assert m[-1] > 1.999 and r['mod_pool_H_ATPase/J_keep'][-1] == pytest.approx(2 * per_unit, rel=1e-3)


def test_targets_fall_as_ATP_falls(fermenter):
    '''Starving (no glucose), with ATP regulation on the pumps, growth and the maintenance: as ATP falls the H+-ATPase's
    goal, growth and the maintenance ATPase fall with it (f_ATP is monotonic in ATP), and ATP takes longer to run out than
    without regulation (below 0.1 mM at ~145 s instead of ~100 s; removing the pumps' copies costs ATP too).'''
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
    assert gone(atp) > 1.3 * gone(off['mod_cytosol/C_ATP'])


def test_maintenance_regulation_avoids_the_ATP_trap_on_half_food(fermenter):
    '''On half the resting medium glucose (0.05 mM) for an hour, the unregulated cell runs its ATP out and cannot use full
    food when it returns: hexokinase needs ATP, so ATP stays near 0. Regulating the pumps, exchangers, flagellum and
    hydrolase by ATP does not help: their copies are removed within minutes (the H+-ATPase to about half), the membrane
    depolarises and ATP runs out sooner. With the maintenance regulated by ATP (protein synthesis, RNA and lipid upkeep,
    repair, the maintenance ATPase) ATP holds above 4 mM and the membrane near -100 mV, while the proteome loses the copies
    its cut synthesis does not replace (about 0.3% in the hour, at k_deg) to the recycled pool; when food returns the cell
    rebuilds them, mostly from those amino acids, and ATP returns to rest (above 9 mM in the hour) and V to -105 mV.'''
    outs = ['mod_cytosol/C_ATP', 'mod_membrane/V', 'mod_pool_H_ATPase/m', 'mod_pool_proteome/m']
    expression = {f'reg_on_ATP_mod_pool_{p}': 1.0 for p in REGULATED}
    fed = _glucose(fermenter.nominal['C_glc_e_mod_membrane_face'])
    (_, off_starve), (_, off_fed) = _phases(fermenter, [(3600.0, _glucose(0.05)), (3600.0, fed)], outs)
    (_, ex_starve), (_, ex_fed) = _phases(fermenter, [(3600.0, {**_glucose(0.05), **expression}), (3600.0, fed)], outs)
    (_, on_starve), (_, on_fed) = _phases(fermenter, [(3600.0, {**_glucose(0.05), **_maintenance_on()}), (3600.0, fed)],
                                          outs)
    assert off_starve['mod_cytosol/C_ATP'][-1] < 1e-3 and off_fed['mod_cytosol/C_ATP'][-1] < 1e-3
    assert ex_fed['mod_cytosol/C_ATP'][-1] < 1e-3
    assert ex_starve['mod_pool_H_ATPase/m'][30] < 0.8 and ex_starve['mod_pool_H_ATPase/m'][-1] < 0.55
    assert on_starve['mod_cytosol/C_ATP'].min() > 4.0 and on_starve['mod_membrane/V'][-1] < -99.0
    assert on_fed['mod_cytosol/C_ATP'][-1] > 9.0
    assert on_fed['mod_membrane/V'][-1] < -104.0
    lost = 1.0 - on_starve['mod_pool_proteome/m'][-1]
    assert 0.002 < lost < 0.007
    assert on_fed['mod_pool_proteome/m'][-1] > on_starve['mod_pool_proteome/m'][-1]


def test_removing_a_large_protein_is_limited_by_the_proteases(fermenter):
    '''Halving the glycolytic enzymes (a sixth of a billion amino acids to remove) asks the proteases for more than their
    capacity (proteolysis_capacity, 5%/h of the cell's protein): they remove Q_prot_max amino acids per second, so the
    enzymes fall linearly, to about 0.56 in an hour, while ATP stays above 8 mM.'''
    nom = fermenter.nominal
    outs = ['mod_pool_Glyc/m', 'mod_protein/s_prot', 'mod_protein/Q_prot_demand_total', 'mod_cytosol/C_ATP']
    (t, r), = _phases(fermenter, [(3600.0, {'m_target_mod_pool_Glyc': 0.5})], outs, dt=60.0)
    served = r['mod_protein/s_prot'] * r['mod_protein/Q_prot_demand_total']
    np.testing.assert_allclose(served, nom['Q_prot_max_mod_protein'], rtol=1e-9)
    m = r['mod_pool_Glyc/m']
    assert 0.53 < m[-1] < 0.6 and np.all(np.diff(m) < 0)
    assert r['mod_cytosol/C_ATP'].min() > 8.0


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
