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
proteases' capacity, and their amino acids are recycled into later building, which then costs less; a cell on half
food whose maintenance is regulated by ATP holds its ATP and recovers when food returns, where the same cell without
regulation, or with only its proteins' targets regulated, runs its ATP out and stays there (hexokinase needs ATP to start
glycolysis again: the ATP trap); and each maintenance cut has its cost: cutting the RNA and lipid upkeep damages the
envelope and loses ribosomes, rebuilt at de novo cost when food returns, and cutting the maintenance ATPase damages
ribosomes and proteins.
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
    return {f'reg_on_{species}{suffix}_mod__{mod}': 1.0 for mod, suffix, _ in MAINTENANCE}


@pytest.fixture(scope='module')
def fermenter(tmp_path_factory):
    v = load_version('generic_adaptive_cell', 'fermenter')
    return checks.ComponentModel(v, work_dir=str(tmp_path_factory.mktemp('fermenter')))


def _glucose(c):
    return {f'C_glc_e_mod__{s}': c for s in CARRIER_INPUTS}


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
    pools = [k[len('reg_on_ATP_mod__'):] for k in fermenter.nominal if k.startswith('reg_on_ATP_mod__pool_')]
    assert len(pools) == 14
    assert all(fermenter.nominal[f'reg_on_{s}_mod__{p}'] == 0.0 for s in SPECIES for p in pools + ['bio'])
    assert all(fermenter.nominal[f'reg_on_{s}{suffix}_mod__{mod}'] == 0.0 for s in SPECIES for mod, suffix, _ in MAINTENANCE)
    outs = ([f'mod/{p}/R_reg' for p in pools] + [f'mod/{p}/m_goal' for p in pools] + [f'mod/{p}/m_target' for p in pools]
            + ['mod/bio/R_reg'] + [f'mod/{mod}/{r}' for mod, _, r in MAINTENANCE] + [f'mod/reg/f_{s}_reg' for s in SPECIES])
    _, r = _phases(fermenter, [(60.0, {})], outs, dt=1.0)[0]
    for p in pools + ['bio']:
        assert np.all(r[f'mod/{p}/R_reg'] == 1.0), p
    for mod, _, name in MAINTENANCE:
        assert np.all(r[f'mod/{mod}/{name}'] == 1.0), (mod, name)
    for p in pools:
        assert np.all(r[f'mod/{p}/m_goal'] == r[f'mod/{p}/m_target']), p   # m_target may be logged once
    for s in SPECIES:
        assert r[f'mod/reg/f_{s}_reg'][0] == 1.0, s


def test_switches_on_change_nothing_at_rest(fermenter):
    '''Every switch on (every species, every pool, growth and the maintenance) leaves the resting cell where it was for
    ten minutes: f = 1 at the resting levels, so no target moves (relative 1e-9 on each output's scale).'''
    outs = ['mod/cytosol/C_ATP', 'mod/membrane/V', 'mod/H/pH_i', 'mod/pool_H_ATPase/m', 'mod/pool_proteome/m',
            'mod/protein/J_keep_total', 'mod/cytosol/J_ATP_maint', 'mod/bio/J_RNA', 'mod/bio/R_repair']
    recycling = ['mod/protein/N_aa', 'mod/protein/J_prot_total', 'mod/protein/J_saved_total', 'mod/protein/J_cat',
                 'mod/protein/J_build_total']
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
    params = {f'reg_on_{species}_mod__pool_H_ATPase': 1.0, f'reg_on_{species}_mod__bio': 1.0,
              'growth_level_mod__bio': 1.0, **_glucose(0.01)}
    outs = ['mod/pool_H_ATPase/m_goal', 'mod/pool_H_ATPase/m_target', f'mod/reg/f_{species}_reg', 'mod/bio/J_growth',
            'mod/bio/throttle', 'mod/cytosol/C_ATP']
    _, r = _phases(fermenter, [(60.0, params)], outs, dt=1.0)[0]
    f = r[f'mod/reg/f_{species}_reg']
    np.testing.assert_allclose(r['mod/pool_H_ATPase/m_goal'], r['mod/pool_H_ATPase/m_target'] * f, rtol=1e-12)
    j_ref = fermenter.nominal['J_ref_mod__bio']
    np.testing.assert_allclose(r['mod/bio/J_growth'], j_ref * r['mod/bio/throttle'] * f, rtol=1e-12, atol=1e-35)
    assert np.all((f > 0) & (f <= 1))


def test_regulation_by_ATP_removes_copies_within_minutes(fermenter):
    '''Starving (no glucose), with ATP regulation on the pumps: as ATP falls the H+-ATPase's goal falls with it and its
    copies follow within minutes (below 0.95 in a minute, about 0.89 when ATP is spent at ~100 s), removed by proteolysis
    at about 1 ATP per amino acid; it all but stops when ATP is gone (the proteases need ATP: below 2% of its peak rate
    under 1 uM ATP), and the amino acids freed go to
    the recycled pool.'''
    outs = ['mod/cytosol/C_ATP', 'mod/pool_H_ATPase/m', 'mod/pool_H_ATPase/m_goal', 'mod/protein/N_aa',
            'mod/pool_H_ATPase/r_prot']
    params = {**_glucose(0.0), 'growth_level_mod__bio': 1.0, **{f'reg_on_ATP_mod__pool_{p}': 1.0 for p in REGULATED}}
    t, r = _phases(fermenter, [(300.0, params)], outs, dt=1.0)[0]
    m = r['mod/pool_H_ATPase/m']
    assert m[60] < 0.95 and 0.85 < m[-1] < 0.92
    assert np.all(np.diff(m) <= 1e-12)
    spent = r['mod/cytosol/C_ATP'] < 1e-3
    assert spent.any() and np.all(r['mod/pool_H_ATPase/r_prot'][spent] < 2e-2 * r['mod/pool_H_ATPase/r_prot'].max())
    assert r['mod/protein/N_aa'][-1] > 1e6


def test_a_lowered_target_is_removed_within_minutes_and_its_amino_acids_recovered(fermenter):
    '''Halving the surface hydrolase: the copies above the target are removed by proteolysis with tau_prot =
    t_half_proteolysis / ln 2 (about 7 min) on top of turnover (k_deg, about 1%/h), so m = 1/2 + 1/2 exp(-(1/tau_prot +
    k_deg) t) while ATP stays near rest (2e-3). Every amino acid of the copies removed goes to the recycled pool, and each
    costs c_deg ATP to remove; the pool's amino acids are conserved (released = in the pool + reused + catabolised), and
    the ATP returned is J_saved + y_aa_cat per amino acid catabolised (relative 2e-3, trapezoid integration).'''
    nom = fermenter.nominal
    k, tau_p = nom['k_deg_mod__pool_Hydrolase'], nom['tau_prot_mod__pool_Hydrolase']
    assert k == pytest.approx(0.01 / 3600, rel=1e-9)
    assert tau_p == pytest.approx(300.0 / np.log(2), rel=1e-9)
    outs = ['mod/pool_Hydrolase/m', 'mod/pool_Hydrolase/r_prot', 'mod/pool_Hydrolase/J_prot', 'mod/protein/N_aa',
            'mod/protein/Q_aa_rel_total', 'mod/protein/Q_aa_use_total', 'mod/protein/r_aa_cat', 'mod/protein/J_cat',
            'mod/protein/J_saved_total', 'mod/protein/J_recycled', 'mod/protein/J_prot_total']
    half = {'s_target_mod__pool_Hydrolase': 0.5 * nom['s_target_mod__pool_Hydrolase']}
    (t, r), = _phases(fermenter, [(1800.0, half)], outs, dt=10.0)
    m = r['mod/pool_Hydrolase/m']
    np.testing.assert_allclose(m, 0.5 + 0.5 * np.exp(-(1 / tau_p + k) * t), rtol=2e-3)
    assert m[90] < 0.6, 'within 15 minutes'
    n_0, L, N_A = nom['n_0_mod__pool_Hydrolase'], nom['L_mod__pool_Hydrolase'], nom['N_A_mod__protein']
    released = np.trapz(r['mod/protein/Q_aa_rel_total'], t)
    assert released == pytest.approx(n_0 * L * (m[0] - m[-1]), rel=2e-3)
    np.testing.assert_allclose(r['mod/pool_Hydrolase/J_prot'],
                               nom['c_deg_mod__pool_Hydrolase'] * n_0 * L * r['mod/pool_Hydrolase/r_prot'] / N_A, rtol=1e-12)
    np.testing.assert_array_equal(r['mod/protein/J_prot_total'], r['mod/pool_Hydrolase/J_prot'])
    reused, catabolised = np.trapz(r['mod/protein/Q_aa_use_total'], t), np.trapz(r['mod/protein/r_aa_cat'], t)
    assert released == pytest.approx(r['mod/protein/N_aa'][-1] + reused + catabolised, rel=2e-3)
    np.testing.assert_allclose(r['mod/protein/J_cat'], nom['y_aa_cat_mod__protein'] * r['mod/protein/r_aa_cat'] / N_A,
                               rtol=1e-12)
    np.testing.assert_allclose(r['mod/protein/J_recycled'], r['mod/protein/J_saved_total'] + r['mod/protein/J_cat'],
                               rtol=1e-12, atol=1e-35)


def test_a_raised_target_costs_ATP_and_less_with_recycled_amino_acids(fermenter):
    '''On food (10 mM glucose), doubling the H+-ATPase costs the ATP to build its copies, (c_syn + c_poly) per amino acid
    from glucose. Doing it while the glycolytic enzymes are cut to 0.8 (their copies removed by proteolysis, at the
    proteases' capacity) builds the same copies with most of their amino acids recycled: the ATP spent building falls by
    more than half, and the build plus what recycling saved is still the full cost (relative 1e-3).'''
    nom = fermenter.nominal
    outs = ['mod/pool_H_ATPase/m', 'mod/protein/J_build_total', 'mod/protein/J_saved_total', 'mod/protein/phi_aa',
            'mod/protein/s_prot']
    double = {**_glucose(10.0), 's_target_mod__pool_H_ATPase': 2.0 * nom['s_target_mod__pool_H_ATPase']}
    (t, plain), = _phases(fermenter, [(1800.0, double)], outs, dt=10.0)
    (_, swap), = _phases(fermenter, [(1800.0, {**double, 'm_target_mod__pool_Glyc': 0.8})], outs, dt=10.0)
    n_0, L, N_A = nom['n_0_mod__pool_H_ATPase'], nom['L_mod__pool_H_ATPase'], nom['N_A_mod__protein']
    full = (nom['c_syn_mod__pool_H_ATPase'] + nom['c_poly_mod__pool_H_ATPase']) * n_0 * L / N_A
    for r in (plain, swap):
        m = r['mod/pool_H_ATPase/m']
        assert m[-1] > 1.98
        assert np.trapz(r['mod/protein/J_build_total'] + r['mod/protein/J_saved_total'], t) == pytest.approx(
            full * (m[-1] - m[0]), rel=1e-3)
    assert np.all(plain['mod/protein/J_saved_total'] == 0.0)
    assert np.trapz(swap['mod/protein/J_build_total'], t) < 0.5 * np.trapz(plain['mod/protein/J_build_total'], t)
    assert swap['mod/protein/phi_aa'].max() > 0.9 and swap['mod/protein/s_prot'][1] < 0.5


def test_holding_an_amount_costs_its_upkeep_every_second(fermenter):
    '''Keeping copies costs their turnover each second: J_keep = n L (c_poly + c_deg + c_recycle) k_deg min(f_ATP, 1) / N_A
    at or below the goal (exact), so a doubled H+-ATPase costs twice the upkeep once built.'''
    nom = fermenter.nominal
    outs = ['mod/pool_H_ATPase/m', 'mod/pool_H_ATPase/J_keep', 'mod/pool_H_ATPase/f_ATP']
    double = {**_glucose(10.0), 's_target_mod__pool_H_ATPase': 2.0 * nom['s_target_mod__pool_H_ATPase']}
    (t, r), = _phases(fermenter, [(3600.0, double)], outs, dt=60.0)
    c = sum(nom[f'{k}_mod__pool_H_ATPase'] for k in ('c_poly', 'c_deg', 'c_recycle'))
    per_unit = nom['n_0_mod__pool_H_ATPase'] * nom['L_mod__pool_H_ATPase'] * c * nom['k_deg_mod__pool_H_ATPase'] / nom[
        'N_A_mod__protein']
    m, f = r['mod/pool_H_ATPase/m'], np.minimum(r['mod/pool_H_ATPase/f_ATP'], 1.0)
    np.testing.assert_allclose(r['mod/pool_H_ATPase/J_keep'], per_unit * m * f, rtol=1e-9)
    assert m[-1] > 1.999 and r['mod/pool_H_ATPase/J_keep'][-1] == pytest.approx(2 * per_unit, rel=1e-3)


def test_targets_fall_as_ATP_falls(fermenter):
    '''Starving (no glucose), with ATP regulation on the pumps, growth and the maintenance: as ATP falls the H+-ATPase's
    goal, growth and the maintenance ATPase fall with it (f_ATP is monotonic in ATP), and ATP takes longer to run out than
    without regulation (below 0.1 mM at ~145 s instead of ~100 s; removing the pumps' copies costs ATP too).'''
    outs = ['mod/cytosol/C_ATP', 'mod/reg/f_ATP_reg', 'mod/pool_H_ATPase/m_goal', 'mod/bio/R_reg', 'mod/cytosol/R_reg']
    base = {**_glucose(0.0), 'growth_level_mod__bio': 1.0}
    reg = {**base, 'reg_on_ATP_mod__bio': 1.0, **{f'reg_on_ATP_mod__pool_{p}': 1.0 for p in REGULATED}, **_maintenance_on()}
    t, on = _phases(fermenter, [(300.0, reg)], outs, dt=1.0)[0]
    _, off = _phases(fermenter, [(300.0, base)], outs, dt=1.0)[0]
    atp = on['mod/cytosol/C_ATP']
    assert np.all(np.diff(on['mod/reg/f_ATP_reg']) * np.diff(atp) >= -1e-15)
    for o in ('mod/pool_H_ATPase/m_goal', 'mod/bio/R_reg', 'mod/cytosol/R_reg'):
        np.testing.assert_allclose(on[o], on['mod/reg/f_ATP_reg'], rtol=1e-12, err_msg=o)
    assert on['mod/reg/f_ATP_reg'][-1] < 0.5 * on['mod/reg/f_ATP_reg'][0]

    def gone(a):
        return t[np.argmax(a < 0.1)]
    assert gone(atp) > 1.3 * gone(off['mod/cytosol/C_ATP'])


def test_maintenance_regulation_avoids_the_ATP_trap_on_half_food(fermenter):
    '''On half the resting medium glucose (0.05 mM) for an hour, the unregulated cell runs its ATP out and cannot use
    full food when it returns: hexokinase needs ATP, so ATP stays near 0. Regulating the pumps, exchangers, flagellum
    and hydrolase by ATP helps: their copies are removed within minutes (the H+-ATPase to about a third in the hour) and
    the membrane depolarises, but ATP holds above 1 mM and recovers on full food. With the maintenance regulated by ATP
    (protein synthesis, RNA and lipid upkeep, repair, the maintenance ATPase) ATP holds higher, above 4 mM, and the
    membrane near -97 mV, while the proteome loses the copies its cut synthesis does not replace (about 0.3% in the
    hour, at k_deg) to the recycled pool, and the cuts cost envelope and ribosome damage and some ribosomes. When the
    resting medium returns the cell rebuilds its proteins, mostly from those amino acids, and repairs the damage, which
    takes most of the medium's small headroom: ATP climbs back (above 6 mM in the hour, above 8.3 mM in four) and V to
    -104 mV, while the envelope damage falls.'''
    outs = ['mod/cytosol/C_ATP', 'mod/membrane/V', 'mod/pool_H_ATPase/m', 'mod/pool_proteome/m']
    expression = {f'reg_on_ATP_mod__pool_{p}': 1.0 for p in REGULATED}
    fed = _glucose(fermenter.nominal['C_glc_e_mod__membrane_face'])
    (_, off_starve), (_, off_fed) = _phases(fermenter, [(3600.0, _glucose(0.05)), (3600.0, fed)], outs)
    (_, ex_starve), (_, ex_fed) = _phases(fermenter, [(3600.0, {**_glucose(0.05), **expression}), (3600.0, fed)], outs)
    outs += ['mod/bio/D_phys']
    (_, on_starve), (_, on_fed), (_, on_later) = _phases(
        fermenter, [(3600.0, {**_glucose(0.05), **_maintenance_on()}), (3600.0, fed), (3 * 3600.0, fed)], outs)
    assert off_starve['mod/cytosol/C_ATP'][-1] < 1e-3 and off_fed['mod/cytosol/C_ATP'][-1] < 1e-3
    assert ex_starve['mod/cytosol/C_ATP'].min() > 1.0 and ex_fed['mod/cytosol/C_ATP'][-1] > 9.0
    assert ex_starve['mod/pool_H_ATPase/m'][30] < 0.8 and ex_starve['mod/pool_H_ATPase/m'][-1] < 0.55
    assert on_starve['mod/cytosol/C_ATP'].min() > 4.0 and on_starve['mod/membrane/V'][-1] < -95.0
    assert on_starve['mod/cytosol/C_ATP'].min() > ex_starve['mod/cytosol/C_ATP'].min(), \
        'maintenance regulation helps more'
    assert on_fed['mod/cytosol/C_ATP'][-1] > 6.0 and on_later['mod/cytosol/C_ATP'][-1] > 8.3
    assert on_later['mod/membrane/V'][-1] < -104.0
    assert on_later['mod/bio/D_phys'][-1] < 0.6 * on_fed['mod/bio/D_phys'][-1]
    lost = 1.0 - on_starve['mod/pool_proteome/m'][-1]
    assert 0.002 < lost < 0.007
    assert on_fed['mod/pool_proteome/m'][-1] > on_starve['mod/pool_proteome/m'][-1]


def test_removing_a_large_protein_is_limited_by_the_proteases(fermenter):
    '''Halving the glycolytic enzymes (a sixth of a billion amino acids to remove) asks the proteases for more than their
    capacity (proteolysis_capacity, 5%/h of the cell's protein): they remove Q_prot_max amino acids per second, so the
    enzymes fall linearly, to about 0.56 in an hour, while ATP stays above 8 mM.'''
    nom = fermenter.nominal
    outs = ['mod/pool_Glyc/m', 'mod/protein/s_prot', 'mod/protein/Q_prot_demand_total', 'mod/cytosol/C_ATP']
    (t, r), = _phases(fermenter, [(3600.0, {'m_target_mod__pool_Glyc': 0.5})], outs, dt=60.0)
    served = r['mod/protein/s_prot'] * r['mod/protein/Q_prot_demand_total']
    np.testing.assert_allclose(served, nom['Q_prot_max_mod__protein'], rtol=1e-9)
    m = r['mod/pool_Glyc/m']
    assert 0.53 < m[-1] < 0.6 and np.all(np.diff(m) < 0)
    assert r['mod/cytosol/C_ATP'].min() > 8.0


def test_a_protein_rises_with_its_mRNA_lifetime_and_within_the_ribosomes_capacity(fermenter):
    '''On food (10 mM glucose), doubling the H+-ATPase (its share of the membrane): a small protein rises with time constant tau = t_half_mRNA / ln 2
    (about 7 min) since its synthesis is well within the ribosomes' capacity: m(tau) = 2 - exp(-1) within 2%. Raising the
    proteome by a fifth asks for more than the ribosomes can make: synthesis spends at most J_cap = build_capacity J_ref
    (about 9 times the resting turnover) and the rise takes hours.'''
    tau = fermenter.nominal['tau_mod__pool_H_ATPase']
    assert tau == pytest.approx(300.0 / np.log(2), rel=1e-9)
    outs = ['mod/pool_H_ATPase/m', 'mod/protein/J_build_total', 'mod/protein/J_cap', 'mod/protein/s_build']
    double = {**_glucose(10.0), 's_target_mod__pool_H_ATPase': 2.0 * fermenter.nominal['s_target_mod__pool_H_ATPase']}
    (t, r), = _phases(fermenter, [(tau, double)], outs, dt=tau / 50)
    assert r['mod/pool_H_ATPase/m'][-1] == pytest.approx(2.0 - np.exp(-1.0), rel=0.02)
    assert np.all(r['mod/protein/s_build'] > 0.999)
    outs = ['mod/pool_proteome/m', 'mod/protein/J_build_total', 'mod/protein/J_cap', 'mod/protein/s_build']
    (t, r), = _phases(fermenter, [(3600.0, {**_glucose(10.0), 'm_target_mod__pool_proteome': 1.2})], outs, dt=60.0)
    assert r['mod/protein/s_build'][1] < 0.5
    assert np.all(r['mod/protein/J_build_total'] <= r['mod/protein/J_cap'] * (1 + 1e-9) + 1e-30)
    assert 1.01 < r['mod/pool_proteome/m'][-1] < 1.15


def test_each_maintenance_cut_costs_damage_or_ribosomes(fermenter):
    '''Three hours on half food (0.05 mM) with the maintenance regulated by ATP, then three hours on food (10 mM): with
    every switch on, the unpaid upkeep damages the envelope (D_phys above 0.02) and loses ribosomes (x_ribo below 0.995),
    and the cut maintenance ATPase damages ribosomes and proteins (D_ribo up by more than 0.025). Leaving the upkeep's
    switches off removes the envelope damage and the ribosome loss; leaving the maintenance ATPase's off removes the
    extra ribosome damage. Protein synthesis can build build_capacity_rest x_ribo (exact). When food returns the
    ribosomes are rebuilt (above 0.998 of rest in three hours, at J_ribo_build > 0), the damage is repaired and ATP
    returns to rest.'''
    nom = fermenter.nominal
    outs = ['mod/cytosol/C_ATP', 'mod/ribosomes/x_ribo', 'mod/ribosomes/J_ribo_build', 'mod/protein/build_capacity',
            'mod/bio/D_phys', 'mod/bio/D_ribo']
    on = _maintenance_on()
    cases = {'all': on, 'upkeep_off': {k: v for k, v in on.items() if '_upkeep_' not in k},
             'maint_off': {k: v for k, v in on.items() if not k.endswith('_mod__cytosol')}}
    runs = {name: _phases(fermenter, [(3 * 3600.0, {**_glucose(0.05), **sw}), (3 * 3600.0, {**_glucose(10.0), **sw})],
                          outs, dt=60.0) for name, sw in cases.items()}
    (_, starve), (_, fed) = runs['all']
    assert starve['mod/bio/D_phys'][-1] > 0.02 and starve['mod/ribosomes/x_ribo'][-1] < 0.995
    d_ribo_rest = nom['D_ribo_init_mod__bio']
    assert starve['mod/bio/D_ribo'][-1] - d_ribo_rest > 0.025
    assert runs['upkeep_off'][0][1]['mod/bio/D_phys'][-1] < 0.005
    assert runs['upkeep_off'][0][1]['mod/ribosomes/x_ribo'][-1] > 0.9995
    assert runs['maint_off'][0][1]['mod/bio/D_ribo'][-1] - d_ribo_rest < 0.015
    for r in (starve, fed):
        np.testing.assert_allclose(r['mod/protein/build_capacity'],
                                   nom['build_capacity_rest_mod__protein'] * r['mod/ribosomes/x_ribo'], rtol=1e-12)
    assert np.all(np.diff(fed['mod/ribosomes/x_ribo']) >= 0) and fed['mod/ribosomes/x_ribo'][-1] > 0.998
    assert fed['mod/ribosomes/J_ribo_build'][1:60].min() > 0
    assert fed['mod/bio/D_phys'][-1] < 0.6 * starve['mod/bio/D_phys'][-1]
    assert fed['mod/bio/D_ribo'][-1] < 0.7 * starve['mod/bio/D_ribo'][-1]
    assert fed['mod/cytosol/C_ATP'][-1] > 9.6


def test_nothing_costs_anything_at_rest(fermenter):
    '''At rest with every maintenance switch on (f = 1): the upkeep is paid in full, no ribosome is lost or rebuilt, and
    neither new damage term moves, for an hour (all at round-off: 1e-12 of their full rates).'''
    nom = fermenter.nominal
    outs = ['mod/bio/x_unpaid', 'mod/ribosomes/R_ribo_loss', 'mod/ribosomes/R_ribo_build', 'mod/ribosomes/x_ribo',
            'mod/bio/R_dmg_lipid', 'mod/bio/R_dmg_qc']
    _, r = _phases(fermenter, [(3600.0, _maintenance_on())], outs, dt=60.0)[0]
    full = {'mod/bio/x_unpaid': 1.0,
            'mod/ribosomes/R_ribo_loss': nom['k_ribo_mod__ribosomes'] * nom['n_ribo_rest_mod__ribosomes'],
            'mod/ribosomes/R_ribo_build': nom['k_ribo_build_mod__ribosomes'] * nom['n_ribo_rest_mod__ribosomes'],
            'mod/bio/R_dmg_lipid': nom['k_lipid_dmg_mod__bio'], 'mod/bio/R_dmg_qc': nom['k_qc_dmg_mod__bio']}
    for o, rate in full.items():
        assert np.all(np.abs(r[o]) <= 1e-12 * rate), o
    np.testing.assert_allclose(r['mod/ribosomes/x_ribo'], 1.0, rtol=0, atol=1e-12)
