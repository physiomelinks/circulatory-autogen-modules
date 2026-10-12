"""
Regulation of protein expression and growth by falling metabolites (cell/expression_regulation ArgusUNPUBLISHED_v01) in
AdaptOrDeity's fermenter (cell/generic_adaptive_cell fermenter): every protein_pool's target and growth are multiplied by
R_reg = prod_s (1 - reg_on_s + reg_on_s f_s), f_s = x^n/(1 + x^n) with x = C/K for ATP, G6P and glucose inside and x =
K_H/H_i for H+; the switches reg_on_s are 0 by default.

    pytest tests/test_expression_regulation.py --include-unreviewed

Checks: the switches off leave the cell exactly as it was; a switch makes a target follow its species; and a starving
cell whose transporters and pumps are regulated by ATP holds its ATP and recovers when food returns, where the same cell
without regulation runs its ATP out and stays there (hexokinase needs ATP to start glycolysis again: the ATP trap).
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
    '''By default every switch is 0: each pool's R_reg and growth's are exactly 1 and every pool goes to its target
    (bitwise), while the signals themselves are at their resting value 16/17 (x = 2: rest is twice K).'''
    pools = [k[len('reg_on_ATP_mod_'):] for k in fermenter.nominal if k.startswith('reg_on_ATP_mod_pool_')]
    assert len(pools) == 14
    assert all(fermenter.nominal[f'reg_on_{s}_mod_{p}'] == 0.0 for s in SPECIES for p in pools + ['bio'])
    outs = ([f'mod_{p}/R_reg' for p in pools] + [f'mod_{p}/m_goal' for p in pools] + [f'mod_{p}/m_target' for p in pools]
            + ['mod_bio/R_reg'] + [f'mod_reg/f_{s}_reg' for s in SPECIES])
    _, r = _phases(fermenter, [(60.0, {})], outs, dt=1.0)[0]
    for p in pools + ['bio']:
        assert np.all(r[f'mod_{p}/R_reg'] == 1.0), p
    for p in pools:
        assert np.all(r[f'mod_{p}/m_goal'] == r[f'mod_{p}/m_target']), p   # m_target may be logged once
    for s in SPECIES:
        assert r[f'mod_reg/f_{s}_reg'][0] == pytest.approx(16 / 17, rel=1e-9), s


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
    '''Starving (no glucose), with ATP regulation on the pumps and growth: as ATP falls the H+-ATPase's goal and growth
    fall with it (f_ATP is monotonic in ATP), and ATP takes longer to run out than without regulation.'''
    outs = ['mod_cytosol/C_ATP', 'mod_reg/f_ATP_reg', 'mod_pool_H_ATPase/m_goal', 'mod_bio/R_reg']
    base = {**_glucose(0.0), 'growth_level_mod_bio': 1.0}
    reg = {**base, 'reg_on_ATP_mod_bio': 1.0, **{f'reg_on_ATP_mod_pool_{p}': 1.0 for p in REGULATED}}
    t, on = _phases(fermenter, [(300.0, reg)], outs, dt=1.0)[0]
    _, off = _phases(fermenter, [(300.0, base)], outs, dt=1.0)[0]
    atp = on['mod_cytosol/C_ATP']
    assert np.all(np.diff(on['mod_reg/f_ATP_reg']) * np.diff(atp) >= -1e-15)
    np.testing.assert_allclose(on['mod_pool_H_ATPase/m_goal'], on['mod_reg/f_ATP_reg'], rtol=1e-12)
    np.testing.assert_allclose(on['mod_bio/R_reg'], on['mod_reg/f_ATP_reg'], rtol=1e-12)
    assert on['mod_reg/f_ATP_reg'][-1] < 0.5 * on['mod_reg/f_ATP_reg'][0]

    def gone(a):
        return t[np.argmax(a < 0.1)]
    assert gone(atp) > gone(off['mod_cytosol/C_ATP'])


def test_regulation_avoids_the_ATP_trap_on_half_food(fermenter):
    '''On half the resting medium glucose (0.05 mM) for an hour, the unregulated cell runs its ATP out (below 0.1 mM in
    about 11 minutes) and cannot use full food when it returns: hexokinase needs ATP, so ATP stays near 0. With ATP
    regulation on its pumps, exchangers, flagellum and hydrolase, those fall to about a third, ATP holds above 4 mM and
    the membrane near -66 mV, and within an hour of full food the cell is back at rest (ATP > 9 mM, V < -100 mV; it
    rebuilds its pumps as fast as the build throttle allows, so not within ten minutes).'''
    outs = ['mod_cytosol/C_ATP', 'mod_membrane/V', 'mod_pool_H_ATPase/m']
    reg = {f'reg_on_ATP_mod_pool_{p}': 1.0 for p in REGULATED}
    fed = _glucose(fermenter.nominal['C_glc_e_mod_membrane_face'])
    (_, off_starve), (_, off_fed) = _phases(fermenter, [(3600.0, _glucose(0.05)), (3600.0, fed)], outs)
    (_, on_starve), (_, on_fed) = _phases(fermenter, [(3600.0, {**_glucose(0.05), **reg}), (3600.0, fed)], outs)
    assert off_starve['mod_cytosol/C_ATP'][-1] < 1e-3 and off_fed['mod_cytosol/C_ATP'][-1] < 1e-3
    assert on_starve['mod_cytosol/C_ATP'].min() > 4.0
    assert on_fed['mod_cytosol/C_ATP'][-1] > 9.0 and on_fed['mod_membrane/V'][-1] < -100.0
    assert on_starve['mod_pool_H_ATPase/m'][-1] < 0.5 < on_fed['mod_pool_H_ATPase/m'][-1]
