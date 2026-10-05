"""
neuron sympathetic: why the somatic AP does not reach the varicosity, and which axon changes fix it.

    venv/bin/python tools/neuron_ap_transmission.py            # -> reviews/neuron_sympathetic_ap_transmission.html
    venv/bin/python tools/neuron_ap_transmission.py --workdir /tmp/apt_work

Simulation only: the committed module files are not changed. The neuron supermodule
(neuron/sympathetic: soma -> axon sympathetic_monolithic_v01 -> varicosity) and the isolated soma
(soma/sympathetic) are generated with the repo harness into a scratch directory. Parameter options
(L_ax, R_i, d_ax, channel_ratio) are runtime parameter changes of that one model. The HH-axon
options are prototypes: a copy of the generated model whose SN_axon component gets a Hodgkin-Huxley
Na/K/leak membrane (HH 1952 rate functions in u = V - V_rest, E_Na +50 mV, E_K V_rest - 12 mV, E_L
chosen for zero net current at V_rest). The multi-compartment estimate is a numpy RC ladder driven
by the recorded somatic voltage (one-way: it ignores the ladder's load on the soma).

Every run passes the full set of option parameters, because a myokit helper keeps parameter values
set by a previous run.
"""
import argparse
import base64
import html
import io
import math
import os
import re
import shutil
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from cam_testing import library, harness  # noqa: E402

OUT = os.path.join(REPO, 'reviews', 'neuron_sympathetic_ap_transmission.html')
CURRENTS = [0.02, 0.05, 0.1, 0.2]          # nA
SIM_TIME, DT = 1.0, 1e-5
L0, RI0, D0, CR0 = 318.3098861837907, 74.02203300817018, 1.0, 0.5
C_SOMA = 0.01 * 5.086427133679043 * 24.285900630052804 ** 2   # pF, soma Cm = c_m k_A d^2
DEFAULTS = {'L_ax_mod_axon': L0, 'R_i_mod_axon': RI0, 'd_ax_mod_axon': D0, 'channel_ratio_mod_axon': CR0}
OUTPUTS = ['mod_soma_membrane/V', 'mod_axon/V', 'mod_varicosity_membrane/V', 'mod_axon/I',
           'mod_varicosity_Ca/Cai', 'mod_varicosity_NE/NE']
COL = {'soma': '#2a78d6', 'axon': '#eb6834', 'var': '#1baf7a', 'iso': '#52514e'}


# ---------------------------------------------------------------- HH axon prototype (scratch copy only)
def _ci(n): return f'<ci>{n}</ci>'
def _cn(v, u='dimensionless'): return f'<cn cellml:units="{u}">{float(v)!r}</cn>'
def _ap(op, *a): return f'<apply><{op}/>' + ''.join(a) + '</apply>'
def _eq(lhs, rhs): return _ap('eq', lhs, rhs)
def _ode(x, rhs): return _ap('eq', _ap('diff', '<bvar><ci>t</ci></bvar>', _ci(x)), rhs)


def hh_axon_component(gNa, gK, gL, phi, Vr=-70.0, ENa=50.0):
    '''SN_axon with an HH Na/K/leak membrane; conductances in microS for the whole compartment.'''
    EK = Vr - 12.0
    am = lambda u: 0.1 * (25 - u) / (math.exp((25 - u) / 10) - 1); bm = lambda u: 4 * math.exp(-u / 18)
    ah = lambda u: 0.07 * math.exp(-u / 20); bh = lambda u: 1 / (math.exp((30 - u) / 10) + 1)
    an = lambda u: 0.01 * (10 - u) / (math.exp((10 - u) / 10) - 1); bn = lambda u: 0.125 * math.exp(-u / 80)
    m0, h0, n0 = am(0) / (am(0) + bm(0)), ah(0) / (ah(0) + bh(0)), an(0) / (an(0) + bn(0))
    EL = Vr + (gNa * m0 ** 3 * h0 * (Vr - ENa) + gK * n0 ** 4 * (Vr - EK)) / gL
    U = _ci('u_hh')
    ex = lambda a: _ap('exp', a)

    def lin(a, b, c):
        d = _ap('minus', _cn(b), U)
        return _ap('divide', _ap('times', _cn(a), d), _ap('minus', ex(_ap('divide', d, _cn(c))), _cn(1)))

    def dec(a, c): return _ap('times', _cn(a), ex(_ap('divide', _ap('minus', U), _cn(c))))

    def gate(x, a, b):
        return _ode(x, _ap('times', _cn(1000.0 * phi, 'per_s'),
                           _ap('minus', _ap('times', a, _ap('minus', _cn(1), _ci(x))), _ap('times', b, _ci(x)))))
    beta_h = _ap('divide', _cn(1), _ap('plus', ex(_ap('divide', _ap('minus', _cn(30), U), _cn(10))), _cn(1)))
    I_ion = _ap('plus',
                _ap('times', _cn(gNa, 'microS'), _ap('power', _ci('m_hh'), _cn(3)), _ci('h_hh'),
                    _ap('minus', _ci('V'), _cn(ENa, 'milliV'))),
                _ap('times', _cn(gK, 'microS'), _ap('power', _ci('n_hh'), _cn(4)), _ap('minus', _ci('V'), _cn(EK, 'milliV'))),
                _ap('times', _cn(gL, 'microS'), _ap('minus', _ci('V'), _cn(EL, 'milliV'))))
    dV = _ap('divide', _ap('times', '<cn cellml:units="milliV_per_kiloV" type="e-notation">1<sep/>6</cn>',
                           _ap('minus', _ci('I'), _ci('I_ion_hh'))), _ci('C'))
    maths = ('<math xmlns="http://www.w3.org/1998/Math/MathML">'
             + _eq(_ci('I'), _ap('divide', _ap('minus', _ci('V_in'), _ci('V')), _ci('R')))
             + _eq(_ci('I_channels'), _ap('times', _ci('I'), _ci('channel_ratio')))    # kept as an output, unused
             + _eq(U, _ap('divide', _ap('minus', _ci('V'), _cn(Vr, 'milliV')), _cn(1, 'milliV')))
             + _eq(_ci('I_ion_hh'), I_ion)
             + gate('m_hh', lin(0.1, 25, 10), dec(4, 18)) + gate('h_hh', dec(0.07, 20), beta_h)
             + gate('n_hh', lin(0.01, 10, 10), dec(0.125, 80)) + _ode('V', dV)
             + _eq(_ci('C'), _ap('times', _ci('c_m'), '<pi/>', _ci('d_ax'), _ci('L_ax'), _cn(0.01, 'picoF_cm2_per_uF_um2')))
             + _eq(_ci('R'), _ap('divide', _ap('times', _cn(4), _ci('R_i'), _ci('L_ax'), _cn(0.01, 'megaOhm_um_per_ohm_cm')),
                                 _ap('times', '<pi/>', _ap('power', _ci('d_ax'), _cn(2)))))
             + '</math>')
    vs = [('t', 'second', 'in'), ('C', 'picoF', 'out'), ('R', 'megaOhm', 'out'), ('c_m', 'uF_per_cm2', 'in'),
          ('d_ax', 'um', 'in'), ('L_ax', 'um', 'in'), ('R_i', 'ohm_cm', 'in'), ('V_in', 'milliV', 'in'),
          ('I', 'nanoA', 'out'), ('V_rest_approx', 'milliV', 'in'), ('channel_ratio', 'dimensionless', 'in'),
          ('I_channels', 'nanoA', 'out'), ('I_ion_hh', 'nanoA', 'out')]
    lines = [f'<variable name="{n}" public_interface="{p}" units="{u}"/>' for n, u, p in vs]
    lines += ['<variable initial_value="V_rest_approx" name="V" public_interface="out" units="milliV"/>',
              '<variable name="u_hh" units="dimensionless"/>',
              f'<variable name="m_hh" units="dimensionless" initial_value="{m0!r}"/>',
              f'<variable name="h_hh" units="dimensionless" initial_value="{h0!r}"/>',
              f'<variable name="n_hh" units="dimensionless" initial_value="{n0!r}"/>']
    return '<component name="SN_axon">\n' + '\n'.join(lines) + '\n' + maths + '\n</component>', EL


def hh_variant(base_dir, dst_dir, **kw):
    if os.path.exists(dst_dir):
        shutil.rmtree(dst_dir)
    shutil.copytree(base_dir, dst_dir)
    name = os.path.basename(base_dir.rstrip('/'))
    path = os.path.join(dst_dir, f'{name}_modules.cellml')
    comp, EL = hh_axon_component(**kw)
    with open(path) as f:
        txt = f.read()
    txt, n = re.subn(r'<component name="SN_axon">.*?</component>', lambda m: comp, txt, count=1, flags=re.S)
    assert n == 1
    with open(path, 'w') as f:
        f.write(txt)
    return os.path.join(dst_dir, f'{name}.cellml'), EL


# ---------------------------------------------------------------- metrics
def up_crossings(t, x, level=0.0):
    i = np.where((x[1:] >= level) & (x[:-1] < level))[0]
    return t[i] + (level - x[i]) / (x[i + 1] - x[i]) * (t[i + 1] - t[i])


def half_width(t, x):
    k = int(np.argmax(x)); base = np.min(x[max(0, k - 2000):k + 1]); half = base + 0.5 * (x[k] - base)
    a = k
    while a > 0 and x[a] > half: a -= 1
    b = k
    while b < len(x) - 1 and x[b] > half: b += 1
    return t[b] - t[a]


def metrics(t, r):
    s, a, v = r['mod_soma_membrane__V'], r['mod_axon__V'], r['mod_varicosity_membrane__V']
    cs, ca, cv = up_crossings(t, s), up_crossings(t, a), up_crossings(t, v)
    m = dict(n_soma=len(cs), n_axon=len(ca), n_var=len(cv), pk_soma=s.max(), pk_axon=a.max(), pk_var=v.max(),
             lat_soma=cs[0] * 1e3 if len(cs) else None,
             delay=(cv[0] - cs[0]) * 1e3 if len(cs) and len(cv) else None,
             I_ax=np.max(np.abs(r['mod_axon__I'])), Cai=r['mod_varicosity_Ca__Cai'].max() * 1e3,
             NE=r['mod_varicosity_NE__NE'].max() * 1e3)
    m['ok'] = m['n_soma'] > 0 and m['n_axon'] == m['n_soma'] and abs(m['pk_axon'] - m['pk_soma']) <= 20
    return m


# ---------------------------------------------------------------- multi-compartment estimate
def ladder_far_end(t, vs, N, R, C, cr=0.0, Vr=-70.0):
    '''Passive RC ladder (N nodes, R/N each, C/N each, sealed end) driven by vs(t); far-end V.'''
    Rk, Ck = R / N * 1e6, C / N * 1e-12 / (1 + cr)
    f = lambda tt: np.interp(tt, t, vs)

    def rhs(tt, y):
        left = np.concatenate(([f(tt)], y[:-1])); right = np.concatenate((y[1:], [y[-1]]))
        return ((left - y) + (right - y)) / Rk / Ck
    sol = solve_ivp(rhs, (t[0], t[-1]), np.full(N, Vr), t_eval=t, method='LSODA', max_step=2e-5, rtol=1e-7, atol=1e-6)
    return sol.y[-1]


# ---------------------------------------------------------------- figures
def png(fig):
    buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=110, bbox_inches='tight', facecolor='white'); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def style(ax):
    ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=0.25, lw=0.6)


def trace_fig(opt, I, t, r, t_iso, v_iso):
    s, a, v = r['mod_soma_membrane__V'], r['mod_axon__V'], r['mod_varicosity_membrane__V']
    cs = up_crossings(t, s)
    t0 = cs[0] if len(cs) else t[np.argmax(s)]
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 2.9), gridspec_kw={'width_ratios': [1.1, 1]})
    for ax, (lo, hi), unit in [(axs[0], (0, SIM_TIME), 1), (axs[1], (t0 - 0.004, t0 + 0.012), 1e3)]:
        k = (t >= lo) & (t <= hi)
        if ax is axs[0]:
            ki = (t_iso >= lo) & (t_iso <= hi)
            ax.plot(t_iso[ki], v_iso[ki], color=COL['iso'], lw=1, alpha=0.55, label='isolated soma')
        ax.plot(t[k] * unit, s[k], color=COL['soma'], lw=1.6, label='soma')
        ax.plot(t[k] * unit, a[k], color=COL['axon'], lw=1.6, label='axon')
        ax.plot(t[k] * unit, v[k], color=COL['var'], lw=1.4, ls=(0, (3, 2)), label='varicosity')
        ax.axhline(0, color='#999', lw=0.6); style(ax); ax.set_ylabel('V (mV)')
    axs[0].set_xlabel('t (s)'); axs[1].set_xlabel('t (ms), first spike')
    axs[0].set_title(f'{opt["label"]}, {I * 1e3:.0f} pA', fontsize=10, loc='left')
    axs[0].legend(fontsize=7.5, frameon=False, loc='upper right', ncol=2)
    fig.tight_layout()
    return png(fig)


# ---------------------------------------------------------------- main
def options():
    A = math.pi * D0 * L0          # um2, axon membrane area (1000 um2)
    dens = lambda ms_cm2: ms_cm2 * 1e-5 * A   # mS/cm2 -> microS for the compartment
    return [
        dict(key='current', label='Current model', change='none (R 300 MOhm, C 10 pF, channel_ratio 0.5)', params={}),
        dict(key='L3', label='L_ax / 3', change='L_ax 318 -> 106 um (R 100 MOhm, C 3.3 pF)', params={'L_ax_mod_axon': L0 / 3}),
        dict(key='L4', label='L_ax / 4', change='L_ax 318 -> 80 um (R 75 MOhm, C 2.5 pF)', params={'L_ax_mod_axon': L0 / 4}),
        dict(key='Ri10', label='R_i / 10', change='R_i 74 -> 7.4 ohm cm (R 30 MOhm, C 10 pF)', params={'R_i_mod_axon': RI0 / 10}),
        dict(key='d5', label='d_ax x 5', change='d_ax 1 -> 5 um (R 12 MOhm, C 50 pF)', params={'d_ax_mod_axon': 5.0}),
        dict(key='cr10', label='channel_ratio 10', change='channel_ratio 0.5 -> 10', params={'channel_ratio_mod_axon': 10.0}),
        dict(key='crneg', label='channel_ratio -0.5', change='sign check: channel_ratio 0.5 -> -0.5', params={'channel_ratio_mod_axon': -0.5}),
        dict(key='hh', label='HH axon (2x Na, 6.3 degC kinetics)',
             change='prototype: HH Na/K/leak in the axon (gNa 240, gK 36, gL 0.3 mS/cm2, phi 1), channel_ratio unused',
             hh=dict(gNa=dens(240), gK=dens(36), gL=dens(0.3), phi=1.0)),
        dict(key='hh22', label='HH axon (3x Na, 22 degC kinetics)',
             change='prototype: HH Na/K/leak (gNa 360, gK 36, gL 0.3 mS/cm2, phi = 3^((22-6.3)/10) = 5.6)',
             hh=dict(gNa=dens(360), gK=dens(36), gL=dens(0.3), phi=3 ** ((22 - 6.3) / 10))),
    ]


def fmt(x, nd=1):
    return '&ndash;' if x is None else f'{x:.{nd}f}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workdir', default=None)
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    work = args.workdir or tempfile.mkdtemp(prefix='neuron_apt_')
    os.makedirs(work, exist_ok=True)

    neuron = library.load_version('neuron', 'sympathetic')
    soma = library.load_version('soma', 'sympathetic')
    spec = dict(neuron.spec, dt=DT, sim_time=SIM_TIME)
    base_path = harness.generate(neuron, os.path.join(work, 'neuron'))
    base_dir = os.path.dirname(base_path)
    soma_path = harness.generate(soma, os.path.join(work, 'soma'))

    # isolated soma
    hs = harness.simulation_helper(soma_path, dict(soma.spec, dt=DT, sim_time=SIM_TIME))
    iso = {}
    for I in CURRENTS:
        t, r = harness.run(hs, ['mod_membrane/V'], params={'parameters/I_in_mod_membrane': I})
        x = r['mod_membrane__V']; c = up_crossings(t, x)
        iso[I] = dict(t=t, v=x, n=len(c), lat=c[0] * 1e3 if len(c) else None, pk=x.max())
        print('isolated soma', I, iso[I]['n'], flush=True)

    h_base = harness.simulation_helper(base_path, spec)
    opts = options()
    results = {}
    for o in opts:
        if 'hh' in o:
            p, o['EL'] = hv = hh_variant(base_dir, os.path.join(work, f'hh_{o["key"]}'), **o['hh'])
            h = harness.simulation_helper(p, spec)
        else:
            h = h_base
        for I in CURRENTS:
            params = {f'parameters/{k}': v for k, v in {**DEFAULTS, **o.get('params', {})}.items()}
            params['parameters/I_in_mod_soma_membrane'] = I
            t, r = harness.run(h, OUTPUTS, params=params)
            results[o['key'], I] = (t, r, metrics(t, r))
            m = results[o['key'], I][2]
            print(o['key'], I, m['n_soma'], m['n_axon'], round(m['pk_soma'], 1), round(m['pk_axon'], 1), m['ok'], flush=True)

    # tau_eff sweep at 50 pA (the single-compartment RC family)
    sweep = {'channel_ratio': [], 'R_i': [], 'L_ax': [], 'd_ax': []}
    for fam, vals in [('channel_ratio', [0, 0.5, 1, 2, 3, 5, 7, 10, 15, 20]),
                      ('R_i', [RI0 / f for f in (1, 2, 3, 5, 7, 10, 15, 20)]),
                      ('L_ax', [L0 / f for f in (1, 1.5, 2, 2.5, 3, 4, 5, 10)]),
                      ('d_ax', [1, 2, 3, 5, 8])]:
        for val in vals:
            pv = dict(DEFAULTS, **{f'{fam}_mod_axon': float(val)})
            Rm = 4 * pv['R_i_mod_axon'] * pv['L_ax_mod_axon'] * 0.01 / (math.pi * pv['d_ax_mod_axon'] ** 2)
            Cp = math.pi * pv['d_ax_mod_axon'] * pv['L_ax_mod_axon'] * 0.01
            tau = Rm * Cp * 1e-3 / (1 + pv['channel_ratio_mod_axon'])     # ms
            params = {f'parameters/{k}': v for k, v in pv.items()}
            params['parameters/I_in_mod_soma_membrane'] = 0.05
            t, r = harness.run(h_base, OUTPUTS, params=params)
            m = metrics(t, r)
            sweep[fam].append((tau, m['pk_axon'] - m['pk_soma'], m['pk_soma']))

    # AP width and multi-compartment estimate, from the current model at 50 pA
    t, r, _ = results['current', 0.05]
    vs = r['mod_soma_membrane__V']
    hw = half_width(t, vs) * 1e3
    c0 = up_crossings(t, vs)[0]; d0 = up_crossings(t, -vs)[0]
    w0 = (d0 - c0) * 1e3          # ms above 0 mV
    w25 = (up_crossings(t, -vs, -25.0)[0] - up_crossings(t, vs, 25.0)[0]) * 1e3   # ms above +25 mV
    k0 = int(np.argmax(vs)); win = slice(max(0, k0 - 1500), k0 + 2500)
    tw, vw = t[win], vs[win]
    R0, C0 = 300.0, 10.0
    ladder = []
    for N, cr in [(1, 0.5), (1, 0.0), (4, 0.0), (16, 0.0), (16, 0.5)]:
        far = ladder_far_end(tw, vw, N, R0, C0, cr)
        ladder.append((N, cr, far.max(), tw, far))
    ladder3 = [(N, ladder_far_end(tw, vw, N, R0 / 3, C0 / 3, 0.5).max()) for N in (1, 16)]

    # ------------------------------------------------ figures
    figs = {}
    for o in opts:
        figs[o['key']] = trace_fig(o, 0.05, *results[o['key'], 0.05][:2], iso[0.05]['t'], iso[0.05]['v'])
    figs['current_200'] = trace_fig(opts[0], 0.2, *results['current', 0.2][:2], iso[0.2]['t'], iso[0.2]['v'])
    for key in ('L3', 'hh'):
        figs[key + '_200'] = trace_fig(next(o for o in opts if o['key'] == key), 0.2, *results[key, 0.2][:2],
                                       iso[0.2]['t'], iso[0.2]['v'])

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for (fam, pts), c, mk in zip(sweep.items(), ['#2a78d6', '#eb6834', '#1baf7a', '#4a3aa7'], 'osD^'):
        pts = sorted(pts)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, ms=5, lw=1.5, color=c, label=f'vary {fam}')
    ax.axhline(-20, color='#e34948', lw=1, ls='--'); ax.text(2.6, -18, 'check: within 20 mV', color='#b3261e', fontsize=8, ha='right')
    ax.axvline(w0, color='#999', lw=0.8, ls=':')
    ax.text(w0 * 0.94, -35, f'soma AP above 0 mV: {w0:.2f} ms', fontsize=8, color='#555', ha='right')
    ax.set_xscale('log'); ax.set_xlabel('tau_eff = R C / (1 + channel_ratio)  (ms)'); ax.set_ylabel('axon peak - soma peak (mV)')
    ax.set_title('Single compartment, 50 pA: attenuation is set by tau_eff alone', fontsize=10, loc='left')
    style(ax); ax.legend(fontsize=8, frameon=False, loc='lower left')
    figs['sweep'] = png(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot((tw - tw[0]) * 1e3, vw, color=COL['soma'], lw=1.8, label='soma (current model, 50 pA)')
    for (N, cr, pk, tt, far), c in zip(ladder, ['#eb6834', '#e87ba4', '#1baf7a', '#4a3aa7', '#008300']):
        ax.plot((tt - tt[0]) * 1e3, far, color=c, lw=1.4, label=f'N={N}, channel_ratio {cr}: peak {pk:.0f} mV')
    style(ax); ax.set_xlabel('t (ms)'); ax.set_ylabel('V (mV)'); ax.legend(fontsize=7.5, frameon=False)
    ax.set_title('Passive RC ladder, R 300 MOhm, C 10 pF in total, far end', fontsize=10, loc='left')
    figs['ladder'] = png(fig)

    # ------------------------------------------------ HTML
    rows = []
    for o in opts:
        for I in CURRENTS:
            m = results[o['key'], I][2]; s = iso[I]
            soma_vs_iso = f'{m["n_soma"]} vs {s["n"]}; {fmt(m["lat_soma"])} vs {fmt(s["lat"])} ms'
            rows.append(
                f'<tr class="{"pass" if m["ok"] else "fail"}"><td>{html.escape(o["label"]) if I == CURRENTS[0] else ""}</td>'
                f'<td>{I * 1e3:.0f}</td><td>{m["n_soma"]}/{m["n_axon"]}/{m["n_var"]}</td>'
                f'<td>{m["pk_soma"]:.1f}</td><td>{m["pk_axon"]:.1f}</td><td>{m["pk_var"]:.1f}</td>'
                f'<td>{fmt(m["delay"], 2)}</td><td>{soma_vs_iso}</td><td>{m["I_ax"]:.2f}</td>'
                f'<td>{m["Cai"]:.2f}</td><td>{m["NE"]:.3f}</td><td>{"pass" if m["ok"] else "fail"}</td></tr>')
    opt_rows = ''.join(f'<tr><td>{html.escape(o["label"])}</td><td>{html.escape(o["change"])}</td></tr>' for o in opts)
    ladder_rows = ''.join(f'<tr><td>{N}</td><td>{cr}</td><td>{pk:.1f}</td></tr>' for N, cr, pk, *_ in ladder)
    cur50 = results['current', 0.05][2]
    L3_50 = results['L3', 0.05][2]
    cr10_50 = results['cr10', 0.05][2]

    def card(key, cap):
        return f'<div class="card"><img alt="{html.escape(cap)}" src="data:image/png;base64,{figs[key]}"><p class="muted">{cap}</p></div>'
    traces = ''.join(card(o['key'], f'{html.escape(o["label"])}: {html.escape(o["change"])}. 50 pA.') for o in opts)
    traces200 = ''.join(card(k, c) for k, c in [('current_200', 'Current model, 200 pA'), ('L3_200', 'L_ax / 3, 200 pA'),
                                                ('hh_200', 'HH axon (2x Na, 6.3 degC kinetics), 200 pA')])
    hh_EL = next(o for o in opts if o['key'] == 'hh')['EL']

    page = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Sympathetic AP transmission</title>
<style>
:root {{ --bg:#ffffff; --fg:#111827; --muted:#4b5563; --card:#f9fafb; --line:#e5e7eb; --accent:#d9480f; --pass:#e7f5ec; --fail:#fdecec; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; --pass:#16301f; --fail:#3a1a1a; }} }}
:root[data-theme="dark"] {{ --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; --pass:#16301f; --fail:#3a1a1a; }}
body {{ background:var(--bg); color:var(--fg); font:15px/1.5 system-ui, -apple-system, Segoe UI, sans-serif; margin:0; padding:24px 16px; }}
main {{ max-width:1180px; margin:0 auto; }}
h1 {{ font-size:1.5rem; margin:0 0 4px; }} h2 {{ font-size:1.15rem; margin:28px 0 8px; border-bottom:1px solid var(--line); padding-bottom:4px; }}
.muted {{ color:var(--muted); font-size:13px; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit, minmax(min(100%, 520px), 1fr)); gap:16px; align-items:start; }}
.card {{ background:var(--card); border:1px solid var(--line); border-radius:8px; padding:12px; overflow-x:auto; }}
.card img {{ width:100%; height:auto; background:#fff; border-radius:4px; display:block; }}
.scroll {{ overflow-x:auto; }}
table {{ border-collapse:collapse; width:100%; font-size:13.5px; }} th, td {{ border-bottom:1px solid var(--line); padding:4px 8px; text-align:left; white-space:nowrap; }}
th {{ color:var(--muted); font-weight:600; }} tr.pass td:last-child {{ background:var(--pass); }} tr.fail td:last-child {{ background:var(--fail); }}
.rec {{ border-left:4px solid var(--accent); padding:8px 12px; background:var(--card); border-radius:4px; }}
code {{ font-size:13px; }}
</style></head><body><main>
<h1>Sympathetic neuron: AP transmission soma &rarr; axon &rarr; varicosity</h1>
<p class="muted">neuron/sympathetic (soma sympathetic, axon sympathetic_monolithic_v01, varicosity sympathetic). Simulation only, no module
file changed. Built by <code>tools/neuron_ap_transmission.py</code>; CVODE (myokit), rtol 1e-8, {SIM_TIME:g} s, output every {DT * 1e3:g} ms;
stimulus <code>I_in_mod_soma_membrane</code> from t = 0. The check is the neuron spec's invariants: equal spike counts (upward 0 mV
crossings, at least one) and axon peak within 20 mV of the soma peak.</p>

<h2>Recommendation</h2>
<div class="rec">
<p><b>Smallest change that passes: shorten the lumped axon, <code>L_ax_axon</code> 318 &rarr; 106 um (L/3)</b>, keeping d_ax, R_i,
channel_ratio. This is a single definitional parameter. It gives R 100 MOhm and C 3.3 pF, so tau_eff = RC/(1+channel_ratio) falls
from 2.0 ms to 0.22 ms. At 50 pA the varicosity peaks at {L3_50["pk_var"]:.1f} mV against the soma's {L3_50["pk_soma"]:.1f} mV, with a
soma&rarr;varicosity delay of {fmt(L3_50["delay"], 2)} ms. Spike counts match at every current, and the soma's spike counts are the
isolated soma's. NE release goes from {cur50["NE"]:.3f} to {L3_50["NE"]:.3f} uM peak. Use L/4 (80 um) if
you want the varicosity closer to the soma: L/4 is within 6 mV, L/3 within 14 mV. The first-spike latency at 50 pA is
{L3_50["lat_soma"]:.1f} ms (L/3) against {iso[0.05]["lat"]:.1f} ms isolated and {cur50["lat_soma"]:.1f} ms now, because the shorter axon loads the soma less.</p>
<p>The axon's length is no longer an anatomical length then. It becomes the electrotonic length of one lumped compartment, which
is honest for a single-node "pass the AP on" axon. <code>channel_ratio</code> ~10 passes too ({cr10_50["pk_var"]:.1f} mV), but it
has no physical meaning at that size: it scales the axial current by 11 on the axon side only, so charge is not conserved. R_i/10
(7.4 ohm cm) passes but is far below any axoplasm value, and d_ax only passes near d ~ 8 um, where C (80 pF) is far larger than the soma's own and the soma's AP shrinks (47 mV at d 5 um).
The HH prototype passes the check, but its axon fires first and drives the soma, and at 22 degC kinetics it oscillates: see option 3.</p>
</div>

<h2>Why it attenuates</h2>
<ul>
<li><b>The varicosity is no load.</b> <code>SN_varicosity_membrane</code> is a pass-through: V = V_in (the axon's V). Its ionic
currents (i_CaN, i_NaCa) drive Ca only and never feed back on V. So the axon node <i>is</i> the varicosity; the two traces are
identical. A real varicosity of d 1.24 um (A 4.84 um<sup>2</sup>) would add only 0.05 pF against the axon's 10 pF, so its load is
negligible either way (option 4: nothing to fix).</li>
<li><b>The axon is a first-order low-pass with unit DC gain.</b> The node has no membrane leak and its far end is sealed. All the
axial current I = (V_soma &minus; V)/R charges the node's own capacitance, scaled by (1 + channel_ratio):
dV/dt = (1 + cr) (V_soma &minus; V)/(R C). That is an RC filter with tau_eff = R C/(1 + cr) = 300 MOhm &times; 10 pF / 1.5 = 2.0 ms.</li>
<li><b>The AP is shorter than tau_eff.</b> The fast part of the somatic spike (above +25 mV) lasts only {w25:.2f} ms;
the trace is above 0 mV for {w0:.2f} ms, mostly on a slow shoulder near 0 mV (half-amplitude width {hw:.2f} ms). The node charges for
less than one time constant, reaches only about {cur50["pk_axon"]:.0f} mV (roughly 1 &minus; exp(&minus;w/tau) of the ~125 mV swing) and never crosses 0 mV. The sweep below shows that every single-compartment change (R_i, L_ax, d_ax, channel_ratio) falls on
one curve of tau_eff. The check needs tau_eff below about 0.3 ms.</li>
<li><b>Load on the soma.</b> The axon is a series R&ndash;C (300 MOhm, C/(1+cr) = 6.7 pF) hung on the soma's {C_SOMA:.0f} pF. During
the AP it draws up to {cur50["I_ax"]:.2f} nA, about {cur50["I_ax"] / 0.05:.0f} times the 50 pA stimulus, but only for the AP's duration. It delays the first spike ({cur50["lat_soma"]:.1f} vs {iso[0.05]["lat"]:.1f} ms isolated at 50 pA), <b>but it does not
cause the single spike:</b> the isolated soma also fires once in 1 s at 20, 50 and 100 pA (and 6 times at 200 pA, the same as in
the neuron). So "the axon drains the soma" is not borne out in this library. The neuron's single spike at 100 pA is the soma's own
behaviour (its M/SK adaptation), which belongs to the soma review.</li>
<li><b>Why the "+1.5" does not help:</b> channel_ratio only shortens tau_eff by 1/(1+cr). A first-order lag can never overshoot its
input, so no linear passive node gives a regenerated AP. It can at most follow the soma closely when tau_eff &lt;&lt; AP width.
A negative channel_ratio slows it further (cr &lt; &minus;1 would make it unstable).</li>
</ul>
<div class="grid">{card('sweep', 'Axon minus soma peak at 50 pA against tau_eff for each parameter family (one model, runtime parameters). The families overlap: tau_eff is the only thing that matters. d_ax gives tau ~ 1/d but raises C, so the soma peak drops (d x 5: soma 47 mV).')}
{card('ladder', 'Option 5 estimate: a passive N-node ladder with the same total R and C, driven by the recorded somatic AP (one-way, no load on the soma). Splitting helps (the distal mode is 4RC/pi^2): 16 nodes reach about +12 mV, still 44 mV short, so no passive split passes at R 300 MOhm, C 10 pF.')}</div>

<h2>Options tested</h2>
<div class="scroll"><table><tr><th>option</th><th>change</th></tr>{opt_rows}</table></div>

<h2>Metrics</h2>
<p class="muted">Counts are soma/axon/varicosity upward 0 mV crossings in {SIM_TIME:g} s. Delay: first somatic to first varicosity
crossing. Soma vs isolated: spike counts and first-spike latency against soma/sympathetic alone (<code>I_in_mod_membrane</code>). I_ax:
peak |axial current|. Cai and NE: varicosity peak cytosolic Ca and peak NE (uM). A negative delay means the axon crossed 0 mV before the soma. NE is <code>mod_varicosity_NE/NE</code>, driven by
Cai above 0.1 uM.</p>
<div class="scroll"><table>
<tr><th>option</th><th>pA</th><th>spikes s/a/v</th><th>soma pk (mV)</th><th>axon pk</th><th>var pk</th><th>delay (ms)</th>
<th>soma vs isolated: n; latency</th><th>I_ax (nA)</th><th>Cai (uM)</th><th>NE (uM)</th><th>check</th></tr>
{''.join(rows)}</table></div>

<h2>Traces at 50 pA</h2>
<p class="muted">Left: the full run, with the isolated soma in grey. Right: 16 ms around the first somatic spike. The varicosity trace
(dashed) lies on the axon's in every option, because the varicosity membrane is a pass-through.</p>
<div class="grid">{traces}</div>

<h2>Traces at 200 pA (repetitive firing)</h2>
<div class="grid">{traces200}</div>

<h2>Notes per option</h2>
<ul>
<li><b>1. Axial resistance.</b> Only the product R C/(1+cr) matters. L_ax scales R and C together (tau ~ L<sup>2</sup>), so it is
the most effective lever and it also lightens the soma's load. R_i alone works but needs R_i ~ 7 ohm cm, about 10 times below
axoplasm (70&ndash;200 ohm cm; Major 1994 gives 170&ndash;340). d_ax alone: tau ~ 1/d but C ~ d, so the soma is loaded (peak 47 mV at
d 5 um) before the axon passes.</li>
<li><b>2. channel_ratio.</b> The sign is consistent: positive values speed up the node. Passing needs channel_ratio ~ 8&ndash;10. Since
I_channels = cr I is added on the axon side only (the soma loses I, not (1+cr) I), this is an unbalanced gain rather than a
membrane current. It is not a physiological fix.</li>
<li><b>3. Regenerative HH current (prototype, not committed).</b> The SN_axon component in a scratch copy of the generated model got
HH Na/K/leak (E_Na +50 mV, E_K &minus;82 mV, E_L {hh_EL:.1f} mV for zero net current at &minus;70 mV). With 2&times; squid Na density and
6.3 degC kinetics it passes the check at every current, but for the wrong reason: the axon fires <i>first</i> (negative delay)
and its AP back-propagates into the soma and triggers it. The soma's firing changes: first spike at 8.7 ms instead of 22.3 ms at
50 pA, 744 ms instead of 297 ms at 20 pA, and 10 spikes instead of 6 at 200 pA. So it fails the "soma close to the isolated soma" criterion.
At the model's 22 degC (phi 5.6) the lumped node only reaches about +8 mV with squid densities. With 3&times; Na it fires on its own,
repetitively, at the axon's resting load: one lumped HH node coupled through 300 MOhm is an oscillator, not a cable. A proper active
axon needs several compartments (see 5) and channel data for sympathetic C-fibres. That is a model change, not a fix.</li>
<li><b>4. Varicosity load.</b> None: the membrane is a pass-through. The varicosity's own geometry is reasonable for a bouton
(d 1.24 um). If a capacitive varicosity were added (C = c_m A = 0.05 pF), it would change nothing at this scale.</li>
<li><b>5. Several compartments.</b> With the same total R and C, a sealed passive ladder's far end is faster than the lumped node
(dominant mode 4RC/pi<sup>2</sup> = 0.41 RC). The peak rises with N but stays more than 40 mV below the soma's:</li>
</ul>
<div class="scroll"><table style="max-width:420px"><tr><th>N</th><th>channel_ratio</th><th>far-end peak (mV)</th></tr>{ladder_rows}</table></div>
<p class="muted">At L/3 (R 100 MOhm, C 3.3 pF, cr 0.5) the far-end peak is {ladder3[0][1]:.1f} mV lumped and {ladder3[1][1]:.1f} mV with 16 nodes,
so the lumped node is a fair stand-in once tau_eff is short. Splitting alone does not rescue a passive 318 um axon; propagation
over that length needs active channels in each compartment.</p>
</main></body></html>
'''
    with open(args.out, 'w') as f:
        f.write(page)
    print('wrote', args.out)


if __name__ == '__main__':
    main()
