"""
i_M Argus2026_v01, parameterisation davis2020_wistar_kinetics: digitise the representative Wistar
M-current deactivation trace of Davis et al. 2020 (Hypertension 76:1915) Fig. 3C into a libcuflynx
obs_data, and build the refit comparison page (current vs calibrated tau_w).

    venv/bin/python tools/i_M_davis2020_kinetics.py obs-data     # digitise -> obs_data.json + source figure
    venv/bin/python -m pytest tests/test_modules.py --component i_M/Argus2026_v01 \
        -k 'calibrate and davis2020_wistar_kinetics' --include-unreviewed   # calibrate (libcuflynx)
    venv/bin/python tools/i_M_davis2020_kinetics.py report       # -> reviews/i_M_Argus2026_v01_tau_w_refit.html
    venv/bin/python tools/i_M_davis2020_kinetics.py gm-report    # -> reviews/i_M_Argus2026_v01_g_M_fit.html (g_M provenance,
                                                                 #    parameterisation davis2020_wistar; adds the thesis crop to its source_figures)

Digitisation (obs-data). The figure's embedded raster (page 4, 1310 x 926 px at 200 ppi) is read
with pdfimages; the Wistar trace (columns 0-300, rows 680-840) is traced column by column as the
median row of its dark pixels. Scales: 21 px = 100 pA (scale bar); time from the step itself,
which Davis 2021 (thesis p. 172) gives as 1 s: onset at column 40, return at 212.5, so
172.5 px = 1 s (the 100 ms scale bar is 17.5-18 px, consistent).

The trace is a raw current (not XE-991-subtracted), so its absolute level and its instantaneous
(ohmic) jumps include leak and other currents. Only the slow relaxations are M current. Each
relaxation is therefore normalised to its own end points, and mapped onto the M gate w with the
module's fixed activation curve (V_mid_w 35 mV, V_den_w 10 mV):

    step to -55 mV:  r(t) = (I(t) - I(0+)) / (I_inf - I(0+)),   w(t) = w_inf(-25) + (w_inf(-55) - w_inf(-25)) r(t)
    back to -25 mV:  r(t) likewise,                             w(t) = w_inf(-55) + (w_inf(-25) - w_inf(-55)) r(t)

with I(0+) and I_inf the start and end of a single exponential fitted to that relaxation (used for
the two end points only; the data points are the digitised samples). The data items are then
series of mod/w, so libcuflynx compares the model gate directly (no custom operation): the fit
constrains the time course (tau_w), not the amplitude (g_M, calibrated in davis2020_wistar).
"""
import argparse
import base64
import io
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
VDIR = os.path.join(REPO, 'modules', 'cell', 'ion_channels', 'i_M', 'versions', 'Argus2026_v01')
PARAMETERISATION = 'davis2020_wistar_kinetics'
IDIR = os.path.join(VDIR, 'parameterisations', PARAMETERISATION)
PDF = os.path.expanduser('~/Zotero/storage/4KCUPWQE/Davis et al. - 2020 - Downregulation of M Current Is Coupled to Membrane.pdf')
PAGE = 4
HTML = os.path.join(REPO, 'reviews', 'i_M_Argus2026_v01_tau_w_refit.html')

# digitisation geometry (pixels of the embedded raster)
ROWS = (680, 840)
DARK = 110
X_ON, X_OFF = 40.0, 212.5        # step onset and return
PX_PER_S = X_OFF - X_ON          # 1 s step
PA_PER_PX = 100.0 / 21.0
X_END = 299                      # last column of the trace
HOLD, STEP = -25.0, -55.0
T_HOLD, T_STEP, T_RETURN = 10.0, 1.0, 0.5
SKIP_STEP, SKIP_RETURN = 1, 3    # samples at the start of each phase inside the jump: std 1 (not fitted)


def w_inf(v, v_mid=35.0, v_den=10.0):
    return 1.0 / (1.0 + np.exp(-(v + v_mid) / v_den))


def tau_w(v, num, tmin, dv=0.0):
    """ms; the module's tau_w (shape constants 2.0, 15, 25 fixed)."""
    return tmin + num / (2.0 * np.exp((v + 35.0 + dv) / 15.0) + np.exp(-(v + 35.0 + dv) / 25.0))


def tau_martin(v):
    """martin2023modelling Methods p. 16: 400/(3.3 exp((V+35)/10) + exp(-(V+35)/20)) ms."""
    return 400.0 / (3.3 * np.exp((v + 35.0) / 10.0) + np.exp(-(v + 35.0) / 20.0))


def raster():
    with tempfile.TemporaryDirectory() as d:
        subprocess.run(['pdfimages', '-png', '-f', str(PAGE), '-l', str(PAGE), PDF, os.path.join(d, 'p')], check=True)
        from PIL import Image
        return Image.open(os.path.join(d, 'p-000.png')).convert('RGB')


def trace(img):
    """{column: median dark row} of the Wistar trace."""
    g = np.asarray(img.convert('L'), dtype=float)
    out = {}
    for x in range(0, X_END + 1):
        ys = np.where(g[ROWS[0]:ROWS[1], x] < DARK)[0]
        if len(ys):
            out[x] = float(np.median(ys)) + ROWS[0]
    return out


def fit_exp(t, y):
    from scipy.optimize import curve_fit
    f = lambda t, a, tau, c: c + a * np.exp(-t / tau)
    p, cov = curve_fit(f, t, y, p0=(y[0] - y[-1], 0.1, y[-1]), maxfev=20000)
    return p, np.sqrt(np.diag(cov)), f


def phase(tr, x0, n, skip):
    """(t grid, current pA at it, exponential fit) of the relaxation starting at column x0."""
    xs = np.array(sorted(tr))
    ys = np.array([tr[x] for x in xs])
    t = np.arange(n) / PX_PER_S
    cur = -np.interp(x0 + np.arange(n), xs, ys) * PA_PER_PX     # pA, up to an offset (rows grow downwards)
    p, se, f = fit_exp(t[skip:], cur[skip:])
    return t, cur, p, se, f


def build_obs_data():
    img = raster()
    tr = trace(img)
    n_step = int(np.floor(T_STEP * PX_PER_S)) + 1 - 1          # samples t = k/172.5 < 1 s
    n_ret = int(min(np.floor(T_RETURN * PX_PER_S), X_END - X_OFF)) + 1
    w25, w55 = w_inf(HOLD), w_inf(STEP)
    items, info = [], {}
    for name, x0, n, skip, wa, wb, sub, label in (
            ('w_deactivation_at_-55mV', X_ON, n_step, SKIP_STEP, w25, w55, 1, 'deactivation at -55 mV'),
            ('w_reactivation_at_-25mV', X_OFF, n_ret, SKIP_RETURN, w55, w25, 2, 'reactivation on return to -25 mV')):
        t, cur, p, se, f = phase(tr, x0, n, skip)
        a, tau, c = p
        i0, iinf = c + a, c
        r = (cur - i0) / (iinf - i0)
        resid = r[skip:] - (f(t[skip:], *p) - i0) / (iinf - i0)
        sd = float(np.std(resid, ddof=3))
        w = wa + (wb - wa) * r
        w[:skip] = wa
        std = np.full(n, abs(wb - wa) * sd)
        std[:skip] = 1.0
        info[name] = {'tau_ms': 1e3 * tau, 'tau_se_ms': 1e3 * se[1], 'amplitude_pA': abs(iinf - i0),
                      'residual_sd_fraction': sd, 'n': n}
        items.append({
            'data_item_name': name,
            'comment': (f'davis2020downregulation Fig. 3C Wistar trace, {label}: digitised (tools/i_M_davis2020_kinetics.py), '
                        f'normalised to the relaxation\'s own end points and mapped onto w (w_inf(-25) = {w25:.4f}, '
                        f'w_inf(-55) = {w55:.4f}). Sample k at k/172.5 s from the start of the sub-experiment (one '
                        f'raster column). Single-exponential fit of the digitised current: tau {1e3 * tau:.0f} +- '
                        f'{1e3 * se[1]:.0f} ms, amplitude {abs(iinf - i0):.0f} pA. std: the residual SD about that '
                        f'fit ({sd:.3f} of the relaxation) in w units; the first {skip} sample(s), inside the '
                        f'voltage jump, have std 1 (not fitted). n = 1 cell, raw current (not XE-991-subtracted), '
                        f'room temperature.'),
            'trace_name_for_plotting': f'w ({label})',
            'data_type': 'series',
            'operation': None,
            'operands': ['mod/w'],
            'unit': 'dimensionless',
            'cost_type': 'gaussian_MLE',
            'weight': 1.0,
            'value': [round(float(v), 6) for v in w],
            'std': [round(float(v), 6) for v in std],
            'obs_dt': 1.0 / PX_PER_S,
            'experiment_idx': 0,
            'subexperiment_idx': sub,
        })
    obs = {
        'obs_data_name': PARAMETERISATION,
        'comment': ('Davis et al. 2020 (Hypertension 76:1915) Fig. 3C, the representative Wistar M-current deactivation '
                    'trace (perforated patch, room temperature; raw current, not XE-991-subtracted; n = 1): held at '
                    '-25 mV, 1 s step to -55 mV (davis2021transcriptomic p. 172), back to -25 mV. Digitised and normalised by '
                    'tools/i_M_davis2020_kinetics.py (see its docstring): each slow relaxation is mapped onto the gate '
                    'w with the module\'s fixed w_inf, so the data constrain tau_w only. Protocol: the channel module '
                    'alone, V clamped: sub-experiment 0 holds -25 mV for 10 s (steady state), 1 steps to -55 mV for '
                    '1 s, 2 returns to -25 mV for 0.5 s (the trace ends ~0.5 s after the return). To redo: edit this '
                    'file (or re-digitise with the tool) and run venv/bin/python -m pytest tests/test_modules.py '
                    '--component i_M/Argus2026_v01 -k \'calibrate and davis2020_wistar_kinetics\' --include-unreviewed'),
        'protocol_info': {
            'pre_times': [0.0],
            'sim_times': [[T_HOLD, T_STEP, T_RETURN]],
            'params_to_change': {'mod/V': [[HOLD, STEP, HOLD]]},
            'experiment_labels': ['davis2020downregulation Fig. 3C: hold -25 mV, 1 s at -55 mV, back to -25 mV'],
        },
        'prediction_items': [],
        'data_items': items,
    }
    os.makedirs(IDIR, exist_ok=True)
    with open(os.path.join(IDIR, f'{PARAMETERISATION}_obs_data.json'), 'w') as f:
        json.dump(obs, f, indent=1, ensure_ascii=False)
    # the source figure: panel C, cropped from the raster at its native 200 ppi
    sdir = os.path.join(IDIR, 'source_figures')
    os.makedirs(sdir, exist_ok=True)
    img.crop((0, 495, 1310, 926)).save(os.path.join(sdir, 'davis2020_fig3c.png'), optimize=True)
    with open(os.path.join(sdir, 'source_figures.json'), 'w') as f:
        json.dump([{'file': 'davis2020_fig3c.png', 'source': 'davis2020downregulation; Fig. 3C',
                    'caption': ('Left: the representative Wistar deactivation trace (hold -25 mV, 1 s step to -55 mV, '
                                'raw current, n = 1), digitised column by column; each slow relaxation normalised and '
                                'mapped onto the gate w (the obs_data series). Scale bar 100 pA / 100 ms.')}],
                  f, indent=2)
    for k, v in info.items():
        print(f'{k}: tau {v["tau_ms"]:.1f} +- {v["tau_se_ms"]:.1f} ms, amplitude {v["amplitude_pA"]:.0f} pA, '
              f'residual SD {v["residual_sd_fraction"]:.3f}, n {v["n"]}')
    return obs, tr



# ---- the refit comparison page ------------------------------------------------------------------

def _chain(helper, names, segments, outputs, starts=None):
    """Runs the helper through ``segments`` [(params {name: value}, n_runs)], continuing from the
    previous state; returns (t, {output: array}) concatenated. ``starts`` (a list) receives each
    segment's first sample index (its first sample is the state at the switch, at the new parameters)."""
    from cam_testing import harness
    import contextlib
    helper.reset_states()
    ts, outs, t0 = [], {o: [] for o in outputs}, 0.0
    n_done = 0
    for params, n_runs in segments:
        if starts is not None:
            starts.append(n_done)
        helper.set_param_vals(list(params), list(params.values()))
        for _ in range(n_runs):
            with contextlib.redirect_stdout(io.StringIO()):
                if helper.run() is False:
                    raise RuntimeError('simulation failed')
            t = np.asarray(helper.get_time(), dtype=float)
            ts.append(t - t[0] + t0)
            t0 = ts[-1][-1]
            n_done += len(t)
            for o in outputs:
                outs[o].append(np.asarray(helper.get_results([[harness.output_name(o)]], flatten=True)[0], dtype=float).ravel())
    return np.concatenate(ts), {o: np.concatenate(v) for o, v in outs.items()}


def _png(fig):
    import matplotlib.pyplot as plt
    b = io.BytesIO()
    fig.savefig(b, format='png', dpi=110, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def w_exact(t, v_seq, num, tmin):
    """w(t) of the protocol (hold -25 to steady state, then the step and the return), exact."""
    out = np.empty_like(t)
    w0 = w_inf(HOLD)
    edges = [0.0, T_STEP, T_STEP + T_RETURN]
    for k, (v, a, b) in enumerate(zip(v_seq, edges[:-1], edges[1:])):
        m = (t >= a) & (t <= b) if k == len(v_seq) - 1 else (t >= a) & (t < b)
        tau = tau_w(v, num, tmin) / 1000.0
        out[m] = w_inf(v) + (w0 - w_inf(v)) * np.exp(-(t[m] - a) / max(tau, 1e-12))
        w0 = w_inf(v) + (w0 - w_inf(v)) * np.exp(-(b - a) / max(tau, 1e-12))
    return out


def build_report():
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cam_testing import checks, harness
    from cam_testing.library import load_version

    with open(os.path.join(IDIR, f'{PARAMETERISATION}_obs_data.json')) as f:
        obs = json.load(f)
    with open(os.path.join(IDIR, f'{PARAMETERISATION}_calibration.json')) as f:
        cal = json.load(f)
    cur = {'tau_w_num': 1200.0, 'tau_w_min': 30.0}
    with open(os.path.join(VDIR, 'parameterisations', 'default', 'default_parameters.csv'), newline='') as f:
        for r in csv.DictReader(f):
            if r['variable_name'] in cur:
                cur[r['variable_name']] = float(r['value'])
    new = {k: float(cal['calibrated_parameters'][k]) for k in ('tau_w_num', 'tau_w_min')}
    fmt = lambda x: '0' if abs(x) < 1e-6 else f'{x:.4g}'
    sets = [('current', cur, '#6b7280'), ('calibrated', new, '#d9480f')]

    # data (w units), by phase, with the samples that carry data
    items = obs['data_items']
    data = []
    for it, off in zip(items, (0.0, T_STEP)):
        t = off + np.arange(len(it['value'])) * it['obs_dt']
        v, sd = np.asarray(it['value']), np.asarray(it['std'])
        keep = sd < 0.5
        data.append((it['data_item_name'], t, v, sd, keep))

    def residuals(p):
        r = []
        for name, t, v, sd, keep in data:
            r.append((w_exact(t[keep], [STEP, HOLD], p[0], p[1]) - v[keep]) / sd[keep])
        return np.concatenate(r)

    def nrmse(p):
        out = {}
        for name, t, v, sd, keep in data:
            m = w_exact(t[keep], [STEP, HOLD], p[0], p[1])
            out[name] = float(np.sqrt(np.mean((m - v[keep]) ** 2)) / np.ptp(v[keep]))
        return out

    # Laplace (Gauss-Newton) standard errors at the calibrated values, from the same residuals
    p0 = np.array([new['tau_w_num'], new['tau_w_min']])
    J = np.empty((residuals(p0).size, 2))
    for j in range(2):
        h = max(1e-3 * abs(p0[j]), 0.5)
        pp, pm = p0.copy(), p0.copy()
        pp[j] += h
        pm[j] = max(pm[j] - h, 0.0)
        J[:, j] = (residuals(pp) - residuals(pm)) / (pp[j] - pm[j])
    try:
        cov = np.linalg.inv(J.T @ J)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        se = [np.nan, np.nan]

    # the module simulated (CVODE) through the protocol at both parameter sets
    imv = load_version('i_M', 'Argus2026_v01')
    cm = checks.ComponentModel(imv)
    by_var = {p.variable_name: p for p in cm.parameters()}
    pn = lambda v: 'parameters/' + harness.parameter_name(by_var[v])
    helper = harness.simulation_helper(cm.model_path, dict(cm.spec, sim_time=0.5, dt=0.001))
    sims = {}
    for label, p, _c in sets:
        base = {pn('tau_w_num'): p['tau_w_num'], pn('tau_w_min'): p['tau_w_min']}
        t, o = _chain(helper, None, [(dict(base, **{pn('V'): HOLD}), 20), (dict(base, **{pn('V'): STEP}), 2),
                                     (dict(base, **{pn('V'): HOLD}), 1)], ['w', 'i_M'])
        sims[label] = (t - 10.0, o)
        exact = w_exact(np.clip(t - 10.0, 0, None)[t >= 10.0], [STEP, HOLD], p['tau_w_num'], p['tau_w_min'])
        sims[label + '_maxerr'] = float(np.max(np.abs(o['w'][t >= 10.0] - exact)))

    # figure 1: the protocol, data and both models (w), and i_M at the default g_M
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True, gridspec_kw={'height_ratios': [3, 2]})
    ax = axes[0]
    for name, t, v, sd, keep in data:
        ax.errorbar(t[keep], v[keep], yerr=sd[keep], fmt='o', ms=2.5, color='#1f2937', ecolor='#9ca3af', elinewidth=0.6,
                    label='davis2020downregulation Fig. 3C, digitised (w units)' if name.startswith('w_de') else None)
    for label, p, c in sets:
        t, o = sims[label]
        ax.plot(t, o['w'], color=c, lw=1.8, label=f'{label}: tau_w_num {fmt(p["tau_w_num"])}, tau_w_min {fmt(p["tau_w_min"])} ms')
    ax.axvspan(0, T_STEP, color='#e5e7eb', zorder=0)
    ax.set_ylabel('M gate w')
    ax.set_xlim(-0.15, T_STEP + T_RETURN)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title('hold -25 mV, 1 s at -55 mV (shaded), back to -25 mV')
    ax = axes[1]
    for label, p, c in sets:
        t, o = sims[label]
        ax.plot(t, 1e3 * o['i_M'], color=c, lw=1.5, label=label)
    ax.set_ylabel('i_M (pA), default g_M')
    ax.set_xlabel('time from the step (s)')
    ax.axvspan(0, T_STEP, color='#e5e7eb', zorder=0)
    ax.legend(fontsize=8)
    fig_trace = _png(fig)

    # figure 2: tau_w(V)
    V = np.linspace(-90, 10, 401)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for label, p, c in sets:
        ax.plot(V, tau_w(V, p['tau_w_num'], p['tau_w_min']), color=c, lw=2, label=f'{label}')
    ax.plot(V, tau_martin(V), color='#2563eb', lw=1.5, ls='--', label='Martin & Pedersen 2023 (as published, 37 degC model)')
    taus = []
    for it in items:
        c = it['comment']
        i = c.index('tau ') + 4
        tau_v = float(c[i:].split(' ')[0])
        tau_e = float(c[i:].split('+- ')[1].split(' ')[0])
        taus.append((tau_v, tau_e))
    ax.errorbar([STEP, HOLD], [t for t, e in taus], yerr=[e for t, e in taus], fmt='s', color='#111827', ms=6, capsize=3,
                label='davis2020downregulation Fig. 3C trace (single-exp. fits, +- SE)')
    ax.set_xlabel('V (mV)')
    ax.set_ylabel('tau_w (ms)')
    ax.set_ylim(0, 560)
    ax.legend(fontsize=8)
    fig_tau = _png(fig)

    # the neuron: soma supermodule, current vs calibrated tau_w
    sv = load_version('soma', 'sympathetic')
    scm = checks.ComponentModel(sv)
    sby = {p.variable_name: p for p in scm.parameters()}
    spn = lambda v: 'parameters/' + harness.parameter_name(sby[v])
    Vout = 'mod_membrane/V'
    rest_helper = harness.simulation_helper(scm.model_path, dict(scm.spec, sim_time=30.0, dt=1e-4))
    step_helper = harness.simulation_helper(scm.model_path, dict(scm.spec, sim_time=1.0, dt=5e-5))
    neuron = {}
    spikes = lambda t, v: int(np.sum((v[:-1] < 0) & (v[1:] >= 0)))
    for label, p, c in sets:
        base = {spn('tau_w_num_i_M'): p['tau_w_num'], spn('tau_w_min_i_M'): p['tau_w_min']}
        t, o = harness.run(rest_helper, [Vout], dict(base, **{spn('I_in_membrane'): 0.0}))
        v = o[harness.output_key(Vout)]
        neuron[label] = {'rest': (t, v, spikes(t, v), int(np.sum((v[:-1] < 0) & (v[1:] >= 0) & (t[1:] >= 20.0))))}   # spikes in total, and in the last 10 s
        for pa in (5, 9, 50):
            t, o = _chain(step_helper, None, [(dict(base, **{spn('I_in_membrane'): 0.0}), 3),
                                               (dict(base, **{spn('I_in_membrane'): pa * 1e-3}), 1)], [Vout])
            v = o[Vout]
            m = t >= 3.0
            neuron[label][pa] = (t[m] - 3.0, v[m], spikes(t[m], v[m]), t[t < 3.0], v[t < 3.0])
    fig, axes = plt.subplots(4, 1, figsize=(7.2, 8.4))
    ax = axes[0]
    for label, p, c in sets:
        t, v, n, n_last = neuron[label]['rest']
        ax.plot(t, v, color=c, lw=0.8, label=f'{label}: {n} spike(s) in 30 s, {n_last} in the last 10 s')
    ax.set_title('rest, 0 pA, 30 s')
    ax.set_ylabel('V (mV)')
    ax.legend(fontsize=8)
    for ax, pa in zip(axes[1:], (5, 9, 50)):
        for label, p, c in sets:
            t, v, n = neuron[label][pa][:3]
            ax.plot(t, v, color=c, lw=0.8, label=f'{label}: {n} spikes')
        ax.set_title(f'{pa} pA step, 1 s (after 3 s at 0 pA)')
        ax.set_ylabel('V (mV)')
        ax.legend(fontsize=8)
    axes[-1].set_xlabel('time (s)')
    fig.tight_layout()
    fig_neuron = _png(fig)

    with open(os.path.join(IDIR, 'source_figures', 'davis2020_fig3c.png'), 'rb') as f:
        src = base64.b64encode(f.read()).decode()

    # g_M if tau_w were the calibrated one: the davis2020_wistar relaxation (1 s, last 10 %) is then complete
    def relax(p):
        tt = np.arange(0, 1.0001, 0.001)
        tau = tau_w(STEP, p['tau_w_num'], p['tau_w_min']) / 1000
        ww = w_inf(STEP) + (w_inf(HOLD) - w_inf(STEP)) * np.exp(-tt / tau)
        n = len(tt)
        return ww[int(0.9 * (n - 1)):int(n - 1)].mean() - ww[0]
    ek = -96.4979
    g_cur = -0.12096 / (relax(cur) * (STEP - ek))
    g_new = -0.12096 / (relax(new) * (STEP - ek))
    n_cur, n_new = nrmse([cur['tau_w_num'], cur['tau_w_min']]), nrmse([new['tau_w_num'], new['tau_w_min']])
    cost_cur = float(np.mean(residuals(np.array([cur['tau_w_num'], cur['tau_w_min']])) ** 2) / 2)
    cost_new = float(np.mean(residuals(p0) ** 2) / 2)

    rows = ''.join(
        f'<tr><td>{k}</td><td>{fmt(cur[k])} ms</td><td>{fmt(new[k])} ms</td><td>{s_:.3g} ms</td></tr>'
        for k, s_ in zip(('tau_w_num', 'tau_w_min'), se))
    tau_rows = ''.join(
        f'<tr><td>{v:g} mV</td><td>{tau_w(v, cur["tau_w_num"], cur["tau_w_min"]):.0f}</td>'
        f'<td>{tau_w(v, new["tau_w_num"], new["tau_w_min"]):.0f}</td><td>{tau_martin(v):.0f}</td><td>{d}</td></tr>'
        for v, d in ((-70, '-'), (-55, f'{taus[0][0]:.0f} +- {taus[0][1]:.0f}'), (-40, '-'),
                     (-25, f'{taus[1][0]:.0f} +- {taus[1][1]:.0f}'), (0, '-')))
    fit_rows = ''.join(
        f'<tr><td>{k}</td><td>{n_cur[k]:.3f}</td><td>{n_new[k]:.3f}</td></tr>' for k in n_cur)
    neu_rows = ''.join(
        f'<tr><td>{lab}</td><td>{neuron["current"][key][2]}</td><td>{neuron["calibrated"][key][2]}</td></tr>'
        for lab, key in (('5 pA, 1 s', 5), ('9 pA, 1 s', 9), ('50 pA, 1 s', 50)))
    r_c, r_n = neuron['current']['rest'], neuron['calibrated']['rest']
    neu_rows = (f'<tr><td>rest 0 pA, 30 s: all (last 10 s)</td><td>{r_c[2]} ({r_c[3]})</td><td>{r_n[2]} ({r_n[3]})</td></tr>'
                + neu_rows)
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>i_M tau_w refit</title>
<style>
:root {{ --bg:#ffffff; --fg:#111827; --muted:#4b5563; --card:#f9fafb; --line:#e5e7eb; --accent:#d9480f; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; }} }}
:root[data-theme="dark"] {{ --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; }}
body {{ background:var(--bg); color:var(--fg); font:15px/1.5 system-ui, -apple-system, Segoe UI, sans-serif; margin:0; padding:24px 16px; }}
main {{ max-width:1180px; margin:0 auto; }}
h1 {{ font-size:1.5rem; margin:0 0 4px; }} h2 {{ font-size:1.15rem; margin:28px 0 8px; border-bottom:1px solid var(--line); padding-bottom:4px; }}
p, li {{ color:var(--fg); }} .muted {{ color:var(--muted); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit, minmax(340px, 1fr)); gap:16px; align-items:start; }}
.card {{ background:var(--card); border:1px solid var(--line); border-radius:8px; padding:12px; }}
.card img {{ width:100%; height:auto; background:#fff; border-radius:4px; display:block; }}
table {{ border-collapse:collapse; width:100%; font-size:14px; }} th, td {{ border-bottom:1px solid var(--line); padding:4px 8px; text-align:left; }}
th {{ color:var(--muted); font-weight:600; }}
.warn {{ border-left:4px solid var(--accent); padding:8px 12px; background:var(--card); }}
code {{ font-size:13px; }}
</style></head><body><main>
<h1>i_M tau_w refit to Davis 2020 Fig. 3C</h1>
<p class="muted">i_M Argus2026_v01, parameterisation <code>{PARAMETERISATION}</code>, libcuflynx {cal['method']} calibration of tau_w_num and tau_w_min ({cal['date']}).
The defaults are <b>not</b> changed; keep or reject on this page. Rebuilt by <code>venv/bin/python tools/i_M_davis2020_kinetics.py report</code>.</p>

<div class="warn"><b>Caveats.</b> The data are one cell (n = 1), the representative trace of Fig. 3C. It is a raw current (not XE-991-subtracted),
so only the slow relaxations are taken as M current, each normalised to its own end points and mapped onto the gate w with the module's
fixed w_inf. The recording was at room temperature, and the trace was digitised from the published raster (about 6 ms per pixel).
The tau_w shape (prefactor 2.0, slopes 15/25 mV) is kept fixed, so tau(-25) &lt; tau(-55) for any parameter values, while the trace has
tau(-25) about 2x tau(-55): <b>no parameter values fit both phases</b>, and the fit is a compromise, which here is a nearly voltage-independent tau_w.</div>

<h2>The source and the fit</h2>
<div class="grid">
<div class="card"><img alt="Davis 2020 Fig. 3C" src="data:image/png;base64,{src}">
<p class="muted">Davis et al. 2020, Hypertension 76:1915, Fig. 3C (cropped from the PDF raster, 200 ppi). Left: the Wistar deactivation trace that was digitised
(scale bar 100 pA, 100 ms; step 1 s). The density panel gives the g_M calibration (parameterisation davis2020_wistar).</p></div>
<div class="card"><img alt="model vs data" src="data:image/png;base64,{fig_trace}">
<p class="muted">The module simulated with CVODE through the same protocol (time from the step). Top: the gate w at the current (grey) and calibrated (orange)
tau_w, over the digitised data in w units (bars: the std used in the fit). Bottom: i_M at the default g_M (0.005479 uS).
The simulated w matches the exact solution to {max(sims['current_maxerr'], sims['calibrated_maxerr']):.1e}.</p></div>
</div>

<h2>Parameters</h2>
<div class="grid"><div class="card"><table>
<tr><th>parameter</th><th>current (default)</th><th>calibrated</th><th>SE (Gauss-Newton)</th></tr>{rows}
</table>
<p class="muted">Range: tau_w_num 0-3000 ms, tau_w_min 0-200 ms. libcuflynx cost (gaussian_MLE, mean per weighted sample) {cal['cost']:.4g}.
CMA-ES gives no uncertainty; the SEs are from the Jacobian of the same std-weighted residuals at the optimum (tau_w_num is at its lower bound, so its SE is one-sided).
Mean 0.5 (residual/std)^2: current {cost_cur:.3g}, calibrated {cost_new:.3g}.</p></div>
<div class="card"><table><tr><th>fit (nrmse of w, samples with data)</th><th>current</th><th>calibrated</th></tr>{fit_rows}</table>
<p class="muted">The calibration test passes at nrmse &lt;= 0.25.</p></div></div>

<h2>tau_w(V)</h2>
<div class="grid"><div class="card"><img alt="tau_w(V)" src="data:image/png;base64,{fig_tau}"></div>
<div class="card"><table><tr><th>V</th><th>current (ms)</th><th>calibrated (ms)</th><th>Martin as published (ms)</th><th>Davis trace (ms)</th></tr>{tau_rows}</table>
<p class="muted">Martin &amp; Pedersen 2023 (Methods p. 16) take tau_w from Zhou et al. 2018, with no temperature factor; its recording temperature is not established.</p></div></div>

<h2>Effect on the neuron (soma sympathetic supermodule)</h2>
<div class="grid"><div class="card"><img alt="neuron" src="data:image/png;base64,{fig_neuron}"></div>
<div class="card"><table><tr><th>run</th><th>spikes, current tau_w</th><th>spikes, calibrated tau_w</th></tr>{neu_rows}</table>
<p class="muted">Default soma parameters (g_M 0.005479 uS) except tau_w_num_i_M and tau_w_min_i_M. Steps are 1 s, after 3 s at 0 pA from the model's initial state.
Spikes are upward crossings of 0 mV; the one spike at rest is the start-up transient (t &lt; 0.1 s) from the model's initial state.</p>
<p><b>Coupling with g_M.</b> The g_M calibration (davis2020_wistar) fits the relaxation within Davis's 1 s step, so it depends on tau_w: with the current tau_w the relaxation is 87 % complete at 1 s, giving g_M {g_cur:.4g} uS (the default). With the calibrated tau_w it is complete, giving g_M {g_new:.4g} uS. If you keep the new tau_w, rerun the davis2020_wistar calibration as well.</p></div></div>

<h2>Redo</h2>
<p><code>venv/bin/python tools/i_M_davis2020_kinetics.py obs-data</code> re-digitises the trace. Alternatively, edit
<code>parameterisations/{PARAMETERISATION}/{PARAMETERISATION}_obs_data.json</code> or <code>_params_for_id.csv</code> directly. Then run
<code>venv/bin/python -m pytest tests/test_modules.py --component i_M/Argus2026_v01 -k 'calibrate and {PARAMETERISATION}' --include-unreviewed</code>
and <code>venv/bin/python tools/i_M_davis2020_kinetics.py report</code>.</p>
</main></body></html>
"""
    os.makedirs(os.path.dirname(HTML), exist_ok=True)
    with open(HTML, 'w') as f:
        f.write(html)
    print('wrote', HTML)
    print('current', cur, 'calibrated', new, 'se', se, 'nrmse', n_cur, n_new)
    print('neuron spikes', {l: {k: neuron[l][k][2] for k in neuron[l]} for l in neuron})
    print('g_M current/new tau', g_cur, g_new)


# ---- the g_M provenance page --------------------------------------------------------------------

THESIS = os.path.expanduser('~/Documents/books/Physiology/Davis_2021_A_transcriptomic_and.pdf')
GM_PARAMETERISATION = 'davis2020_wistar'
GM_IDIR = os.path.join(VDIR, 'parameterisations', GM_PARAMETERISATION)
GM_HTML = os.path.join(REPO, 'reviews', 'i_M_Argus2026_v01_g_M_fit.html')
# davis2020downregulation Fig. 3C, Wistar M-current densities (pA/pF), digitised from the PDF (pdftoppm -r 600; the dot
# centres against the axis ticks 0/-2/-4/-6, +- about 0.05 pA/pF; six isolated dots confirmed by blob
# detection on the 200 ppi raster). Their median is -4.01 (Davis: -4.032).
WISTAR_DENSITIES = [-2.43, -2.38, -2.18, -4.01, -3.99, -3.96, -4.23, -4.35, -4.82, -4.87, -5.64]
# the density panel's y axis in the raster: pixel rows of 0 and -6 pA/pF
DENSITY_Y0, DENSITY_Y6 = 538.5, 883.0
CM_PF = 30.0
E_K = -96.4979


def thesis_crop():
    """The Davis 2021 thesis passage defining the measure (pp. 172-173 = PDF pages 173-174), at 200 dpi."""
    from PIL import Image
    s = 200 / 72.0
    parts = []
    with tempfile.TemporaryDirectory() as d:
        for page, (y0, y1) in ((173, (634, 768)), (174, (64, 200))):
            subprocess.run(['pdftoppm', '-r', '200', '-png', '-f', str(page), '-l', str(page), THESIS,
                            os.path.join(d, 't')], check=True, stderr=subprocess.DEVNULL)
            f = next(x for x in sorted(os.listdir(d)) if x.startswith('t') and x.endswith('.png'))
            im = Image.open(os.path.join(d, f)).convert('RGB')
            os.remove(os.path.join(d, f))
            parts.append(im.crop((int(75 * s), int(y0 * s), int(525 * s), int(y1 * s))))
    gap = 24
    out = Image.new('RGB', (max(p.width for p in parts), sum(p.height for p in parts) + gap), 'white')
    out.paste(parts[0], (0, 0))
    out.paste(parts[1], (0, parts[0].height + gap))
    return out


def fig3c_marked(img):
    """Fig. 3C crop with the digitised Wistar points, their median and IQR drawn beside the dot plot."""
    from PIL import ImageDraw
    crop = img.crop((0, 495, 1310, 926)).copy()
    dr = ImageDraw.Draw(crop)
    ypx = lambda v: DENSITY_Y0 + (DENSITY_Y6 - DENSITY_Y0) * (v / -6.0) - 495
    q1, med, q3 = np.percentile(WISTAR_DENSITIES, [25, 50, 75])
    x0 = 1012
    dr.rectangle((x0, ypx(q3), x0 + 26, ypx(q1)), outline=(217, 72, 15), width=2)
    dr.line((x0 - 4, ypx(med), x0 + 30, ypx(med)), fill=(217, 72, 15), width=3)
    dr.line((x0 + 30, ypx(-4.032), x0 + 46, ypx(-4.032)), fill=(37, 99, 235), width=3)
    for v in WISTAR_DENSITIES:
        y = ypx(v)
        dr.ellipse((x0 + 8, y - 3, x0 + 16, y + 3), outline=(217, 72, 15), width=2)
    return crop, (q1, med, q3)


def build_gm_report():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cam_testing import checks, harness
    from cam_testing.library import load_version

    with open(os.path.join(GM_IDIR, f'{GM_PARAMETERISATION}_obs_data.json')) as f:
        obs = json.load(f)
    with open(os.path.join(GM_IDIR, f'{GM_PARAMETERISATION}_calibration.json')) as f:
        cal = json.load(f)
    item = obs['data_items'][0]
    rho = float(cal['calibrated_parameters']['rho_M'])
    run = cal['runs'][0]
    zitem = run['constant_items'][item['data_item_name']]

    img = raster()
    marked, (q1, med, q3) = fig3c_marked(img)
    thesis = thesis_crop()
    sdir = os.path.join(GM_IDIR, 'source_figures')
    os.makedirs(sdir, exist_ok=True)
    thesis.save(os.path.join(sdir, 'davis2021thesis_p172_173.png'), optimize=True)
    with open(os.path.join(sdir, 'source_figures.json')) as f:
        figs = json.load(f)
    figs = [x for x in figs if x['file'] != 'davis2021thesis_p172_173.png'] + [{
        'file': 'davis2021thesis_p172_173.png', 'source': 'davis2021transcriptomic; pp. 172-173 (Ch. 4)',
        'caption': ('The definition of the Fig. 3C / Fig. 4.10C measure: XE-991-subtracted deactivation, held at -25 mV, '
                    '1 s hyperpolarisation to -55 mV, "the difference in current following the initial ohmic change, and the '
                    'steady-state current towards the end of the step" (the obs_data operation mean_in_range_minus_initial).')}]
    with open(os.path.join(sdir, 'source_figures.json'), 'w') as f:
        json.dump(figs, f, indent=2)

    def b64(im):
        b = io.BytesIO()
        im.save(b, format='PNG', optimize=True)
        return base64.b64encode(b.getvalue()).decode()

    # simulate the module through the protocol (CVODE), three g_M
    imv = load_version('i_M', 'Argus2026_v01')
    cm = checks.ComponentModel(imv)
    by_var = {p.variable_name: p for p in cm.parameters()}
    pn = lambda v: 'parameters/' + harness.parameter_name(by_var[v])
    A, gamma = float(by_var['A_mem'].value), float(by_var['gamma_M'].value)
    g_of = lambda r: A * r * gamma * 1e-6
    cases = [('calibrated (relaxation reading)', rho, '#d9480f'),
             ('steady-state reading (first fit)', 0.0025561 / (A * gamma * 1e-6), '#2563eb'),
             ('old default (0.00263 uS)', 0.00263 / (A * gamma * 1e-6), '#6b7280')]
    helper = harness.simulation_helper(cm.model_path, dict(cm.spec, sim_time=0.5, dt=0.001))
    res = []
    for label, r, c in cases:
        starts = []
        t, o = _chain(helper, None, [({pn('rho_M'): r, pn('V'): HOLD}, 20), ({pn('rho_M'): r, pn('V'): STEP}, 2)],
                      ['i_M', 'w'], starts)
        t = t - 10.0
        i = 1e3 * o['i_M']                                     # pA
        step = np.arange(len(t)) >= starts[1]
        # each 0.5 s run repeats the previous run's last time point as its first: drop the repeat
        # inside the step (the second run's first sample), keeping the switch sample
        n_run = (len(t) - starts[1]) // 2
        step[starts[1] + n_run] = False
        ts, xs = t[step], i[step]
        n = len(xs)
        a, b = int(0.9 * (n - 1)), int(1.0 * (n - 1))
        init, end = xs[0], xs[a:b].mean()
        ss_hold = i[starts[1] - 1]
        res.append({'label': label, 'rho': r, 'g': g_of(r), 'c': c, 't': t, 'i': i, 'w': o['w'], 'init': init, 'end': end,
                    'relax': end - init, 'ssdiff': end - ss_hold, 'a': ts[a], 'b': ts[b - 1]})

    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.6), sharex=True, gridspec_kw={'height_ratios': [3, 2]})
    ax = axes[0]
    for r in res:
        ax.plot(r['t'], r['i'], color=r['c'], lw=1.8, label=f"{r['label']}: g_M {r['g']:.4g} uS")
    r = res[0]
    ax.plot([0], [r['init']], 'o', color='#111827', ms=6, zorder=5)
    ax.annotate(f"first sample of the step\n(after the ohmic jump) {r['init']:.1f} pA", (0, r['init']), (0.12, r['init'] + 40),
                fontsize=8, arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.hlines(r['end'], r['a'], r['b'], color='#111827', lw=3, zorder=5)
    ax.annotate(f"mean of the last 10 %\n{r['end']:.1f} pA", (r['b'], r['end']), (0.55, r['end'] + 110), fontsize=8,
                arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.annotate('', (0.97, r['end']), (0.97, r['init']), arrowprops=dict(arrowstyle='<->', lw=1.2, color='#111827'))
    ax.text(0.955, (r['init'] + r['end']) / 2 - 10, f"relaxation\n{r['relax']:.1f} pA\n(data {1e3 * item['value']:.1f})",
            fontsize=8, ha='right', va='center')
    ax.axvspan(0, 1.0, color='#e5e7eb', zorder=0)
    ax.set_ylabel('i_M (pA)')
    ax.set_xlim(-0.3, 1.02)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title('held at -25 mV (10 s), 1 s step to -55 mV (shaded), E_K -96.5 mV')
    ax = axes[1]
    ax.plot(res[0]['t'], res[0]['w'], color='#111827', lw=1.6)
    ax.axhline(w_inf(STEP), color='#9ca3af', ls=':', lw=1)
    ax.text(0.5, w_inf(STEP) + 0.02, 'w_inf(-55 mV)', fontsize=8, color='#4b5563')
    ax.axvspan(0, 1.0, color='#e5e7eb', zorder=0)
    ax.set_ylabel('M gate w')
    ax.set_xlabel('time from the step (s)')
    fig_sim = _png(fig)

    rows_sim = ''.join(
        f"<tr><td>{r['label']}</td><td>{r['g']:.5g}</td><td>{r['rho']:.5g}</td><td>{r['relax']:.1f}</td>"
        f"<td>{r['ssdiff']:.1f}</td></tr>" for r in res)
    # analytic check and the tau_w dependence
    def relax_per_g(num, tmin):
        tt = np.arange(0, 1.0001, 0.001)
        tau = tau_w(STEP, num, tmin) / 1000
        ww = w_inf(STEP) + (w_inf(HOLD) - w_inf(STEP)) * np.exp(-tt / tau)
        n = len(tt)
        return (ww[int(0.9 * (n - 1)):int(n - 1)].mean() - ww[0]) * (STEP - E_K)
    g_an = item['value'] / relax_per_g(1200, 30)
    try:
        with open(os.path.join(IDIR, f'{PARAMETERISATION}_calibration.json')) as f:
            kin = json.load(f)['calibrated_parameters']
        g_kin = item['value'] / relax_per_g(float(kin['tau_w_num']), float(kin['tau_w_min']))
        kin_txt = (f"with the refitted tau_w of parameterisation {PARAMETERISATION} (tau_w_num {float(kin['tau_w_num']):.3g}, tau_w_min "
                   f"{float(kin['tau_w_min']):.4g} ms) the relaxation is complete within the step and g_M would be "
                   f"<b>{g_kin:.4g} uS</b>")
    except FileNotFoundError:
        kin_txt = 'the tau_w refit is not available'
    ss_per_g = w_inf(STEP) * (STEP - E_K) - w_inf(HOLD) * (HOLD - E_K)
    frac = (relax_per_g(1200, 30) / ((w_inf(STEP) - w_inf(HOLD)) * (STEP - E_K)))

    pi = obs['protocol_info']
    proto = (f"<tr><td>0</td><td>{pi['params_to_change']['mod/V'][0][0]:g} mV</td><td>{pi['sim_times'][0][0]:g} s</td><td>steady holding state (w = w_inf(-25 mV) = {w_inf(HOLD):.4f})</td></tr>"
             f"<tr><td>1</td><td>{pi['params_to_change']['mod/V'][0][1]:g} mV</td><td>{pi['sim_times'][0][1]:g} s</td><td>the deactivation step; the data item is read here</td></tr>")
    item_show = {k: v for k, v in item.items() if k != 'comment'}
    import html as _h
    q = lambda x: _h.escape(str(x))
    sem = np.std(WISTAR_DENSITIES, ddof=1) / np.sqrt(len(WISTAR_DENSITIES))
    body = f"""
<h1>i_M g_M: where it comes from</h1>
<p class="muted">i_M Argus2026_v01, parameterisation <code>{GM_PARAMETERISATION}</code>. g_M = A_mem rho_M gamma_M, with rho_M calibrated by libcuflynx to one
obs_data item. Rebuilt by <code>venv/bin/python tools/i_M_davis2020_kinetics.py gm-report</code>.</p>

<h2>1. The source</h2>
<div class="grid">
<div class="card"><img alt="Davis 2020 Fig. 3C" src="data:image/png;base64,{b64(marked)}">
<p class="muted">Davis et al. 2020, Hypertension 76:1915, Fig. 3C (the PDF raster, 200 ppi). Left: representative raw deactivation traces.
Right: the XE-991-sensitive M-current density per cell. Drawn on it, in orange: the 11 Wistar points as digitised (open circles), their median
({med:.2f} pA/pF, bar) and IQR ({q1:.2f} to {q3:.2f}, box). In blue: Davis's stated median -4.032 pA/pF.
The plotted black bars match the mean +- SEM of the points ({np.mean(WISTAR_DENSITIES):.2f} +- {sem:.2f}), not the median +- IQR the legend states.</p></div>
<div class="card"><img alt="Davis 2021 thesis pp. 172-173" src="data:image/png;base64,{b64(thesis)}">
<p class="muted">Davis 2021 DPhil thesis, Ch. 4, pp. 172-173 (the same data, Fig. 4.10C). It defines the measure: held at -25 mV, a 1 s step
to -55 mV, "the difference in current following the initial ohmic change, and the steady-state current towards the end of the step",
of the 10 uM XE-991-sensitive current. So the measure is the <b>deactivation relaxation</b>, not the steady-state difference between -25 and -55 mV.</p></div>
</div>

<h2>2. The obs_data</h2>
<p><code>parameterisations/{GM_PARAMETERISATION}/{GM_PARAMETERISATION}_obs_data.json</code>: the channel module alone, its boundary conditions set in the parameterisation
parameters (E_K -96.4979 mV: Nernst at Ko 4.5 / Ki 200 mM, 295.15 K; A_mem 3000 um^2; T 295.15 K). V is changed by <code>params_to_change</code>:</p>
<div class="grid"><div class="card"><table><tr><th>sub-experiment</th><th>V (mod/V)</th><th>duration</th><th></th></tr>{proto}</table></div>
<div class="card"><pre>{q(json.dumps(item_show, indent=1))}</pre></div></div>
<ul>
<li><b>Operation.</b> <code>mean_in_range_minus_initial(i_M, start_frac 0.9, end_frac 1.0)</code> on sub-experiment 1. It takes the mean of i_M over the last 10 %
of the step minus its first sample. The gate w is continuous and V jumps, so the first sample is the current just after the ohmic change. This is an existing libcuflynx operation.</li>
<li><b>Value.</b> -4.032 pA/pF (Wistar median, n = 11) x Cm {CM_PF:g} pF = -120.96 pA = <b>-0.12096 nA</b>.
The sign is negative because the outward current falls during the step.</li>
<li><b>std.</b> The SEM of the digitised points, {sem:.2f} pA/pF, x {CM_PF:g} pF = <b>0.0102 nA</b>.
Davis report a median +- IQR, with no SEM.</li>
</ul>
<div class="warn"><b>Assumptions.</b>
(1) <b>Cm = 30 pF</b> (the SN model's soma capacitance): Davis report densities only, so g_M scales with this choice.
(2) The perforated-patch <b>liquid junction potential of 24.3 mV is not corrected</b>: Davis 2019 bioRxiv Methods say it was not corrected,
so the true potentials may differ from the nominal -25 and -55 mV. The protocol uses the nominal values, as the source does.
(3) <b>Protocol voltages.</b> The bioRxiv preprint's legend says -20 to -50 mV; the published paper and the thesis say -25 to -55 mV. The published values are used.
(4) Only the channel is modelled. Davis's XE-991 subtraction removes the other currents, so the module alone is the right comparison.</div>

<h2>3. The simulation</h2>
<div class="grid"><div class="card"><img alt="simulation" src="data:image/png;base64,{fig_sim}"></div>
<div class="card"><table><tr><th>g_M</th><th>uS</th><th>rho_M (/um^2)</th><th>relaxation reading (pA)</th><th>end of step minus holding current (pA)</th></tr>{rows_sim}</table>
<p class="muted">CVODE simulation of the module, through the obs_data protocol. Data: -121.0 pA.
At a given g_M, the relaxation is {res[0]['relax'] / res[0]['ssdiff']:.2f} of the total change from the holding current to the end of the step. (The first fit read the data as the change between the two true steady states, w_inf(-55)(-55 - E_K) - w_inf(-25)(-25 - E_K) = {ss_per_g:.1f} mV x g_M, which gives -121 pA at 0.002556 uS by construction; within the 1 s step it reaches the value in the last column.)
The steady-state difference also includes the instantaneous ohmic drop of the M current (w(-25) times the 30 mV step), and Davis exclude that drop.
Reading the data as a steady-state difference (the first fit, 0.002556 uS) therefore gives a g_M {res[0]['g'] / res[1]['g']:.2f}x smaller than the relaxation reading.
The old default 0.00263 uS was a hand calculation of the same steady-state kind with E_K -93.7 mV.</p></div></div>

<h2>4. The fit</h2>
<div class="grid"><div class="card"><table>
<tr><th></th><th></th></tr>
<tr><td>method</td><td>libcuflynx {q(cal['method'])}, rho_M in [0.01, 1] /um^2, start 0.3</td></tr>
<tr><td>calibrated rho_M</td><td>{rho:.7g} /um^2</td></tr>
<tr><td>g_M = A rho gamma (A {A:g} um^2, gamma {gamma:g} pS)</td><td><b>{g_of(rho):.5g} uS</b></td></tr>
<tr><td>model item / data</td><td>{zitem['model']:.6g} / {zitem['data']:.6g} nA ({zitem['z']:.2g} std)</td></tr>
<tr><td>cost (gaussian_MLE)</td><td>{cal['cost']:.3g}</td></tr>
<tr><td>analytic check (exact w(t), same window)</td><td>{g_an:.5g} uS</td></tr>
<tr><td>date</td><td>{q(cal['date'])}</td></tr>
</table></div>
<div class="card"><p><b>Dependence on tau_w.</b> The reading is taken within Davis's 1 s step, so it depends on how far the relaxation has gone.
With the present tau_w (466 ms at -55 mV) it is {100 * frac:.0f} % complete, and g_M = <b>{g_an:.4g} uS</b>; {kin_txt}.
If tau_w changes, recalibrate this parameterisation.</p>
<p><b>Redo.</b> Edit <code>parameterisations/{GM_PARAMETERISATION}/{GM_PARAMETERISATION}_obs_data.json</code> (value, std, Cm, protocol), then:</p>
<pre>venv/bin/python -m pytest tests/test_modules.py --component i_M/Argus2026_v01 \\
    -k 'calibrate and {GM_PARAMETERISATION}' --include-unreviewed
venv/bin/python -m cam_testing.calibrate apply i_M/Argus2026_v01 {GM_PARAMETERISATION} --dry-run
venv/bin/python -m cam_testing.calibrate apply i_M/Argus2026_v01 {GM_PARAMETERISATION}
venv/bin/python tools/i_M_davis2020_kinetics.py gm-report</pre>
<p class="muted"><code>apply</code> writes rho_M to the i_M default, soma sympathetic (rho_M_i_M), the monolithic soma (rho_M) and SN_simple / SN_simple_flat.</p></div></div>
"""
    os.makedirs(os.path.dirname(GM_HTML), exist_ok=True)
    with open(GM_HTML, 'w') as f:
        f.write(HTML_PAGE.format(title='i_M g_M fit', body=body))
    print('wrote', GM_HTML)
    for r in res:
        print(r['label'], r['g'], r['relax'], r['ssdiff'])
    print('analytic', g_an, 'median/IQR', med, q1, q3)


HTML_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
:root {{ --bg:#ffffff; --fg:#111827; --muted:#4b5563; --card:#f9fafb; --line:#e5e7eb; --accent:#d9480f; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; }} }}
:root[data-theme="dark"] {{ --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; }}
body {{ background:var(--bg); color:var(--fg); font:15px/1.5 system-ui, -apple-system, Segoe UI, sans-serif; margin:0; padding:24px 16px; }}
main {{ max-width:1180px; margin:0 auto; }}
h1 {{ font-size:1.5rem; margin:0 0 4px; }} h2 {{ font-size:1.15rem; margin:28px 0 8px; border-bottom:1px solid var(--line); padding-bottom:4px; }}
.muted {{ color:var(--muted); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit, minmax(340px, 1fr)); gap:16px; align-items:start; }}
.card {{ background:var(--card); border:1px solid var(--line); border-radius:8px; padding:12px; overflow-x:auto; }}
.card img {{ width:100%; height:auto; background:#fff; border-radius:4px; display:block; }}
table {{ border-collapse:collapse; width:100%; font-size:14px; }} th, td {{ border-bottom:1px solid var(--line); padding:4px 8px; text-align:left; }}
th {{ color:var(--muted); font-weight:600; }}
.warn {{ border-left:4px solid var(--accent); padding:8px 12px; background:var(--card); }}
pre {{ font-size:12px; white-space:pre-wrap; word-break:break-word; margin:0; }} code {{ font-size:13px; }}
</style></head><body><main>
{body}
</main></body></html>
"""

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('what', choices=['obs-data', 'report', 'gm-report'])
    a = ap.parse_args()
    if a.what == 'obs-data':
        build_obs_data()
    elif a.what == 'gm-report':
        build_gm_report()
    else:
        build_report()
