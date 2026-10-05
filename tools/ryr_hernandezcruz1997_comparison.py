"""
RyR HernandezCruz1997_v01 vs Argus2026_v01 in the soma supermodule (soma sympathetic).

    venv/bin/python tools/ryr_hernandezcruz1997_comparison.py    # -> reviews/RyR_HernandezCruz1997_comparison.html

The soma is not changed. The script copies modules/ to a temporary directory, switches the RyR
submodule of soma/sympathetic to HernandezCruz1997_v01 there (and drops the soma instance's
Argus-only RyR rows, so the RyR takes its own default instance), and points the library and
libcuflynx at the copy to generate the second model. The current model is generated from the
repository itself.

Protocols (each model): 30 s at 0 pA from the model's initial state (rest); from the end of the
rest, a 1 s step of 5, 9 or 50 pA followed by 2 s at 0 pA. Caffeine (HernandezCruz1997, Eq. A16):
the activating-site on-rate kon_a_RyR switched from 10 to 2500 /(mM s) (Kd 200 -> 0.8 uM) at the
end of the rest, 30 s at 0 pA; then the same steps from the end of that caffeine exposure.

CICR measure: the RyR release during the step and the 2 s after it, beyond the resting release
(integral of j_RyR - j_RyR_rest), divided by the Ca that entered through L-type channels in the
same window (integral of j_Ca_L): the release gain.
"""
import base64
import contextlib
import csv
import functools
import html
import io
import json
import os
import shutil
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
OUT = os.path.join(REPO, 'reviews', 'RyR_HernandezCruz1997_comparison.html')
SOMA_REL = os.path.join('cell', 'neuron', 'soma', 'versions', 'sympathetic')
NEW = 'HernandezCruz1997_v01'
DT = 1e-4
REST_S = 30
LONG_REST_S = 180
STEPS_PA = (5, 9, 50)
STEP_S, POST_S = 1, 2
ENTRY_MIN = 1e-18      # mol (1 amol): below it the gain is not computed (no Ca entry)
KON_CAFF = 2500.0          # /(mM s): HernandezCruz1997 k_on+Caff 2.5 /uM/s

V = 'mod_membrane/V'
CAI = 'mod_Ca/Cai'
CAB = 'mod_Ca/Ca_bulk'
CAER = 'mod_Ca/Ca_ER'
JRYR = 'mod_RyR/j_RyR'
JCAL = 'mod_i_CaL/j_Ca_L'
HC_GATES = ['mod_RyR/a_RyR', 'mod_RyR/i_RyR', 'mod_RyR/l_RyR', 'mod_RyR/P_RyR']
COLORS = {'current': '#2a78d6', 'HC1997': '#eb6834', 'HC1997 + caffeine': '#1baf7a'}


def copy_library(tmp):
    """modules/ without generated files, with soma/sympathetic's RyR switched to the new version."""
    dst = os.path.join(tmp, 'modules')
    shutil.copytree(os.path.join(REPO, 'modules'), dst,
                    ignore=shutil.ignore_patterns('results', 'plots', 'risk', '*.omex', '*.html'))
    sdir = os.path.join(dst, SOMA_REL)
    cfg_path = os.path.join(sdir, 'soma_sympathetic_modules_config.json')
    with open(cfg_path) as f:
        cfg = json.load(f)
    n = 0
    for s in cfg[0]['submodules']:
        if s['name'] == 'RyR':
            s['module_subtype'] = NEW
            n += 1
    assert n == 1, 'RyR submodule not found'
    with open(cfg_path, 'w') as f:
        json.dump(cfg, f, indent=2)
    csv_path = os.path.join(sdir, 'instances', 'default', 'default_parameters.csv')
    with open(csv_path, newline='') as f:
        rows = list(csv.reader(f))
    keep = [r for r in rows if not (r and r[0].endswith('_RyR') and r[0] != 'variable_name')]
    dropped = [r[0] for r in rows if r not in keep]
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerows(keep)
    return dst, dropped


def use_library(path):
    from cam_testing import harness, library
    library.MODULES_DIR = path
    harness.MODULES_DIR = path
    for name in dir(library):
        obj = getattr(library, name)
        if isinstance(obj, functools._lru_cache_wrapper):
            obj.cache_clear()


class Soma:
    def __init__(self, label, work_dir):
        from cam_testing import checks, harness
        from cam_testing.library import load_version
        sv = load_version('soma', 'sympathetic')
        self.cm = checks.ComponentModel(sv, work_dir=work_dir)
        by = {p.variable_name: p for p in self.cm.parameters()}
        self.pn = lambda v: 'parameters/' + harness.parameter_name(by[v])
        self.has = lambda v: v in by
        self.h = harness.simulation_helper(self.cm.model_path, dict(self.cm.spec, sim_time=1.0, dt=DT))
        self.label = label
        self.ryr = [n.split('.', 1)[1] for n in self.h.get_all_variable_names() if n.startswith('mod_RyR_module.')]
        self.outputs = [V, CAI, CAB, CAER, JRYR, JCAL] + [g for g in HC_GATES if g.split('/')[1] in self.ryr]

    def segment(self, params, seconds):
        from cam_testing import harness
        self.h.set_param_vals([self.pn(k) for k in params], list(params.values()), change_states=False)
        ts, outs = [], {o: [] for o in self.outputs}
        for _ in range(seconds):
            with contextlib.redirect_stdout(io.StringIO()):
                if self.h.run() is False:
                    raise RuntimeError(f'{self.label}: simulation failed')
            t = np.asarray(self.h.get_time(), dtype=float)
            ts.append(t - t[0] + (ts[-1][-1] + DT if ts else 0.0))
            for o in self.outputs:
                outs[o].append(np.asarray(self.h.get_results([[harness.output_name(o)]], flatten=True)[0], dtype=float).ravel())
        return np.concatenate(ts), {o: np.concatenate(v) for o, v in outs.items()}

    def state(self):
        return list(self.h.simulation.state())

    def set_state(self, s):
        self.h.simulation.set_state(s)


def spikes(v, mask=None):
    up = (v[:-1] < 0) & (v[1:] >= 0)
    if mask is not None:
        up &= mask[1:]
    return int(np.sum(up))


def steps(soma, start_state, base, j_rest):
    res = {}
    for pa in STEPS_PA:
        soma.set_state(start_state)
        t1, o1 = soma.segment(dict(base, I_in_membrane=pa * 1e-3), STEP_S)
        t2, o2 = soma.segment(dict(base, I_in_membrane=0.0), POST_S)
        t = np.concatenate([t1, t2 + t1[-1] + DT])
        o = {k: np.concatenate([o1[k], o2[k]]) for k in o1}
        rel = np.trapezoid(o[JRYR] - j_rest, t)
        ent = np.trapezoid(np.abs(o[JCAL]), t)
        res[pa] = dict(t=t, o=o, spikes=spikes(o[V], t < STEP_S), gain=rel / ent if ent > ENTRY_MIN else float('nan'),
                       p_ratio=(float(np.max(o['mod_RyR/P_RyR']) / o['mod_RyR/P_RyR'][0]) if 'mod_RyR/P_RyR' in o else float('nan')),
                       release=rel, entry=ent, j_peak=float(np.max(o[JRYR])), cai_peak=float(np.max(o[CAI])),
                       er0=float(o[CAER][0]), er_min=float(np.min(o[CAER])), er_end=float(o[CAER][-1]))
    return res


def simulate(soma, caffeine):
    base = {}
    t, o = soma.segment(dict(I_in_membrane=0.0), REST_S)
    rest = dict(t=t, o=o, spikes=spikes(o[V]), spikes_last10=spikes(o[V], t >= REST_S - 10))
    st = soma.state()
    j_rest = float(np.mean(o[JRYR][t >= REST_S - 1]))
    rest['j_rest'] = j_rest
    out = {'rest': rest, 'steps': steps(soma, st, base, j_rest)}
    # a longer rest (continuing from 30 s to LONG_REST_S), to see where the rest goes
    soma.set_state(st)
    tl, ol = soma.segment(dict(I_in_membrane=0.0), LONG_REST_S - REST_S)
    out['long_rest'] = dict(spikes=spikes(ol[V]), spikes_last30=spikes(ol[V], tl >= tl[-1] - 30), v=float(ol[V][-1]),
                            cai=float(ol[CAI][-1]), er=float(ol[CAER][-1]), j=float(ol[JRYR][-1]))
    if caffeine:
        soma.set_state(st)
        caff = dict(kon_a_RyR_RyR=KON_CAFF)
        tc, oc = soma.segment(dict(caff, I_in_membrane=0.0), REST_S)
        out['caffeine'] = dict(t=tc, o=oc, spikes=spikes(oc[V]), er_min=float(np.min(oc[CAER])),
                               j_peak=float(np.max(oc[JRYR])), cai_peak=float(np.max(oc[CAI])),
                               t_peak=float(tc[np.argmax(oc[JRYR])]))
        stc = soma.state()
        jc = float(np.mean(oc[JRYR][tc >= REST_S - 1]))
        out['caffeine_steps'] = steps(soma, stc, caff, jc)
        out['caffeine']['j_end'] = jc
    return out


# ---- page --------------------------------------------------------------------------------------

def _png(fig):
    import matplotlib.pyplot as plt
    b = io.BytesIO()
    fig.savefig(b, format='png', dpi=100, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def trace_fig(series, title, t_mark=None):
    """series: [(label, t, o)]; rows V, Cai (shell), Ca_ER, j_RyR."""
    import matplotlib.pyplot as plt
    rows = [(V, 'V (mV)', 1), (CAI, 'Cai shell (uM)', 1e3), (CAER, 'Ca_ER (mM)', 1), (JRYR, 'j_RyR (amol/s)', 1e18)]
    fig, axes = plt.subplots(len(rows), 1, figsize=(7.4, 7.6), sharex=True)
    for ax, (k, yl, sc) in zip(axes, rows):
        for label, t, o in series:
            ax.plot(t, o[k] * sc, color=COLORS[label], lw=0.8, label=label)
        ax.set_ylabel(yl, fontsize=9)
        ax.grid(alpha=0.25)
        ax.tick_params(labelsize=8)
        if t_mark is not None:
            ax.axvspan(*t_mark, color='#d4d4d4', alpha=0.35, lw=0)
    axes[0].set_title(title, fontsize=10)
    axes[0].legend(fontsize=8, loc='upper right')
    axes[-1].set_xlabel('time (s)')
    fig.tight_layout()
    return _png(fig)


def gates_fig(series, title):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(7.4, 2.8))
    styles = {'a_RyR': '-', 'i_RyR': '--', 'l_RyR': ':', 'P_RyR': '-.'}
    for label, t, o in series:
        for g in HC_GATES:
            if g in o:
                n = g.split('/')[1]
                y = o[g] * (1e3 if n == 'P_RyR' else 1)
                ax.plot(t, y, styles[n], color=COLORS[label], lw=1, label=f'{label}: {n}' + (' x1000' if n == 'P_RyR' else ''))
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('occupancy')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    return _png(fig)


def fmt(x, p=3):
    return 'n/a' if x != x else f'{x:.{p}g}'


# Reading of the results (written after running this script on 2026-10-05; the numbers above it are recomputed on each run)
FINDINGS = [
    'With the current RyR, each spontaneous spike at rest coincides with a store release: Ca_ER climbs to about 0.55 mM '
    '(past Ca_ER_RyR_mid 0.5 mM), then j_RyR jumps to 40,000-85,000 amol/s and Ca_ER falls to about 0.07 mM. The 0.2 Hz rest '
    'firing of the soma appears to be paced by these luminal-triggered releases. With HC1997 there are no such events: the soma '
    'is silent over the 30 s rest (one spike between 30 and 180 s while it slowly depolarises) and settles near -68 mV.',
    'Without caffeine, HC1997 gives no CICR from influx. At 5 and 9 pA there is no spike and no extra release. At 50 pA the single '
    'spike raises the shell Ca to about 3 uM, P_RyR rises less than 2-fold and j_RyR peaks at about 160 amol/s. The release follows '
    'Cai smoothly, with no regenerative surge: the activating site (Kd 200 uM) barely responds to a few uM, and the inactivating site '
    '(Kd 0.63 uM) occupies about half the channels at the soma\'s resting Ca. This matches the paper\'s finding that influx alone does not '
    'trigger CICR in rat SCG neurons. The release gain of about 1.9 is graded leak, not CICR. It is large because the soma\'s ER holds about 0.4 mM: '
    'the luminal site is about 76 % occupied (l^7 about 0.15), versus 37 % (l^7 about 9e-4) at the paper\'s resting 75 uM, '
    'so the luminal gate is about 170 times more open than in the source model.',
    'Caffeine (kon_a 2500 /(mM s)) gives the paper\'s two components. A transient release (TR) peaks at about 850 amol/s within 60 ms and '
    'empties the store from 0.41 to about 0.12 mM in about 0.5 s. A persistent release (PR) then decays over tens of seconds as the stores run down '
    '(cf. Fig. 1 and Fig. 10). The soma does not fire: K(Ca) hyperpolarises it. After 30 s in caffeine the store is depleted (0.1 mM), and the steps '
    'release about 2.6 times the L-type entry, again without a regenerative surge.',
    'The current model does not fire during the 50 pA step either. The step starts 1.3 s after a release-driven spike, and the current RyR keeps the '
    'shell Ca near 1 uM at rest (j_RyR about 220 amol/s), so K(Ca) holds the membrane down. HC1997 lowers the resting shell Ca to 0.36-0.43 uM, '
    'still above the measured 102 nM (Wanaverbecq2003).',
    'Caveats. (1) The flux magnitude of HC1997 (rho_RyR 5403 /um^2) is a conversion of the paper\'s RyR_max on its cytosolic shell, not a '
    'measured density. (2) The soma\'s ER Ca (about 0.4 mM) is about 5 times the paper\'s 75 uM, which opens the luminal gate far more than in the '
    'source model. (3) HC1997\'s rest is not at steady state at 30 s (V is still drifting), so the steps start from a slowly changing state. '
    '(4) The rest of the soma (SERCA, buffers, K(Ca)) was tuned with the current RyR in place, so any switch would need recalibration through the soma\'s '
    'calibration workflow.',
]


def build_page(R, dropped):
    cur, hc = R['current'], R['HC1997']
    figs = []
    figs.append(('Rest, 0 pA, 30 s', trace_fig([('current', cur['rest']['t'], cur['rest']['o']),
                                                  ('HC1997', hc['rest']['t'], hc['rest']['o'])], 'rest, 0 pA, 30 s')))
    figs.append(('', gates_fig([('HC1997', hc['rest']['t'], hc['rest']['o'])], 'HC1997 site occupancies at rest')))
    for pa in STEPS_PA:
        s = [('current', cur['steps'][pa]['t'], cur['steps'][pa]['o']), ('HC1997', hc['steps'][pa]['t'], hc['steps'][pa]['o']),
             ('HC1997 + caffeine', hc['caffeine_steps'][pa]['t'], hc['caffeine_steps'][pa]['o'])]
        figs.append((f'{pa} pA step', trace_fig(s, f'{pa} pA for 1 s (shaded), then 2 s at 0 pA; from the end of the rest '
                                                  f'(caffeine: from the end of 30 s in caffeine)', (0, STEP_S))))
        figs.append(('', gates_fig(s[1:], f'{pa} pA: HC1997 site occupancies')))
    c = hc['caffeine']
    figs.append(('Caffeine at rest', trace_fig([('HC1997 + caffeine', c['t'], c['o'])],
                                               'caffeine applied at t = 0 after the 30 s rest (kon_a 10 -> 2500 /(mM s)), 0 pA')))
    figs.append(('', gates_fig([('HC1997 + caffeine', c['t'], c['o'])], 'caffeine: HC1997 site occupancies')))

    def last(o, k, t, w=1.0):
        return float(np.mean(o[k][t >= t[-1] - w]))

    rest_rows = []
    for label, r in (('current', cur['rest']), ('HC1997', hc['rest'])):
        t, o = r['t'], r['o']
        rest_rows.append(f"<tr><td>{label}</td><td>{r['spikes']}</td><td>{r['spikes_last10']}</td>"
                         f"<td>{fmt(last(o, V, t))}</td><td>{fmt(last(o, CAI, t) * 1e6)}</td><td>{fmt(last(o, CAER, t))}</td>"
                         f"<td>{fmt(r['j_rest'] * 1e18)}</td></tr>")
    step_rows = []
    for pa in STEPS_PA:
        for label, s in (('current', cur['steps'][pa]), ('HC1997', hc['steps'][pa]), ('HC1997 + caffeine', hc['caffeine_steps'][pa])):
            step_rows.append(f"<tr><td>{pa}</td><td>{label}</td><td>{s['spikes']}</td><td>{fmt(s['cai_peak'] * 1e3)}</td>"
                             f"<td>{fmt(s['er0'])} / {fmt(s['er_min'])} / {fmt(s['er_end'])}</td><td>{fmt(s['j_peak'] * 1e18)}</td>"
                             f"<td>{fmt(s['release'] * 1e18)}</td><td>{fmt(s['entry'] * 1e18)}</td><td>{fmt(s['gain'])}</td><td>{fmt(s['p_ratio'])}</td></tr>")
    summary = R['summary']
    body = []
    body.append('<h1>RyR: Hern&aacute;ndez-Cruz 1997 vs Argus2026_v01 in the soma</h1>')
    body.append('<p class="sub">soma <code>sympathetic</code> with its RyR submodule switched to <code>HernandezCruz1997_v01</code> '
                'in a temporary copy of the library (the soma is unchanged in the repository). Generated by '
                '<code>tools/ryr_hernandezcruz1997_comparison.py</code>.</p>')
    body.append('<h2>Summary (computed)</h2><ul>' + ''.join(f'<li>{s}</li>' for s in summary) + '</ul>')
    body.append('<h2>Reading</h2><ul>' + ''.join(f'<li>{html.escape(s)}</li>' for s in FINDINGS) + '</ul>')
    body.append('<h2>Setup</h2><ul>'
                '<li><b>current</b>: RyR Argus2026_v01 with the soma instance\'s values.</li>'
                '<li><b>HC1997</b>: RyR HernandezCruz1997_v01 at its default instance (activating Kd 200 uM, inactivating 0.63 uM, '
                'luminal 130 uM with N = 7; k_RyR = A_ER rho k_single, rho 5403 /um^2, k_single 2.591 um^3/s; T_ref 295.65 K). '
                'The soma instance rows dropped in the copy (Argus-only parameters): ' + ', '.join(f'<code>{html.escape(d)}</code>' for d in dropped) + '.</li>'
                '<li><b>HC1997 + caffeine</b>: kon_a_RyR = 2500 /(mM s) (k_on+Caff, Kd 0.8 uM), the source\'s representation of 10 mM caffeine.</li>'
                f'<li>Rest: 30 s at 0 pA from the model\'s initial state. Steps: from the end of the rest (caffeine: from the end of 30 s in caffeine), '
                f'{STEP_S} s at 5, 9, 50 pA then {POST_S} s at 0 pA. Spikes: upward crossings of 0 mV during the step.</li>'
                '<li>Release gain = (integral of j_RyR - resting j_RyR) / (integral of |j_Ca_L|) over the step and the 2 s after it: '
                'Ca released by the RyR per Ca entered through L-type channels (n/a when less than 1 amol entered, i.e. no spike). '
                'A gain of order 1 or more means the release amplifies the influx; whether it is regenerative CICR shows in P_RyR '
                '(a regenerative release keeps rising after the trigger) and in the occupancy plots.</li></ul>')
    body.append('<h2>Rest</h2><table><tr><th>model</th><th>spikes in 30 s</th><th>in last 10 s</th><th>V (mV)</th>'
                '<th>Cai shell (nM)</th><th>Ca_ER (mM)</th><th>j_RyR (amol/s)</th></tr>' + ''.join(rest_rows) + '</table>'
                '<p class="note">Values over the last second of the rest.</p>')
    body.append('<h2>Steps</h2><table><tr><th>pA</th><th>model</th><th>spikes</th><th>peak Cai shell (uM)</th>'
                '<th>Ca_ER start / min / end (mM)</th><th>peak j_RyR (amol/s)</th><th>RyR release above rest (amol)</th>'
                '<th>L-type entry (amol)</th><th>release gain</th><th>peak P_RyR / P_RyR at step onset (HC1997)</th></tr>' + ''.join(step_rows) + '</table>')
    body.append(f"<h2>Caffeine at rest</h2><p>HC1997 with caffeine applied after the 30 s rest: peak j_RyR {fmt(c['j_peak'] * 1e18)} amol/s "
                f"at {fmt(c['t_peak'])} s, peak Cai shell {fmt(c['cai_peak'] * 1e3)} uM, Ca_ER minimum {fmt(c['er_min'])} mM, "
                f"{c['spikes']} spike(s) in 30 s; j_RyR at the end {fmt(c['j_end'] * 1e18)} amol/s.</p>")
    body.append('<h2>Traces</h2>')
    for h, png in figs:
        if h:
            body.append(f'<h3>{html.escape(h)}</h3>')
        body.append(f'<img alt="{html.escape(h or "occupancies")}" src="data:image/png;base64,{png}">')
    css = ('body{font-family:system-ui,sans-serif;max-width:980px;margin:0 auto;padding:16px;color:#111827;background:#fff;line-height:1.45}'
           'table{border-collapse:collapse;font-size:13px;margin:8px 0;display:block;overflow-x:auto}'
           'td,th{border:1px solid #d1d5db;padding:3px 7px;text-align:left}th{background:#f3f4f6}'
           'img{max-width:100%;display:block;margin:6px 0 14px}.sub,.note{color:#4b5563;font-size:14px}code{font-size:90%}')
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>RyR HC1997 comparison</title><style>{css}</style></head><body>{"".join(body)}</body></html>')
    with open(OUT, 'w') as f:
        f.write(page)


def summarise(R):
    cur, hc = R['current'], R['HC1997']
    s = []
    s.append(f"Rest (0 pA, 30 s): current {cur['rest']['spikes']} spike(s), HC1997 {hc['rest']['spikes']}; resting Ca_ER "
             f"{fmt(float(cur['rest']['o'][CAER][-1]))} mM (current) vs {fmt(float(hc['rest']['o'][CAER][-1]))} mM (HC1997); "
             f"resting j_RyR {fmt(cur['rest']['j_rest'] * 1e18)} vs {fmt(hc['rest']['j_rest'] * 1e18)} amol/s.")
    for pa in STEPS_PA:
        a, b, c = cur['steps'][pa], hc['steps'][pa], hc['caffeine_steps'][pa]
        s.append(f"{pa} pA: spikes current {a['spikes']}, HC1997 {b['spikes']}, HC1997 + caffeine {c['spikes']}; release gain "
                 f"{fmt(a['gain'])} / {fmt(b['gain'])} / {fmt(c['gain'])}; Ca_ER minimum {fmt(a['er_min'])} / {fmt(b['er_min'])} / {fmt(c['er_min'])} mM.")
    for label, r in (('current', cur), ('HC1997', hc)):
        lr = r['long_rest']
        s.append(f"Long rest ({label}, 0 pA continued to {LONG_REST_S} s): {lr['spikes']} spike(s) in 30-{LONG_REST_S} s "
                 f"({lr['spikes_last30']} in the last 30 s); at {LONG_REST_S} s V {fmt(lr['v'])} mV, Cai shell {fmt(lr['cai'] * 1e6)} nM, "
                 f"Ca_ER {fmt(lr['er'])} mM, j_RyR {fmt(lr['j'] * 1e18)} amol/s.")
    c = hc['caffeine']
    s.append(f"Caffeine at rest (kon_a 2500 /(mM s)): release peaks at {fmt(c['j_peak'] * 1e18)} amol/s after {fmt(c['t_peak'])} s, "
             f"Ca_ER falls to {fmt(c['er_min'])} mM, Cai shell peaks at {fmt(c['cai_peak'] * 1e3)} uM, {c['spikes']} spike(s).")
    return s


def main():
    tmp = tempfile.mkdtemp(prefix='ryr_hc97_')
    try:
        R = {}
        R['current'] = simulate(Soma('current', os.path.join(tmp, 'work_current')), caffeine=False)
        lib, dropped = copy_library(tmp)
        use_library(lib)
        R['HC1997'] = simulate(Soma('HC1997', os.path.join(tmp, 'work_hc')), caffeine=True)
        R['summary'] = summarise(R)
        for line in R['summary']:
            print(html.unescape(line))
        build_page(R, dropped)
        print('wrote', os.path.relpath(OUT, REPO))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()
