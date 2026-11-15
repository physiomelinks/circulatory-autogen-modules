"""
The cardiac clocks' phase computed from the beat phase s (docs/heart/cardiac_clock_phase.md).

Every heart and chamber version with a chamber phase clock (a component variable chi_afloor, chi_vfloor
or chi_floor) has it in one of two forms:
  integrated  chi is a state with the piecewise rate of the original law (start window, wait at the integer)
  phase       chi_floor = P(s), the closed form, and chi = floor(s - t_s) + P(s), computed from the beat phase s
The tests:
  - the closed form P(s) (written here as the document writes it) equals the original law integrated
    exactly, for constant periods and for discrete period changes at any point of the beat, and the output
    chi equals the integrated state, across the onsets;
  - each version's own clock equations, whichever form, give P(s) and that chi;
  - the converter between the forms reproduces each version's file (round trip);
  - each version family simulated in both forms gives the same trajectories (test_forms_agree);
  - the atrial onset does not move with the solver tolerance or a 1e-12 change of the initial state
    (test_atrial_onset_ignores_round_off). The integrated form fails this, by round-off, in about one run
    in ten (the document's section 6); it is expected to fail while a version still integrates chi.
"""
import contextlib
import glob
import io
import math
import os
import re
import shutil

import numpy as np
import pytest
from lxml import etree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, 'modules')
MATHML = 'http://www.w3.org/1998/Math/MathML'
M = '{%s}' % MATHML

# one entry per clock: its variables, and the original law's comparison at 1/4 and 1/2 (strict: > for relaxation
# and rest, the flat hearts' ventricle; otherwise >=)
CLOCKS = {
    'chi_a': dict(floor='chi_afloor', final='chi_afloor_final', ts='t_astart_norm', ts_raw='t_astart', tc='T_ac', tr='T_ar',
                  phi='phi_a', strict=False, eff=('T_ac_wCont', 'T_ar_wCont', 'T_wCont')),
    'chi_v': dict(floor='chi_vfloor', final='chi_vfloor_final', ts='t_vstart_norm', ts_raw='t_vstart', tc='T_vc', tr='T_vr',
                  phi='phi_v', strict=True, eff=('T_vc_wCont', 'T_vr_wCont', 'T_wCont')),
    'chi': dict(floor='chi_floor', final='chi_floor_final', ts='t_start_norm', ts_raw='t_start', tc='T_c', tr='T_r',
                phi='phi', strict=False, eff=('T_c_eff', 'T_r_eff', 'T_eff')),
}
EPS_1, EPS_2 = 0.07, 0.02          # the original law's start window and dead band (fractions of a beat / of a cycle)


# ---- the maths, as the document writes it ---------------------------------------------------------------

def P(s, T, Tc, Tr, ts):
    """The closed form: the wrapped phase as a function of the beat phase s (ts = t_start / T)."""
    s = np.asarray(s, float)
    p = s - ts
    phi = p - np.floor(p)
    x = phi * T
    out = np.where(x <= Tc, 0.25 * x / Tc,
                   np.where(x < Tc + Tr, 0.25 + 0.25 * (x - Tc) / Tr, 0.5 + 0.5 * (x - Tc - Tr) / (T - Tc - Tr)))
    return np.where(s < ts, 0.0, out)


def chi_out(s, T, Tc, Tr, ts):
    """The output chi: 0 before the first onset, then floor(s - ts) + P(s)."""
    s = np.asarray(s, float)
    return np.where(s < ts, 0.0, np.floor(s - ts) + P(s, T, Tc, Tr, ts))


def original_rate(T, Tc, Tr, ts, strict):
    """d chi/dt of the original law, as a function of (psi, mt, T_w): the pieces in order, first match wins."""
    def rate(psi, m, Tw):
        Tcw, Trw = Tc * Tw / T, Tr * Tw / T
        if m >= ts and m <= ts + EPS_1 and psi <= 0.25:
            return 0.25 / Tcw
        if ts + EPS_1 > 1.0 and m <= ts + EPS_1 - 1.0 and psi <= 0.25:
            return 0.25 / Tcw
        if psi > EPS_2 and psi <= 0.25:
            return 0.25 / Tcw
        if (psi > 0.25 if strict else psi >= 0.25) and psi < 0.5:
            return 0.25 / Trw
        if (psi > 0.5 if strict else psi >= 0.5):
            return 0.5 / (Tw - Tcw - Trw)
        return 0.0
    return rate


def integrate(rate, schedule, t_end, sample_t, chi0=0.0, s0=0.0):
    """
    The law d chi/dt = rate(chi - floor(chi), s - floor(s), T_w), ds/dt = 1/T_w, integrated exactly: the rate is
    constant between switches, so each step is exact and each switch is located by bisection (to 1e-15 s).
    schedule: [(t_i, T_i)], T_w = T_i from t_i on (t_0 = 0). Returns (s, chi) at the times sample_t.
    """
    t, s, chi = 0.0, s0, chi0
    times = sorted(set([float(x) for x in sample_t] + [t_ for t_, _ in schedule if 0 < t_ < t_end] + [t_end]))
    out_s, out_chi, k = {}, {}, 0
    Tw_at = lambda tt: [T_ for t_, T_ in schedule if t_ <= tt + 1e-15][-1]

    def r_at(tt, ss, cc, Tw):
        return rate(cc - math.floor(cc), ss - math.floor(ss), Tw)

    for target in times:
        while t < target:
            Tw = Tw_at(t)
            h = min(0.01, target - t)
            r0 = r_at(t, s, chi, Tw)
            r1 = r_at(t + h, s + h / Tw, chi + r0 * h, Tw)
            if r1 != r0:
                lo, hi = 0.0, h                 # the rate switches inside (0, h]: find where
                while hi - lo > 1e-15:
                    mid = 0.5 * (lo + hi)
                    if r_at(t + mid, s + mid / Tw, chi + r0 * mid, Tw) == r0:
                        lo = mid
                    else:
                        hi = mid
                h = hi
                # the strict comparisons give a rate of 0 at the single point psi = 1/2 (docs section 1): the bisection
                # lands on it exactly, which a solver never does, so step 1e-12 s through it
                if r_at(t + h, s + h / Tw, chi + r0 * h, Tw) == 0.0 and r0 > 0.0 and \
                        abs((chi + r0 * h) - math.floor(chi + r0 * h) - 0.5) < 1e-12:
                    h += 1e-12
            t, s, chi = t + h, s + h / Tw, chi + r0 * h
        out_s[target], out_chi[target] = s, chi
    return np.array([out_s[float(x)] for x in sample_t]), np.array([out_chi[float(x)] for x in sample_t])


# the clocks the library's versions use (T = 1 s): atrium onset at mt 0.8, ventricle at 0
TIMINGS = {'atrium': dict(T=1.0, Tc=0.17, Tr=0.17, ts=0.8), 'ventricle': dict(T=1.0, Tc=0.34, Tr=0.15, ts=0.0),
           'chamber_lv': dict(T=1.0, Tc=0.3, Tr=0.15, ts=0.0), 'other': dict(T=0.9, Tc=0.2, Tr=0.25, ts=0.35)}

# discrete changes of the period: a step 1.0 -> 0.6 s or 0.6 -> 1.2 s at a point of the beat (beat k, phase mt),
# a ramp (as steps of 10 ms) and a period held per beat (docs: section 6, "Discrete period changes")
PLACES = {'wait': (0, 0.5), 'wait_pre_onset': (0, 0.799), 'rest_pre_onset': (1, 0.799), 'window': (1, 0.83),
          'a_contraction': (1, 0.9), 'a_relaxation': (2, 0.05), 'a_rest': (1, 0.5), 'v_systole': (1, 0.2)}


def period_scenarios():
    sc = {'constant': [(0.0, 1.0)]}
    for d, (ta, tb) in {'down': (1.0, 0.6), 'up': (0.6, 1.2)}.items():
        for p, (k, mt) in PLACES.items():
            sc[f'{d}_{p}'] = [(0.0, ta), ((k + mt) * ta, tb)]
    sc['ramp_down'] = [(0.0, 1.0)] + [(1.85 + 0.01 * i, 1.0 - 0.4 * (i + 1) / 30) for i in range(30)]
    seq, t = [], 0.0
    for T_ in (1.0, 0.6, 1.2, 0.7, 1.1, 0.6, 1.0):     # held for one beat each
        seq.append((t, T_)); t += T_
    sc['held_sequence'] = seq
    return sc


SCENARIOS = period_scenarios()


@pytest.mark.parametrize('timing', sorted(TIMINGS))
@pytest.mark.parametrize('strict', [False, True], ids=['ge', 'strict'])
@pytest.mark.parametrize('scenario', sorted(SCENARIOS))
def test_closed_form_is_the_integrated_law(timing, strict, scenario):
    """docs section 3: P(s) is the original law integrated from chi = 0, s = 0, whatever the period does;
    section 4: the output chi is the integrated state, also just before the first onset and at each onset."""
    c = TIMINGS[timing]
    t_end = 5.0
    sample = np.linspace(0.0, t_end, 2001)
    s, chi = integrate(original_rate(c['T'], c['Tc'], c['Tr'], c['ts'], strict), SCENARIOS[scenario], t_end, sample)
    psi = chi - np.floor(chi)
    ref = P(s, c['T'], c['Tc'], c['Tr'], c['ts'])
    d_psi = (psi - ref + 0.5) % 1.0 - 0.5                 # psi wraps: compare modulo 1
    assert np.abs(d_psi).max() < 1e-9, (np.abs(d_psi).max(), s[np.argmax(np.abs(d_psi))])
    assert np.abs(chi - chi_out(s, c['T'], c['Tc'], c['Tr'], c['ts'])).max() < 1e-9


@pytest.mark.parametrize('timing', ['atrium', 'other'])
def test_chi_output_at_the_onsets(timing):
    """The output chi equals the integrated state at s just below and at each onset k + ts, and just before the
    first one (docs section 4)."""
    c = TIMINGS[timing]
    T = c['T']
    rate = original_rate(T, c['Tc'], c['Tr'], c['ts'], False)
    onsets = [(k + c['ts']) * T for k in range(4)]
    sample = sorted([max(o - d, 0.0) for o in onsets for d in (1e-9, 1e-6, 0.0)] + [o + 1e-6 for o in onsets])
    s, chi = integrate(rate, [(0.0, T)], 4.5 * T, sample)
    assert np.abs(chi - chi_out(s, T, c['Tc'], c['Tr'], c['ts'])).max() < 1e-8
    first = sample.index(onsets[0] - 1e-6)
    assert chi[first] == 0.0 and chi_out(s[first], T, c['Tc'], c['Tr'], c['ts']) == 0.0


def test_absolute_durations_break_the_closed_form():
    """docs section 5: with absolute durations (T_c, T_r in seconds whatever the period) chi is not a function of s:
    a period step moves psi off P(s)."""
    c = TIMINGS['atrium']

    def absolute(psi, m, Tw):
        Tc, Tr, ts = c['Tc'], c['Tr'], c['ts']
        if m >= ts and m <= ts + EPS_1 and psi <= 0.25:
            return 0.25 / Tc
        if psi > EPS_2 and psi <= 0.25:
            return 0.25 / Tc
        if 0.25 <= psi < 0.5:
            return 0.25 / Tr
        if psi >= 0.5:
            return 0.5 / (Tw - Tc - Tr)
        return 0.0
    sample = np.linspace(0, 5.0, 1001)
    s, chi = integrate(absolute, SCENARIOS['down_a_rest'], 5.0, sample)
    d = (chi - np.floor(chi) - P(s, c['T'], c['Tc'], c['Tr'], c['ts']) + 0.5) % 1.0 - 0.5
    assert np.abs(d).max() > 0.05


# ---- reading a version's clock ------------------------------------------------------------------------------

def _loc(e):
    return etree.QName(e).localname


def clock_files():
    """{path: [clock names]} for every module CellML file with a chamber phase clock."""
    out = {}
    for f in sorted(glob.glob(os.path.join(MODULES, '**', '*_modules.cellml'), recursive=True)):
        text = open(f).read()
        found = [c for c, v in CLOCKS.items() if re.search(r'name="%s"' % v['floor'], text)]
        if found:
            out[os.path.relpath(f, ROOT)] = found
    return out


CLOCK_FILES = clock_files()


def component_of(root, var):
    for comp in root.iter('{*}component'):
        if any(v.get('name') == var for v in comp.iter('{*}variable')):
            return comp
    raise KeyError(var)


def equations(comp):
    """{lhs variable: rhs element} for the algebraic equations and {state: rate element} for the ODEs."""
    alg, ode = {}, {}
    for math_el in comp.iter(M + 'math'):
        for eq in math_el:
            if not isinstance(eq.tag, str) or _loc(eq) != 'apply' or _loc(eq[0]) != 'eq':
                continue
            lhs, rhs = eq[1], eq[2]
            if _loc(lhs) == 'ci':
                alg[lhs.text.strip()] = rhs
            elif _loc(lhs) == 'apply' and _loc(lhs[0]) == 'diff':
                ode[[c for c in lhs if _loc(c) == 'ci'][0].text.strip()] = rhs
    return alg, ode


def form(text, clock):
    """'integrated' (chi is a state) or 'phase' (chi_floor computed from s) for one clock of a file."""
    root = etree.fromstring(text.encode())
    _alg, ode = equations(component_of(root, CLOCKS[clock]['floor']))
    return 'integrated' if clock in ode else 'phase'


def integrates_chi(path):
    text = open(os.path.join(ROOT, path)).read()
    return any(form(text, c) == 'integrated' for c in CLOCK_FILES[path])


OPS = {'plus': '+', 'minus': '-', 'times': '*', 'divide': '/', 'lt': '<', 'leq': '<=', 'gt': '>', 'geq': '>=',
       'eq': '==', 'and': ' and ', 'or': ' or '}


def py(e):
    """A Python expression of a MathML element (scalars), names as they are."""
    t = _loc(e)
    if t == 'ci':
        return e.text.strip()
    if t == 'cn':
        return repr(float(e.text.strip()))
    if t == 'pi':
        return repr(math.pi)
    if t == 'piecewise':
        expr = 'math.nan'
        for c in reversed(list(e)):
            if _loc(c) == 'otherwise':
                expr = py(c[0])
            else:
                expr = f'(({py(c[0])}) if ({py(c[1])}) else ({expr}))'
        return expr
    if t == 'apply':
        op, args = _loc(e[0]), [a for a in e[1:] if _loc(a) not in ('bvar',)]
        if op == 'minus' and len(args) == 1:
            return f'(-({py(args[0])}))'
        if op in OPS:
            return '(' + OPS[op].join(f'({py(a)})' for a in args) + ')'
        if op in ('floor', 'cos', 'sin', 'exp', 'tanh'):
            return f'math.{op}({py(args[0])})'
        if op in ('abs', 'ln'):
            return f'{dict(abs="abs", ln="math.log")[op]}({py(args[0])})'
        if op == 'root' and not any(_loc(a) == 'degree' for a in e[1:]):
            return f'math.sqrt({py(args[0])})'
        if op in ('max', 'min'):
            return f'{op}(' + ', '.join(py(a) for a in args) + ')'
        if op == 'power':
            return f'(({py(args[0])}) ** ({py(args[1])}))'
    raise NotImplementedError(etree.tostring(e)[:200])


def clock_functions(text, clock, values):
    """
    The clock of one component as Python functions: its rate (integrated form; a function of chi, s, T_w) or its
    chi_floor and chi (phase form; functions of s, T_w). ``values``: its inputs (T, T_c, T_r, t_start; T_w is set
    per call as T_wCont or T_eff). Every other variable the clock needs is computed from the component's equations.
    """
    v = CLOCKS[clock]
    root = etree.fromstring(text.encode())
    alg, ode = equations(component_of(root, v['floor']))
    code = {}
    for k, e in alg.items():      # equations the clock does not need may use operators py() does not know
        try:
            code[k] = py(e)
        except NotImplementedError:
            code[k] = 'math.nan'
    comp = component_of(root, v['floor'])
    for var in comp.iter('{*}variable'):       # constants given as initial values (eps_1, eps_2)
        iv = var.get('initial_value')
        if iv and var.get('name') not in code and var.get('name') not in ode and re.match(r'^[-0-9.e+]+$', iv):
            code[var.get('name')] = repr(float(iv))

    def make(target):
        order, seen = [], set()

        def visit(name):
            if name in seen or name not in code:
                return
            seen.add(name)
            for dep in re.findall(r'[A-Za-z_][A-Za-z0-9_]*', code[name]):
                if dep not in ('math', 'and', 'or', 'if', 'else', 'max', 'min', 'float', 'nan'):
                    visit(dep)
            order.append(name)
        for dep in re.findall(r'[A-Za-z_][A-Za-z0-9_]*', target):
            visit(dep)
        body = ''.join(f'    {n} = {code[n]}\n' for n in order if n not in ('mt', 'T_wCont', 'T_eff') and n not in values)
        src = (f'def f({clock}, s, Tw, {", ".join(values)}):\n    mt = s - math.floor(s)\n'
               f'    {v["eff"][2]} = Tw\n    T_eff = Tw\n    T_wCont = Tw\n{body}    return {target}\n')
        env = {'math': math}
        exec(src, env)
        return lambda chi, s, Tw: env['f'](chi, s, Tw, *values.values())
    if clock in ode:
        return {'rate': make(py(ode[clock]))}
    return {'floor': make(v['floor']), 'chi': make(clock)}


def _values(text, clock, c):
    """Inputs for clock_functions: the component's own names for T, T_c, T_r and t_start."""
    v = CLOCKS[clock]
    return {'T': c['T'], v['tc']: c['Tc'], v['tr']: c['Tr'], v['ts_raw']: c['ts'] * c['T']}


@pytest.mark.parametrize('path,clock', [(p, c) for p, cs in CLOCK_FILES.items() for c in cs])
@pytest.mark.parametrize('scenario', ['constant', 'down_window', 'up_a_relaxation'])
def test_version_clock_gives_the_closed_form(path, clock, scenario):
    """Each version's own clock equations, in whichever form it has, give P(s) and the output chi."""
    text = open(os.path.join(ROOT, path)).read()
    if scenario != 'constant' and not re.search(r'name="(T_wCont|T_eff)"', text):
        pytest.skip('a constant-period version: its T_w is T')
    c = TIMINGS['atrium' if clock == 'chi_a' else 'ventricle' if clock == 'chi_v' else 'other']
    fns = clock_functions(text, clock, _values(text, clock, c))
    t_end = 4.0
    sample = np.linspace(0.0, t_end, 801)
    if 'rate' in fns:
        # the rate reads chi only through chi - floor(chi), and s through mt
        s, chi = integrate(fns['rate'], SCENARIOS[scenario], t_end, sample)
        floor_ = chi - np.floor(chi)
    else:
        s, _ = integrate(lambda psi, m, Tw: 0.0, SCENARIOS[scenario], t_end, sample)
        floor_ = np.array([fns['floor'](0.0, x, 1.0) for x in s])
        chi = np.array([fns['chi'](0.0, x, 1.0) for x in s])
    ref = P(s, c['T'], c['Tc'], c['Tr'], c['ts'])
    assert np.abs((floor_ - ref + 0.5) % 1.0 - 0.5).max() < 1e-9
    assert np.abs(chi - chi_out(s, c['T'], c['Tc'], c['Tr'], c['ts'])).max() < 1e-9


# ---- converting between the forms ---------------------------------------------------------------------------

def _cn(x):
    return f'<cn cellml:units="dimensionless">{x}</cn>'


def phase_mathml(clock):
    v = CLOCKS[clock]
    ts, tc, tr, phi, flr = v['ts'], v['tc'], v['tr'], v['phi'], v['floor']
    p = f'<apply><minus/><ci>s</ci><ci>{ts}</ci></apply>'
    phiT = f'<apply><times/><ci>{phi}</ci><ci>T</ci></apply>'
    before = f'<apply><lt/><ci>s</ci><ci>{ts}</ci></apply>'
    tctr = f'<apply><plus/><ci>{tc}</ci><ci>{tr}</ci></apply>'
    e_phi = f'<apply><eq/><ci>{phi}</ci><apply><minus/>{p}<apply><floor/>{p}</apply></apply></apply>'
    e_flr = (f'<apply><eq/><ci>{flr}</ci><piecewise>'
             f'<piece>{_cn("0.0")}{before}</piece>'
             f'<piece><apply><divide/><apply><times/>{_cn("0.25")}{phiT}</apply><ci>{tc}</ci></apply>'
             f'<apply><leq/>{phiT}<ci>{tc}</ci></apply></piece>'
             f'<piece><apply><plus/>{_cn("0.25")}<apply><divide/><apply><times/>{_cn("0.25")}'
             f'<apply><minus/>{phiT}<ci>{tc}</ci></apply></apply><ci>{tr}</ci></apply></apply>'
             f'<apply><lt/>{phiT}{tctr}</apply></piece>'
             f'<otherwise><apply><plus/>{_cn("0.5")}<apply><divide/><apply><times/>{_cn("0.5")}'
             f'<apply><minus/>{phiT}{tctr}</apply></apply><apply><minus/><ci>T</ci>{tctr}</apply></apply></apply></otherwise>'
             f'</piecewise></apply>')
    e_chi = (f'<apply><eq/><ci>{clock}</ci><piecewise><piece>{_cn("0.0")}{before}</piece>'
             f'<otherwise><apply><plus/><apply><floor/>{p}</apply><ci>{flr}</ci></apply></otherwise></piecewise></apply>')
    return [e_phi, e_flr], [e_chi]


def integrated_mathml(clock, eff):
    """The original law: chi_floor = chi - floor(chi), and d chi/dt with the start window and the wait."""
    v = CLOCKS[clock]
    ts, flr = v['ts'], v['floor']
    tcw, trw, tw = eff
    ge = 'gt' if v['strict'] else 'geq'
    rc = f'<apply><divide/>{_cn("0.25")}<ci>{tcw}</ci></apply>'
    le_q = f'<apply><leq/><ci>{flr}</ci>{_cn("0.25")}</apply>'
    tse = f'<apply><plus/><ci>{ts}</ci><ci>eps_1</ci></apply>'
    e_flr = f'<apply><eq/><ci>{flr}</ci><apply><minus/><ci>{clock}</ci><apply><floor/><ci>{clock}</ci></apply></apply></apply>'
    e_ode = (f'<apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>{clock}</ci></apply><piecewise>'
             f'<piece>{rc}<apply><and/><apply><geq/><ci>mt</ci><ci>{ts}</ci></apply><apply><leq/><ci>mt</ci>{tse}</apply>{le_q}</apply></piece>'
             f'<piece>{rc}<apply><and/><apply><gt/>{tse}{_cn("1.0")}</apply>'
             f'<apply><leq/><ci>mt</ci><apply><minus/>{tse}{_cn("1.0")}</apply></apply>{le_q}</apply></piece>'
             f'<piece>{rc}<apply><and/><apply><gt/><ci>{flr}</ci><ci>eps_2</ci></apply>{le_q}</apply></piece>'
             f'<piece><apply><divide/>{_cn("0.25")}<ci>{trw}</ci></apply><apply><and/><apply><{ge}/><ci>{flr}</ci>{_cn("0.25")}</apply>'
             f'<apply><lt/><ci>{flr}</ci>{_cn("0.5")}</apply></apply></piece>'
             f'<piece><apply><divide/>{_cn("0.5")}<apply><minus/><apply><minus/><ci>{tw}</ci><ci>{tcw}</ci></apply><ci>{trw}</ci></apply></apply>'
             f'<apply><{ge}/><ci>{flr}</ci>{_cn("0.5")}</apply></piece>'
             f'<otherwise><cn cellml:units="per_s">0.0</cn></otherwise></piecewise></apply>')
    return [e_flr], [e_ode]


_TAG = re.compile(r'<(/?)([A-Za-z:]+)[^>]*?(/?)>')


def _math_children(text):
    spans = []
    for m in re.finditer(r'<math\b[^>]*>', text):
        pos, depth, start = m.end(), 0, None
        while True:
            t = _TAG.search(text, pos)
            closing, selfclose = t.group(1), t.group(3)
            if closing and depth == 0:
                break
            if not closing and not selfclose:
                if depth == 0:
                    start = t.start()
                depth += 1
            elif closing:
                depth -= 1
                if depth == 0:
                    spans.append((start, t.end()))
            elif depth == 0:
                spans.append((t.start(), t.end()))
            pos = t.end()
    return spans


def _kind(block):
    el = etree.fromstring(block.replace('cellml:units', 'units'))
    kids = list(el)
    if _loc(el) != 'apply' or _loc(kids[0]) != 'eq':
        return None
    if _loc(kids[1]) == 'ci':
        return ('alg', kids[1].text.strip())
    if _loc(kids[1]) == 'apply' and _loc(kids[1][0]) == 'diff':
        return ('ode', [c for c in kids[1] if _loc(c) == 'ci'][0].text.strip())
    return None


def _pretty(blocks, indent, unit):
    ns = f'xmlns="{MATHML}" xmlns:cellml="http://www.cellml.org/cellml/1.1#"'
    out = []
    for b in blocks:
        el = etree.fromstring(f'<math {ns}>{b}</math>')[0]
        etree.indent(el, space=unit)
        s = re.sub(r' xmlns(:cellml)?="[^"]*"', '', etree.tostring(el).decode())
        out.append('\n'.join(indent + line for line in s.split('\n')))
    return '\n'.join(out).lstrip()


def convert(text, to):
    """The file with every clock in the form ``to`` ('phase' or 'integrated'); text edits only."""
    unit = '  ' if '\n  <component' in text or '\n  <units' in text else '    '
    clocks = [c for c, v in CLOCKS.items() if re.search(r'name="%s"' % v['floor'], text)]
    for clock in clocks:
        v = CLOCKS[clock]
        if form(text, clock) == to:
            continue
        spans = _math_children(text)
        kinds = [_kind(text[a:b]) for a, b in spans]
        if to == 'phase':
            first, second = phase_mathml(clock)
            targets = {('alg', v['floor']): first, ('ode', clock): second}
        else:
            eff = next(e for e in [v['eff'], (v['tc'], v['tr'], 'T')] if re.search(r'name="%s"' % e[0], text))
            first, second = integrated_mathml(clock, eff)
            targets = {('alg', v['floor']): first, ('alg', clock): second}
        drop = {('alg', v['phi'])} if to == 'integrated' else set()
        edits = []
        for (a, b), k in zip(spans, kinds):
            ind = text[text.rfind('\n', 0, a) + 1:a]
            if k in targets:
                edits.append((a, b, _pretty(targets[k], ind, unit)))
            elif k in drop:
                edits.append((text.rfind('\n', 0, a), b, ''))
        for a, b, new in sorted(edits, reverse=True):
            text = text[:a] + new + text[b:]
        decl = re.compile(r'(<variable\b[^>]*\bname="%s"[^>]*?)(\s*/>)' % clock)
        if to == 'phase':
            text = re.sub(r'(<variable\b[^>]*?)\s+initial_value="[^"]*"([^>]*\bname="%s"[^>]*/>)' % clock, r'\1\2', text, count=1)
            text = re.sub(r'(<variable\b[^>]*\bname="%s"[^>]*?)\s+initial_value="[^"]*"' % clock, r'\1', text, count=1)
            m = re.search(r'\n([ \t]*)<variable\b[^>]*\bname="%s"[^>]*/>' % v['final'], text)
            text = text[:m.end()] + '\n' + m.group(1) + f'<variable name="{v["phi"]}" units="dimensionless"/>' + text[m.end():]
        else:
            text = decl.sub(lambda m: m.group(1) + ' initial_value="0.0"' + m.group(2), text, count=1)
            text = re.sub(r'\n[ \t]*<variable\b[^>]*\bname="%s"[^>]*/>' % v['phi'], '', text)
    if to == 'integrated' and not re.search(r'name="eps_1"', text):
        m = re.search(r'\n([ \t]*)<variable\b[^>]*\bname="%s"[^>]*/>' % CLOCKS[clocks[0]]['final'], text)
        add = ''.join(f'\n{m.group(1)}<variable name="{n}" units="dimensionless" initial_value="{x}"/>' for n, x in (('eps_1', EPS_1), ('eps_2', EPS_2)))
        text = text[:m.end()] + add + text[m.end():]
    if to == 'phase':
        for eps in ('eps_1', 'eps_2'):
            if not re.search(r'<ci>\s*%s\s*</ci>' % eps, text):
                text = re.sub(r'\n[ \t]*<variable\b[^>]*\bname="%s"[^>]*/>' % eps, '', text)
    return text


def math_signature(text):
    """The file's equations and variable declarations, order-free: for comparing two forms of one file."""
    root = etree.fromstring(text.encode())
    eqs, decls = [], []
    for comp in root.iter('{*}component'):
        for var in comp.iter('{*}variable'):
            a = dict(var.attrib)
            if 'initial_value' in a:
                a['initial_value'] = repr(float(a['initial_value'])) if re.match(r'^[-0-9.e+]+$', a['initial_value']) else a['initial_value']
            decls.append((comp.get('name'), tuple(sorted(a.items()))))
        for math_el in comp.iter(M + 'math'):
            for eq in math_el:
                if not isinstance(eq.tag, str):
                    continue
                s = re.sub(r'\s+', '', etree.tostring(eq).decode())
                s = re.sub(r'<cn([^>]*)>([^<]*)</cn>', lambda m: f'<cn{m.group(1)}>{float(m.group(2))!r}</cn>', s)
                eqs.append((comp.get('name'), re.sub(r' xmlns(:\w+)?="[^"]*"', '', s)))
    return sorted(eqs), sorted(decls)


@pytest.mark.parametrize('path', sorted(CLOCK_FILES))
def test_converter_round_trip(path):
    """Converting a version's file to the other form and back gives its equations again, so the simulations below
    compare the version with exactly the original law (or exactly the closed form)."""
    text = open(os.path.join(ROOT, path)).read()
    here = 'integrated' if integrates_chi(path) else 'phase'
    other = 'phase' if here == 'integrated' else 'integrated'
    there = convert(text, other)
    assert all(form(there, c) == other for c in CLOCK_FILES[path])
    assert math_signature(convert(there, here)) == math_signature(text)


# ---- simulations of the version families --------------------------------------------------------------------

def _version(version_id):
    from cam_testing import library
    for v in library.all_versions():
        if v.id == version_id:
            return v
    pytest.skip(f'{version_id} is not in the library')


@contextlib.contextmanager
def _quiet():
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def generate_forms(version_id, work_dir):
    """The version's test network generated as it is, and a copy with every clock in the other form."""
    from cam_testing import harness
    pytest.importorskip('libcuflynx')
    v = _version(version_id)
    with _quiet():
        path = harness.generate(v, os.path.join(work_dir, 'as_is'))
    other_dir = os.path.join(work_dir, 'other')
    shutil.copytree(os.path.dirname(path), other_dir)
    forms_ = set()
    for f in glob.glob(os.path.join(other_dir, '*_modules.cellml')):
        text = open(f).read()
        clocks = [c for c, vv in CLOCKS.items() if re.search(r'name="%s"' % vv['floor'], text)]
        if not clocks:
            continue
        here = form(text, clocks[0])
        forms_.add(here)
        own = [v.get('name') for v in component_of(etree.fromstring(text.encode()), CLOCKS[clocks[0]]['floor']).iter('{*}variable')]
        if here == 'integrated' and 's' not in own:
            pytest.skip('the chamber has no beat phase s input (its cardiac_clock_port carries mt and T_eff only): '
                        'the phase form needs the port change of the module PR')
        open(f, 'w').write(convert(text, 'phase' if here == 'integrated' else 'integrated'))
    assert len(forms_) == 1, forms_
    as_is = forms_.pop()
    other = os.path.join(other_dir, os.path.basename(path))
    return {as_is: path, ('phase' if as_is == 'integrated' else 'integrated'): other}


def myokit_model(path):
    from libcuflynx.solver_wrappers import get_simulation_helper
    with _quiet():
        h = get_simulation_helper(model_path=path, solver='CVODE_myokit', model_type='cellml', dt=0.01, sim_time=1.0,
                                  solver_info={}, pre_time=0.0)
    return h.model


def simulate(model, t_end, tol, max_step=1e-3, steps=(), state_factor=None, log_interval=1e-3, breaks=()):
    """Myokit CVODE; ``steps``: [(t, {constant qname: value})] applied during the run (a period change); ``breaks``:
    times at which the run stops and CVODE restarts, so no step spans them. Logged on one grid of ``log_interval``."""
    import myokit
    sim = myokit.Simulation(model)
    sim.set_tolerance(abs_tol=tol, rel_tol=tol)
    sim.set_max_step_size(max_step)
    if state_factor is not None:
        sim.set_state([x * state_factor for x in model.initial_values(as_floats=True)])
    names = [s.qname() for s in model.states()] + [model.time().qname()]
    clocks = [x.qname() for x in model.variables(deep=True) if x.name() in ('e_a', 'e_v', 'e', 'chi_a', 'chi_v', 'chi')]
    grid = np.round(np.arange(int(round(t_end / log_interval)) + 1) * log_interval, 12)
    snap = lambda x: float(grid[np.argmin(np.abs(grid - x))]) if np.abs(grid - x).min() < 1e-9 else float(x)
    change_at = {snap(t_): c for t_, c in steps}
    ends = sorted(set(list(change_at) + [snap(b) for b in breaks if 0 < b < t_end] + [float(t_end)]))
    out, t = {k: [] for k in names + clocks}, 0.0
    for t_next in ends:
        last = t_next == t_end
        times = grid[(grid >= t) & ((grid <= t_next) if last else (grid < t_next))]
        # Myokit logs times in [start, end): the last run goes half an interval past t_end to log t_end
        d = sim.run(t_next - t + (log_interval / 2 if last else 0.0), log=names + clocks, log_times=times)
        for k in out:
            out[k].append(np.asarray(d[k]))
        for q, val in change_at.get(t_next, {}).items():
            sim.set_constant(q, val)
        t = t_next
    return {k: np.concatenate(v) for k, v in out.items()}


def clock_switch_times(model, t, s):
    """The times at which an integrated chi law switches piece: each onset k + t_s', and t_s' + T_c/T and
    t_s' + (T_c + T_r)/T after it. The constants are read from ``model``'s closed form (chi_floor =
    piecewise(s < t_s', 0, phi T <= T_c, ..., phi T < T_c + T_r, ...)); s(t) is the beat phase logged on the times
    t (piecewise linear in t, with its kinks at the period steps, which are on the grid)."""
    out = []
    for x in model.variables(deep=True):
        if x.name() not in [v['floor'] for v in CLOCKS.values()]:
            continue
        c = [list(e) for e in x.rhs().conditions()]
        ts, T, tc, tctr = c[0][1].eval(), list(c[1][0])[1].eval(), c[1][1].eval(), c[2][1].eval()
        for k in range(int(np.floor(s[-1] - ts)) + 2):
            for frac in (0.0, tc / T, tctr / T):
                if s[0] <= k + ts + frac <= s[-1]:
                    out.append(float(np.interp(k + ts + frac, s, t)))
    return sorted(set(out))


def _period_input(model):
    """The qname of the test network's Delta_T parameter (the controlled hearts' period change), or None."""
    c = [v.qname() for v in model.variables(const=True, deep=True) if re.search(r'(^|\.)Delta_T(_|$)', v.name())]
    return c[0] if c else None


# version families: the flat constant-period heart, the flat controlled heart (with period steps), a chamber in its
# cardiac clock harness, and the controlled heart supermodule (chambers reading the clock's s, with period steps)
FAMILIES = {'flat_constant': 'heart__vp', 'flat_controlled': 'heart__vp_wCont', 'chamber': 'chamber__Liang2009NoVisco_v02',
            'supermodule_controlled': 'heart__controlled_Argus2026_v01'}
# period steps (t, Delta_T) from T = 1 s: 1.0 -> 0.65 s in atrial contraction (s = 1.9), 0.65 -> 1.2 s in atrial
# relaxation (s = 3.05), 1.2 -> 1.0 s in ventricular systole (s = 4.2)
PERIOD_STEPS = [(1.9, -0.35), (2.65, 0.2), (4.03, 0.0)]
FORMS_TOL = 1e-6      # max |difference| / (max |value| + 1e-12) over each state and the activations


@pytest.mark.parametrize('family', sorted(FAMILIES))
def test_forms_agree(family, tmp_path):
    """Each version family simulated in both forms gives the same trajectories. CVODE at rtol = atol = 1e-12, and the
    integrated form restarted at each of its switches (known from the closed form): stepping across a switch, its own
    error depends on the SUNDIALS build (1.5e-5 in zeta_miv_pre on CI's, 1e-7 on a native 7.7.0) and would dominate the
    comparison; the closed form has no switch in a state and moves by 3e-12 between tolerances 1e-10 and 1e-12."""
    paths = generate_forms(FAMILIES[family], str(tmp_path))
    models = {k: myokit_model(p) for k, p in paths.items()}
    period = _period_input(models['integrated'])
    steps = [(t, {period: d}) for t, d in PERIOD_STEPS] if period and 'controlled' in family else []
    runs = {'phase': simulate(models['phase'], 5.5, 1e-12, steps=steps)}
    m = models['phase']
    s_q = [x.qname() for x in m.states() if x.name() == 's'][0]
    breaks = clock_switch_times(m, runs['phase'][m.time().qname()], runs['phase'][s_q])
    runs['integrated'] = simulate(models['integrated'], 5.5, 1e-12, steps=steps, breaks=breaks)
    a, b = runs['integrated'], runs['phase']
    common = sorted(k for k in set(a) & set(b) if not k.endswith('time') and a[k].shape == b[k].shape)
    assert any(k.endswith(('.e_a', '.e')) for k in common)
    worst = {}
    for k in common:
        x, y = a[k], b[k]
        if re.search(r'\.chi(_[av])?$', k):
            x, y = x - np.floor(x + 1e-9), y - np.floor(y + 1e-9)   # the phase, modulo the integer it just crossed
            x, y = np.minimum(x, 1 - x), np.minimum(y, 1 - y)
        worst[k] = float(np.abs(x - y).max() / (np.abs(x).max() + 1e-12))
    k = max(worst, key=worst.get)
    print(f'{family}: worst {k} {worst[k]:.2e}; ' + ', '.join(f'{n} {worst[n]:.1e}' for n in worst if re.search(r'\.(e_a|e_v|e|q_la|q_lv|q)$', n)))
    assert worst[k] < FORMS_TOL, (k, worst[k])


def _onset_state(model, s_q, at):
    """The model's state just before the beat phase s reaches ``at`` (the first atrial onset), by CVODE at 1e-10."""
    import myokit
    sim = myokit.Simulation(model)
    sim.set_tolerance(1e-12, 1e-10)
    sim.set_max_step_size(1e-3)
    d = sim.run(0.5, log=[s_q, model.time().qname()])
    T_w = 0.5 / d[s_q][-1]
    sim.reset()
    sim.run(at * T_w - 1e-7, log=myokit.LOG_NONE)
    return sim.state()


ONSET_VERSIONS = {'flat_controlled': 'heart__vp_wCont', 'supermodule_controlled': 'heart__controlled_Argus2026_v01'}


@pytest.mark.parametrize('family', sorted(ONSET_VERSIONS))
def test_clock_rates_have_no_jump_at_the_onset(family, tmp_path, request):
    """
    The mechanism of the onset shift (docs section 6): just before the first atrial onset the integrated chi_a waits at
    exactly 0, and a change of it by round-off (-1e-12) moves its rate from the wait (0) to the rest piece, so a
    finite-difference Jacobian taken there is ~1e12 and Newton cannot correct chi_a. Here: no state's rate may move by
    more than 1e3 x the change when one clock state (chi, chi_a, chi_v or s) changes by 1e-12.
    """
    v = _version(ONSET_VERSIONS[family])
    from cam_testing import harness
    pytest.importorskip('libcuflynx')
    with _quiet():
        path = harness.generate(v, str(tmp_path))
    text = ''.join(open(f).read() for f in glob.glob(os.path.join(os.path.dirname(path), '*_modules.cellml')))
    if re.search(r'<diff/>\s*<bvar>\s*<ci>\s*t\s*</ci>\s*</bvar>\s*<ci>\s*chi(_a|_v)?\s*</ci>', text):
        request.applymarker(pytest.mark.xfail(reason='the version still integrates chi: its rate jumps at the wait '
                                                     '(docs/heart/cardiac_clock_phase.md, section 6)', strict=True))
    model = myokit_model(path)
    states = [x.qname() for x in model.states()]
    s_q = [q for q in states if q.endswith('.s')][0]
    x0 = _onset_state(model, s_q, 0.8)
    r0 = np.array(model.evaluate_derivatives(x0))
    worst = (0.0, '')
    for j, q in enumerate(states):
        if not re.search(r'\.(s|chi|chi_a|chi_v)$', q):
            continue
        for sign in (-1.0, 1.0):
            x = list(x0)
            x[j] += sign * 1e-12
            slope = np.abs(np.array(model.evaluate_derivatives(x)) - r0).max() / 1e-12
            worst = max(worst, (slope, f'{q} {sign:+.0f}e-12'))
    assert worst[0] < 1e3, worst


# the sensitivity sweep: tolerances x {as is, initial states x (1 + 1e-12)}; each atrial onset must stay within 1 ms of
# s = k + t_s (the integrated form's bad runs are 35-76 ms off; docs section 6)
SWEEP_TOLS = np.logspace(-10, -7, 13)


def atrial_onsets(t, s, e, ts=0.8, Tc_frac=0.17, level=1e-4):
    """(found - expected) in ms for each atrial onset: e crosses ``level`` after s passes k + ts; expected from s and
    the exact rise e = (1 - cos(4 pi psi))/2, psi = (s - k - ts) / (4 Tc/T)."""
    out = []
    for k in range(int(np.floor(s[-1] - ts)) + 1):
        if np.interp(k + ts + 0.05, s, t) > t[-1]:
            break
        ds = np.arccos(1 - 2 * level) / (4 * np.pi) * 4 * Tc_frac          # s from the onset to e = level
        expected = np.interp(k + ts + ds, s, t)
        i = np.searchsorted(t, np.interp(k + ts - 0.05, s, t))
        j = i + int(np.argmax(e[i:] > level))
        found = t[j - 1] + (level - e[j - 1]) / (e[j] - e[j - 1]) * (t[j] - t[j - 1])
        out.append((found - expected) * 1e3)
    return out


@pytest.mark.parametrize('family', sorted(ONSET_VERSIONS))
def test_atrial_onset_ignores_round_off(family, tmp_path, request):
    """
    The atrial onset does not move with the solver tolerance (1e-10 to 1e-7) or a 1e-12 change of the initial state.
    With the integrated chi, round-off moves it 35-76 ms in about 1 run in 10 (docs section 6): how often depends on
    the platform's SUNDIALS build and the solver (it showed in libOpenCOR, Myokit on a generic SUNDIALS 7.7.0 build
    and a C replica, not in Myokit on a -march=native build), so the integrated form is an expected failure here that
    may pass.
    """
    v = _version(ONSET_VERSIONS[family])
    from cam_testing import harness
    pytest.importorskip('libcuflynx')
    with _quiet():
        path = harness.generate(v, str(tmp_path))
    text = ''.join(open(f).read() for f in glob.glob(os.path.join(os.path.dirname(path), '*_modules.cellml')))
    if re.search(r'<diff/>\s*<bvar>\s*<ci>\s*t\s*</ci>\s*</bvar>\s*<ci>\s*chi(_a|_v)?\s*</ci>', text):
        request.applymarker(pytest.mark.xfail(reason='the version still integrates chi: its atrial onset moves by '
                                                     'round-off in some runs (docs/heart/cardiac_clock_phase.md, '
                                                     'section 6)', strict=False))
    model = myokit_model(path)
    s_q = [x.qname() for x in model.states() if x.name() == 's'][0]
    e_q = [x.qname() for x in model.variables(deep=True) if x.name() == 'e_a' or x.qname().endswith('ra_module.e')][0]
    bad = []
    for tol in SWEEP_TOLS:
        for factor in (None, 1 + 1e-12):
            r = simulate(model, 2.0, float(tol), state_factor=factor, log_interval=1e-4)
            errs = atrial_onsets(r[model.time().qname()], r[s_q], r[e_q])
            assert len(errs) == 2
            if max(abs(x) for x in errs) > 1.0:
                bad.append((f'{tol:.3g}', factor, [round(x, 2) for x in errs]))
    assert not bad, bad
