"""
Failure-risk analysis over the joint parameter space (on request, not part of the test run).

The one-at-a-time verified ranges only describe a cross through the nominal point. Here each
component's parameters are sampled jointly with a scrambled Sobol sequence over the same box
the sweep uses (bc_sweep.ranges, or bc_sweep.factors x nominal; log scale where the spec says
so). Every sample is run with the component's solver settings and counts as a failure when the
run errors, an output is non-finite, or an invariant is violated.

Outputs, per component (modules/<name>/results/<component>/risk_analysis.json + samples .npz):
  - overall failure probability with a 95% interval
  - per-parameter failure-risk curves P(fail | parameter) with 95% intervals
  - a first-order importance of each parameter for failure: the correlation ratio eta^2, the
    fraction of the failure indicator's variance explained by that parameter alone
  - a corner plot of pairwise failure rates (the most important parameters)

    python -m cam_testing.risk --module Lotka_Volterra [--samples 1024] [--component ID]
"""
import argparse
import json
import math
import os

import numpy as np

from cam_testing import checks, harness, plots
from cam_testing.library import load_module, module_names

DEFAULT_SAMPLES = 512
N_BINS = 8


def wilson(k, n, z=1.96):
    '''95% Wilson score interval for k failures out of n.'''
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def parameter_box(cm):
    '''[(variable_name, parameter_name, lo, hi, log)] for every swept parameter.'''
    sweep_spec = cm.spec.get('bc_sweep') or {}
    ranges = sweep_spec.get('ranges') or {}
    factors = [float(f) for f in sweep_spec.get('factors', [0.5, 2.0])]
    exclude = set(sweep_spec.get('exclude') or [])
    box = []
    for p in cm.component.parameters():
        if p.is_todo or p.variable_name in exclude:
            continue
        pname = harness.parameter_name(p)
        r = ranges.get(p.variable_name)
        if isinstance(r, dict):
            box.append((p.variable_name, pname, float(r['min']), float(r['max']), r.get('scale') == 'log'))
        elif r is not None:
            box.append((p.variable_name, pname, float(r[0]), float(r[1]), False))
        else:
            nominal = cm.nominal[pname]
            if nominal == 0:
                continue   # no scale to vary it by; give it a range in the spec
            lo, hi = sorted((nominal * min(factors), nominal * max(factors)))
            box.append((p.variable_name, pname, lo, hi, lo > 0))
    return box


def sample_box(box, n, seed=0):
    from scipy.stats import qmc
    m = int(math.ceil(math.log2(max(n, 2))))
    u = qmc.Sobol(d=len(box), scramble=True, seed=seed).random_base2(m)[:n]
    x = np.empty_like(u)
    for j, (_, _, lo, hi, log) in enumerate(box):
        x[:, j] = np.exp(np.log(lo) + u[:, j] * (np.log(hi) - np.log(lo))) if log else lo + u[:, j] * (hi - lo)
    return x


def _bin_edges(lo, hi, log, n_bins=N_BINS):
    return np.geomspace(lo, hi, n_bins + 1) if log else np.linspace(lo, hi, n_bins + 1)


def marginal_risk(values, failed, lo, hi, log):
    edges = _bin_edges(lo, hi, log)
    idx = np.clip(np.searchsorted(edges, values, side='right') - 1, 0, len(edges) - 2)
    bins = []
    for b in range(len(edges) - 1):
        m = idx == b
        k, n = int(failed[m].sum()), int(m.sum())
        bins.append({'lo': float(edges[b]), 'hi': float(edges[b + 1]), 'n': n, 'failures': k,
                     'risk': k / n if n else float('nan'), 'ci': wilson(k, n)})
    return bins


def correlation_ratio(bins, overall):
    '''eta^2 = between-bin variance of the failure indicator / its total variance.'''
    n_tot = sum(b['n'] for b in bins)
    var = overall * (1 - overall)
    if var == 0 or n_tot == 0:
        return 0.0
    between = sum(b['n'] * (b['risk'] - overall) ** 2 for b in bins if b['n']) / n_tot
    return float(between / var)


def _mode(spec, message):
    '''The invariant's description (or expression) that a failure message refers to.'''
    for inv in spec.get('invariants') or []:
        expr = inv['expr'] if isinstance(inv, dict) else inv
        if message.startswith(expr):
            return (inv.get('description') if isinstance(inv, dict) else None) or expr
    return message


def analyse(cm, n_samples=DEFAULT_SAMPLES, seed=0, corner_max=6):
    component = cm.component
    box = parameter_box(cm)
    if not box:
        return None
    x = sample_box(box, n_samples, seed)
    failed = np.zeros(len(x), dtype=bool)
    reasons, modes = [], []
    for i, row in enumerate(x):
        params = {pname: float(v) for (_, pname, *_), v in zip(box, row)}
        try:
            t, outputs = cm.run(params)
            bad = checks._non_finite(outputs)
            with np.errstate(all='ignore'):
                inv = checks._check_invariants(cm.spec, t, outputs, cm.param_env(params))
            if bad or inv:
                failed[i] = True
                reasons.append((i, ('non-finite ' + ', '.join(bad)) if bad else inv[0]))
                modes.append('non-finite output' if bad else _mode(cm.spec, inv[0]))
        except Exception as e:
            failed[i] = True
            reasons.append((i, f'{type(e).__name__}: {str(e)[:120]}'))
            modes.append('simulation error')

    k, n = int(failed.sum()), len(failed)
    overall = k / n
    per_param = {}
    for j, (var, pname, lo, hi, log) in enumerate(box):
        bins = marginal_risk(x[:, j], failed, lo, hi, log)
        per_param[var] = {'range': [lo, hi], 'log': log, 'bins': bins,
                          'importance': correlation_ratio(bins, overall)}
    ranked = sorted(per_param, key=lambda v: -per_param[v]['importance'])

    res_dir = os.path.join(component.module.results_dir, component.id)
    os.makedirs(res_dir, exist_ok=True)
    np.savez_compressed(os.path.join(res_dir, 'risk_samples.npz'), x=x, failed=failed,
                        names=np.array([b[0] for b in box]))

    figs = []
    if k:
        corner_vars = ranked[:corner_max]
        figs.append(plots.plot_risk_corner(checks.plot_path(component, 'risk_corner'), x, failed,
                                           [b[0] for b in box], {v: per_param[v] for v in corner_vars},
                                           corner_vars, f'{component.label}: failure risk'))
    figs.append(plots.plot_risk_marginals(checks.plot_path(component, 'risk_marginals'), per_param, ranked,
                                          overall, f'{component.label}: failure risk by parameter'))
    result = {
        'n_samples': n, 'n_failures': k, 'failure_probability': overall, 'ci': wilson(k, n),
        'upper_bound_if_none_failed': (3.0 / n) if k == 0 else None,
        'seed': seed, 'box': [{'variable': b[0], 'min': b[2], 'max': b[3], 'log': b[4]} for b in box],
        'failure_modes': {m: modes.count(m) for m in sorted(set(modes), key=lambda m: -modes.count(m))},
        'ranking': [{'variable': v, 'importance': per_param[v]['importance']} for v in ranked],
        'per_parameter': per_param,
        'example_failures': [{'sample': dict(zip([b[0] for b in box], map(float, x[i]))), 'reason': r}
                             for i, r in reasons[:10]],
        'plots': [os.path.relpath(f, component.module.dir) for f in figs],
    }
    with open(os.path.join(res_dir, 'risk_analysis.json'), 'w') as f:
        json.dump(result, f, indent=1, default=checks._json_default)
    return result


def load(component):
    path = os.path.join(component.module.results_dir, component.id, 'risk_analysis.json')
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def local_risk(component, values, k=25):
    '''
    Failure fraction among the k stored samples nearest to ``values`` ({variable: value}),
    distance measured in the sampling box's unit coordinates (log where sampled in log).
    '''
    res_dir = os.path.join(component.module.results_dir, component.id)
    path = os.path.join(res_dir, 'risk_samples.npz')
    risk = load(component)
    if risk is None or not os.path.isfile(path):
        return None
    data = np.load(path)
    names = list(data['names'])
    box = {b['variable']: b for b in risk['box']}

    def unit(col, v, name):
        b = box[name]
        if b['log']:
            return (np.log(v) - np.log(b['min'])) / (np.log(b['max']) - np.log(b['min']))
        return (v - b['min']) / (b['max'] - b['min'])
    use = [j for j, nme in enumerate(names) if nme in values]
    if not use:
        return None
    X = np.column_stack([unit(j, data['x'][:, j], names[j]) for j in use])
    q = np.array([unit(j, values[names[j]], names[j]) for j in use])
    nearest = np.argsort(np.sum((X - q) ** 2, axis=1))[:k]
    kf = int(data['failed'][nearest].sum())
    outside = bool(np.any((q < 0) | (q > 1)))
    return {'risk': kf / len(nearest), 'ci': wilson(kf, len(nearest)), 'k': len(nearest), 'outside_box': outside}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--module', action='append', default=[])
    parser.add_argument('--component', action='append', default=[])
    parser.add_argument('--samples', type=int, default=DEFAULT_SAMPLES)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args(argv)
    for name in args.module or module_names():
        module = load_module(name)
        for component in module.components():
            if args.component and component.id not in args.component:
                continue
            if component.spec.get('skip'):
                continue
            cm = checks.ComponentModel(component)
            try:
                r = analyse(cm, args.samples, args.seed)
            except harness.MissingParameters as e:
                print(f'{name}/{component.id}: skipped ({e})')
                continue
            if r is None:
                print(f'{name}/{component.id}: no parameters to sample')
                continue
            lo, hi = r['ci']
            top = ', '.join(f"{x['variable']} ({x['importance']:.2f})" for x in r['ranking'][:3])
            print(f"{name}/{component.id}: P(fail) = {r['failure_probability']:.3f} [{lo:.3f}, {hi:.3f}] "
                  f"over {r['n_samples']} samples; most influential: {top}")


if __name__ == '__main__':
    main()
