"""
Failure-risk analysis over the joint parameter space (on request, not part of the test run).

The one-at-a-time verified ranges only describe a cross through the nominal point. Here each
component's parameters are sampled jointly with a scrambled Sobol sequence over the same box
the sweep uses (bc_sweep.ranges, or bc_sweep.factors x nominal; log scale where the spec says
so). Every sample is run with the component's solver settings and counts as a failure when the
run errors, an output is non-finite, or an invariant is violated.

Outputs, per version, committed in versions/<version>/risk/ so they are available on any clone:
  <module_type>_<version>_risk.json          summaries + a fitted failure classifier (coefficients only)
  <module_type>_<version>_risk_samples.npz   the samples and whether each failed
  - overall failure probability with a 95% interval
  - per-parameter failure-risk curves P(fail | parameter) with 95% intervals
  - a first-order importance of each parameter for failure: the correlation ratio eta^2, the
    fraction of the failure indicator's variance explained by that parameter alone
  - a corner plot of pairwise failure rates (the most important parameters)
  - a logistic-regression classifier on quadratic features of the unit-box coordinates, giving a
    smooth P(fail | parameters) anywhere in the box; its 5-fold cross-validated AUC is stored

    python -m cam_testing.risk --module Lotka_Volterra [--samples 1024] [--component Lotka_Volterra/nn]
"""
import argparse
import json
import math
import os

import numpy as np

from cam_testing import checks, harness, plots
from cam_testing.library import all_versions

DEFAULT_SAMPLES = 512
DISCRETE = {}   # id(ComponentModel) -> {variable: allowed values} for discrete sweep ranges
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
    discrete = DISCRETE.setdefault(id(cm), {})
    for p in cm.component.parameters():
        if p.is_todo or p.variable_name in exclude:
            continue
        pname = harness.parameter_name(p)
        r = ranges.get(p.variable_name)
        if isinstance(r, dict) and 'values' in r:
            vals = sorted(float(v) for v in r['values'])
            box.append((p.variable_name, pname, vals[0], vals[-1], False))
            discrete[p.variable_name] = vals
        elif isinstance(r, dict):
            box.append((p.variable_name, pname, float(r['min']), float(r['max']), r.get('scale') == 'log'))
        elif r is not None:
            box.append((p.variable_name, pname, float(r[0]), float(r[1]), False))
        else:
            nominal = cm.nominal[pname]
            if nominal == 0:
                continue   # no scale to vary it by; give it a range in the spec
            lo, hi = sorted((nominal * min(factors), nominal * max(factors)))
            box.append((p.variable_name, pname, lo, hi, lo > 0))
    # clipped to each parameter's valid range, as the sweep is (bc_sweep.bounds, gate initial values)
    bounds = checks.parameter_bounds(cm.component)
    clipped = []
    for var, pname, lo, hi, log in box:
        if var in bounds:
            lo, hi = max(lo, bounds[var][0]), min(hi, bounds[var][1])
            if lo >= hi:
                continue
            log = log and lo > 0
        clipped.append((var, pname, lo, hi, log))
    return clipped


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
    for j, b in enumerate(box):                      # discrete parameters take their allowed values
        vals = DISCRETE.get(id(cm), {}).get(b[0])
        if vals:
            x[:, j] = np.asarray(vals)[np.clip(np.floor((x[:, j] - b[2]) / max(b[3] - b[2], 1e-300) * len(vals)),
                                                0, len(vals) - 1).astype(int)]
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

    os.makedirs(risk_dir(component), exist_ok=True)
    np.savez_compressed(samples_path(component), x=x, failed=failed, names=np.array([b[0] for b in box]))

    result = {
        'n_samples': n, 'n_failures': k, 'failure_probability': overall, 'ci': wilson(k, n),
        'upper_bound_if_none_failed': (3.0 / n) if k == 0 else None,
        'seed': seed, 'box': [{'variable': b[0], 'min': b[2], 'max': b[3], 'log': b[4]} for b in box],
        'failure_modes': {m: modes.count(m) for m in sorted(set(modes), key=lambda m: -modes.count(m))},
        'ranking': [{'variable': v, 'importance': per_param[v]['importance']} for v in ranked],
        'per_parameter': per_param,
        'example_failures': [{'sample': dict(zip([b[0] for b in box], map(float, x[i]))), 'reason': r}
                             for i, r in reasons[:10]],
        'classifier': fit_classifier(unit_coords(x, box), failed),
    }
    result['plots'] = make_plots(component, result, x, failed)
    with open(result_path(component), 'w') as f:
        json.dump(result, f, indent=1, default=checks._json_default)
    return result


def risk_dir(component):
    return component.risk_dir


def result_path(component):
    return os.path.join(risk_dir(component), f'{component.stem}_risk.json')


def samples_path(component):
    return os.path.join(risk_dir(component), f'{component.stem}_risk_samples.npz')


def unit_coords(x, box):
    """Samples in the box's unit coordinates (log scale where sampled in log)."""
    u = np.empty_like(np.asarray(x, dtype=float))
    for j, b in enumerate(box):
        lo, hi, log = (b[2], b[3], b[4]) if isinstance(b, tuple) else (b['min'], b['max'], b['log'])
        col = np.asarray(x, dtype=float)[:, j]
        u[:, j] = (np.log(col) - np.log(lo)) / (np.log(hi) - np.log(lo)) if log else (col - lo) / (hi - lo)
    return u


def _quadratic_powers(d):
    powers = [[0] * d]
    for i in range(d):
        p = [0] * d; p[i] = 1; powers.append(p)
    for i in range(d):
        for j in range(i, d):
            p = [0] * d; p[i] += 1; p[j] += 1; powers.append(p)
    return powers


def _features(u, powers):
    u = np.atleast_2d(u)
    return np.column_stack([np.prod(u ** np.array(p), axis=1) for p in powers])


def fit_classifier(u, failed):
    """
    Logistic regression on quadratic features of the unit coordinates. Stored as plain
    coefficients so it can be evaluated anywhere without the fitting library.
    """
    failed = np.asarray(failed, dtype=int)
    d = u.shape[1]
    if failed.min() == failed.max():
        return {'kind': 'constant', 'p_fail': float(failed.mean()),
                'note': 'every sample ' + ('failed' if failed[0] else 'passed')}
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    powers = _quadratic_powers(d)
    X = _features(u, powers)[:, 1:]          # intercept fitted separately
    model = LogisticRegression(C=10.0, max_iter=5000)
    n_splits = int(min(5, failed.sum(), (1 - failed).sum()))
    auc = None
    if n_splits >= 2:
        auc = float(np.mean(cross_val_score(model, X, failed, scoring='roc_auc',
                                            cv=StratifiedKFold(n_splits, shuffle=True, random_state=0))))
    model.fit(X, failed)
    return {'kind': 'logistic_quadratic', 'powers': powers[1:], 'coef': model.coef_[0].tolist(),
            'intercept': float(model.intercept_[0]), 'cv_auc': auc}


def predict(classifier, u):
    """P(fail) at unit coordinates u (n x d) from a stored classifier."""
    u = np.atleast_2d(u)
    if classifier['kind'] == 'constant':
        return np.full(len(u), classifier['p_fail'])
    z = classifier['intercept'] + _features(u, classifier['powers']) @ np.asarray(classifier['coef'])
    return 1.0 / (1.0 + np.exp(-z))


def make_plots(component, result, x=None, failed=None):
    """(Re)draws the risk plots into plots/ from the stored results and samples."""
    if x is None:
        if not os.path.isfile(samples_path(component)):
            return []
        data = np.load(samples_path(component))
        x, failed = data['x'], data['failed'].astype(bool)
    names = [b['variable'] for b in result['box']]
    per_param = result['per_parameter']
    ranked = [r['variable'] for r in result['ranking']]
    figs = []
    if result['n_failures']:
        corner = ranked[:6]
        figs.append(plots.plot_risk_corner(checks.plot_path(component, 'risk_corner'), x, failed, names,
                                           {v: per_param[v] for v in corner}, corner,
                                           f'{component.label}: failure risk'))
    figs.append(plots.plot_risk_marginals(checks.plot_path(component, 'risk_marginals'), per_param, ranked,
                                          result['failure_probability'],
                                          f'{component.label}: failure risk by parameter'))
    return [os.path.relpath(f, component.dir) for f in figs]


def load(component):
    if not os.path.isfile(result_path(component)):
        return None
    with open(result_path(component)) as f:
        return json.load(f)


def local_risk(component, values, k=25):
    """
    P(fail) at a parameter set ({variable: value}, missing ones at nominal): from the stored
    classifier, with the failure fraction of the k nearest stored samples alongside.
    """
    risk = load(component)
    if risk is None:
        return None
    box = risk['box']
    nominal = {p.variable_name: p.float_value for p in component.parameters() if not p.is_todo}
    q = np.array([[values.get(b['variable'], nominal.get(b['variable'], np.nan)) for b in box]], dtype=float)
    if np.isnan(q).any():
        return None
    uq = unit_coords(q, box)[0]
    out = {'outside_box': bool(np.any((uq < 0) | (uq > 1)))}
    if risk.get('classifier'):
        out['model'] = float(predict(risk['classifier'], uq)[0])
        out['cv_auc'] = risk['classifier'].get('cv_auc')
    if os.path.isfile(samples_path(component)):
        data = np.load(samples_path(component))
        U = unit_coords(data['x'], box)
        nearest = np.argsort(np.sum((U - uq) ** 2, axis=1))[:k]
        kf = int(data['failed'][nearest].sum())
        out.update({'risk': kf / len(nearest), 'ci': wilson(kf, len(nearest)), 'k': len(nearest)})
    else:
        out.update({'risk': out.get('model', float('nan')), 'ci': (float('nan'), float('nan')), 'k': 0})
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--module', action='append', default=[], help='module_type (with those nested in it) or path under modules/ (repeatable)')
    parser.add_argument('--component', action='append', default=[], help='version, <module_type>/<version> (repeatable)')
    parser.add_argument('--samples', type=int, default=DEFAULT_SAMPLES)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--reviewed-only', action='store_true', help='only versions whose spec has reviewed: true')
    args = parser.parse_args(argv)
    for component in all_versions(args.module):
        if args.reviewed_only and not component.reviewed:
            continue
        if args.component and component.key not in args.component and component.id not in args.component:
            continue
        if component.spec.get('skip') or component.is_supermodule:
            continue
        cm = checks.ComponentModel(component)
        try:
            r = analyse(cm, args.samples, args.seed)
        except harness.MissingParameters as e:
            print(f'{component.key}: skipped ({e})')
            continue
        if r is None:
            print(f'{component.key}: no parameters to sample')
            continue
        lo, hi = r['ci']
        top = ', '.join(f"{x['variable']} ({x['importance']:.2f})" for x in r['ranking'][:3])
        auc = r['classifier'].get('cv_auc')
        print(f"{component.key}: P(fail) = {r['failure_probability']:.3f} [{lo:.3f}, {hi:.3f}] "
              f"over {r['n_samples']} samples; most influential: {top}"
              + (f"; classifier CV AUC {auc:.2f}" if auc is not None else ''))


if __name__ == '__main__':
    main()
