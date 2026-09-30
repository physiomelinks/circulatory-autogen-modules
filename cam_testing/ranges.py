"""
Parameter ranges a module has been tested and validated over, and a check of a system
model's parameter values against them.

    tested range     the range each parameter is varied over by verification_test_BC and the
                     joint risk analysis: bc_sweep.ranges in <name>_tests.yaml, or
                     bc_sweep.factors x the nominal value
    validated spread the values at which the component matched data in a passing validation
                     test (published parameter intervals and calibrated values)

    python -m cam_testing.ranges check --vessel-array sys_vessel_array.csv --parameters sys_parameters.csv
        lists each value as inside / outside the tested range and the validated spread, with
        the estimated failure risk near the values (from `make risk`)
"""
import argparse
import csv
import sys

from cam_testing import checks
from cam_testing.library import load_module, module_names


def tested_range(component, param):
    """(lo, hi, log) the parameter is varied over in the sweep and the risk analysis, or None."""
    sweep = component.spec.get('bc_sweep') or {}
    if param.variable_name in (sweep.get('exclude') or []) or param.is_todo:
        return None
    r = (sweep.get('ranges') or {}).get(param.variable_name)
    if isinstance(r, dict):
        return float(r['min']), float(r['max']), r.get('scale') == 'log'
    if r is not None:
        return float(r[0]), float(r[1]), False
    nominal = param.float_value
    if nominal == 0:
        return None
    factors = [float(f) for f in sweep.get('factors', [0.5, 2.0])]
    lo, hi = sorted((nominal * min(factors), nominal * max(factors)))
    return lo, hi, lo > 0


def validated_spread(component):
    """{variable: [min, max]} over the passing validation tests (needs their results)."""
    spans = {}

    def add(var, lo, hi):
        spans[var] = [min(spans[var][0], lo), max(spans[var][1], hi)] if var in spans else [lo, hi]

    for test, kind in (('validation_test_baseline', 'baseline'), ('validation_test_calibrate', 'calibrate')):
        r = checks.load(component, test)
        if r is None or r.status != checks.PASSED:
            continue
        spec = (component.spec.get('validation') or {}).get(kind) or {}
        for var, (lo, hi) in (spec.get('parameter_ranges') or {}).items():
            add(var, float(lo), float(hi))
        for var, val in (r.metrics.get('validated_values') or {}).items():
            for x in (val if isinstance(val, list) else [val]):
                add(var, float(x), float(x))
    return spans


def _module_index():
    index = {}
    for name in module_names():
        for component in load_module(name).components():
            index[(component.vessel_type, component.BC_type)] = component
    return index


def check(vessel_array_path, parameters_path, out=sys.stdout):
    """Returns the number of values outside the tested range."""
    from cam_testing import risk
    from cam_testing import vessel_array
    vessels = [{'name': r['name'], 'vessel_type': r.get('module_type'), 'BC_type': r.get('module_subtype')}
               for r in vessel_array.read_records(vessel_array_path)]
    with open(parameters_path) as f:
        values = {r['variable_name'].strip(): r['value'].strip()
                  for r in csv.DictReader(f, skipinitialspace=True) if r.get('variable_name')}
    index = _module_index()
    rows, outside, unknown, risks = [], 0, set(), []
    for vessel in vessels:
        component = index.get((vessel['vessel_type'], vessel['BC_type']))
        if component is None:
            unknown.add((vessel['vessel_type'], vessel['BC_type']))
            continue
        spread = validated_spread(component)
        mine = {}
        for p in component.parameters():
            name = p.variable_name if p.vessel_type == 'global' else f"{p.variable_name}_{vessel['name']}"
            try:
                val = float(values[name])
            except (KeyError, ValueError):
                continue
            mine[p.variable_name] = val
            tr = tested_range(component, p)
            va = spread.get(p.variable_name)
            in_t = None if tr is None else tr[0] <= val <= tr[1]
            in_v = None if va is None else va[0] <= val <= va[1]
            outside += in_t is False
            rows.append((name, val, component.module.name, tr, in_t, va, in_v))
        lr = risk.local_risk(component, mine)
        if lr is not None:
            risks.append((vessel['name'], component.module.name, lr))

    def span(r):
        return '—' if r is None else f'[{r[0]:.4g}, {r[1]:.4g}]'

    def mark(ok):
        return 'no range' if ok is None else ('inside' if ok else 'OUTSIDE')

    print(f'{"parameter":32s} {"value":>12s}  {"module":24s} {"tested range":24s} {"":9s} {"validated":24s}', file=out)
    for name, val, mod, tr, in_t, va, in_v in sorted(rows, key=lambda r: (r[4] is not False, r[0])):
        print(f'{name:32s} {val:12.4g}  {mod:24s} {span(tr):24s} {mark(in_t):9s} {span(va):24s} {mark(in_v)}', file=out)
    if risks:
        print('\nestimated failure risk at these values (fitted classifier; nearest stored samples):', file=out)
        for vname, mod, lr in sorted(risks, key=lambda r: -r[2]['risk']):
            note = ' (outside the sampled box: extrapolated)' if lr['outside_box'] else ''
            model = (f"model {lr['model']:.2f}" + (f" (CV AUC {lr['cv_auc']:.2f})" if lr.get('cv_auc') else '')
                     if 'model' in lr else 'no model')
            print(f"  {vname:24s} {mod:24s} {model};  {lr['k']} nearest samples: {lr['risk']:.2f} "
                  f"[{lr['ci'][0]:.2f}, {lr['ci'][1]:.2f}]{note}", file=out)
    else:
        print('\nno failure-risk samples found: run `make risk MODULE=<name>` for the modules used', file=out)
    if unknown:
        print(f'\nnot in this library: {sorted(unknown)}', file=out)
    print(f'\n{outside} value(s) outside a tested range; {len(rows)} checked', file=out)
    return outside


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='cmd', required=True)
    ch = sub.add_parser('check', help="check a system model's parameters against module ranges")
    ch.add_argument('--vessel-array', required=True)
    ch.add_argument('--parameters', required=True)
    args = parser.parse_args(argv)
    return 1 if check(args.vessel_array, args.parameters) else 0


if __name__ == '__main__':
    sys.exit(main())
