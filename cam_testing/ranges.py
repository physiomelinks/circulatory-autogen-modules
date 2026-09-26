"""
Verified and validated parameter ranges.

Each row of <name>_parameters.csv carries four range columns:

    verified_min, verified_max    the span of values around the nominal for which the module
                                  ran, stayed finite and kept its invariants in
                                  verification_test_BC (one parameter varied at a time)
    validated_min, validated_max  the span of values at which the module matched data in a
                                  passing validation test (published parameter intervals and
                                  calibrated values)

Commands:

    python -m cam_testing.ranges update [--module NAME]
        write the ranges from the latest test results into the parameters files (commit them)

    python -m cam_testing.ranges check --vessel-array sys_vessel_array.csv --parameters sys_parameters.csv
        check a system model's parameter values against the ranges of the modules it uses
"""
import argparse
import csv
import os
import sys

from cam_testing import checks
from cam_testing.library import PARAMETER_COLUMNS, load_module, module_names

RANGE_COLUMNS = ['verified_min', 'verified_max', 'validated_min', 'validated_max']


def verified_ranges_from_sweep(runs, nominal):
    '''
    runs: [(value, passed)] for one parameter. The verified range is the contiguous run of
    passing values that contains the nominal value (so a failure between two passing values
    ends the range).
    '''
    pts = sorted(runs)
    values = [v for v, _ in pts]
    if not pts:
        return None
    # index of the value closest to nominal
    i0 = min(range(len(values)), key=lambda i: abs(values[i] - nominal))
    if not pts[i0][1]:
        return None
    lo = hi = i0
    while lo - 1 >= 0 and pts[lo - 1][1]:
        lo -= 1
    while hi + 1 < len(pts) and pts[hi + 1][1]:
        hi += 1
    return [values[lo], values[hi]]


def validated_ranges(component):
    '''{variable_name: [min, max]} over every passing validation test of the component.'''
    spans = {}

    def add(var, lo, hi):
        if var in spans:
            spans[var] = [min(spans[var][0], lo), max(spans[var][1], hi)]
        else:
            spans[var] = [lo, hi]

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


def _fmt(x):
    return '' if x is None else f'{x:.6g}'


def update(name):
    module = load_module(name)
    by_key = {}
    for component in module.components():
        bc = checks.load(component, 'verification_test_BC')
        verified = (bc.metrics.get('verified_ranges') or {}) if bc else {}
        validated = validated_ranges(component)
        for p in component.parameters():
            key = (p.vessel_type, p.BC_type, p.variable_name)
            vr = verified.get(p.variable_name)
            va = validated.get(p.variable_name)
            if p.vessel_type == 'global':
                # a global constant is shared: keep the intersection of what each component verified
                prev = by_key.get(key)
                if prev and prev[0] is not None and vr is not None:
                    vr = [max(prev[0][0], vr[0]), min(prev[0][1], vr[1])]
                elif prev and vr is None:
                    vr = prev[0]
            by_key[key] = (vr, va)
    rows = []
    for p in module.parameters:
        vr, va = by_key.get((p.vessel_type, p.BC_type, p.variable_name), (None, None))
        row = {k: getattr(p, k) for k in PARAMETER_COLUMNS if k not in RANGE_COLUMNS}
        row.update({'verified_min': _fmt(vr[0]) if vr else p.verified_min,
                    'verified_max': _fmt(vr[1]) if vr else p.verified_max,
                    'validated_min': _fmt(va[0]) if va else p.validated_min,
                    'validated_max': _fmt(va[1]) if va else p.validated_max})
        rows.append(row)
    path = os.path.join(module.dir, f'{name}_parameters.csv')
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=PARAMETER_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return path


def recorded_range(param, which):
    lo, hi = getattr(param, f'{which}_min'), getattr(param, f'{which}_max')
    if lo in ('', None) or hi in ('', None):
        return None
    return float(lo), float(hi)


def _module_index():
    '''(vessel_type, BC_type) -> (module, component) across the library.'''
    index = {}
    for name in module_names():
        for component in load_module(name).components():
            index[(component.vessel_type, component.BC_type)] = component
    return index


def check(vessel_array_path, parameters_path, out=sys.stdout):
    '''
    Compares a system model's parameter values with the verified and validated ranges of
    the modules its vessels use. Returns the number of values outside a verified range.
    '''
    with open(vessel_array_path) as f:
        vessels = [{k.strip(): (v or '').strip() for k, v in r.items()} for r in csv.DictReader(f, skipinitialspace=True)]
    with open(parameters_path) as f:
        values = {r['variable_name'].strip(): r['value'].strip()
                  for r in csv.DictReader(f, skipinitialspace=True) if r.get('variable_name')}
    index = _module_index()
    from cam_testing import risk
    rows, outside, unknown_modules, risks = [], 0, set(), []
    for vessel in vessels:
        component = index.get((vessel['vessel_type'], vessel['BC_type']))
        if component is None:
            unknown_modules.add((vessel['vessel_type'], vessel['BC_type']))
            continue
        vals = {}
        for p in component.parameters():
            name = p.variable_name if p.vessel_type == 'global' else f"{p.variable_name}_{vessel['name']}"
            try:
                vals[p.variable_name] = float(values[name])
            except (KeyError, ValueError):
                pass
        lr = risk.local_risk(component, vals)
        if lr is not None:
            risks.append((vessel['name'], component.module.name, lr))
        for p in component.parameters():
            name = p.variable_name if p.vessel_type == 'global' else f"{p.variable_name}_{vessel['name']}"
            if name not in values:
                continue
            try:
                val = float(values[name])
            except ValueError:
                continue
            ver, vad = recorded_range(p, 'verified'), recorded_range(p, 'validated')
            in_ver = None if ver is None else ver[0] <= val <= ver[1]
            in_vad = None if vad is None else vad[0] <= val <= vad[1]
            if in_ver is False:
                outside += 1
            rows.append((name, val, component.module.name, ver, in_ver, vad, in_vad))

    def span(r):
        return '—' if r is None else f'[{r[0]:.4g}, {r[1]:.4g}]'

    def mark(ok):
        return 'no range' if ok is None else ('inside' if ok else 'OUTSIDE')

    print(f'{"parameter":32s} {"value":>12s}  {"module":24s} {"verified":24s} {"":9s} {"validated":24s}', file=out)
    for name, val, mod, ver, in_ver, vad, in_vad in sorted(rows, key=lambda r: (r[4] is not False, r[0])):
        print(f'{name:32s} {val:12.4g}  {mod:24s} {span(ver):24s} {mark(in_ver):9s} {span(vad):24s} {mark(in_vad)}', file=out)
    if risks:
        print('\nestimated failure risk near these values (nearest samples of `make risk`):', file=out)
        for vname, mod, lr in sorted(risks, key=lambda r: -r[2]['risk']):
            note = ' (outside the sampled box: extrapolated)' if lr['outside_box'] else ''
            print(f"  {vname:24s} {mod:24s} {lr['risk']:.2f}  [{lr['ci'][0]:.2f}, {lr['ci'][1]:.2f}] from {lr['k']} samples{note}", file=out)
    if unknown_modules:
        print(f'\nnot in this library: {sorted(unknown_modules)}', file=out)
    print(f'\n{outside} value(s) outside a verified range; {len(rows)} checked', file=out)
    return outside


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='cmd', required=True)
    up = sub.add_parser('update', help='write ranges from test results into <name>_parameters.csv')
    up.add_argument('--module', action='append', default=[])
    ch = sub.add_parser('check', help="check a system model's parameters against module ranges")
    ch.add_argument('--vessel-array', required=True)
    ch.add_argument('--parameters', required=True)
    args = parser.parse_args(argv)
    if args.cmd == 'update':
        for name in args.module or module_names():
            print(f'updated {update(name)}')
        return 0
    return 1 if check(args.vessel_array, args.parameters) else 0


if __name__ == '__main__':
    sys.exit(main())
