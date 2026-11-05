"""
Compares the maths of module versions in two module libraries, pair by pair: the equations each
version's CellML component defines, compared symbolically with parameters as symbols, so two
libraries can parameterise a module differently and still be checked for the same equations.

    python -m cam_testing.compare_math --pairs counterparts.yaml --a-root CAM --b-root OTHER \\
        [--markdown report.md] [--json report.json] [--strict]

The pairs file (YAML) lists the counterparts, how the second library's names map to the first's, and
the known divergences with their reasons:

    pairs:
      - a: cell/glycolysis/detailed_ArgusUNPUBLISHED_v01      # under <a-root>/modules
        b: cell/glycolysis/coupled_ArgusUNPUBLISHED_v01       # under <b-root>/modules
        rename: {C_ATP_c: C_ATP}                              # b's name -> a's name
        divergences:
          - variable: v_PK
            reason: the cardiomyocyte's PK is not activated by F16BP

"cam" and "here" are accepted for "a" and "b" (a library comparing itself with circulatory-autogen-modules).

An equation is identified by the variable it defines (d/dt x for an ODE). A difference is an equation
only one version has, or one whose right-hand sides differ algebraically. Each difference is explained
(listed in divergences), unexplained, or a listed divergence is stale (no longer differs). Parameter
values (the default parameterisations) are compared for information only. The report never fails
unless --strict is given, which exits 1 on an unexplained difference or a stale divergence.
"""
import argparse
import csv
import glob
import json
import os
import sys
import xml.etree.ElementTree as ET

import sympy
import yaml

from cam_testing.mathml import CELLML_NS, MATHML_NS, _local, equation_target

_FUNCTIONS = {
    'sin': sympy.sin, 'cos': sympy.cos, 'tan': sympy.tan, 'sec': sympy.sec, 'csc': sympy.csc, 'cot': sympy.cot,
    'sinh': sympy.sinh, 'cosh': sympy.cosh, 'tanh': sympy.tanh, 'coth': sympy.coth,
    'arcsin': sympy.asin, 'arccos': sympy.acos, 'arctan': sympy.atan, 'arcsinh': sympy.asinh,
    'arccosh': sympy.acosh, 'arctanh': sympy.atanh, 'exp': sympy.exp, 'ln': sympy.log, 'abs': sympy.Abs,
    'floor': sympy.floor, 'ceiling': sympy.ceiling,
}
_RELATIONS = {'eq': sympy.Eq, 'neq': sympy.Ne, 'lt': sympy.Lt, 'gt': sympy.Gt, 'leq': sympy.Le, 'geq': sympy.Ge}
_CONSTANTS = {'pi': sympy.pi, 'exponentiale': sympy.E, 'true': sympy.true, 'false': sympy.false,
              'infinity': sympy.oo, 'notanumber': sympy.nan}


def _number(text):
    '''An exact rational for a CellML number, so 2.0 and 2 compare equal.'''
    return sympy.Rational(text.strip())


def to_sympy(el, rename=None):
    '''A content MathML element as a sympy expression; variable names pass through rename (old -> new).'''
    rename = rename or {}
    tag = _local(el.tag)
    if tag == 'ci':
        name = (el.text or '').strip()
        return sympy.Symbol(rename.get(name, name))
    if tag == 'cn':
        if el.get('type') == 'e-notation':
            parts = [el.text or ''] + [(s.tail or '') for s in el]
            return _number(f'{parts[0].strip()}e{parts[1].strip()}')
        return _number(el.text or '')
    if tag in _CONSTANTS:
        return _CONSTANTS[tag]
    if tag == 'piecewise':
        pieces = []
        for child in el:
            parts = list(child)
            if _local(child.tag) == 'piece':
                pieces.append((to_sympy(parts[0], rename), to_sympy(parts[1], rename)))
            elif _local(child.tag) == 'otherwise':
                pieces.append((to_sympy(parts[0], rename), True))
        return sympy.Piecewise(*pieces)
    if tag != 'apply':
        raise ValueError(f'unsupported MathML element <{tag}>')
    children = list(el)
    op = _local(children[0].tag)
    qualifiers = {_local(c.tag): c for c in children[1:] if _local(c.tag) in ('bvar', 'degree', 'logbase')}
    args = [to_sympy(c, rename) for c in children[1:] if _local(c.tag) not in ('bvar', 'degree', 'logbase')]
    if op == 'plus':
        return sympy.Add(*args)
    if op == 'minus':
        return -args[0] if len(args) == 1 else args[0] - args[1]
    if op == 'times':
        return sympy.Mul(*args)
    if op == 'divide':
        return args[0] / args[1]
    if op == 'power':
        return args[0] ** args[1]
    if op == 'root':
        degree = to_sympy(list(qualifiers['degree'])[0], rename) if 'degree' in qualifiers else 2
        return args[0] ** (1 / sympy.sympify(degree))
    if op == 'log':
        base = to_sympy(list(qualifiers['logbase'])[0], rename) if 'logbase' in qualifiers else 10
        return sympy.log(args[0], base)
    if op in ('min', 'max'):
        return (sympy.Min if op == 'min' else sympy.Max)(*args)
    if op in _FUNCTIONS:
        return _FUNCTIONS[op](*args)
    if op in _RELATIONS:
        return _RELATIONS[op](*args)
    if op in ('and', 'or', 'xor'):
        return {'and': sympy.And, 'or': sympy.Or, 'xor': sympy.Xor}[op](*args)
    if op == 'not':
        return sympy.Not(args[0])
    raise ValueError(f'unsupported MathML operator <{op}>')


def version_dir(root, path):
    '''<root>/modules/<module_type>/versions/<version> for a path "module_type/version"
    (module_type may be nested, e.g. cell/ion_channels/i_leak_Ca).'''
    module_type, version = path.rstrip('/').rsplit('/', 1)
    return os.path.join(root, 'modules', module_type, 'versions', version)


def version_equations(vdir, rename=None):
    '''{key: expression} for a version's CellML component: key "x" for x = f (expression f),
    "d/dt x" for an ODE, "0 = f #i" for an implicit equation (expression lhs - rhs).'''
    rename = rename or {}
    configs = glob.glob(os.path.join(vdir, '*_modules_config.json'))
    if not configs:
        raise FileNotFoundError(f'no module version (no *_modules_config.json) in {vdir}')
    with open(configs[0]) as f:
        config = json.load(f)[0]
    root = ET.parse(os.path.join(vdir, config['component_file'])).getroot()
    out = {}
    for comp in root.iter(f'{{{CELLML_NS}}}component'):
        if comp.get('name') != config['component_type']:
            continue
        implicit = 0
        for math in comp.iter(f'{{{MATHML_NS}}}math'):
            for eq in math:
                if _local(eq.tag) != 'apply':
                    continue
                target, kind = equation_target(eq)
                lhs, rhs = list(eq)[1], list(eq)[2]
                if target is None:
                    out[f'0 = f #{implicit}'] = to_sympy(lhs, rename) - to_sympy(rhs, rename)
                    implicit += 1
                    continue
                name = rename.get(target, target)
                out[f'd/dt {name}' if kind == 'ode' else name] = to_sympy(rhs, rename)
    return out


def default_parameters(vdir, rename=None):
    '''{name: value} from the version's default parameterisation (CSV rows with a numeric value).'''
    rename = rename or {}
    values = {}
    for path in glob.glob(os.path.join(vdir, 'parameterisations', 'default', '*_parameters.csv')):
        with open(path) as f:
            for row in csv.DictReader(f):
                try:
                    values[rename.get(row['variable_name'], row['variable_name'])] = float(row['value'])
                except (KeyError, TypeError, ValueError):
                    continue
    return values


def same(a, b):
    '''Whether two expressions are equal for every value of their symbols.'''
    if a == b:
        return True
    difference = a - b
    if sympy.expand(difference) == 0:
        return True
    return sympy.simplify(difference) == 0


def compare_pair(a_root, b_root, pair):
    '''The differences between one pair of versions: {"a", "b", "differences": [...], "parameters": [...]}.'''
    pair = {'a': pair.get('a', pair.get('cam')), 'b': pair.get('b', pair.get('here')), **pair}
    rename = pair.get('rename') or {}
    a_dir, b_dir = version_dir(a_root, pair['a']), version_dir(b_root, pair['b'])
    a_eqs, b_eqs = version_equations(a_dir), version_equations(b_dir, rename)
    listed = {d['variable']: d.get('reason', '') for d in pair.get('divergences') or []}
    differences = []
    for key in sorted(set(a_eqs) | set(b_eqs)):
        if key not in b_eqs:
            kind = f"only in {pair['a']}"
        elif key not in a_eqs:
            kind = f"only in {pair['b']}"
        elif not same(a_eqs[key], b_eqs[key]):
            kind = 'differs'
        else:
            continue
        variable = key[len('d/dt '):] if key.startswith('d/dt ') else key
        reason = listed.get(key, listed.get(variable))
        differences.append({'equation': key, 'kind': kind, 'status': 'explained' if reason else 'unexplained',
                            'reason': reason or '', 'a': sympy.sstr(a_eqs[key]) if key in a_eqs else None,
                            'b': sympy.sstr(b_eqs[key]) if key in b_eqs else None})
    found = {d['equation'] for d in differences} | {d['equation'][len('d/dt '):] for d in differences
                                                     if d['equation'].startswith('d/dt ')}
    for variable, reason in listed.items():
        if variable not in found:
            differences.append({'equation': variable, 'kind': 'listed, no longer differs', 'status': 'stale',
                                'reason': reason, 'a': None, 'b': None})
    a_par, b_par = default_parameters(a_dir), default_parameters(b_dir, rename)
    parameters = [{'name': n, 'a': a_par[n], 'b': b_par[n]} for n in sorted(set(a_par) & set(b_par))
                  if not (a_par[n] == b_par[n] or abs(a_par[n] - b_par[n]) <= 1e-12 * max(abs(a_par[n]), abs(b_par[n])))]
    return {'a': pair['a'], 'b': pair['b'], 'differences': differences, 'parameters': parameters}


def compare(pairs_file, a_root, b_root):
    with open(pairs_file) as f:
        pairs = yaml.safe_load(f)['pairs']
    return [compare_pair(a_root, b_root, pair) for pair in pairs]


def markdown(results, a_label='a', b_label='b'):
    lines = ['# Maths comparison', '']
    for r in results:
        counts = {s: sum(d['status'] == s for d in r['differences']) for s in ('unexplained', 'explained', 'stale')}
        lines += [f"## `{r['a']}` ({a_label}) vs `{r['b']}` ({b_label})", '',
                  f"{counts['unexplained']} unexplained, {counts['explained']} explained, {counts['stale']} stale; "
                  f"{len(r['parameters'])} parameter values differ (information only).", '']
        if r['differences']:
            lines += ['| Equation | Status | What | Reason |', '|---|---|---|---|']
            for d in r['differences']:
                what = d['kind']
                if d['kind'] == 'differs':
                    what = f"{a_label}: `{d['a']}`<br>{b_label}: `{d['b']}`"
                lines.append(f"| `{d['equation']}` | {d['status']} | {what} | {d['reason']} |")
            lines.append('')
        if r['parameters']:
            lines += ['<details><summary>Parameter values that differ</summary>', '',
                      f'| Parameter | {a_label} | {b_label} |', '|---|---|---|']
            lines += [f"| `{p['name']}` | {p['a']:.6g} | {p['b']:.6g} |" for p in r['parameters']]
            lines += ['', '</details>', '']
    return '\n'.join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.strip().split('\n\n')[0])
    p.add_argument('--pairs', required=True, help='YAML file of counterpart versions')
    p.add_argument('--a-root', required=True, help='root of the first library (has modules/)')
    p.add_argument('--b-root', required=True, help='root of the second library (has modules/)')
    p.add_argument('--a-label', default='a')
    p.add_argument('--b-label', default='b')
    p.add_argument('--markdown', help='write the report as Markdown here (default: stdout)')
    p.add_argument('--json', help='write the results as JSON here')
    p.add_argument('--strict', action='store_true', help='exit 1 on an unexplained difference or a stale divergence')
    args = p.parse_args(argv)
    results = compare(args.pairs, args.a_root, args.b_root)
    report = markdown(results, args.a_label, args.b_label)
    if args.markdown:
        with open(args.markdown, 'w') as f:
            f.write(report + '\n')
    else:
        print(report)
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(results, f, indent=2)
    bad = sum(d['status'] in ('unexplained', 'stale') for r in results for d in r['differences'])
    return 1 if args.strict and bad else 0


if __name__ == '__main__':
    sys.exit(main())
