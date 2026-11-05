"""Tiny infix -> CellML 1.1 content-MathML converter, for writing module files in the library's style.

Infix forms: + - * / ** with variables and numbers (dimensionless), exp(x), ln(x), abs(x), pi, and
u(value, units) for a number carrying units (a conversion factor such as u(1e-9, A_per_nanoA)), so the
equations stay unit-consistent.
"""
import ast
from xml.sax.saxutils import escape

OPS = {ast.Add: 'plus', ast.Sub: 'minus', ast.Mult: 'times', ast.Div: 'divide', ast.Pow: 'power'}
FUNCS = {'exp': 'exp', 'ln': 'ln', 'abs': 'abs'}


def _cn(value, units):
    """A number in the library's form: plain decimals, or type="e-notation" with <sep/> for exponents."""
    text = repr(float(value))
    if 'e' in text:
        mantissa, exponent = text.split('e')
        return f'<cn cellml:units="{units}" type="e-notation">{mantissa}<sep/>{int(exponent)}</cn>'
    return f'<cn cellml:units="{units}">{text}</cn>'


def _num(v):
    return _cn(v, 'dimensionless')


def _expr(node, ind):
    pad = ' ' * ind
    if isinstance(node, ast.BinOp):
        # flatten chains of + and * like the library's hand-written files
        op = OPS[type(node.op)]
        args = _flatten(node, type(node.op)) if op in ('plus', 'times') else [node.left, node.right]
        inner = '\n'.join(_expr(a, ind + 4) for a in args)
        return f'{pad}<apply>\n{pad}    <{op}/>\n{inner}\n{pad}</apply>'
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return f'{pad}<apply>\n{pad}    <minus/>\n{_expr(node.operand, ind + 4)}\n{pad}</apply>'
    if isinstance(node, ast.Name):
        if node.id == 'pi':
            return f'{pad}<pi/>'
        return f'{pad}<ci>{node.id}</ci>'
    if isinstance(node, ast.Constant):
        return pad + _num(float(node.value))
    if isinstance(node, ast.Call) and node.func.id == 'u':
        value, units = node.args
        return pad + _cn(value.value, units.id)
    if isinstance(node, ast.Call) and node.func.id in FUNCS:
        (arg,) = node.args
        return f'{pad}<apply>\n{pad}    <{FUNCS[node.func.id]}/>\n{_expr(arg, ind + 4)}\n{pad}</apply>'
    raise ValueError(f'unsupported expression node {ast.dump(node)}')


def _flatten(node, op):
    if isinstance(node, ast.BinOp) and isinstance(node.op, op):
        return _flatten(node.left, op) + _flatten(node.right, op)
    return [node]


def equation(lhs, rhs, ind=12):
    pad = ' ' * ind
    body = _expr(ast.parse(rhs, mode='eval').body, ind + 4)
    if lhs.startswith('d/dt '):
        var = lhs.split()[1]
        left = (f'{pad}    <apply>\n{pad}        <diff/>\n{pad}        <bvar>\n{pad}            <ci>t</ci>\n'
                f'{pad}        </bvar>\n{pad}        <ci>{var}</ci>\n{pad}    </apply>')
    else:
        left = f'{pad}    <ci>{lhs}</ci>'
    return f'{pad}<apply>\n{pad}    <eq/>\n{left}\n{body}\n{pad}</apply>'


def component_model(model_name, component, variables, equations):
    """variables: [(name, units, interface, initial_value or None, comment or None)]"""
    lines = ["<?xml version='1.0' encoding='UTF-8'?>",
             f'<model name="{model_name}" xmlns="http://www.cellml.org/cellml/1.1#" '
             'xmlns:cellml="http://www.cellml.org/cellml/1.1#">',
             f'    <component name="{component}">']
    for name, units, iface, init, comment in variables:
        if comment:
            lines.append(f'        <!-- {escape(comment)} -->')
        init_attr = f'initial_value="{init}" ' if init is not None else ''
        lines.append(f'        <variable {init_attr}name="{name}" public_interface="{iface}" units="{units}"/>')
    lines.append('        <math xmlns="http://www.w3.org/1998/Math/MathML">')
    for lhs, rhs in equations:
        lines.append(equation(lhs, rhs))
    lines += ['        </math>', '    </component>', '</model>', '']
    return '\n'.join(lines)


def units_model(defs):
    """defs: {name: [(units, prefix, exponent, multiplier)]}"""
    lines = ["<?xml version='1.0' encoding='UTF-8'?>",
             '<model name="Units" xmlns="http://www.cellml.org/cellml/1.1#" xmlns:cellml="http://www.cellml.org/cellml/1.1#">']
    for name, parts in defs.items():
        lines.append(f'    <units name="{name}">')
        for u, prefix, exp, mult in parts:
            attrs = ''.join([f' exponent="{exp:g}"' if exp != 1 else '', f' multiplier="{mult:g}"' if mult != 1 else '',
                             f' prefix="{prefix}"' if prefix else ''])
            lines.append(f'        <unit{attrs} units="{u}"/>')
        lines.append('    </units>')
    lines += ['</model>', '']
    return '\n'.join(lines)
