'''A small builder for the CellML 1.1 MathML the lumped constitutive modules are written in.'''


def ci(name):
    return f'<ci>{name}</ci>'


def cn(value, units='dimensionless'):
    return f'<cn cellml:units="{units}">{value}</cn>'


def apply(op, *args):
    return f'<apply><{op}/>' + ''.join(args) + '</apply>'


def plus(*args):
    return apply('plus', *args) if len(args) > 1 else args[0]


def minus(a, b=None):
    return apply('minus', a) if b is None else apply('minus', a, b)


def times(*args):
    return apply('times', *args) if len(args) > 1 else args[0]


def divide(a, b):
    return apply('divide', a, b)


def power(a, b):
    return apply('power', a, b)


def eq(lhs, rhs):
    return apply('eq', lhs, rhs)


def ode(state, rhs, t='t'):
    return f'<apply><eq/><apply><diff/><bvar>{ci(t)}</bvar>{ci(state)}</apply>{rhs}</apply>'


def component(name, variables, equations):
    '''variables: [(name, units, interface, initial_value or None)].'''
    lines = [f'    <component name="{name}">']
    for var, units, interface, init in variables:
        init_attr = f' initial_value="{init}"' if init is not None else ''
        lines.append(f'        <variable{init_attr} name="{var}" public_interface="{interface}" units="{units}"/>')
    lines.append('        <math xmlns="http://www.w3.org/1998/Math/MathML">')
    for e in equations:
        lines.append('            ' + e)
    lines.append('        </math>')
    lines.append('    </component>')
    return '\n'.join(lines)


def model(name, components):
    return ("<?xml version='1.0' encoding='UTF-8'?>\n"
            f'<model name="{name}" xmlns="http://www.cellml.org/cellml/1.1#" '
            'xmlns:cellml="http://www.cellml.org/cellml/1.1#">\n'
            + '\n'.join(components) + '\n</model>\n')
