"""
Content MathML (as used in CellML 1.1) -> LaTeX, for rendering module equations with KaTeX.

Browsers only render Presentation MathML, so CellML equations are converted here. Anything
not understood is rendered as \\text{?tag} and reported, so gaps are visible rather than
silently wrong.
"""
import re
import xml.etree.ElementTree as ET

MATHML_NS = 'http://www.w3.org/1998/Math/MathML'
CELLML_NS = 'http://www.cellml.org/cellml/1.1#'

# Binding strength, for deciding where brackets are needed.
PREC = {'eq': 1, 'neq': 1, 'lt': 1, 'gt': 1, 'leq': 1, 'geq': 1, 'and': 0, 'or': 0, 'xor': 0,
        'plus': 2, 'minus': 2, 'times': 3, 'divide': 3, 'neg': 4, 'power': 5}
ATOM = 10

RELATIONS = {'eq': '=', 'neq': r'\neq', 'lt': '<', 'gt': '>', 'leq': r'\leq', 'geq': r'\geq'}
LOGIC = {'and': r'\land', 'or': r'\lor', 'xor': r'\oplus'}
FUNCTIONS = {
    'sin': r'\sin', 'cos': r'\cos', 'tan': r'\tan', 'sec': r'\sec', 'csc': r'\csc', 'cot': r'\cot',
    'sinh': r'\sinh', 'cosh': r'\cosh', 'tanh': r'\tanh', 'coth': r'\coth',
    'sech': r'\operatorname{sech}', 'csch': r'\operatorname{csch}',
    'arcsin': r'\arcsin', 'arccos': r'\arccos', 'arctan': r'\arctan',
    'arcsinh': r'\operatorname{arcsinh}', 'arccosh': r'\operatorname{arccosh}',
    'arctanh': r'\operatorname{arctanh}', 'ln': r'\ln', 'exp': r'\exp',
    'min': r'\min', 'max': r'\max', 'rem': r'\operatorname{rem}',
}
CONSTANTS = {'pi': r'\pi', 'exponentiale': 'e', 'true': r'\mathrm{true}', 'false': r'\mathrm{false}',
             'infinity': r'\infty', 'notanumber': r'\mathrm{NaN}'}


def _local(tag):
    return tag.split('}', 1)[-1]


_GREEK = {'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa',
          'lambda', 'mu', 'nu', 'xi', 'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi',
          'omega', 'Gamma', 'Delta', 'Theta', 'Lambda', 'Xi', 'Pi', 'Sigma', 'Phi', 'Psi', 'Omega'}


def variable(name):
    '''u_ra_A -> u_{\\mathrm{ra\\_A}}; alpha -> \\alpha; long names kept upright.'''
    name = name.strip()
    head, _, sub = name.partition('_')
    head_tex = '\\' + head if head in _GREEK else (head if len(head) == 1 else r'\mathrm{' + head.replace('_', r'\_') + '}')
    if not sub:
        return head_tex
    return head_tex + '_{' + r'\mathrm{' + sub.replace('_', r'\_') + '}}'


def number(text):
    text = text.strip()
    m = re.fullmatch(r'([-+]?[0-9.]+)[eE]([-+]?[0-9]+)', text)
    if m:
        return f'{m.group(1)} \\times 10^{{{int(m.group(2))}}}'
    return text


class Converter(object):
    def __init__(self):
        self.unsupported = set()

    def convert(self, el):
        return self._node(el)[0]

    def _wrap(self, tex_prec, min_prec):
        tex, prec = tex_prec
        return f'\\left({tex}\\right)' if prec < min_prec else tex

    def _node(self, el):
        '''Returns (latex, precedence).'''
        tag = _local(el.tag)
        if tag == 'ci':
            return variable(el.text or ''), ATOM
        if tag == 'cn':
            if el.get('type') == 'e-notation':
                parts = [el.text or ''] + [(s.tail or '') for s in el]
                return number(f'{parts[0].strip()}e{parts[1].strip()}'), ATOM
            return number(el.text or ''), ATOM
        if tag in CONSTANTS:
            return CONSTANTS[tag], ATOM
        if tag == 'apply':
            return self._apply(el)
        if tag == 'piecewise':
            return self._piecewise(el)
        if tag in ('math', 'semantics'):
            return ' \\\\ '.join(self.convert(c) for c in el), 0
        self.unsupported.add(tag)
        return r'\text{?' + tag + '}', ATOM

    def _piecewise(self, el):
        rows = []
        for child in el:
            kind = _local(child.tag)
            parts = list(child)
            if kind == 'piece':
                rows.append(f'{self.convert(parts[0])} & \\text{{if }} {self.convert(parts[1])}')
            elif kind == 'otherwise':
                rows.append(f'{self.convert(parts[0])} & \\text{{otherwise}}')
        return r'\begin{cases} ' + r' \\ '.join(rows) + r' \end{cases}', ATOM

    def _apply(self, el):
        children = list(el)
        op = _local(children[0].tag)
        qualifiers = [c for c in children[1:] if _local(c.tag) in ('bvar', 'degree', 'logbase')]
        args = [c for c in children[1:] if c not in qualifiers]

        if op in RELATIONS:
            sides = [self._wrap(self._node(a), 2) for a in args]
            return f' {RELATIONS[op]} '.join(sides), PREC[op]
        if op in LOGIC:
            return f' {LOGIC[op]} '.join(self._wrap(self._node(a), 1) for a in args), PREC[op]
        if op == 'not':
            return r'\lnot ' + self._wrap(self._node(args[0]), ATOM), 4
        if op == 'plus':
            terms = [self._wrap(self._node(args[0]), 2)]
            for a in args[1:]:
                t = self._wrap(self._node(a), 2)
                terms.append(t if t.startswith('-') else '+ ' + t)
            return ' '.join(terms), 2
        if op == 'minus':
            if len(args) == 1:
                return '-' + self._wrap(self._node(args[0]), 3), PREC['neg']
            return f'{self._wrap(self._node(args[0]), 2)} - {self._wrap(self._node(args[1]), 3)}', 2
        if op == 'times':
            return r' \, '.join(self._wrap(self._node(a), 3) for a in args), 3
        if op == 'divide':
            return r'\frac{' + self.convert(args[0]) + '}{' + self.convert(args[1]) + '}', ATOM
        if op == 'power':
            base = self._wrap(self._node(args[0]), ATOM)
            return base + '^{' + self.convert(args[1]) + '}', PREC['power']
        if op == 'root':
            degree = [q for q in qualifiers if _local(q.tag) == 'degree']
            inner = self.convert(args[0])
            if degree:
                return r'\sqrt[' + self.convert(list(degree[0])[0]) + ']{' + inner + '}', ATOM
            return r'\sqrt{' + inner + '}', ATOM
        if op == 'abs':
            return r'\left|' + self.convert(args[0]) + r'\right|', ATOM
        if op in ('floor', 'ceiling'):
            left, right = (r'\lfloor', r'\rfloor') if op == 'floor' else (r'\lceil', r'\rceil')
            return f'\\left{left} {self.convert(args[0])} \\right{right}', ATOM
        if op == 'log':
            base = [q for q in qualifiers if _local(q.tag) == 'logbase']
            sub = '_{' + self.convert(list(base[0])[0]) + '}' if base else '_{10}'
            return r'\log' + sub + r'\left(' + self.convert(args[0]) + r'\right)', ATOM
        if op == 'diff':
            bvar = [q for q in qualifiers if _local(q.tag) == 'bvar'][0]
            var = self.convert(list(bvar)[0])
            return r'\frac{d' + self._wrap(self._node(args[0]), ATOM) + '}{d' + var + '}', ATOM
        if op == 'factorial':
            return self._wrap(self._node(args[0]), ATOM) + '!', ATOM
        if op in FUNCTIONS:
            inner = ', '.join(self.convert(a) for a in args)
            return FUNCTIONS[op] + r'\left(' + inner + r'\right)', ATOM
        self.unsupported.add(op)
        inner = ', '.join(self.convert(a) for a in args)
        return r'\operatorname{?' + op + r'}\left(' + inner + r'\right)', ATOM


def component_equations(cellml_path, component_name):
    '''
    The equations of one component as a list of LaTeX strings, plus the set of MathML tags
    that could not be converted.
    '''
    tree = ET.parse(cellml_path)
    root = tree.getroot()
    converter = Converter()
    equations = []
    for comp in root.iter(f'{{{CELLML_NS}}}component'):
        if comp.get('name') != component_name:
            continue
        for math in comp.iter(f'{{{MATHML_NS}}}math'):
            for eq in math:
                if _local(eq.tag) == 'apply':
                    equations.append(converter.convert(eq))
        break
    return equations, converter.unsupported


def equation_target(eq):
    '''(variable, kind) an equation (a top-level MathML apply) defines: ("x", "ode") for d(x)/dt = ...,
    ("y", "algebraic") for y = ..., (None, None) otherwise (e.g. 0 = f(x)).'''
    children = list(eq)
    if len(children) < 2 or _local(children[0].tag) != 'eq':
        return None, None
    lhs = children[1]
    if _local(lhs.tag) == 'ci':
        return (lhs.text or '').strip(), 'algebraic'
    if _local(lhs.tag) == 'apply' and len(lhs) and _local(lhs[0].tag) == 'diff':
        ci = [c for c in lhs if _local(c.tag) == 'ci']
        if ci:
            return (ci[0].text or '').strip(), 'ode'
    return None, None


def component_equation_targets(cellml_path, component_name):
    '''What each equation of component_equations() defines, in the same order: [(variable, kind)].'''
    root = ET.parse(cellml_path).getroot()
    for comp in root.iter(f'{{{CELLML_NS}}}component'):
        if comp.get('name') != component_name:
            continue
        return [equation_target(eq) for math in comp.iter(f'{{{MATHML_NS}}}math') for eq in math
                if _local(eq.tag) == 'apply']
    return []


def component_variables(cellml_path, component_name):
    '''[(name, units, initial_value or None, interface)] for one component.'''
    root = ET.parse(cellml_path).getroot()
    for comp in root.iter(f'{{{CELLML_NS}}}component'):
        if comp.get('name') == component_name:
            return [(v.get('name'), v.get('units'), v.get('initial_value'),
                     v.get('public_interface') or '')
                    for v in comp.findall(f'{{{CELLML_NS}}}variable')]
    return []


def component_names(cellml_path):
    root = ET.parse(cellml_path).getroot()
    return [c.get('name') for c in root.iter(f'{{{CELLML_NS}}}component')]
