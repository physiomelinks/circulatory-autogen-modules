"""The heart_simple_wCont_smoothfloor period floor: T_wCont = (x + T_min + sqrt((x - T_min)^2 + eps^2))/2.

The equation is read from the module's CellML and evaluated, so the test checks the file itself.
"""
import math
import os
import xml.etree.ElementTree as ET

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CELLML = os.path.join(ROOT, 'modules', 'heart', 'versions', 'vp_wCont_smoothfloor',
                      'heart_vp_wCont_smoothfloor_modules.cellml')
MATHML = '{http://www.w3.org/1998/Math/MathML}'
T_MIN, EPS = 0.25, 0.02

pytestmark = pytest.mark.skipif(not os.path.exists(CELLML), reason='heart vp_wCont_smoothfloor is not in this library')

pytestmark = pytest.mark.skipif(not os.path.exists(CELLML), reason='heart vp_wCont_smoothfloor is not in this library')


def _evaluate(node, env):
    tag = node.tag.replace(MATHML, '')
    if tag == 'ci':
        return env[node.text.strip()]
    if tag == 'cn':
        return float(node.text)
    op, args = node[0].tag.replace(MATHML, ''), [_evaluate(n, env) for n in node[1:]]
    if op == 'plus':
        return sum(args)
    if op == 'minus':
        return args[0] - args[1] if len(args) == 2 else -args[0]
    if op == 'divide':
        return args[0] / args[1]
    if op == 'power':
        return args[0] ** args[1]
    if op == 'root':
        return math.sqrt(args[0])
    raise AssertionError(f'unsupported operator {op}')


def _floor_equation():
    for apply in ET.parse(CELLML).iter(MATHML + 'apply'):
        if apply[0].tag == MATHML + 'eq' and apply[1].tag == MATHML + 'ci' and apply[1].text.strip() == 'T_wCont':
            return apply[2]
    raise AssertionError('no T_wCont equation')


def _t_wcont(x):
    return _evaluate(_floor_equation(), {'T': x, 'Delta_T': 0.0, 'T_min': T_MIN, 'eps_T_floor': EPS})


def test_floor_is_the_square_root_smooth_max():
    for x in (0.05, 0.25, 0.4, 0.85, 1.2):
        assert _t_wcont(x) == pytest.approx((x + T_MIN + math.sqrt((x - T_MIN) ** 2 + EPS ** 2)) / 2, abs=1e-15)


@pytest.mark.parametrize('x, excess_ms', [(0.85, 0.17), (0.6, 0.29), (0.4, 0.66), (0.25, 10.0), (0.2, 1.93), (0.15, 0.99)])
def test_excess_over_hard_max(x, excess_ms):
    assert (_t_wcont(x) - max(x, T_MIN)) * 1e3 == pytest.approx(excess_ms, abs=0.01)


def test_matches_hard_max_above_0p6_s():
    for i in range(0, 241):
        x = 0.6 + 0.005 * i
        assert 0 < _t_wcont(x) - max(x, T_MIN) < 0.3e-3 + 1e-9


def test_converges_to_t_min_from_above():
    values = [_t_wcont(x) for x in (0.2, 0.1, 0.0, -0.5, -2.0, -10.0, -1000.0)]
    assert all(v > T_MIN for v in values)
    assert values == sorted(values, reverse=True)
    assert values[-1] - T_MIN < 1e-6
