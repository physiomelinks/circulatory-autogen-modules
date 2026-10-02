"""
Fixed-step integration of a generated model, for the timestep-convergence verification.

CVODE picks its own steps, so shrinking its maximum step shows solver-tolerance effects
rather than a clean order of accuracy. Here the model's right-hand side (from the same
Myokit model libcuflynx simulates) is compiled into a Python function and integrated with a
classical fixed-step scheme, so the error against a fine reference must shrink as
O(dt^order) when the module is implemented consistently.
"""
import math

import numpy as np
from myokit.formats.python import PythonExpressionWriter

SCHEMES = {
    # name: (order, butcher tableau (a, b, c))
    'euler': (1, ([[]], [1.0], [0.0])),
    'heun': (2, ([[], [1.0]], [0.5, 0.5], [0.0, 1.0])),
    'rk4': (4, ([[], [0.5], [0.0, 0.5], [0.0, 0.0, 1.0]], [1 / 6, 1 / 3, 1 / 3, 1 / 6], [0.0, 0.5, 0.5, 1.0])),
}


def _ident(var):
    return 'v_' + var.qname().replace('.', '__')


def compile_rhs(model):
    '''
    Returns (rhs, y0, state_qnames, observe) for a Myokit model:
      rhs(t, y) -> dy/dt as a numpy array
      observe(t, y, qnames) -> values of any variables at (t, y)
    '''
    model = model.clone()
    states = list(model.states())
    state_qnames = [s.qname() for s in states]
    time_var = model.time()

    writer = PythonExpressionWriter()
    writer.set_lhs_function(lambda lhs: ('d_' if lhs.is_derivative() else '') + _ident(lhs.var()))

    lines = []
    lines.append(f'    {_ident(time_var)} = t')
    for i, s in enumerate(states):
        lines.append(f'    {_ident(s)} = y[{i}]')
    ordered = [eq for comp in model.solvable_order().values() for eq in comp]
    for eq in ordered:
        var = eq.lhs.var()
        if var is time_var or var.binding() is not None:
            continue
        lines.append(f'    {writer.eq(eq)}')
    body = '\n'.join(lines)

    names = {'math': math, 'np': np}
    names.update({k: getattr(math, k) for k in dir(math) if not k.startswith('_')})
    rhs_src = ('def rhs(t, y):\n' + body + '\n'
               '    return np.array([' + ', '.join('d_' + _ident(s) for s in states) + '])\n')
    obs_src = ('def observe(t, y, wanted):\n' + body + '\n'
               '    loc = locals()\n'
               '    return [loc[w] for w in wanted]\n')
    exec(compile(rhs_src, '<cam_rhs>', 'exec'), names)
    exec(compile(obs_src, '<cam_observe>', 'exec'), names)
    rhs, observe_raw = names['rhs'], names['observe']

    def observe(t, y, qnames):
        return observe_raw(t, y, ['v_' + q.replace('.', '__') for q in qnames])

    y0 = np.array(model.initial_values(as_floats=True), dtype=float)
    return rhs, y0, state_qnames, observe


def integrate(rhs, y0, t_end, dt, scheme='rk4', t0=0.0):
    '''Fixed-step explicit Runge-Kutta; returns (t, Y) with Y[i] the state at t[i].'''
    _order, (a, b, c) = SCHEMES[scheme]
    n = int(round((t_end - t0) / dt))
    if not math.isclose(n * dt, t_end - t0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError(f'dt={dt} does not divide the interval {t_end - t0}')
    t = t0 + dt * np.arange(n + 1)
    Y = np.empty((n + 1, len(y0)))
    Y[0] = y = np.array(y0, dtype=float)
    stages = len(b)
    for i in range(n):
        ti = t[i]
        k = []
        for s in range(stages):
            ys = y + dt * sum(a[s][j] * k[j] for j in range(len(a[s]))) if s else y
            k.append(rhs(ti + c[s] * dt, ys))
        y = y + dt * sum(b[s] * k[s] for s in range(stages))
        Y[i + 1] = y
    return t, Y
