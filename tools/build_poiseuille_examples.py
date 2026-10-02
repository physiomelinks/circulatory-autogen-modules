"""
Builds the Poiseuille transport example system models (modules/system/poiseuille/): straight square
channels of poiseuille_transport volumes and faces between an inlet and an outlet reservoir.

  square_channel_advection   steady pressure drop: advection and diffusion of a solute from the inlet
                             reservoir (C = 1 mM) to the outlet (C = 0); checked against the exact
                             Poiseuille flow of a square duct and the exact steady 1D solution
  square_channel_diffusion   no pressure drop: pure diffusion, a linear steady profile
  square_channel_pulsatile   sinusoidal inlet pressure whose amplitude reverses the flow; checked for
                             fluid and solute conservation, the bounds of the concentration, and (in
                             tests/test_poiseuille_channels.py) against a direct scipy solve

The channel is N cells of length dx; face k joins cell k-1 to cell k (face 0 joins the inlet reservoir
to cell 0, face N cell N-1 to the outlet), with d_up = d_down = dx/2 inside and 0 on the reservoir side.
Rewrites the vessel array and parameters, and writes the spec only if there isn't one.

    python tools/build_poiseuille_examples.py
"""
import csv
import os

import numpy as np
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'modules', 'system', 'poiseuille')

SIDE = 2.0e-5        # m, square cross-section side
LENGTH = 2.0e-4      # m
N = 20               # cells
MU = 1.0e-3          # Pa s, water at 20 C
SIGMA = 1.0e-9       # m2/s, a small solute in water
C_V = 1.0e-17        # m3/Pa per cell: compliance, a regularisation (hydraulic time constant ~1e-3 s)


def square_duct_factor(terms=200):
    '''Q = a^4 dP / (12 mu L) (1 - 192/pi^5 sum_{n odd} tanh(n pi/2)/n^5) = a^4 dP / (K mu L); returns K.'''
    n = np.arange(1, 2 * terms, 2, dtype=float)
    return 12.0 / (1.0 - 192.0 / np.pi ** 5 * np.sum(np.tanh(n * np.pi / 2) / n ** 5))


K_SQUARE = square_duct_factor()


def geometry(n=N, length=LENGTH, side=SIDE):
    dx = length / n
    faces = [(0.0 if k == 0 else dx / 2, 0.0 if k == n else dx / 2) for k in range(n + 1)]
    return dx, side ** 2, faces


def exact_steady(P_in, P_out, C_in, C_out, n=N, length=LENGTH, side=SIDE, mu=MU, sigma=SIGMA):
    '''Exact steady flow, cell-centre concentrations and solute flux of the straight channel.'''
    dx, A, _ = geometry(n, length, side)
    Q = A ** 2 * (P_in - P_out) / (K_SQUARE * mu * length)
    x = (np.arange(n) + 0.5) * dx
    Pe = Q * length / (sigma * A)
    if abs(Pe) < 1e-12:
        C = C_in + (C_out - C_in) * x / length
        J = sigma * A * (C_in - C_out) / length
    else:
        C = C_in + (C_out - C_in) * np.expm1(Pe * x / length) / np.expm1(Pe)
        J = Q * (C_in * np.exp(Pe) - C_out) / np.expm1(Pe)
    return Q, C, J


def bernoulli(x):
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-4
    safe = np.where(small, 1.0, x)
    out = np.where(safe > 0, safe * np.exp(-np.abs(safe)) / -np.expm1(-np.abs(safe)), safe / np.expm1(safe))
    return np.where(small, 1 - x / 2 + x ** 2 / 12, out)


def scipy_solve(case, t, n=N, length=LENGTH, side=SIDE, mu=MU, sigma=SIGMA, c_v=C_V):
    '''Direct solve of the same finite-volume equations, independent of libcuflynx, for checking the model.'''
    from scipy.integrate import solve_ivp
    dx, A, faces = geometry(n, length, side)
    L = np.array([du + dd for du, dd in faces])
    G = A ** 2 / (K_SQUARE * mu * L)
    V_P = A * dx

    def p_in(tt):
        return case['P_in'] + case['P_amp'] * np.sin(2 * np.pi * case['f_P'] * tt)

    def rhs(tt, y):
        P, C = y[:n], y[n:]
        Pf = np.concatenate([[p_in(tt)], P, [case['P_out']]])
        Cf = np.concatenate([[case['C_in']], C, [case['C_out']]])
        dP = Pf[:-1] - Pf[1:]
        Q = G * dP
        D = sigma * A / L
        Pe = A * dP / (K_SQUARE * mu * sigma)
        J = D * (bernoulli(-Pe) * Cf[:-1] - bernoulli(Pe) * Cf[1:])
        qnet, jnet = Q[:-1] - Q[1:], J[:-1] - J[1:]
        V = V_P + c_v * P
        return np.concatenate([qnet / c_v, (jnet - C * qnet) / V])

    y0 = np.concatenate([np.full(n, case['P0']), np.full(n, case['C0'])])
    sol = solve_ivp(rhs, (t[0], t[-1]), y0, t_eval=t, method='BDF', rtol=1e-10, atol=1e-14)
    return sol.y[:n], sol.y[n:]


CASES = {
    'square_channel_advection': dict(
        description='Advection and diffusion of a solute along a straight 20 um square channel, 200 um long (20 cells), '
                    'under a steady 1 Pa pressure drop (cell Peclet number 14 over the channel).',
        P_in=1.0, P_amp=0.0, f_P=1.0, P_out=0.0, C_in=1.0, C_out=0.0, P0=0.0, C0=0.0, sim_time=40.0, dt=0.05),
    'square_channel_diffusion': dict(
        description='Pure diffusion along a straight 20 um square channel, 200 um long (20 cells): no pressure drop.',
        P_in=0.0, P_amp=0.0, f_P=1.0, P_out=0.0, C_in=1.0, C_out=0.0, P0=0.0, C0=0.0, sim_time=200.0, dt=0.5),
    'square_channel_pulsatile': dict(
        description='A straight 20 um square channel, 200 um long (20 cells), with inlet pressure 1 + 2 sin(2 pi t) Pa: '
                    'the flow reverses every cycle.',
        P_in=1.0, P_amp=2.0, f_P=1.0, P_out=0.0, C_in=1.0, C_out=0.0, P0=0.0, C0=0.0, sim_time=10.0, dt=0.01),
}


def cells(n=N):
    return [f'c{k}' for k in range(n)]


def write(model, case):
    d = os.path.join(OUT, model)
    os.makedirs(d, exist_ok=True)
    dx, A, faces = geometry()
    names = cells()
    rows = [['inlet', 'nn', 'poiseuille_transport_inlet', '', 'f0']]
    params = [['sigma_solute', 'm2_per_s', SIGMA, 'example: small solute in water'],
              ['mu_visc', 'Js_per_m3', MU, 'example: water at 20 C'],
              ['P_mean_inlet', 'J_per_m3', case['P_in'], 'example: inlet mean pressure'],
              ['P_amp_inlet', 'J_per_m3', case['P_amp'], 'example: inlet pressure amplitude'],
              ['f_P_inlet', 'Hz', case['f_P'], 'example: inlet pressure frequency'],
              ['C_res_inlet', 'millimolar', case['C_in'], 'example: inlet concentration'],
              ['P_mean_outlet', 'J_per_m3', case['P_out'], 'example: outlet pressure'],
              ['P_amp_outlet', 'J_per_m3', 0.0, 'example: constant outlet pressure'],
              ['f_P_outlet', 'Hz', 1.0, 'example (unused with P_amp = 0)'],
              ['C_res_outlet', 'millimolar', case['C_out'], 'example: outlet concentration']]
    for k, c in enumerate(names):
        rows.append([c, 'nn', 'poiseuille_transport_volume', f'f{k}', f'f{k + 1}'])
        params += [[f'C_init_{c}', 'millimolar', case['C0'], 'example: initial concentration'],
                   [f'P_init_{c}', 'J_per_m3', case['P0'], 'example: initial pressure'],
                   [f'V_P_{c}', 'm3', A * dx, 'example mesh: cell volume a^2 dx'],
                   [f'C_V_{c}', 'm6_per_J', C_V, 'example: compliance (regularisation)']]
    for k, (du, dd) in enumerate(faces):
        up = 'inlet' if k == 0 else names[k - 1]
        down = 'outlet' if k == N else names[k]
        rows.append([f'f{k}', 'nn', 'poiseuille_transport_face', up, down])
        params += [[f'A_f_f{k}', 'm2', A, 'example mesh: square cross-section a^2'],
                   [f'd_up_f{k}', 'metre', du, f'example mesh: {up} to the face'],
                   [f'd_down_f{k}', 'metre', dd, f'example mesh: face to {down}'],
                   [f'open_flag_f{k}', 'dimensionless', 1, 'open face'],
                   [f'K_shape_f{k}', 'dimensionless', K_SQUARE, 'square duct Poiseuille factor (series solution)'],
                   [f'k_TA_f{k}', 'dimensionless', 0, 'no Taylor-Aris dispersion (exact 1D solution)']]
    rows.append(['outlet', 'nn', 'poiseuille_transport_outlet', f'f{N}', ''])
    with open(os.path.join(d, f'{model}_vessel_array.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['name', 'BC_type', 'vessel_type', 'inp_vessels', 'out_vessels'])
        w.writerows(rows)
    with open(os.path.join(d, f'{model}_parameters.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['variable_name', 'units', 'value', 'data_reference'])
        w.writerows(params)

    spec_path = os.path.join(d, f'{model}_system.yaml')
    if os.path.isfile(spec_path):
        return
    cs = ', '.join(f'{c}__C' for c in names)
    vs = ', '.join(f'{c}__V' for c in names)
    amount = f'np.sum(np.stack([{vs}])*np.stack([{cs}]), axis=0)'
    volume = f'np.sum(np.stack([{vs}]), axis=0)'
    invariants = [
        {'description': 'solute conservation: the change in sum V C equals the net amount the reservoirs gave up '
                        '(relative 1e-6 of the largest term)',
         'expr': f'np.max(np.abs(({amount} - {amount}[0]) + inlet__N_received + outlet__N_received)) <= '
                 f'1e-6*max(np.max(np.abs(inlet__N_received)), np.max(np.abs({amount})), 1e-30)'},
        {'description': 'fluid conservation: the change in sum V equals the net volume the reservoirs gave up (relative 1e-6)',
         'expr': f'np.max(np.abs(({volume} - {volume}[0]) + inlet__V_received + outlet__V_received)) <= '
                 f'1e-6*max(np.max(np.abs(inlet__V_received)), np.max(np.abs({volume} - {volume}[0])), 1e-30)'},
        {'description': 'maximum principle: every concentration stays within [min, max] of the reservoir and initial values',
         'expr': f'np.min(np.stack([{cs}])) >= {min(case["C_in"], case["C_out"], case["C0"])} - 1e-9 and '
                 f'np.max(np.stack([{cs}])) <= {max(case["C_in"], case["C_out"], case["C0"])} + 1e-9'},
    ]
    if case['P_amp'] == 0:
        Q, C, J = exact_steady(case['P_in'], case['P_out'], case['C_in'], case['C_out'])
        qs = ', '.join(f'f{k}__Q' for k in range(N + 1))
        invariants += [
            {'description': f'steady flow in every face = exact square-duct Poiseuille flow a^4 dP/(12 mu L)(1 - 192/pi^5 '
                            f'sum tanh(n pi/2)/n^5) = {Q:.6g} m3/s (relative 1e-6 at t_end)',
             'expr': f'np.max(np.abs(np.stack([{qs}])[:, -1] - {Q!r})) <= 1e-6*abs({Q!r}) + 1e-24'},
            {'description': 'steady concentrations at the cell centres = exact 1D advection-diffusion solution '
                            '(1e-6 mM at t_end)',
             'expr': f'np.max(np.abs(np.stack([{cs}])[:, -1] - np.array({[float(c) for c in C]!r}))) <= 1e-6'},
            {'description': f'steady solute flux into the outlet = exact {J:.6g} mol/s (relative 1e-6 at t_end)',
             'expr': f'abs(f{N}__J[-1] - {J!r}) <= 1e-6*abs({J!r})'},
        ]
    spec = {
        'model': model, 'category': 'poiseuille', 'reviewed': False,
        'description': case['description'],
        'sim_time': case['sim_time'], 'pre_time': 0.0, 'dt': case['dt'],
        'solver_info': {'rtol': 1e-10, 'atol': 1e-16, 'MaximumStep': 0.01},
        'equivalence': {'status': 'not_applicable',
                        'reason': 'an example built from this library (tools/build_poiseuille_examples.py); there is no '
                                  'circulatory_autogen original'},
        'invariants': invariants,
        'notes': ['Built by tools/build_poiseuille_examples.py: re-run it after changing the channel.'],
    }
    with open(spec_path, 'w') as f:
        yaml.safe_dump(spec, f, sort_keys=False, width=120)


def main():
    for model, case in CASES.items():
        write(model, case)


if __name__ == '__main__':
    main()
