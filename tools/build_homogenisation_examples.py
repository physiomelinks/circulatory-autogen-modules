"""
Builds the homogenisation system-model pairs (system_models/homogenisation/): for each pair, a fine model
and a coarse model of the same physical object, so that the coarse level's derived parameters can be
checked against the finer level (tests/test_systems.py::system_homogenisation_test).

Every volume is diffusion_volume's tissue_diffusion_volume (Michaelis-Menten O2 consumption) and every face
tissue_diffusion_face. The surroundings are a well-mixed bath: one volume of 1 m^3 with no consumption,
so its concentration is fixed to within ~1e-13 mM over a run. Only in-plane faces exchange; the top and
bottom of the 2D sheet (depth DEPTH) are insulated.

  level 0 -> 1  cell_disc_fv      a disc cell (radius R) meshed into N_RINGS concentric finite volumes
                cell_disc_lumped  the same cell as one well-mixed volume. Its face to the bath uses the
                                  derived effective distance d_eff = R/4: for uniform consumption q the
                                  steady profile is C(r) = C_s - q (R^2 - r^2)/(4 sigma), whose area mean
                                  C_s - q R^2/(8 sigma) a single volume reproduces with J = sigma A (C_s - C)/d_eff,
                                  A = 2 pi R DEPTH, J = q pi R^2 DEPTH  =>  d_eff = R/4.
  level 1 -> 2  tissue_square_cells    an N_SIDE x N_SIDE patch of square cells (one volume per cell, faces at
                                       half-widths: the cell level), exposed to the bath on all four sides
                tissue_square_lumped   the same patch as one tissue volume. The continuum gives d_eff = 4 m L
                                       (m = 0.0351442, the mean of u for -lap(u) = 1 on the unit square with u = 0
                                       on its boundary), but the 8 x 8 cell level is discrete and differs from the
                                       continuum by ~6%, so d_eff is calibrated to the cell level instead.

Upscaling rule: use the derived formula where it reproduces the finer level's mean deficit (uniform
consumption, linear problem) to 1%; otherwise calibrate to the finer level. Both values are recorded in the
coarse spec's homogenisation.upscaling.

Each model's spec gets invariants for mass conservation (consumption and bath exchange accounted, i.e.
zero mass loss) and the second law (every face's dissipation J ln(C_up/C_down) >= 0). The first law
(enthalpy) is covered by mass conservation here: ideal diffusion of a dilute solute releases no heat and
no thermal domain is modelled, so enthalpy is conserved exactly when mass is. The coarse spec of each pair
gets a `homogenisation:` section. Rewrites module arrays, parameters and specs.

    python tools/build_homogenisation_examples.py
"""
import csv
import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing import module_array  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'system_models', 'homogenisation')
SIGMA = 2.41e-9          # m2/s, tissue O2 diffusivity (Secomb 2020, Table 1), as build_diffusion_examples
DEPTH = 2.0e-5           # m, depth of the 2D sheet
C_BATH = 0.2             # mM, air-saturated water at ~37 C
BATH_V = 1.0             # m3: the bath's concentration is effectively fixed
M_DEFAULT = -0.0558      # mM/s, tissue_diffusion_volume default instance (CelayaAlcala2021)
C50 = 0.0133             # mM, default instance
R_CELL = 1.0e-5          # m, disc cell radius (20 um cell)
N_RINGS = 32
L_CELL = 2.0e-5          # m, square cell side
N_SIDE = 8               # cells per side of the tissue patch
SQUARE_MEAN = 0.0351442  # mean of u, -lap u = 1 on the unit square, u = 0 on the boundary

# Parameters swept in both models of a pair: each name is a glob over the parameters file.
SWEEP_PARAMETERS = {
    'bath_O2': {'pattern': 'bath/C_P_init', 'units': 'millimolar'},
    'M': {'pattern': 'c*/M', 'units': 'millimolar_per_s'},
}
SWEEP_POINTS = [
    {'bath_O2': 0.2, 'M': M_DEFAULT},
    {'bath_O2': 0.3, 'M': M_DEFAULT},
    {'bath_O2': 0.05, 'M': M_DEFAULT},
    {'bath_O2': 0.02, 'M': M_DEFAULT},
    {'bath_O2': 0.2, 'M': -0.02},
    {'bath_O2': 0.2, 'M': -0.15},
]


# ---------------------------------------------------------------- meshes
def disc_fv():
    '''Concentric rings: ring 0 is the central disc (node at the centre), ring i>0 has its node mid-ring.'''
    dr = R_CELL / N_RINGS
    edges = [i * dr for i in range(N_RINGS + 1)]
    nodes = [0.0] + [0.5 * (edges[i] + edges[i + 1]) for i in range(1, N_RINGS)]
    cells = {f'c{i}': {'V': np.pi * (edges[i + 1] ** 2 - edges[i] ** 2) * DEPTH} for i in range(N_RINGS)}
    faces = [(f'c{i}', f'c{i + 1}', 2 * np.pi * edges[i + 1] * DEPTH, edges[i + 1] - nodes[i], nodes[i + 1] - edges[i + 1])
             for i in range(N_RINGS - 1)]
    exposed = [(f'c{N_RINGS - 1}', 2 * np.pi * R_CELL * DEPTH, R_CELL - nodes[-1])]
    return cells, faces, exposed


def disc_lumped():
    cells = {'c0': {'V': np.pi * R_CELL ** 2 * DEPTH}}
    return cells, [], [('c0', 2 * np.pi * R_CELL * DEPTH, R_CELL / 4)]


def square_cells():
    cells = {f'c{i}_{j}': {'V': L_CELL * L_CELL * DEPTH} for i in range(N_SIDE) for j in range(N_SIDE)}
    faces, exposed = [], []
    for i in range(N_SIDE):
        for j in range(N_SIDE):
            if i + 1 < N_SIDE:
                faces.append((f'c{i}_{j}', f'c{i + 1}_{j}', L_CELL * DEPTH, L_CELL / 2, L_CELL / 2))
            if j + 1 < N_SIDE:
                faces.append((f'c{i}_{j}', f'c{i}_{j + 1}', L_CELL * DEPTH, L_CELL / 2, L_CELL / 2))
            sides = (i == 0) + (i == N_SIDE - 1) + (j == 0) + (j == N_SIDE - 1)
            if sides:
                exposed.append((f'c{i}_{j}', sides * L_CELL * DEPTH, L_CELL / 2))
    return cells, faces, exposed


def square_lumped(d_eff=None):
    side = N_SIDE * L_CELL
    return {'c0': {'V': side * side * DEPTH}}, [], [('c0', 4 * side * DEPTH, d_eff or 4 * SQUARE_MEAN * side)]


# ---------------------------------------------------------------- numerical check of the derivations
def steady_linear(cells, faces, exposed, q=1.0, sigma=1.0, cs=0.0):
    '''Steady state with uniform consumption q and the bath at cs = 0: returns the volume-weighted mean
    deficit below the bath (positive). Solves sum_faces g (c_i - c_j) + g_e c_i = q V_i for the deficits c.'''
    names = list(cells)
    idx = {n: k for k, n in enumerate(names)}
    n = len(names)
    A = np.zeros((n, n))
    b = np.array([-q * cells[k]['V'] for k in names])
    for a, c, area, da, dc in faces:
        g = sigma * area / (da + dc)
        A[idx[a], idx[a]] -= g
        A[idx[c], idx[c]] -= g
        A[idx[a], idx[c]] += g
        A[idx[c], idx[a]] += g
    for a, area, d in exposed:
        g = sigma * area / d
        A[idx[a], idx[a]] -= g
        b[idx[a]] -= g * cs
    c = np.linalg.solve(A, b)
    v = np.array([cells[k]['V'] for k in names])
    return float(np.dot(v, c) / v.sum())


def calibrated_distance(fine, coarse):
    '''Effective distance that makes the coarse volume reproduce the fine level's mean deficit exactly
    (uniform consumption, linear problem): deficit = q V d / (sigma A).'''
    cells, _, exposed = coarse
    (V,) = [c['V'] for c in cells.values()]
    (area,) = [e[1] for e in exposed]
    return steady_linear(*fine) * area / V


def upscale(name, fine, derived):
    '''Derived where the formula reproduces the finer level to 1%, calibrated to the finer level where not.'''
    calibrated = calibrated_distance(fine, derived)
    d_derived = derived[2][0][2]
    rel = abs(d_derived - calibrated) / calibrated
    use = 'derived' if rel < 0.01 else 'calibrated'
    print(f'{name}: d_eff derived {d_derived:.6g} m, calibrated to the fine level {calibrated:.6g} m '
          f'(differ by {rel:.2%}) -> using {use}')
    return {'derived_m': d_derived, 'calibrated_m': calibrated, 'relative_difference': rel, 'used': use,
            'value_m': d_derived if use == 'derived' else calibrated}


# ---------------------------------------------------------------- writing
def write(model, description, cells, faces, exposed, sim_time, homogenisation=None):
    d = os.path.join(OUT, model)
    os.makedirs(d, exist_ok=True)
    rows = []
    params = [['sigma_diff', 'm2_per_s', SIGMA, 'secomb2020mass; Table 1 tissue O2 diffusivity 2410 um^2/s']]
    all_faces = [(a, c, area, da, dc, f'f{k}') for k, (a, c, area, da, dc) in enumerate(faces)]
    all_faces += [('bath', a, area, 0.0, dist, f'e{k}') for k, (a, area, dist) in enumerate(exposed)]

    def add_volume(name, V, c0, m, ref):
        inp = [f[5] for f in all_faces if f[1] == name]
        out = [f[5] for f in all_faces if f[0] == name]
        rows.append([name, 'Fang2008_v01', 'tissue_diffusion_volume', ' '.join(inp), ' '.join(out), 'default'])
        params.extend([[f'{name}/C_P_init', 'millimolar', c0, ref],
                       [f'{name}/V_P', 'm3', V, ref],
                       [f'{name}/flux_c', 'mol_per_s', 0, 'no capillary supply'],
                       [f'{name}/M', 'millimolar_per_s', m, ref],
                       [f'{name}/C50', 'millimolar', C50, 'CelayaAlcala2021 via the default instance']])

    add_volume('bath', BATH_V, C_BATH, 0, 'homogenisation example: well-mixed bath, no consumption')
    for name, c in cells.items():
        # Cells start empty, so each run also contains the transient (an initial-condition step).
        add_volume(name, c['V'], 0, M_DEFAULT, 'homogenisation example: geometry and default consumption')
    for a, c, area, da, dc, fname in all_faces:
        kind = 'derived effective distance (see tools/build_homogenisation_examples.py)' \
            if a == 'bath' and model.endswith('lumped') else 'mesh geometry'
        rows.append([fname, 'Fang2008_v01', 'tissue_diffusion_face', a, c, 'default'])
        params += [[f'{fname}/A_f', 'm2', area, 'mesh geometry: face area'],
                   [f'{fname}/d_up', 'metre', da, kind],
                   [f'{fname}/d_down', 'metre', dc, kind],
                   [f'{fname}/open_flag', 'dimensionless', 1, 'open face']]
    module_array.write_records(os.path.join(d, f'{model}_module_array.json'), module_array.from_rows(rows))
    with open(os.path.join(d, f'{model}_parameters.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['variable_name', 'units', 'value', 'data_reference'])
        w.writerows(params)

    vols = ' + '.join(f'{c["V"]!r}*{n}__C_P' for n, c in cells.items())
    cons = ' + '.join(f'{c["V"]!r}*{n}__M_reduced' for n, c in cells.items())
    supply = ' + '.join(f'{f[5]}__J' for f in all_faces if f[0] == 'bath')
    integral = lambda e: f'np.concatenate([[0.0], np.cumsum(0.5*np.diff(t)*(({e})[1:] + ({e})[:-1]))])'
    invariants = [
        {'description': 'mass conservation (zero loss): the amount in the cells changes only by bath supply and '
                        'consumption, sum V C(t) - sum V C(0) = int (supply + sum V M_reduced) dt (relative 1e-6)',
         'expr': f'np.max(np.abs(({vols}) - ({vols})[0] - {integral(f"{supply} + {cons}")})) '
                 f'<= 1e-6*max(np.max(np.abs({vols})), np.max(np.abs({integral(supply)})))'},
        {'description': 'second law: every face dissipates, J ln(C_up/C_down) >= 0 (using the connected volumes, so '
                        'it also checks the couplings)',
         'expr': ' and '.join(f'np.all({f[5]}__J*(np.log(np.maximum({f[0]}__C_P, 1e-300)) - '
                              f'np.log(np.maximum({f[1]}__C_P, 1e-300))) >= -1e-30)' for f in all_faces)},
    ]
    spec = {
        'model': model, 'category': 'homogenisation', 'reviewed': False,
        'description': description,
        'sim_time': sim_time, 'pre_time': 0.0, 'dt': 0.01,
        'solver_info': {'rtol': 1e-10, 'atol': 1e-16},
        'equivalence': {'status': 'not_applicable',
                        'reason': 'a homogenisation example built from this library '
                                  '(tools/build_homogenisation_examples.py); there is no circulatory_autogen original'},
        'invariants': invariants,
        'notes': ['Built by tools/build_homogenisation_examples.py: re-run it after changing the geometry.',
                  'First law: ideal diffusion releases no heat and no thermal domain is modelled, so enthalpy is '
                  'conserved exactly when mass is (the mass invariant).'],
    }
    if homogenisation:
        spec['homogenisation'] = homogenisation
    with open(os.path.join(d, f'{model}_system.yaml'), 'w') as f:
        yaml.safe_dump(spec, f, sort_keys=False, width=120)


def homogenisation_spec(fine_model, level, fine_cells, derivation, upscaling):
    total_v = sum(c['V'] for c in fine_cells.values())
    return {
        'fine': f'homogenisation/{fine_model}',
        'level': level,
        'derivation': derivation,
        'upscaling': {k: float(v) if isinstance(v, (int, float, np.floating)) else v for k, v in upscaling.items()},
        'quantities': [
            {'name': 'mean O2 concentration (mM)',
             'fine': {'weighted_mean': [[f'{n}/C_P', c['V'] / total_v] for n, c in fine_cells.items()]},
             'coarse': {'weighted_mean': [['c0/C_P', 1.0]]}},
            {'name': 'total uptake from the bath (mol/s)',
             'fine': {'sum_matching': 'e*/J'},
             'coarse': {'sum_matching': 'e*/J'}},
        ],
        'steady_state_tol': 0.01,
        'transient_nrmse_tol': 0.05,
        'sweep': {'parameters': SWEEP_PARAMETERS, 'points': SWEEP_POINTS},
    }


def main():
    disc_up = upscale('disc (level 0 -> 1)', disc_fv(), disc_lumped())
    square_up = upscale('square (level 1 -> 2)', square_cells(), square_lumped())
    disc_cells = disc_fv()[0]
    write('cell_disc_fv', f'Level 0: a {2e6 * R_CELL:g} um disc cell meshed into {N_RINGS} concentric finite volumes '
          f'({1e6 * DEPTH:g} um deep), consuming O2 in a bath.', *disc_fv(), sim_time=1.0)
    write('cell_disc_lumped', 'Level 1: the same disc cell as one well-mixed volume; face to the bath at the derived '
          'effective distance R/4.', *disc_lumped(), sim_time=1.0,
          homogenisation=homogenisation_spec('cell_disc_fv', '0->1', disc_cells,
                                             'd_eff = R/4 from the steady profile of uniform consumption in a disc',
                                             disc_up))
    square = square_cells()[0]
    side_um = N_SIDE * L_CELL * 1e6
    write('tissue_square_cells', f'Level 1: a {side_um:g} um square patch of {N_SIDE}x{N_SIDE} cells (one volume per '
          'cell) consuming O2, bath on all four sides.', *square_cells(), sim_time=60.0)
    write('tissue_square_lumped', 'Level 2: the same patch as one homogenised tissue volume; face to the bath at an '
          f'effective distance {square_up["used"]} from the cell level (continuum formula 4 m L, m = 0.0351442).',
          *square_lumped(square_up['value_m']), sim_time=60.0,
          homogenisation=homogenisation_spec('tissue_square_cells', '1->2', square,
                                             'continuum: d_eff = 4 m L, m = mean of the unit-square Poisson solution; '
                                             'the discrete cell level differs, so d_eff is calibrated to it', square_up))


if __name__ == '__main__':
    main()
