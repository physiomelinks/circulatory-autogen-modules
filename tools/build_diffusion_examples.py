"""
Builds the diffusion example system models (system_models/diffusion/) from simple meshes of
diffusion_volume's tissue_diffusion_volume cells and tissue_diffusion_face faces:

  diffusion_2d_triangles  a 2D mesh of triangles (a 2 x 2 grid of squares, each cut along a diagonal),
                          20 um deep: non-square cells, faces of different areas and orientations
  diffusion_3d_boxes      a 2 x 2 x 2 block of box cells of different sizes, faces in x, y and z

Each cell is a volume (V_P); each shared face joins two cells with its area A_f and the distances
d_up, d_down from the two cell centres (centroids) to the face. The domains are closed, with
different initial concentrations and no sources or consumption, so the total amount sum V C is
conserved and every cell tends to the volume-weighted mean. Rewrites the vessel array and
parameters, and writes the spec only if there isn't one.

    python tools/build_diffusion_examples.py
"""
import csv
import itertools
import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing import vessel_array  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'system_models', 'diffusion')
SIGMA = 2.41e-9         # m2/s, tissue O2 diffusivity (Secomb 2020, Table 1)
DEPTH = 2.0e-5          # m, depth of the 2D mesh


def triangles_2d():
    h = 2.0e-5
    nodes = {(i, j): np.array([i * h, j * h]) for i in range(3) for j in range(3)}
    tris = []
    for i, j in itertools.product(range(2), range(2)):
        a, b, c, d = (i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)
        tris += [(a, b, c), (a, c, d)] if (i + j) % 2 == 0 else [(a, b, d), (b, c, d)]
    cells = {}
    for k, tri in enumerate(tris):
        p = [nodes[n] for n in tri]
        area = 0.5 * abs((p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]))
        cells[f't{k}'] = {'nodes': tri, 'centre': sum(p) / 3, 'V': area * DEPTH}
    faces = []
    names = list(cells)
    for a, b in itertools.combinations(names, 2):
        shared = set(cells[a]['nodes']) & set(cells[b]['nodes'])
        if len(shared) != 2:
            continue
        p, q = (nodes[n] for n in sorted(shared))
        edge = q - p
        normal = np.array([edge[1], -edge[0]]) / np.linalg.norm(edge)
        dist = lambda c: abs(float(np.dot(c - p, normal)))   # centroid to the edge's line
        faces.append((a, b, float(np.linalg.norm(edge)) * DEPTH, dist(cells[a]['centre']), dist(cells[b]['centre'])))
    return cells, faces


def boxes_3d():
    widths = {'x': [1.5e-5, 2.5e-5], 'y': [2.0e-5, 3.0e-5], 'z': [1.0e-5, 2.0e-5]}
    cells = {}
    for i, j, k in itertools.product(range(2), range(2), range(2)):
        cells[f'b{i}{j}{k}'] = {'ijk': (i, j, k), 'V': widths['x'][i] * widths['y'][j] * widths['z'][k]}
    faces = []
    for (a, ca), (b, cb) in itertools.combinations(cells.items(), 2):
        d = [abs(x - y) for x, y in zip(ca['ijk'], cb['ijk'])]
        if sum(d) != 1:
            continue
        axis = d.index(1)
        others = [ax for n, ax in enumerate('xyz') if n != axis]
        area = np.prod([widths[ax][ca['ijk']['xyz'.index(ax)]] for ax in others])
        w = widths['xyz'[axis]]
        faces.append((a, b, float(area), w[ca['ijk'][axis]] / 2, w[cb['ijk'][axis]] / 2))
    return cells, faces


def initial(cells):
    rng = np.random.default_rng(7)
    return {n: float(v) for n, v in zip(cells, rng.uniform(0.0, 0.12, len(cells)))}


def numpy_solve(cells, faces, C0, t):
    '''Direct finite-volume solve of the same mesh, for checking the generated model.'''
    from scipy.integrate import solve_ivp
    names = list(cells)
    idx = {n: k for k, n in enumerate(names)}
    V = np.array([cells[n]['V'] for n in names])

    def rhs(_, c):
        d = np.zeros_like(c)
        for a, b, A, da, db in faces:
            J = SIGMA * A / (da + db) * (c[idx[a]] - c[idx[b]])
            d[idx[b]] += J
            d[idx[a]] -= J
        return d / V
    sol = solve_ivp(rhs, (t[0], t[-1]), [C0[n] for n in names], t_eval=t, method='LSODA', rtol=1e-12, atol=1e-18)
    return {n: sol.y[idx[n]] for n in names}


def write(model, description, cells, faces):
    d = os.path.join(OUT, model)
    os.makedirs(d, exist_ok=True)
    C0 = initial(cells)
    rows, params = [], [['sigma_diff', 'm2_per_s', SIGMA, 'Secomb2020; Table 1 tissue O2 diffusivity 2410 um^2/s']]
    for n, c in cells.items():
        inp = [f'f{k}' for k, f in enumerate(faces) if f[1] == n]
        out = [f'f{k}' for k, f in enumerate(faces) if f[0] == n]
        rows.append([n, 'nn', 'tissue_diffusion_volume', ' '.join(inp), ' '.join(out), 'default'])
        params += [[f'C_P_init_{n}', 'millimolar', round(C0[n], 6), 'example: random initial concentration'],
                   [f'V_P_{n}', 'm3', c['V'], 'example mesh: cell volume'],
                   [f'flux_c_{n}', 'mol_per_s', 0, 'example: closed domain, no capillary supply'],
                   [f'M_{n}', 'millimolar_per_s', 0, 'example: no consumption'],
                   [f'C50_{n}', 'millimolar', 0.0133, 'example (unused with M = 0)']]
    for k, (a, b, A, da, db) in enumerate(faces):
        rows.append([f'f{k}', 'nn', 'tissue_diffusion_face', a, b, 'default'])
        params += [[f'A_f_f{k}', 'm2', A, 'example mesh: face area'],
                   [f'd_up_f{k}', 'metre', da, f'example mesh: centre of {a} to the face'],
                   [f'd_down_f{k}', 'metre', db, f'example mesh: centre of {b} to the face'],
                   [f'open_flag_f{k}', 'dimensionless', 1, 'open face']]
    vessel_array.write_records(os.path.join(d, f'{model}_vessel_array.json'), vessel_array.from_rows(rows))
    with open(os.path.join(d, f'{model}_parameters.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['variable_name', 'units', 'value', 'data_reference'])
        w.writerows(params)

    spec_path = os.path.join(d, f'{model}_system.yaml')
    if os.path.isfile(spec_path):
        return
    total = ' + '.join(f'{c["V"]!r}*{n}__C_P' for n, c in cells.items())
    vtot = sum(c['V'] for c in cells.values())
    mean = sum(c['V'] * round(C0[n], 6) for n, c in cells.items()) / vtot
    cps = ', '.join(f'{n}__C_P' for n in cells)
    spec = {
        'model': model, 'category': 'diffusion', 'reviewed': False,
        'description': description,
        'sim_time': 2.0, 'pre_time': 0.0, 'dt': 0.01,
        'solver_info': {'rtol': 1e-10, 'atol': 1e-16},
        'equivalence': {'status': 'not_applicable',
                        'reason': 'an example built from this library (tools/build_diffusion_examples.py); there is no '
                                  'circulatory_autogen original'},
        'invariants': [
            {'description': 'the total amount sum V_P C_P is conserved in the closed domain (relative 1e-9)',
             'expr': f'np.ptp({total}) <= 1e-9*np.max(np.abs({total}))'},
            {'description': 'maximum principle: no cell goes above the initial maximum or below the initial minimum',
             'expr': f'np.max(np.stack([{cps}])) <= np.max(np.stack([{cps}])[:, 0]) + 1e-12 and '
                     f'np.min(np.stack([{cps}])) >= np.min(np.stack([{cps}])[:, 0]) - 1e-12'},
            {'description': f'every cell tends to the volume-weighted mean {mean:.6g} mM (within 1e-6 mM at t_end)',
             'expr': f'np.max(np.abs(np.stack([{cps}])[:, -1] - {mean!r})) <= 1e-6'},
        ],
        'notes': ['Built by tools/build_diffusion_examples.py: re-run it after changing the mesh.'],
    }
    with open(spec_path, 'w') as f:
        yaml.safe_dump(spec, f, sort_keys=False, width=120)


def main():
    write('diffusion_2d_triangles', 'Finite-volume diffusion on a 2D triangle mesh (8 triangles, 20 um deep): non-square '
          'cells joined by faces of different areas and orientations.', *triangles_2d())
    write('diffusion_3d_boxes', 'Finite-volume diffusion in 3D: a 2 x 2 x 2 block of box cells of different sizes, joined '
          'by faces in x, y and z.', *boxes_3d())


if __name__ == '__main__':
    main()
