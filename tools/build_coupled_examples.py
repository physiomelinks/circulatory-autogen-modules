"""
Builds the coupled example system models (system_models/coupled/): 0D modules exchanging a
solute with a 3D box of tissue, modelled either by a finite-volume grid of CellML cells
(tissue_diffusion_volume + tissue_diffusion_face) or by the FEniCS model tissue_diffusion_FEniCS
(an external model: module_format external_api, coupled to the generated C++ by
libcuflynx.coupling). The two describe the same box on the same grid, so they can be compared.

  microvasc_O2_FV        two capillaries (capillary pp + GE_capillary) exchanging O2 with a
  microvasc_O2_FEniCS    cube of tissue (N x N x N cells of 20 um, N = 3: 60 um)
  SN_NE_original         periodic stimulus -> soma -> axon -> varicosity (NE an internal state)
  SN_NE_single_volume    the same with varicosity NEexchange_v01 and one well-mixed extracellular
                         volume equal to the varicosity volume: exactly SN_NE_original
  SN_NE_FV               ... NEexchange_v01 releasing NE into a cube of extracellular space
  SN_NE_FEniCS           (N x N x N cells of 1 um, the varicosity volume)

The coupled 0D module k exchanges with the grid cell the FEniCS model gives neighbour k of N:
x index floor((k + 1/2) n / N), the middle cell in y and z. Rewrites the vessel arrays and
parameters, and writes each spec only if there isn't one.

    python tools/build_coupled_examples.py
"""
import csv
import json
import os

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'system_models', 'coupled')
N = 3   # grid cells per direction (--n): the CellML grid has N^3 cells and 3 N^2 (N - 1) faces
H_O2, H_NE = 2e-5, 1e-6   # cell sizes: the tissue_diffusion_volume default cell, and the varicosity volume

O2 = {'L': None, 'sigma': 2.41e-9, 'C_init': 0.0597, 'M': -0.0558, 'C50': 0.0133, 'instance': 'default'}
NE = {'L': None, 'sigma': 3.0e-10, 'C_init': 0.0, 'M': 0.0, 'C50': 1.0, 'instance': 'NE_extracellular'}


def rec(name, module_type, subtype, inp=(), out=(), instance='default'):
    return {'name': name, 'module_type': module_type, 'module_subtype': subtype, 'instance': instance,
            'inp_instances': list(inp), 'out_instances': list(out)}


def coupled_cells(n_regions):
    return [(int((k + 0.5) * N / n_regions), N // 2, N // 2) for k in range(n_regions)]


def fv_grid(sources, props):
    """Records and parameters of an N^3 grid of tissue_diffusion_volume cells joined by
    tissue_diffusion_face faces; sources[k] (a 0D module) exchanges with coupled_cells()[k]."""
    h = props['L'] / N
    cell = lambda i, j, k: f'c{i}{j}{k}'
    targets = dict(zip(coupled_cells(len(sources)), sources))
    inputs = {cell(i, j, k): [] for i in range(N) for j in range(N) for k in range(N)}
    outputs = {c: [] for c in inputs}
    faces = []
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for axis, (di, dj, dk) in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1))):
                    if i + di < N and j + dj < N and k + dk < N:
                        up, down = cell(i, j, k), cell(i + di, j + dj, k + dk)
                        f = f'f{"xyz"[axis]}{i}{j}{k}'
                        faces.append(rec(f, 'tissue_diffusion_face', 'Fang2008_v01', [up], [down]))
                        outputs[up].append(f)
                        inputs[down].append(f)
    records, params = [], []
    for (i, j, k), c in ((key, cell(*key)) for key in sorted({(i, j, k) for i in range(N) for j in range(N)
                                                                 for k in range(N)})):
        src = targets.get((i, j, k))
        records.append(rec(c, 'tissue_diffusion_volume', 'Fang2008_v01', ([src] if src else []) + inputs[c], outputs[c]))
        params += [(f'V_P_{c}', 'm3', h ** 3, f'grid cell, {h * 1e6:g} um cube'),
                   (f'C_P_init_{c}', 'millimolar', props['C_init'], 'initial concentration (as the FEniCS instance)'),
                   (f'M_{c}', 'millimolar_per_s', props['M'], 'as the FEniCS instance'),
                   (f'C50_{c}', 'millimolar', props['C50'], 'as the FEniCS instance')]
        if src is None:
            params.append((f'flux_c_{c}', 'mol_per_s', 0.0, 'no source: the cell exchanges only with its neighbours'))
    for f in faces:
        params += [(f"A_f_{f['name']}", 'm2', h * h, 'face of the uniform grid'),
                   (f"d_up_{f['name']}", 'metre', h / 2, 'centre-to-face distance'),
                   (f"d_down_{f['name']}", 'metre', h / 2, 'centre-to-face distance'),
                   (f"open_flag_{f['name']}", 'dimensionless', 1, 'open face')]
    params.append(('sigma_diff', 'm2_per_s', props['sigma'], 'as the FEniCS instance'))
    return records + faces, params, [cell(*c) for c in coupled_cells(len(sources))]


def fenics_row(sources, props):
    """The FEniCS module, on the same box and grid as fv_grid (overriding its instance's size)."""
    params = [(f'{k}_tissue', 'metre', props['L'], 'the box of the CellML grid variant') for k in ('Lx', 'Ly', 'Lz')]
    params += [(f'{k}_tissue', 'dimensionless', N, 'the grid of the CellML grid variant')
               for k in ('grid_nx', 'grid_ny', 'grid_nz')]
    return [rec('tissue', 'tissue_diffusion_FEniCS', 'box_v01', sources, [], props['instance'])], params


# --- the 0D sides --------------------------------------------------------------------------------

# GE_capillary constants of CA_user/capillary_network (Spencer O2/CO2 dissociation curves, solubilities)
GE_CONSTANTS = [('exp_threshold', 'millimolar', 0.5), ('S_O2', 'mol_per_m3_per_Pa', 9.75975975975976e-06),
                ('S_CO2', 'mol_per_m3_per_Pa', 0.00022522522522522523), ('k_1', 'mmHg', 14.99), ('k_2', 'mmHg', 194.4),
                ('a_1', 'dimensionless', 0.3836), ('a_2', 'dimensionless', 1.819), ('alpha_1', 'per_mmHg', 0.03198),
                ('alpha_2', 'per_mmHg', 0.05591), ('beta_1', 'per_mmHg', 0.008275), ('beta_2', 'per_mmHg', 0.03255),
                ('C_O2_sat', 'mol_per_m3', 9.0), ('C_CO2_sat', 'mol_per_m3', 86.11)]


def capillaries(tissue_names):
    """Two capillaries (capillary pp) with their O2 exchange (GE_capillary), each sending O2 into
    tissue_names[k]. Values from CA_user/capillary_network (capillary_6, capillary_GE_6)."""
    records, params = [], []
    for k, target in enumerate(tissue_names):
        cap, ge = f'cap_{k}', f'capillary_GE_{k}'
        records += [rec(cap, 'capillary', 'pp', [], [ge]), rec(ge, 'GE_capillary', 'Albanese2016_v01', [cap], [target])]
        params += [(f'u_in_{cap}', 'J_per_m3', 2000.0, 'capillary inlet ~15 mmHg (chosen for the example)'),
                   (f'u_out_{cap}', 'J_per_m3', 666.6, 'capillary_network u_out_capillary_18'),
                   (f'u_ext_{cap}', 'J_per_m3', 0.0, 'capillary_network'),
                   (f'E_{cap}', 'J_per_m3', 400000.0, 'capillary_network E_capillary_6'),
                   (f'l_{cap}', 'metre', 0.0002, 'capillary_network l_capillary_6'),
                   (f'r_{cap}', 'metre', 2e-06, 'capillary_network r_capillary_6'),
                   (f'q_C_init_{cap}', 'm3', 0.0, 'capillary_network q_C_init_capillary_6'),
                   (f'C_O2_a_{ge}', 'mol_per_m3', 8.75, 'capillary_network C_O2_a_capillary_GE_0'),
                   (f'C_CO2_a_{ge}', 'mol_per_m3', 21.43, 'capillary_network C_CO2_a_capillary_GE_0'),
                   (f'C_O2_c_init_{ge}', 'mol_per_m3', 3.0, 'capillary_network C_O2_c_init_capillary_GE_6'),
                   (f'perm_O2_{ge}', 'm_per_s', 4e-06, 'capillary_network perm_O2_capillary_GE_6'),
                   (f'ub_O2_t_{ge}', 'mol_per_m3', O2['C_init'], 'starting tissue O2 (= C_P_init of the tissue)')]
        # the rest of GE_capillary's constants (its default instance has none yet)
        params += [(f'{n}_{ge}', u, v, f'capillary_network {n}_capillary_GE_6') for n, u, v in GE_CONSTANTS]
    return records, params


def neuron(varicosity_subtype, targets):
    """Periodic stimulus -> soma -> axon -> varicosity (library defaults)."""
    records = [rec('stim', 'i_stim_periodic', 'ArgusUNPUBLISHED_v01', [], ['soma_SN']),
               rec('soma_SN', 'soma', 'sympathetic_monolithic_v01', ['stim'], ['axon_SN']),
               rec('axon_SN', 'axon', 'sympathetic_monolithic_v01', ['soma_SN'], ['var_SN']),
               rec('var_SN', 'varicosity', varicosity_subtype, ['axon_SN'], targets)]
    params = [('t_start_stim', 'second', 0.01, 'first pulse after 10 ms at rest'),
              ('t_end_stim', 'second', 1000.0, 'pulses throughout'),
              ('period_stim', 'second', 0.05, '20 Hz train (a high sympathetic firing rate)'),
              ('pulse_duration_stim', 'second', 0.002, '2 ms current pulse'),
              ('i_amplitude_stim', 'nanoA', -1.0, 'suprathreshold: -1 nA (the soma takes I_in as an outward current, so a depolarising pulse is negative) fires one action potential per pulse'),
              ('stim_flag_stim', 'dimensionless', 1, 'pulse train on')]
    return records, params


def single_volume(sources):
    """One well-mixed extracellular volume per source, of the varicosity's volume (Vol 1e-18 m3)."""
    records, params = [], []
    for k, src in enumerate(sources):
        c = f'ecs_{k}'
        records.append(rec(c, 'tissue_diffusion_volume', 'Fang2008_v01', [src], []))
        params += [(f'V_P_{c}', 'm3', 1e-18, 'the varicosity volume Vol (k_V d^3 at the defaults)'),
                   (f'C_P_init_{c}', 'millimolar', 0.0, 'NE_init of the varicosity'),
                   (f'M_{c}', 'millimolar_per_s', 0.0, 'no clearance but NET (in the varicosity)'),
                   (f'C50_{c}', 'millimolar', 1.0, 'unused with M = 0')]
    return records, params, [r['name'] for r in records]


def write(model, records, params, spec):
    d = os.path.join(OUT, model)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f'{model}_module_array.json'), 'w') as f:
        f.write('[\n' + ',\n'.join(' ' + json.dumps(r) for r in records) + '\n]\n')
    with open(os.path.join(d, f'{model}_parameters.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['variable_name', 'units', 'value', 'data_reference'])
        w.writerows(params)
    spec_path = os.path.join(d, f'{model}_system.yaml')
    if not os.path.isfile(spec_path):
        with open(spec_path, 'w') as f:
            yaml.safe_dump({'model': model, 'category': 'coupled', 'reviewed': False, **spec}, f, sort_keys=False,
                           width=120)


COUPLED_NOTE = ('Built by tools/build_coupled_examples.py. Generated as C++ and run, with the FEniCS-coupled '
                'variants, by tests/test_coupled_systems.py (outputs and run times compared).')
SKIP = 'a coupled example: generated as C++ and run by tests/test_coupled_systems.py, not by the CellML system tests'


def main(argv=None):
    import argparse
    global N
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--n', type=int, default=N, help=f'grid cells per direction (default {N})')
    N = parser.parse_args(argv).n
    O2['L'], NE['L'] = N * H_O2, N * H_NE
    sources = ['capillary_GE_0', 'capillary_GE_1']
    o2_common = {'sim_time': 10.0, 'pre_time': 0.0, 'dt': 0.1, 'solver_info': {'rtol': 1e-8, 'atol': 1e-12},
                 'equivalence': {'status': 'not_applicable', 'reason': 'an example built from this library'},
                 'skip': SKIP}
    grid, grid_params, cells = fv_grid(sources, O2)
    caps, cap_params = capillaries(cells)
    write('microvasc_O2_FV', caps + grid, cap_params + grid_params + [
        ('saturation_cap', 'dimensionless', 0.9999, 'pulmonary_GE default (numerical stability)')],
        dict(o2_common, description='Two capillaries exchanging O2 with a cube of tissue (20 um cells) modelled as a '
             'finite-volume grid of CellML cells.', notes=[COUPLED_NOTE]))
    fe, fe_params = fenics_row(sources, O2)
    caps, cap_params = capillaries(['tissue', 'tissue'])
    write('microvasc_O2_FEniCS', caps + fe, cap_params + fe_params + [
        ('saturation_cap', 'dimensionless', 0.9999, 'pulmonary_GE default (numerical stability)'),
        ('sigma_diff', 'm2_per_s', O2['sigma'], 'tissue O2 diffusivity (secomb2020mass)')],
        dict(o2_common, description='Two capillaries exchanging O2 with a cube of tissue (20 um cells) modelled in FEniCSx '
             '(tissue_diffusion_FEniCS).', notes=[COUPLED_NOTE]))

    ne_common = {'sim_time': 0.2, 'pre_time': 0.0, 'dt': 1e-4, 'solver_info': {'rtol': 1e-8, 'atol': 1e-12},
                 'equivalence': {'status': 'not_applicable', 'reason': 'an example built from this library'},
                 'skip': SKIP}
    nrn, nrn_params = neuron('sympathetic_monolithic_v01', [])
    write('SN_NE_original', nrn, nrn_params,
          dict(ne_common, description='Periodic stimulus -> soma -> axon -> varicosity, NE an internal state of the '
               'varicosity: the reference for SN_NE_single_volume.', notes=[COUPLED_NOTE]))
    vol, vol_params, names = single_volume(['var_SN'])
    nrn, nrn_params = neuron('NEexchange_v01', names)
    write('SN_NE_single_volume', nrn + vol, nrn_params + vol_params,
          dict(ne_common, description='The varicosity (NEexchange_v01) releasing NE into one well-mixed volume equal '
               'to its own: exactly SN_NE_original.', notes=[COUPLED_NOTE]))
    grid, grid_params, cells = fv_grid(['var_SN'], NE)
    nrn, nrn_params = neuron('NEexchange_v01', cells)
    write('SN_NE_FV', nrn + grid, nrn_params + grid_params,
          dict(ne_common, description='The varicosity releasing NE into a cube of extracellular space, a '
               'finite-volume grid of CellML cells (1 um, the varicosity volume).', notes=[COUPLED_NOTE]))
    fe, fe_params = fenics_row(['var_SN'], NE)
    nrn, nrn_params = neuron('NEexchange_v01', ['tissue'])
    write('SN_NE_FEniCS', nrn + fe, nrn_params + fe_params + [
        ('sigma_diff', 'm2_per_s', NE['sigma'], 'effective NE diffusivity (as the NE_extracellular instance)')],
        dict(ne_common, description='The varicosity releasing NE into a cube of extracellular space (1 um cells) modelled in '
             'FEniCSx (tissue_diffusion_FEniCS, instance NE_extracellular).', notes=[COUPLED_NOTE]))


if __name__ == '__main__':
    main()
