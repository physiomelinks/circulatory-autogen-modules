"""Builds the space between the cells of a colony (AdaptOrDeity level 4, a prokaryote biofilm): finite volumes of
extracellular water holding the species the cells exchange and the matrix they make, and the faces between them.

    python tools/build_extracellular_space.py

Module versions (equations as infix text, tools/cellml_writer.py):

  transport/extracellular_volume/biofilm_ArgusUNPUBLISHED_v01   glucose, lactate, H+ (buffered), O2 and matrix in one
                                                                volume of the space between cells
  transport/extracellular_face/biofilm_ArgusUNPUBLISHED_v01     diffusion of the four species between two volumes (or a
                                                                volume and open water), slowed by the matrix on each side

Values come from modules/cell/generic_adaptive_cell/parameters.json (the cell family's parameters and their sources);
the matrix itself is made by the cell (energy_metabolism adaptive_ArgusUNPUBLISHED_v01, matrix_export port, built by
tools/build_adaptive_cell.py).

Geometry of the default instances: cells of radius R_cell on a triangular lattice, colony_gap apart, h_cell deep (the
game's sheet). One volume is the water of one lattice triangle less the three 60 degree sectors of its cells; a face is
the gap between two neighbouring cells, shared by the two triangles on either side of their link.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import module_writer  # noqa: E402
from module_writer import BIB, REPO, register_bib, write_version  # noqa: E402
from matrix_references import MATRIX_BIB  # noqa: E402

VERSION = 'biofilm_ArgusUNPUBLISHED_v01'
PARAMS = json.load(open(os.path.join(REPO, 'modules', 'cell', 'generic_adaptive_cell', 'parameters.json')))['parameters']
SPECIES = ['glc', 'lac', 'H', 'O2']
LN10 = math.log(10.0)

module_writer.UNITS.update({
    'm2_per_s': [('metre', None, 2, 1), ('second', None, -1, 1)],
    'm3_per_s': [('metre', None, 3, 1), ('second', None, -1, 1)],
    'm3_per_mol': [('metre', None, 3, 1), ('mole', None, -1, 1)],
})
register_bib(MATRIX_BIB)
register_bib({
    'argus2026unpublished': dict(
        type='unpublished', author='Argus, Finbar', title='Unpublished modules of the circulatory-autogen module library',
        year='2026', note='Placeholder citation (required_citations, 2026-10-06): replace with the real reference when it '
                          'exists: the owner\'s unpublished or modified models, no publication'),
})


def val(key):
    return PARAMS[key]['value']


def row(name, key, how=None):
    """A CSV row whose value is the parameters.json entry (its source and status), or derived from it (how)."""
    e = PARAMS[key]
    if how is None:
        return (name, float(e['value']), f"{e['source']}; {e['note']}", e['status'] == 'sourced')
    raise ValueError(how)


def derived(name, value, how, *keys):
    ref = 'derived: ' + how + ' (' + '; '.join(f"{k}: {PARAMS[k]['source']}, {PARAMS[k]['status']}" for k in keys) + ')'
    return (name, float(value), ref, all(PARAMS[k]['status'] in ('sourced', 'definitional') for k in keys))


# ---- geometry of the default instances (um, then SI)
R_CELL, H_CELL, GAP = val('R_cell'), val('h_cell'), val('colony_gap')
SPACING = 2 * R_CELL + GAP
AREA_UM2 = math.sqrt(3) / 4 * SPACING ** 2 - math.pi * R_CELL ** 2 / 2   # triangle less three 60 degree sectors
V_E = AREA_UM2 * H_CELL * 1e-18                                         # m3
A_F = GAP * H_CELL * 1e-12                                              # m2: the gap between two cells, h deep
D_HALF = SPACING / (2 * math.sqrt(3)) * 1e-6                            # m: triangle centroid to the shared edge
GEOM = (f'geometry: cells of radius {R_CELL:g} um, {GAP:g} um apart on a triangular lattice, {H_CELL:g} um deep '
        f'(R_cell, h_cell, colony_gap)')

# ---- seawater and the matrix
H_SEA = 1000.0 * 10 ** -val('pH_o')                                     # mM free H+
B_H = val('beta_seawater') / (LN10 * H_SEA)                             # bound/free: d(bound)/d(free) of a linear buffer
V_MX = val('M_glucose_unit') * val('polysaccharide_specific_volume')    # m3 of polymer per mol glucose units
C_MX_MAX = val('matrix_gel_density_max') / val('M_glucose_unit')        # mol/m3 = mM glucose units
D_H_EFF = val('D_H_water') + B_H * val('D_buffer_water')                # free protons plus those the buffer carries

# =================================================================== the volume
vol_vars = [('t', 'second', 'in', None, None)]
for s in SPECIES:
    vol_vars += [(f'C_{s}_init', 'millimolar', 'in', None, None),
                 (f'C_{s}', 'millimolar', 'out', f'C_{s}_init',
                  'free H+ (the bound share is b_H C_H)' if s == 'H' else None),
                 (f'J_{s}_in', 'mol_per_s', 'in', None,
                  'summed fluxes of the faces on the entrance and exit extracellular_face_ports' if s == 'glc' else None),
                 (f'J_{s}_out', 'mol_per_s', 'in', None, None),
                 (f'j_{s}_cell', 'mol_per_s', 'in', None,
                  'summed exchange with the cells around the volume, positive into it' if s == 'glc' else None)]
vol_vars += [
    ('C_mx_init', 'millimolar', 'in', None, None),
    ('C_mx', 'millimolar', 'out', 'C_mx_init', 'matrix (exopolysaccharide) as mM of glucose units; it does not diffuse'),
    ('j_mx_cell', 'mol_per_s', 'in', None, 'matrix the cells deliver into the volume (their matrix_export, summed)'),
    ('V_e', 'm3', 'in', None, 'water volume between the cells'),
    ('b_H', 'dimensionless', 'in', None, 'buffered H+ per free H+ (a linear buffer: capacity V_e (1 + b_H))'),
    ('c_std', 'millimolar', 'in', None, None),
    ('v_mx', 'm3_per_mol', 'in', None, 'polymer volume per mol of glucose units (partial specific volume x molar mass)'),
    ('k_mx', 'per_s', 'in', None, 'breakdown rate of the matrix'),
    ('f_mx_glc', 'dimensionless', 'in', None, 'share of the broken-down matrix released as glucose'),
    ('C_mx_max', 'millimolar', 'in', None, 'densest matrix gel the volume holds'),
    ('k_ovf', 'per_s', 'in', None, 'rate at which matrix above C_mx_max leaves for open water'),
    ('phi', 'dimensionless', 'out', None, 'polymer volume fraction'),
    ('f_hind', 'dimensionless', 'out', None, 'diffusion in the matrix over diffusion in water (Mackie and Meares 1955)'),
    ('r_mx', 'millimolar_per_s', 'out', None, 'matrix broken down'),
    ('r_ovf', 'millimolar_per_s', 'out', None, 'matrix squeezed out above the densest gel'),
    ('J_mx_lost', 'mol_per_s', 'out', None, 'glucose units of matrix leaving the colony (not released as glucose, '
                                            'or squeezed out)'),
    ('N_mx_lost', 'mole', 'out', 0.0, 'glucose units of matrix that have left the colony since the start (the '
                                      'running total of J_mx_lost, for carbon accounting)'),
    ('pH', 'dimensionless', 'out', None, None),
]
vol_eqs = [
    ('phi', 'v_mx * C_mx'),
    ('f_hind', '((1.0 - phi) / (1.0 + phi)) ** 2.0'),
    ('r_mx', 'k_mx * C_mx'),
    ('r_ovf', 'k_ovf * C_mx * (C_mx / C_mx_max) ** 8.0'),
    ('J_mx_lost', 'V_e * ((1.0 - f_mx_glc) * r_mx + r_ovf)'),
    ('pH', '-ln(C_H / c_std) / 2.302585092994046'),
    ('d/dt C_glc', '(J_glc_in + J_glc_out + j_glc_cell) / V_e + f_mx_glc * r_mx'),
    ('d/dt C_lac', '(J_lac_in + J_lac_out + j_lac_cell) / V_e'),
    ('d/dt C_H', '(J_H_in + J_H_out + j_H_cell) / (V_e * (1.0 + b_H))'),
    ('d/dt C_O2', '(J_O2_in + J_O2_out + j_O2_cell) / V_e'),
    ('d/dt C_mx', 'j_mx_cell / V_e - r_mx - r_ovf'),
    ('d/dt N_mx_lost', 'J_mx_lost'),
]
vol_roles = {v[0]: 'constant' for v in vol_vars if v[0] != 't'}
vol_roles.update({f'C_{s}': 'variable' for s in SPECIES + ['mx']})
vol_roles.update({f'J_{s}_{d}': 'variable' for s in SPECIES for d in ('in', 'out')})
vol_roles.update({f'j_{s}_cell': 'boundary_condition' for s in SPECIES + ['mx']})
vol_roles.update({n: 'variable' for n in ('phi', 'f_hind', 'r_mx', 'r_ovf', 'J_mx_lost', 'N_mx_lost', 'pH')})
shared = [f'C_{s}' for s in SPECIES] + ['f_hind']
vol_ports = {
    'entrance_ports': [{'port_type': 'extracellular_face_port', 'variables': shared + [f'J_{s}_in' for s in SPECIES],
                        'multi_port': ['True'] * 5 + ['sum'] * 4}],
    'exit_ports': [{'port_type': 'extracellular_face_port', 'variables': shared + [f'J_{s}_out' for s in SPECIES],
                    'multi_port': ['True'] * 5 + ['sum'] * 4}],
    'general_ports': [
        # the cells around the volume, as one exchange (a single port, not a summed multiport, so that unconnected
        # its fluxes stay parameters a host such as the game sets): they read its concentrations and add their
        # exchange and matrix
        {'port_type': 'extracellular_exchange_port',
         'variables': [f'C_{s}' for s in SPECIES] + [f'j_{s}_cell' for s in SPECIES + ['mx']]},
    ],
}
SEA = {'glc': (val('glc_e'), 'glc_e'), 'lac': (0.0, None), 'H': (H_SEA, 'pH_o'), 'O2': (val('C_O2_proterozoic'),
                                                                                        'C_O2_proterozoic')}
vol_params = []
for s in SPECIES:
    v, key = SEA[s]
    vol_params.append(derived(f'C_{s}_init', v, 'free H+ = 10^-pH_o M', 'pH_o') if s == 'H' else
                      (row(f'C_{s}_init', key) if key else (f'C_{s}_init', 0.0, 'definitional; no lactate in the sea',
                                                            False)))
vol_params += [(f'J_{s}_{d}', 0.0, 'no face connected', False) for s in SPECIES for d in ('in', 'out')]
vol_params += [(f'j_{s}_cell', 0.0, 'no cell connected', False) for s in SPECIES + ['mx']]
vol_params += [
    ('C_mx_init', 0.0, 'definitional; no matrix until the cells make it', False),
    ('V_e', V_E, f'{GEOM}: triangle of side {SPACING:g} um less three 60 degree sectors, x depth', False),
    derived('b_H', B_H, 'beta_seawater / (ln 10 x free H+ at pH_o): bound per free H+ of a linear buffer', 'beta_seawater',
            'pH_o'),
    ('c_std', 1000.0, 'definitional; standard concentration 1 M (for pH)', True),
    derived('v_mx', V_MX, 'M_glucose_unit x polysaccharide_specific_volume', 'M_glucose_unit',
            'polysaccharide_specific_volume'),
    row('k_mx', 'matrix_decay_rate'), row('f_mx_glc', 'matrix_glucose_return'),
    derived('C_mx_max', C_MX_MAX, 'matrix_gel_density_max / M_glucose_unit', 'matrix_gel_density_max', 'M_glucose_unit'),
    row('k_ovf', 'matrix_overflow_rate'),
]

# harness: the volume between two open-water faces, with cells exchanging and delivering matrix
FACES = ['fN', 'fS']
harness_params = []
for f in FACES:
    side = 'up' if f == 'fN' else 'down'
    harness_params += [[f'A_f_{f}', 'm2', A_F, 'harness: one cell-cell gap, h deep'],
                       [f'd_up_{f}', 'metre', D_HALF, 'harness: centroid to the gap'],
                       [f'd_down_{f}', 'metre', 0.0 if f == 'fS' else D_HALF, 'harness: the open side'],
                       [f'open_flag_{f}', 'dimensionless', 1, 'harness: open face'],
                       [f'f_{side}_{f}', 'dimensionless', 1.0, 'harness: open water, no matrix']]
    for s in SPECIES:
        harness_params.append([f'C_{s}_{side}_{f}', 'millimolar', SEA[s][0] * (1.5 if f == 'fN' else 1.0) +
                               (0.01 if s == 'lac' and f == 'fN' else 0.0), 'harness: open water'])
face_sum = {s: ' + '.join(['fN__J_%s_down' % s, 'fS__J_%s_up' % s]) for s in SPECIES}


def cap(s):
    return 'V_e*(1 + b_H)' if s == 'H' else 'V_e'


vol_invariants = []
for s in SPECIES:
    src = f'({face_sum[s]}) + j_{s}_cell' + (' + V_e*f_mx_glc*k_mx*C_mx' if s == 'glc' else '')
    moved = f'np.trapz(np.abs(fN__J_{s}_down) + np.abs(fS__J_{s}_up) + abs(j_{s}_cell), t)'
    vol_invariants.append({
        'description': f'{s} balance: {cap(s)} dC_{s} = integral of the face fluxes and the cells\' exchange'
                       + (' plus the glucose released by matrix breakdown' if s == 'glc' else '')
                       + ' (trapezoid; 1e-3 of the amount moved)',
        'expr': f'abs({cap(s)}*(C_{s}[-1] - C_{s}[0]) - np.trapz({src}, t)) <= 1e-3*{moved} + 1e-30'})
GLC_MOVED = 'np.trapz(np.abs(fN__J_glc_down) + np.abs(fS__J_glc_up) + abs(j_glc_cell) + j_mx_cell + J_mx_lost, t)'
vol_invariants += [
    {'description': 'matrix balance: V_e dC_mx = integral of the matrix delivered less breakdown and overflow (trapezoid; '
                    '1e-3 of the amount moved)',
     'expr': 'abs(V_e*(C_mx[-1] - C_mx[0]) - np.trapz(j_mx_cell - V_e*(r_mx + r_ovf), t)) <= '
             '1e-3*np.trapz(j_mx_cell + V_e*(r_mx + r_ovf), t) + 1e-30'},
    {'description': 'carbon (glucose units): V_e (C_glc + C_mx) changes by the face and cell glucose fluxes and the matrix '
                    'delivered, less J_mx_lost (trapezoid; 1e-3 of the amount moved)',
     'expr': f'abs(V_e*((C_glc[-1] + C_mx[-1]) - (C_glc[0] + C_mx[0])) - np.trapz(({face_sum["glc"]}) + j_glc_cell + '
             f'j_mx_cell - J_mx_lost, t)) <= 1e-3*{GLC_MOVED} + 1e-30'},
    {'description': 'carbon held: V_e (C_glc + C_mx) + N_mx_lost changes only by the face and cell glucose fluxes and '
                    'the matrix delivered (trapezoid; 1e-3 of the amount moved)',
     'expr': f'abs(V_e*((C_glc[-1] + C_mx[-1]) - (C_glc[0] + C_mx[0])) + N_mx_lost[-1] - N_mx_lost[0] - '
             f'np.trapz(({face_sum["glc"]}) + j_glc_cell + j_mx_cell, t)) <= 1e-3*{GLC_MOVED} + 1e-30'},
    {'description': 'definitions (exact): phi = v_mx C_mx, f_hind = ((1 - phi)/(1 + phi))^2, pH = -log10(C_H/c_std), '
                    'J_mx_lost = V_e ((1 - f_mx_glc) k_mx C_mx + k_ovf C_mx (C_mx/C_mx_max)^8)',
     'expr': 'np.allclose(phi, v_mx*C_mx, rtol=1e-12, atol=0) and np.allclose(f_hind, ((1 - phi)/(1 + phi))**2, '
             'rtol=1e-12, atol=0) and np.allclose(pH, -np.log10(C_H/c_std), rtol=1e-12, atol=0) and '
             'np.allclose(J_mx_lost, V_e*((1 - f_mx_glc)*k_mx*C_mx + k_ovf*C_mx*(C_mx/C_mx_max)**8), rtol=1e-12, atol=0)'},
    {'applies': 'run',
     'description': 'steady state of the lactate and O2 by t_end (no source in the volume): the face fluxes balance '
                    'the cells\' exchange (1e-6 of the largest term)',
     'expr': f'abs(({face_sum["lac"]})[-1] + j_lac_cell) <= 1e-6*max(np.max(np.abs(fN__J_lac_down)), abs(j_lac_cell), '
             f'1e-30) and abs(({face_sum["O2"]})[-1] + j_O2_cell) <= 1e-6*max(np.max(np.abs(fN__J_O2_down)), '
             f'abs(j_O2_cell), 1e-30)'},
    {'description': 'with no cell exchange the concentrations stay between the open-water values and their start '
                    '(maximum principle; glucose may rise by the matrix it releases)',
     'expr': '(j_lac_cell != 0) or (np.all(C_lac >= min(C_lac[0], fN__C_lac_up[0], fS__C_lac_down[0]) - 1e-12) and '
             'np.all(C_lac <= max(C_lac[0], fN__C_lac_up[0], fS__C_lac_down[0]) + 1e-12))'},
    {'description': 'the matrix stays between 0 and its densest gel plus 15% (overflow holds the ceiling)',
     'expr': 'np.all(C_mx >= -1e-12) and np.all(C_mx <= 1.15*max(C_mx_max, C_mx[0]))'},
    {'description': 'second law: each species diffuses down its gradient across both faces, J (C_upstream - C_downstream) '
                    '>= 0',
     'expr': ' and '.join(f'np.all(fN__J_{s}_down*(fN__C_{s}_up - C_{s}) >= -1e-30) and np.all(fS__J_{s}_up*(fS__C_{s}_down - '
                          f'C_{s}) >= -1e-30)' for s in SPECIES)},
]
vol_verif = {
    'sim_time': 0.1, 'pre_time': 0.0, 'dt': 0.0001, 'solver': 'CVODE_myokit',
    'solver_info': {'rtol': 1e-10, 'atol': 1e-16, 'MaximumStep': 0.0001},
    'outputs': [f'C_{s}' for s in SPECIES] + ['C_mx', 'N_mx_lost', 'phi', 'f_hind', 'r_mx', 'r_ovf', 'J_mx_lost', 'pH', 'V_e'] +
               [f'fN/J_{s}_down' for s in SPECIES] + [f'fS/J_{s}_up' for s in SPECIES] +
               [f'fN/C_{s}_up' for s in SPECIES] + [f'fS/C_{s}_down' for s in SPECIES],
    # a cell takes glucose and releases lactate and acid into the volume and delivers matrix; the matrix starts above
    # the densest gel so the overflow runs
    'run_parameters': {'j_glc_cell': -2e-20, 'j_lac_cell': 4e-20, 'j_H_cell': 4e-20, 'j_O2_cell': 0.0,
                       'j_mx_cell': 1e-20, 'C_mx_init': 1.1 * C_MX_MAX},
    'harness': {
        'module_array': [['mod', VERSION, 'extracellular_volume', 'fN', 'fS', 'default'],
                         ['fN', VERSION, 'extracellular_face', '', 'mod', 'default'],
                         ['fS', VERSION, 'extracellular_face', 'mod', '', 'default']],
        'parameters': harness_params,
    },
    'invariants': vol_invariants,
    'bc_sweep': {
        'factors': [0.5, 0.75, 1.0, 1.5, 2.0],
        'ranges': {'C_mx_init': {'min': 0.0, 'max': 1.1 * C_MX_MAX, 'points': 4},
                   'j_glc_cell': {'min': -4e-20, 'max': 0.0, 'points': 3},
                   'j_lac_cell': {'min': 0.0, 'max': 8e-20, 'points': 3},
                   'j_mx_cell': {'min': 0.0, 'max': 4.5e-20, 'points': 3},
                   'V_e': {'min': 0.3 * V_E, 'max': 3 * V_E, 'points': 3, 'scale': 'log'},
                   'f_mx_glc': {'min': 0.0, 'max': 1.0, 'points': 3},
                   'k_mx': {'min': 0.0, 'max': 1e-4, 'points': 3}},
        'exclude': [f'J_{s}_{d}' for s in SPECIES for d in ('in', 'out')],
        'rationale': 'the cells\' exchange over the range a 2 um fermenter gives (resting glucose flux 4.5e-19 mol/s over '
                     'its 6 volumes, up to 10x when fed), the matrix from none to above the densest gel, the volume from '
                     'tight to loose packing; J_in and J_out are the summed face fluxes',
    },
    'timestep': {'scheme': 'rk4', 'dts': [4e-5, 2e-5, 1e-5, 5e-6], 't_end': 0.004, 'min_order': 0.8, 'tol': 0.0001,
                 'cvode_tol': 0.001},
    'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-06, 'atol': 1e-12}]},
    'validation': {'default': {
        'baseline': {'status': 'not_applicable',
                     'reason': 'one volume of the space between cells has no measured counterpart on its own; verified '
                               'by the balances, the steady state and the definitions. Biofilm diffusion data (Stewart '
                               '2003: De/Daq 0.25 for organic solutes, 0.6 for gases, cells and matrix together) apply to '
                               'an assembled colony.'},
        'calibrate': {'note': 'not_applicable: as baseline'}}},
}
vol_notes = (
    'One finite volume of the water between the cells of a colony (AdaptOrDeity level 4: a prokaryote biofilm, unsealed): '
    'glucose, lactate, H+ and O2 as concentrations, and the matrix (exopolysaccharide, EPS) as mM of glucose units. '
    'dC_s/dt = (J_s_in + J_s_out + j_s_cell)/V_e for each species: J_s_in and J_s_out sum the extracellular faces on the '
    'entrance and exit ports (list-form multi_port, as tissue_diffusion_volume), j_s_cell is the exchange with the cells around it '
    '(extracellular_exchange_port, a single port so that unconnected the exchange stays a parameter its host sets; '
    'positive into the volume). H+ is buffered by the seawater (a linear buffer, b_H bound '
    'per free H+, so its capacity is V_e (1 + b_H) and its pH = -log10(C_H/c_std)); the faces carry the buffer\'s share of '
    'the protons in their D_H. Na+, K+ and Ca2+ stay at seawater values and are not carried. The matrix does not diffuse: '
    'dC_mx/dt = j_mx_cell/V_e - k_mx C_mx - k_ovf C_mx (C_mx/C_mx_max)^8. The cells deliver it (their matrix_export, '
    'energy_metabolism adaptive_ArgusUNPUBLISHED_v01); the colony\'s hydrolases break it down at k_mx, releasing f_mx_glc of '
    'it as glucose in the volume (EPS is biodegradable, a reserve in starvation: Zhang and Bishop 2003); matrix above the '
    'densest gel C_mx_max (biofilm EPS of order 10-100 g/L, Flemming and Wingender 2010) is squeezed out into open water. '
    'J_mx_lost is what leaves the colony, and N_mx_lost its running total. The polymer volume fraction phi = v_mx C_mx slows diffusion by f_hind = '
    '((1 - phi)/(1 + phi))^2 (Mackie and Meares 1955), shared with the faces; valid for phi < 1, and the ceiling keeps phi '
    'about 0.03. Default instance: the gap between three cells of radius 1 um, 0.3 um apart on a triangular lattice, '
    '2 um deep (1.44 um3), filled with Proterozoic seawater. Tested between two open-water faces with a cell taking '
    'glucose and releasing lactic acid and matrix: the balance of each species, of the matrix and of carbon, the steady '
    'state, the definitions, the ceiling and the second law across the faces.')
write_version('transport', 'extracellular_volume', 'extracellular_volume', vol_vars, vol_eqs,
              {'required_citations': ['argus2026unpublished', 'mackie1955diffusion'],
               **vol_ports, 'roles': vol_roles},
              vol_params, vol_verif, {'reviewed': False, 'notes': vol_notes},
              ['argus2026unpublished', 'amsden1998solute', 'flemming2010biofilm', 'mackie1955diffusion',
               'stewart2003diffusion', 'zhang2003biodegradability'], VERSION)

# =================================================================== the face
face_vars = [('t', 'second', 'in', None, None)]
for s in SPECIES:
    face_vars += [(f'C_{s}_up', 'millimolar', 'in', None, None), (f'C_{s}_down', 'millimolar', 'in', None, None),
                  (f'D_{s}', 'm2_per_s', 'in', None, 'diffusivity in water' + (' (free H+ plus the buffer\'s share)'
                                                                                if s == 'H' else '')),
                  (f'G_{s}', 'm3_per_s', 'out', None, 'conductance through the two half-cells in series'),
                  (f'J_{s}', 'mol_per_s', 'out', None, 'into the downstream volume'),
                  (f'J_{s}_up', 'mol_per_s', 'out', None, None), (f'J_{s}_down', 'mol_per_s', 'out', None, None)]
face_vars += [
    ('f_up', 'dimensionless', 'in', None, 'matrix hindrance of the upstream volume (its f_hind; 1 in open water)'),
    ('f_down', 'dimensionless', 'in', None, None),
    ('A_f', 'm2', 'in', None, 'face area: the gap between two cells times the depth'),
    ('d_up', 'metre', 'in', None, 'upstream centre to the face'), ('d_down', 'metre', 'in', None, None),
    ('open_flag', 'dimensionless', 'in', None, None),
]
face_eqs = []
for s in SPECIES:
    face_eqs += [(f'G_{s}', f'open_flag * A_f * D_{s} * f_up * f_down / (d_up * f_down + d_down * f_up)'),
                 (f'J_{s}', f'G_{s} * (C_{s}_up - C_{s}_down)'), (f'J_{s}_down', f'J_{s}'), (f'J_{s}_up', f'-J_{s}')]
face_roles = {v[0]: 'constant' for v in face_vars if v[0] != 't'}
face_roles.update({n: 'boundary_condition' for s in SPECIES for n in (f'C_{s}_up', f'C_{s}_down')})
face_roles.update({'f_up': 'boundary_condition', 'f_down': 'boundary_condition'})
face_roles.update({n: 'variable' for s in SPECIES for n in (f'G_{s}', f'J_{s}', f'J_{s}_up', f'J_{s}_down')})
face_ports = {
    'entrance_ports': [{'port_type': 'extracellular_face_port',
                        'variables': [f'C_{s}_up' for s in SPECIES] + ['f_up'] + [f'J_{s}_up' for s in SPECIES]}],
    'exit_ports': [{'port_type': 'extracellular_face_port',
                    'variables': [f'C_{s}_down' for s in SPECIES] + ['f_down'] + [f'J_{s}_down' for s in SPECIES]}],
    'general_ports': [],
}
face_params = []
for s in SPECIES:
    v, key = SEA[s]
    for side in ('up', 'down'):
        name = f'C_{s}_{side}'
        face_params.append(derived(name, v, 'free H+ = 10^-pH_o M: open water when this side is unconnected', 'pH_o')
                           if s == 'H' else (row(name, key) if key else (name, 0.0, 'definitional; no lactate in the sea',
                                                                         False)))
face_params += [row('D_glc', 'D_glc_water'), row('D_lac', 'D_lac_water'),
                derived('D_H', D_H_EFF, 'D_H_water + b_H x D_buffer_water: free protons plus those the seawater buffer '
                                        'carries (b_H = beta_seawater/(ln 10 H+))', 'D_H_water', 'D_buffer_water',
                        'beta_seawater', 'pH_o'),
                row('D_O2', 'D_O2_water'),
                ('f_up', 1.0, 'definitional; no matrix (open water) when unconnected', False),
                ('f_down', 1.0, 'definitional; no matrix (open water) when unconnected', False),
                ('A_f', A_F, f'{GEOM}: gap x depth', False),
                ('d_up', D_HALF, f'{GEOM}: lattice-triangle centroid to the gap, s/(2 sqrt 3)', False),
                ('d_down', D_HALF, f'{GEOM}: lattice-triangle centroid to the gap, s/(2 sqrt 3)', False),
                ('open_flag', 1.0, 'definitional; open face', False)]

face_harness = []
for v, c0, mx in (('va', 1.0, 0.0), ('vb', 0.2, 0.5)):
    face_harness += [[f'C_glc_init_{v}', 'millimolar', 0.1 * c0, 'harness: two volumes at different concentrations'],
                     [f'C_lac_init_{v}', 'millimolar', 0.5 * c0, 'harness'],
                     [f'C_H_init_{v}', 'millimolar', H_SEA * (1 + c0), 'harness'],
                     [f'C_O2_init_{v}', 'millimolar', val('C_O2_proterozoic') * c0, 'harness'],
                     [f'C_mx_init_{v}', 'millimolar', mx * C_MX_MAX, 'harness: matrix in one volume'],
                     [f'k_mx_{v}', 'per_s', 0.0, 'harness: matrix held, so the hindrance is constant'],
                     [f'V_e_{v}', 'm3', V_E * (1.0 if v == 'va' else 1.5), 'harness: two volumes of different size']]
face_invariants = [
    {'description': 'G_s = open_flag A_f D_s f_up f_down/(d_up f_down + d_down f_up) (two half-cells in series), J_s = '
                    'G_s (C_up - C_down) with C_up, C_down and f_up, f_down from the connected volumes, J_down = J, J_up = '
                    '-J (exact)',
     'expr': ' and '.join(f'np.allclose(G_{s}, open_flag*A_f*D_{s}*va__f_hind*vb__f_hind/(d_up*vb__f_hind + '
                          f'd_down*va__f_hind), rtol=1e-12, atol=0) and np.allclose(J_{s}, G_{s}*(va__C_{s} - vb__C_{s}), '
                          f'rtol=1e-9, atol=1e-35) and np.array_equal(J_{s}_down, J_{s}) and np.array_equal(J_{s}_up, '
                          f'-J_{s})' for s in SPECIES)},
    {'description': 'mass is conserved across the face for every species: V_a C_a + V_b C_b constant (H+ with its buffer, '
                    'V (1 + b_H) C_H; relative 1e-9)',
     'expr': ' and '.join(f'np.ptp({w}va__V_e*va__C_{s} + {w}vb__V_e*vb__C_{s}) <= 1e-9*np.max(np.abs({w}va__V_e*va__C_{s}'
                          f' + {w}vb__V_e*vb__C_{s}))' for s in SPECIES
                          for w in ['(1 + va__b_H)*' if s == 'H' else ''])},
    {'description': 'exact solution with the matrix held: C_a - C_b = (C_a0 - C_b0) exp(-G (1/K_a + 1/K_b) t), K the '
                    'capacity (V, or V (1 + b_H) for H+) (relative 1e-6 of the initial difference)',
     'expr': ' and '.join(f'np.max(np.abs((va__C_{s} - vb__C_{s}) - (va__C_{s}[0] - vb__C_{s}[0])*np.exp(-G_{s}[0]*'
                          f'(1/({k}va__V_e[0]) + 1/({k}vb__V_e[0]))*t))) <= 1e-6*abs(va__C_{s}[0] - vb__C_{s}[0]) + 1e-18'
                          for s in SPECIES for k in ['(1 + va__b_H)*' if s == 'H' else ''])},
    {'description': 'second law: each species flows down its gradient, J_s (C_up - C_down) >= 0, so the entropy '
                    'produced, sum J_s R ln(C_up/C_down) / T, is never negative',
     'expr': ' and '.join(f'np.all(J_{s}*(va__C_{s} - vb__C_{s}) >= 0)' for s in SPECIES)},
    {'description': 'the matrix slows diffusion: with matrix on one side G_s is below the open-water conductance '
                    'A_f D_s/(d_up + d_down)',
     'expr': ' and '.join(f'np.all(G_{s} <= open_flag*A_f*D_{s}/(d_up + d_down)*(1 + 1e-12))' for s in SPECIES)},
]
face_verif = {
    'sim_time': 0.05, 'pre_time': 0.0, 'dt': 0.0001, 'solver': 'CVODE_myokit',
    'solver_info': {'rtol': 1e-11, 'atol': 1e-18, 'MaximumStep': 0.0005},
    'outputs': [n for s in SPECIES for n in (f'G_{s}', f'J_{s}', f'J_{s}_up', f'J_{s}_down', f'va/C_{s}', f'vb/C_{s}')] +
               ['va/V_e', 'vb/V_e', 'va/f_hind', 'vb/f_hind', 'va/b_H', 'D_glc', 'D_lac', 'D_H', 'D_O2'],
    'harness': {
        'module_array': [['va', VERSION, 'extracellular_volume', '', 'mod', 'default'],
                         ['mod', VERSION, 'extracellular_face', 'va', 'vb', 'default'],
                         ['vb', VERSION, 'extracellular_volume', 'mod', '', 'default']],
        'parameters': face_harness,
    },
    'invariants': face_invariants,
    'bc_sweep': {
        'factors': [0.5, 0.75, 1.0, 1.5, 2.0],
        'ranges': {'A_f': {'min': 0.2 * A_F, 'max': 5 * A_F, 'points': 3, 'scale': 'log'},
                   'd_up': {'min': 0.5 * D_HALF, 'max': 2 * D_HALF, 'points': 3},
                   'd_down': {'min': 0.5 * D_HALF, 'max': 2 * D_HALF, 'points': 3},
                   'open_flag': {'min': 0, 'max': 1, 'points': 2},
                   'D_glc': {'min': 1e-10, 'max': 1e-9, 'points': 3, 'scale': 'log'},
                   'C_mx_init': {'min': 0.0, 'max': C_MX_MAX, 'points': 3}},
        'exclude': [f'C_{s}_{d}' for s in SPECIES for d in ('up', 'down')] + ['f_up', 'f_down'],
        'rationale': 'the face geometry from tight to loose packing, the diffusivity over the species\' range and the '
                     'matrix from none to the densest gel; concentrations and hindrance come from the volumes',
    },
    'timestep': {'scheme': 'rk4', 'dts': [4e-5, 2e-5, 1e-5, 5e-6], 't_end': 0.004, 'min_order': 0.8, 'tol': 0.0001,
                 'cvode_tol': 0.001},
    'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-06, 'atol': 1e-12}]},
    'validation': {'default': {
        'baseline': {'status': 'not_applicable',
                     'reason': 'a face is the flux law between two volumes; verified against the exact two-volume '
                               'solution, conservation and the second law. Biofilm diffusion data apply to an assembled '
                               'colony.'},
        'calibrate': {'note': 'not_applicable: as baseline'}}},
}
face_notes = (
    'Diffusion of glucose, lactate, H+ and O2 through the gap between two neighbouring cells, from one volume of the '
    'space between cells to the next (or to open water: with one side unconnected its concentrations are fixed open-water '
    'values and its f is 1). Each half of the path is slowed by its volume\'s matrix (f_hind, Mackie and Meares 1955), the '
    'two halves in series: G_s = open_flag A_f D_s / (d_up/f_up + d_down/f_down), written as A_f D_s f_up f_down/(d_up '
    'f_down + d_down f_up) so a half of length 0 (open water at the mouth) is allowed. J_s = G_s (C_up - C_down) into the '
    'downstream volume (J_down = J, J_up = -J), so the pair conserves mass. D_H is the protons\' effective diffusivity in '
    'buffered seawater (free H+ plus the buffer\'s share, D_H_water + b_H D_buffer_water), matching the volume\'s buffered '
    'capacity. Diffusivities in water at 20 C: Stewart 2003 (glucose, O2) and limiting ionic conductivities (lactate, '
    'H+). Tested between two closed volumes of different sizes, one holding matrix: the flux law, conservation of every '
    'species, the exact two-volume solution and the second law.')
write_version('transport', 'extracellular_face', 'extracellular_face', face_vars, face_eqs,
              {'required_citations': ['argus2026unpublished', 'mackie1955diffusion'],
               **face_ports, 'roles': face_roles},
              face_params, face_verif, {'reviewed': False, 'notes': face_notes},
              ['argus2026unpublished', 'mackie1955diffusion', 'stewart2003diffusion'], VERSION)

print(f'volume {V_E * 1e18:.4g} um3, face {A_F * 1e12:.3g} um2 over {2 * D_HALF * 1e6:.4g} um; b_H {B_H:.4g}, D_H '
      f'{D_H_EFF:.4g} m2/s; v_mx {V_MX:.4g} m3/mol, C_mx_max {C_MX_MAX:.4g} mM (phi {V_MX * C_MX_MAX:.3g}, f_hind '
      f'{((1 - V_MX * C_MX_MAX) / (1 + V_MX * C_MX_MAX)) ** 2:.3g})')
