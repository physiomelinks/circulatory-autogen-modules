"""Writes the mitochondrial energy-metabolism version, the heat modules (thermal_volume, heat_conduction_face) and
membrane_carrier_face, with their calibration, from equations written as infix text (tools/cellml_writer.py).

    python tools/build_cell_modules.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from module_writer import SOLVER, UNITS, register_bib, write_version  # noqa: E402

VERSION = 'ArgusUNPUBLISHED_v01'
V_CELL = math.pi * 1e-5 ** 2 * 2e-5     # m3: 20 um diameter disc cell, 20 um deep (the game's starting cell)

# ------------------------------------------------------------------ references
BIB = {
    'argus2026unpublished': dict(type='unpublished', author='Argus, Finbar',
                                 title='Unpublished modules of the circulatory-autogen module library', year='2026',
                                 note='Placeholder citation (required_citations, 2026-10-06): replace with the real '
                                      'reference when it exists: the owner\'s unpublished or modified models, no publication'),
    'gnaiger1990anaerobic': dict(type='article', author='Gnaiger, Erich and Kemp, Richard B.',
                                 title='Anaerobic metabolism in aerobic mammalian cells: information from the ratio of '
                                       'calorimetric heat flux and respirometric oxygen flux',
                                 journal='Biochimica et Biophysica Acta (BBA) - Bioenergetics', year='1990',
                                 volume='1016', number='3', pages='328--332',
                                 note='Oxycaloric equivalent -470 kJ/mol O2 for aerobic catabolism of glucose; net '
                                      'glucose to lactate in bicarbonate buffer -63 kJ/mol lactate'),
    'hinkle2005p': dict(type='article', author='Hinkle, Peter C.',
                             title='P/O ratios of mitochondrial oxidative phosphorylation',
                             journal='Biochimica et Biophysica Acta (BBA) - Bioenergetics', year='2005',
                             volume='1706', number='1-2', pages='1--11', doi='10.1016/j.bbabio.2004.09.004',
                             note='P/O about 2.5 (NADH) and 1.5 (succinate); about 30 ATP per glucose oxidised'),
    'gajewski1986thermodynamics': dict(type='article', author='Gajewski, Ellen and Steckler, Donald K. and Goldberg, Robert N.',
                                 title='Thermodynamics of the hydrolysis of adenosine 5\'-triphosphate to adenosine '
                                       '5\'-diphosphate',
                                 journal='Journal of Biological Chemistry', year='1986', volume='261', number='27',
                                 pages='12733--12737',
                                 note='Delta H = -20.5 +/- 0.4 kJ/mol at 25 C for ATP4- + H2O = ADP3- + HPO4 2- + H+'),
    'traut1994physiological': dict(type='article', author='Traut, Thomas W.',
                                   title='Physiological concentrations of purines and pyrimidines',
                                   journal='Molecular and Cellular Biochemistry', year='1994', volume='140',
                                   number='1', pages='1--22', doi='10.1007/BF00928361',
                                   note='Mean ATP in mammalian cells 3152 +/- 1698 uM'),
    'mueckler2013slc': dict(type='article', author='Mueckler, Mike and Thorens, Bernard',
                             title='The SLC2 (GLUT) family of membrane transporters',
                             journal='Molecular Aspects of Medicine', year='2013', volume='34', number='2-3',
                             pages='121--138', doi='10.1016/j.mam.2012.07.001',
                             note='GLUT1 Km for glucose influx about 2 mM'),
    'halestrap1999proton': dict(type='article', author='Halestrap, Andrew P. and Price, Nigel T.',
                                title='The proton-linked monocarboxylate transporter (MCT) family: structure, function '
                                      'and regulation',
                                journal='Biochemical Journal', year='1999', volume='343', number='2', pages='281--299',
                                note='MCT1 Km for L-lactate 3-5 mM'),
    'celayaalcala2021simulation': dict(type='article',
                                       author='Celaya-Alcala, Jose T. and Lee, Grace V. and Smith, Amy F. and Li, Bohan '
                                              'and Sakad{\\v{z}}i{\\\'c}, Sava and Boas, David A. and Secomb, Timothy W.',
                                       title='Simulation of oxygen transport and estimation of tissue perfusion in '
                                             'extensive microvascular networks: Application to cerebral cortex',
                                       journal='Journal of Cerebral Blood Flow \\& Metabolism', year='2021', volume='41',
                                       number='3', pages='656--669', doi='10.1177/0271678X20927100',
                                       note='Table 1: Michaelis constant P0 = 10.5 mmHg, maximal O2 demand M0 = 7.5 cm3 '
                                            'O2/(100 cm3 min); mean tissue PO2 about 47 mmHg under reference conditions'),
    'fang2008oxygen': dict(type='article',
                           author='Fang, Qianqian and Sakad{\\v{z}}i{\\\'c}, Sava and Ruvinskaya, Lana and Devor, Anna '
                                  'and Dale, Anders M. and Boas, David A.',
                           title='Oxygen advection and diffusion in a three dimensional vascular anatomical network',
                           journal='Optics Express', year='2008', volume='16', number='22', pages='17530--17541',
                           doi='10.1364/OE.16.017530',
                           note='Bunsen solubility alpha = 1.27e-15 umol/(um3 mmHg) = 1.27e-3 mol/(m3 mmHg)'),
}
register_bib(BIB)



# ================================================================== energy_metabolism
O2_DEFAULT = 0.0597          # mM: 47 mmHg x 1.27 uM/mmHg, as the library's tissue O2 modules
O2_TARGET = 0.0558 * O2_DEFAULT / (O2_DEFAULT + 0.0133)   # mM/s: CelayaAlcala2021 demand at that PO2
EM = dict(T=310.15, T_ref=310.15, dT_Q10=10.0, Q10=2.0, V_c=V_CELL, C_O2=O2_DEFAULT, C_glc=1.0, C_lac=1.0,
          C_CO2=1.2, ATP_init=3.152, A_tot=3.5, K_O2=0.0133, K_glc=0.05, K_ADP=0.35, K_i_O2=0.0133, V_gly=0.05,
          K_ATP=0.15, n_ATP_ox=30.0, n_ATP_gly=2.0, h_O2_ox=470e3, h_lac_gly=63e3, h_ATP=20.5e3)


def em_rates(p, ATP):
    ADP = p['A_tot'] - ATP
    f_T = p['Q10'] ** ((p['T'] - p['T_ref']) / p['dT_Q10'])
    s_glc = p['C_glc'] / (p['K_glc'] + p['C_glc'])
    s_adp = ADP / (p['K_ADP'] + ADP)
    r_ox = p['V_ox'] * f_T * p['C_O2'] / (p['K_O2'] + p['C_O2']) * s_glc * s_adp
    r_gly = p['V_gly'] * f_T * s_glc * s_adp * p['K_i_O2'] / (p['K_i_O2'] + p['C_O2'])
    r_dem = p['k_dem'] * f_T * ATP / (p['K_ATP'] + ATP)
    return r_ox, r_gly, r_dem


# Calibrate V_ox and k_dem: at the default operating point the steady state has ATP = ATP_init (Traut 1994) and
# O2 consumption 6 r_ox = the library's tissue demand (CelayaAlcala2021 M0 C/(C + P0)).
ADP0 = EM['A_tot'] - EM['ATP_init']
s_glc0 = EM['C_glc'] / (EM['K_glc'] + EM['C_glc'])
s_adp0 = ADP0 / (EM['K_ADP'] + ADP0)
EM['V_ox'] = O2_TARGET / 6 / (EM['C_O2'] / (EM['K_O2'] + EM['C_O2']) * s_glc0 * s_adp0)
EM['k_dem'] = 1.0  # placeholder for em_rates
r_ox0, r_gly0, _ = em_rates(EM, EM['ATP_init'])
EM['k_dem'] = (EM['n_ATP_ox'] * r_ox0 + EM['n_ATP_gly'] * r_gly0) / (EM['ATP_init'] / (EM['K_ATP'] + EM['ATP_init']))
r_ox0, r_gly0, r_dem0 = em_rates(EM, EM['ATP_init'])
Q_MET0 = EM['V_c'] * (r_ox0 * (6 * EM['h_O2_ox']) + r_gly0 * (2 * EM['h_lac_gly']))
print(f'energy_metabolism: V_ox {EM["V_ox"]:.6g} mM/s, k_dem {EM["k_dem"]:.6g} mM/s; default O2 use '
      f'{6 * r_ox0:.6g} mM/s, glycolytic ATP share {EM["n_ATP_gly"] * r_gly0 / r_dem0:.3%}, heat {Q_MET0:.4g} W/cell '
      f'({Q_MET0 / EM["V_c"]:.4g} W/m3)')

em_vars = [
    ('t', 'second', 'in', None, None),
    ('T', 'kelvin', 'in', None, 'temperature of the volume (heat_source_port)'),
    ('T_ref', 'kelvin', 'in', None, None), ('dT_Q10', 'kelvin', 'in', None, None), ('Q10', 'dimensionless', 'in', None, None),
    ('V_c', 'm3', 'in', None, 'volume this metabolism occupies: a cell, part of a cell, or a tissue volume'),
    ('C_O2', 'millimolar', 'in', None, 'concentrations in the volume\'s pools (exchange ports)'),
    ('C_glc', 'millimolar', 'in', None, None), ('C_lac', 'millimolar', 'in', None, None), ('C_CO2', 'millimolar', 'in', None, None),
    ('ATP_init', 'millimolar', 'in', None, None),
    ('ATP', 'millimolar', 'out', 'ATP_init', 'adenine nucleotides: ATP + ADP = A_tot (AMP and Pi not modelled)'),
    ('A_tot', 'millimolar', 'in', None, None), ('ADP', 'millimolar', 'out', None, None),
    ('f_T', 'dimensionless', 'out', None, 'Q10 temperature factor on every rate'),
    ('V_ox', 'millimolar_per_s', 'in', None, None), ('K_O2', 'millimolar', 'in', None, None),
    ('K_glc', 'millimolar', 'in', None, None), ('K_ADP', 'millimolar', 'in', None, None),
    ('V_gly', 'millimolar_per_s', 'in', None, None), ('K_i_O2', 'millimolar', 'in', None, None),
    ('k_dem', 'millimolar_per_s', 'in', None, None), ('K_ATP', 'millimolar', 'in', None, None),
    ('n_ATP_ox', 'dimensionless', 'in', None, None), ('n_ATP_gly', 'dimensionless', 'in', None, None),
    ('h_O2_ox', 'J_per_mol', 'in', None, 'enthalpies released (positive = heat out), see notes'),
    ('h_lac_gly', 'J_per_mol', 'in', None, None), ('h_ATP', 'J_per_mol', 'in', None, None),
    ('r_ox', 'millimolar_per_s', 'out', None, 'glucose oxidised (aerobic, coupled to ATP synthesis)'),
    ('r_gly', 'millimolar_per_s', 'out', None, 'glucose to lactate (anaerobic glycolysis, coupled to ATP synthesis)'),
    ('r_dem', 'millimolar_per_s', 'out', None, 'ATP hydrolysed by the cell\'s work (maintenance)'),
    ('J_ATP', 'millimolar_per_s', 'out', None, 'ATP synthesis rate'),
    ('flux_O2', 'mol_per_s', 'out', None, 'sources into the pools (negative = taken up)'),
    ('flux_glc', 'mol_per_s', 'out', None, None), ('flux_lac', 'mol_per_s', 'out', None, None),
    ('flux_CO2', 'mol_per_s', 'out', None, None),
    ('Q_met', 'watt', 'out', None, 'heat released into the volume (heat_source_port)'),
]
em_eqs = [
    ('ADP', 'A_tot - ATP'),
    ('f_T', 'Q10 ** ((T - T_ref) / dT_Q10)'),
    ('r_ox', 'V_ox * f_T * (C_O2 / (K_O2 + C_O2)) * (C_glc / (K_glc + C_glc)) * (ADP / (K_ADP + ADP))'),
    ('r_gly', 'V_gly * f_T * (C_glc / (K_glc + C_glc)) * (ADP / (K_ADP + ADP)) * (K_i_O2 / (K_i_O2 + C_O2))'),
    ('r_dem', 'k_dem * f_T * (ATP / (K_ATP + ATP))'),
    ('J_ATP', 'n_ATP_ox * r_ox + n_ATP_gly * r_gly'),
    ('d/dt ATP', 'J_ATP - r_dem'),
    ('flux_O2', '-6.0 * r_ox * V_c'),
    ('flux_CO2', '6.0 * r_ox * V_c'),
    ('flux_glc', '-(r_ox + r_gly) * V_c'),
    ('flux_lac', '2.0 * r_gly * V_c'),
    ('Q_met', 'V_c * (r_ox * (6.0 * h_O2_ox - n_ATP_ox * h_ATP) + r_gly * (2.0 * h_lac_gly - n_ATP_gly * h_ATP) '
              '+ r_dem * h_ATP)'),
]
em_roles = {n: 'constant' for n, *_ in em_vars if n != 't'}
em_roles.update({n: 'boundary_condition' for n in ('T', 'C_O2', 'C_glc', 'C_lac', 'C_CO2')})
em_roles.update({n: 'variable' for n in ('ATP', 'ADP', 'f_T', 'r_ox', 'r_gly', 'r_dem', 'J_ATP', 'flux_O2', 'flux_glc',
                                          'flux_lac', 'flux_CO2', 'Q_met')})
em_roles.update({n: 'global_constant' for n in ('T_ref', 'dT_Q10')})
em_roles.pop('t', None)

PH = 'unreferenced_placeholder'
em_params = [
    ('T', EM['T'], 'operating point: 37 C', False),
    ('T_ref', EM['T_ref'], 'reference temperature of the rate constants (37 C)', False),
    ('dT_Q10', EM['dT_Q10'], 'definition of Q10 (rate ratio per 10 K)', True),
    ('Q10', EM['Q10'], f'{PH}; metabolic Q10 is typically 2-3; 2 chosen pending a source', False),
    ('V_c', V_CELL, 'geometry: 20 um diameter disc cell, 20 um deep (the game\'s starting cell)', False),
    ('C_O2', O2_DEFAULT, 'celayaalcala2021simulation; mean tissue PO2 47 mmHg x 1.27 uM/mmHg (fang2008oxygen)', True),
    ('C_glc', EM['C_glc'], f'{PH}; intracellular glucose of order 1 mM (below plasma ~5 mM)', False),
    ('C_lac', EM['C_lac'], f'{PH}; resting lactate of order 1 mM', False),
    ('C_CO2', EM['C_CO2'], f'{PH}; dissolved CO2 at 40 mmHg x 0.03 mM/mmHg (textbook solubility)', False),
    ('ATP_init', EM['ATP_init'], 'traut1994physiological; mean ATP in mammalian cells 3152 uM', True),
    ('A_tot', EM['A_tot'], f'{PH}; ATP + ADP with ATP/ADP about 9', False),
    ('K_O2', EM['K_O2'], 'celayaalcala2021simulation; P0 = 10.5 mmHg x 1.27 uM/mmHg, the library\'s tissue value', True),
    ('K_glc', EM['K_glc'], f'{PH}; hexokinase I affinity for glucose, order 0.03-0.1 mM', False),
    ('K_ADP', EM['K_ADP'], f'{PH}; apparent ADP affinity of respiration in tissue, order 0.3 mM', False),
    ('V_ox', EM['V_ox'], 'calibrated: at the default point the steady state has ATP = ATP_init and O2 use '
                         'M0 C/(C + P0) = 0.0456 mM/s (celayaalcala2021simulation)', False),
    ('V_gly', EM['V_gly'], f'{PH}; anaerobic capacity: 2 V_gly is about 40% of resting ATP demand', False),
    ('K_i_O2', EM['K_i_O2'], f'{PH}; Pasteur effect: glycolysis half-inhibited at K_O2', False),
    ('k_dem', EM['k_dem'], 'calibrated with V_ox (same conditions)', False),
    ('K_ATP', EM['K_ATP'], f'{PH}; ATPase affinity for ATP, order 0.1-0.5 mM', False),
    ('n_ATP_ox', EM['n_ATP_ox'], 'hinkle2005p; about 30 ATP per glucose oxidised', True),
    ('n_ATP_gly', EM['n_ATP_gly'], 'stoichiometry of glycolysis: 2 ATP per glucose to 2 lactate', True),
    ('h_O2_ox', EM['h_O2_ox'], 'gnaiger1990anaerobic; oxycaloric equivalent -470 kJ/mol O2 (glucose oxidation)', True),
    ('h_lac_gly', EM['h_lac_gly'], 'gnaiger1990anaerobic; glucose to lactate in bicarbonate buffer -63 kJ/mol lactate', True),
    ('h_ATP', EM['h_ATP'], 'gajewski1986thermodynamics; Delta H of ATP hydrolysis -20.5 kJ/mol at 25 C', True),
]

h_glc = 6 * EM['h_O2_ox']
h_lac = h_glc / 2 - EM['h_lac_gly']
em_invariants = [
    {'description': 'carbon is conserved: 6 flux_glc + 3 flux_lac + flux_CO2 = 0 (glucose C6, lactate C3) (exact)',
     'expr': 'np.allclose(6*flux_glc + 3*flux_lac + flux_CO2, 0, rtol=0, atol=1e-12*np.max(np.abs(flux_CO2)))'},
    {'description': 'respiratory quotient 1 for glucose: flux_CO2 = -flux_O2 (exact)',
     'expr': 'np.allclose(flux_CO2, -flux_O2, rtol=1e-12, atol=0)'},
    {'description': 'adenine nucleotides conserved: ATP + ADP = A_tot (exact)',
     'expr': 'np.allclose(ATP + ADP, A_tot, rtol=1e-12, atol=0)'},
    {'description': 'first law: the enthalpy taken in with glucose minus that leaving with lactate equals the heat '
                    'released plus the enthalpy stored in ATP: V h_ATP dATP + Q_met dt + (flux_glc h_glc + flux_lac h_lac) dt '
                    f'= 0 integrated over the run, with h_glc = 6 h_O2_ox = {h_glc:g} J/mol and h_lac = h_glc/2 - h_lac_gly '
                    f'= {h_lac:g} J/mol, relative to CO2, H2O and ADP + Pi (relative 1e-6 of the heat released)',
     'expr': 'abs(V_c*h_ATP*(ATP[-1] - ATP[0]) + np.trapz(Q_met + flux_glc*6*h_O2_ox + flux_lac*(3*h_O2_ox - '
             'h_lac_gly), t)) <= 1e-6*np.trapz(np.abs(Q_met), t)'},
    {'description': 'second law: every reaction runs forward (irreversible, far from equilibrium), r_ox, r_gly, '
                    'r_dem >= 0',
     'expr': 'np.all(r_ox >= 0) and np.all(r_gly >= 0) and np.all(r_dem >= 0)'},
    {'description': f'from depleted ATP (run_parameters ATP_init = 2 mM) ATP recovers to the calibrated steady state '
                    f'{EM["ATP_init"]} mM by the end of the run (relative 1e-6)',
     'expr': f'abs(ATP[-1] - {EM["ATP_init"]!r}) <= 1e-6*{EM["ATP_init"]!r}', 'applies': 'run'},
]
em_notes = (
    'For cells with mitochondria (a later stage of AdaptOrDeity; the starting fermenter uses adaptive_ArgusUNPUBLISHED_v01). '
    'Energy metabolism of a well-mixed volume (a cell, part of a cell or a tissue volume): the local kinetics only. '
    'It reads the concentrations of the volume\'s O2, glucose, lactate and CO2 pools and its temperature, holds ATP '
    '(ATP + ADP = A_tot), and returns sources into the pools (mol/s, negative = uptake) and the heat it releases (W). '
    'Transport between volumes is done by face modules (diffusion faces, membrane carriers, heat faces), so the same '
    'version serves every level of detail. Three lumped reactions, each with Michaelis-Menten dependences and a Q10 '
    'temperature factor f_T: glucose oxidation r_ox (6 O2 in, 6 CO2 out, n_ATP_ox ATP made; respiratory control by '
    'ADP), anaerobic glycolysis r_gly (glucose to 2 lactate, 2 ATP; inhibited by O2: the Pasteur effect) and ATP '
    'demand r_dem. Heat is the enthalpy each reaction releases (first law; the free energy that only becomes entropy '
    'is not heat): oxidation 6 h_O2_ox - n_ATP_ox h_ATP, glycolysis 2 h_lac_gly - n_ATP_gly h_ATP, hydrolysis h_ATP per '
    'mol, so at steady ATP the heat is 470 kJ per mol O2 used plus 63 kJ per mol lactate made (Gnaiger and Kemp 1990). '
    'V_ox and k_dem are calibrated so the default point is a steady state at ATP = 3.152 mM using O2 at the library\'s '
    'tissue rate. Ports: general O2/glucose/lactate/CO2 exchange ports (C, flux) to the pools, and heat_source_port '
    '(T, Q_met) to a thermal_volume. Not modelled yet: pH/bicarbonate, AMP and adenylate kinase, '
    'phosphocreatine, the temperature dependence of the enthalpies (Delta Cp).')

write_version(
    'cell', 'energy_metabolism', 'energy_metabolism_mitochondrial', em_vars, em_eqs,
    {'required_citations': ['argus2026unpublished'],
     'entrance_ports': [], 'exit_ports': [],
     'general_ports': [{'port_type': 'O2_exchange_port', 'variables': ['C_O2', 'flux_O2']},
                       {'port_type': 'glucose_exchange_port', 'variables': ['C_glc', 'flux_glc']},
                       {'port_type': 'lactate_exchange_port', 'variables': ['C_lac', 'flux_lac']},
                       {'port_type': 'CO2_exchange_port', 'variables': ['C_CO2', 'flux_CO2']},
                       {'port_type': 'heat_source_port', 'variables': ['T', 'Q_met']}],
     'roles': em_roles},
    [(n, EM[n] if n in EM else v, r, s) for n, v, r, s in em_params],
    {'sim_time': 60.0, 'pre_time': 0.0, 'dt': 0.01, **SOLVER,
     'outputs': ['ATP', 'ADP', 'f_T', 'r_ox', 'r_gly', 'r_dem', 'J_ATP', 'flux_O2', 'flux_glc', 'flux_lac', 'flux_CO2',
                 'Q_met'],
     'run_parameters': {'ATP_init': 2.0},
     'invariants': em_invariants,
     'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
    {'reviewed': False, 'notes': em_notes},
    ['argus2026unpublished', 'gnaiger1990anaerobic', 'hinkle2005p', 'gajewski1986thermodynamics',
     'traut1994physiological', 'celayaalcala2021simulation', 'fang2008oxygen'],
    'mitochondrial_ArgusUNPUBLISHED_v01')

# ================================================================== thermal_volume
th_vars = [
    ('t', 'second', 'in', None, None),
    ('T_init', 'kelvin', 'in', None, None),
    ('T', 'kelvin', 'out', 'T_init', 'temperature of the volume, well mixed'),
    ('rho', 'kg_per_m3', 'in', None, None), ('c_p', 'J_per_kg_K', 'in', None, None),
    ('V_th', 'm3', 'in', None, None),
    ('Q_in', 'watt', 'in', None, 'heat in through the faces on the entrance and exit heat_face_ports (summed)'),
    ('Q_out', 'watt', 'in', None, None),
    ('Q_src', 'watt', 'in', None, 'heat released inside the volume (heat_source_port, e.g. energy_metabolism Q_met)'),
    ('C_th', 'J_per_K', 'out', None, 'heat capacity'),
]
UNITS['J_per_K'] = [('joule', None, 1, 1), ('kelvin', None, -1, 1)]
th_eqs = [('C_th', 'rho * c_p * V_th'), ('d/dt T', '(Q_in + Q_out + Q_src) / C_th')]
th_roles = {'T_init': 'constant', 'T': 'variable', 'rho': 'constant', 'c_p': 'constant', 'V_th': 'constant',
            'Q_in': 'variable', 'Q_out': 'variable', 'Q_src': 'boundary_condition', 'C_th': 'variable'}
th_params = [
    ('T_init', 310.15, 'operating point: 37 C', False),
    ('rho', 1050.0, f'{PH}; soft tissue density of order 1050 kg/m3 (to check against the IT\'IS tissue database)', False),
    ('c_p', 3600.0, f'{PH}; soft tissue specific heat of order 3600 J/(kg K) (to check against the IT\'IS tissue '
                    'database)', False),
    ('V_th', V_CELL, 'geometry: 20 um diameter disc cell, 20 um deep', False),
    ('Q_src', Q_MET0, 'energy_metabolism Q_met at its default operating point', False),
]
write_version(
    'transport', 'thermal_volume', 'thermal_volume', th_vars, th_eqs,
    {'required_citations': ['argus2026unpublished'],
     'entrance_ports': [{'port_type': 'heat_face_port', 'variables': ['T', 'Q_in'], 'multi_port': ['True', 'sum']}],
     'exit_ports': [{'port_type': 'heat_face_port', 'variables': ['T', 'Q_out'], 'multi_port': ['True', 'sum']}],
     'general_ports': [{'port_type': 'heat_source_port', 'variables': ['T', 'Q_src']}],
     'roles': th_roles},
    th_params,
    {'sim_time': 10.0, 'pre_time': 0.0, 'dt': 0.01, **SOLVER, 'outputs': ['T', 'C_th'], 'run_parameters': {},
     'invariants': [{'description': 'with no faces, all the source heat is stored: C_th (T - T_init) = Q_src t (exact '
                                    'to 1e-9 relative)',
                     'expr': 'np.allclose(C_th*(T - T_init), Q_src*t, rtol=1e-9, atol=1e-30)'}],
     'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
    {'reviewed': False, 'notes': (
        'The thermal counterpart of tissue_diffusion_volume: one well-mixed finite volume of heat capacity '
        'C_th = rho c_p V_th, dT/dt = (Q_in + Q_out + Q_src)/C_th. Q_in and Q_out are the summed heat flows of the '
        'heat_conduction_faces on its entrance and exit heat_face_ports (0 when none); Q_src is the heat released '
        'inside it (energy_metabolism Q_met through the heat_source_port). Only in-plane faces carry heat in '
        'AdaptOrDeity; the top and bottom of a 2D sheet are insulated.')},
    ['argus2026unpublished'], VERSION)

# ================================================================== heat_conduction_face
hf_vars = [
    ('t', 'second', 'in', None, None),
    ('k_th', 'W_per_m_K', 'in', None, None), ('A_f', 'm2', 'in', None, None),
    ('d_up', 'metre', 'in', None, None), ('d_down', 'metre', 'in', None, None),
    ('open_flag', 'dimensionless', 'in', None, None),
    ('T_up', 'kelvin', 'in', None, None), ('T_down', 'kelvin', 'in', None, None),
    ('G_th', 'W_per_K', 'out', None, 'thermal conductance'),
    ('Q', 'watt', 'out', None, 'heat flow from the upstream to the downstream volume'),
    ('Q_up', 'watt', 'out', None, None), ('Q_down', 'watt', 'out', None, None),
]
hf_eqs = [('G_th', 'open_flag * k_th * A_f / (d_up + d_down)'), ('Q', 'G_th * (T_up - T_down)'),
          ('Q_up', '-Q'), ('Q_down', 'Q')]
hf_roles = {'k_th': 'constant', 'A_f': 'constant', 'd_up': 'constant', 'd_down': 'constant', 'open_flag': 'constant',
            'T_up': 'boundary_condition', 'T_down': 'boundary_condition', 'G_th': 'variable', 'Q': 'variable',
            'Q_up': 'variable', 'Q_down': 'variable'}
write_version(
    'transport', 'heat_conduction_face', 'heat_conduction_face', hf_vars, hf_eqs,
    {'required_citations': ['argus2026unpublished'],
     'entrance_ports': [{'port_type': 'heat_face_port', 'variables': ['T_up', 'Q_up']}],
     'exit_ports': [{'port_type': 'heat_face_port', 'variables': ['T_down', 'Q_down']}],
     'general_ports': [], 'roles': hf_roles},
    [('k_th', 0.5, f'{PH}; soft tissue thermal conductivity of order 0.5 W/(m K) (water 0.6)', False),
     ('A_f', 4e-10, 'geometry: 20 um x 20 um face', False),
     ('d_up', 1e-5, 'geometry: centre to face, half a 20 um cell', False),
     ('d_down', 1e-5, 'geometry: centre to face, half a 20 um cell', False),
     ('open_flag', 1.0, 'open face', False),
     ('T_up', 310.15, 'operating point', False), ('T_down', 310.0, 'operating point: 0.15 K cooler', False)],
    {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01, **SOLVER, 'outputs': ['G_th', 'Q', 'Q_up', 'Q_down'],
     'run_parameters': {},
     'invariants': [{'description': 'Fourier conduction: Q = open_flag k_th A_f/(d_up + d_down) (T_up - T_down), '
                                    'Q_down = Q, Q_up = -Q (energy conserved across the face) (exact)',
                     'expr': 'np.allclose(Q, open_flag*k_th*A_f/(d_up + d_down)*(T_up - T_down), rtol=1e-12, atol=0) '
                             'and np.array_equal(Q_down, Q) and np.array_equal(Q_up, -Q)'},
                    {'description': 'second law: heat flows from hot to cold, Q (T_up - T_down) >= 0',
                     'expr': 'np.all(Q*(T_up - T_down) >= 0)'}],
     'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
    {'reviewed': False, 'notes': (
        'Heat conduction between two thermal volumes (Fourier), the thermal counterpart of tissue_diffusion_face: '
        'G_th = open_flag k_th A_f/(d_up + d_down), Q = G_th (T_up - T_down) into the downstream volume (Q_down = Q) '
        'and -Q into the upstream one (Q_up = -Q). Conduction moves heat without creating any; its entropy '
        'production Q (1/T_down - 1/T_up) is >= 0. With one side unconnected, that side\'s temperature is a fixed '
        'boundary value.')},
    ['argus2026unpublished'], VERSION)

# ================================================================== membrane_carrier_face
mc_vars = [
    ('t', 'second', 'in', None, None),
    ('V_max', 'mol_per_s', 'in', None, 'maximal flux of the carriers in this membrane area'),
    ('K_m', 'millimolar', 'in', None, None), ('open_flag', 'dimensionless', 'in', None, None),
    ('C_up', 'millimolar', 'in', None, 'upstream side (by convention extracellular)'),
    ('C_down', 'millimolar', 'in', None, 'downstream side (by convention intracellular)'),
    ('J', 'mol_per_s', 'out', None, 'net flux from upstream to downstream'),
    ('J_up', 'mol_per_s', 'out', None, None), ('J_down', 'mol_per_s', 'out', None, None),
]
mc_eqs = [('J', 'open_flag * V_max * (C_up / (K_m + C_up) - C_down / (K_m + C_down))'),
          ('J_up', '-J'), ('J_down', 'J')]
mc_roles = {'V_max': 'constant', 'K_m': 'constant', 'open_flag': 'constant', 'C_up': 'boundary_condition',
            'C_down': 'boundary_condition', 'J': 'variable', 'J_up': 'variable', 'J_down': 'variable'}
glut = [('V_max', 0.1 * V_CELL, f'{PH}; 0.1 mM/s of the 20 um cell\'s volume: more than its maximal glycolytic '
                                'flux; transporter density is a game trait', False),
        ('K_m', 2.0, 'mueckler2013slc; GLUT1 Km for glucose influx about 2 mM', True),
        ('open_flag', 1.0, 'open', False),
        ('C_up', 5.0, f'{PH}; extracellular glucose of order plasma 5 mM', False),
        ('C_down', 1.0, f'{PH}; intracellular glucose of order 1 mM', False)]
mct = [('V_max', 0.2 * V_CELL, f'{PH}; 0.2 mM/s of the 20 um cell\'s volume', False),
       ('K_m', 4.0, 'halestrap1999proton; MCT1 Km for L-lactate 3-5 mM (midpoint)', True),
       ('open_flag', 1.0, 'open', False),
       ('C_up', 1.0, f'{PH}; extracellular lactate of order 1 mM', False),
       ('C_down', 2.0, f'{PH}; intracellular lactate raised by glycolysis', False)]
write_version(
    'transport', 'membrane_carrier_face', 'membrane_carrier_face', mc_vars, mc_eqs,
    {'required_citations': ['argus2026unpublished'],
     'entrance_ports': [{'port_type': 'diffusion_face_port', 'variables': ['C_up', 'J_up']}],
     'exit_ports': [{'port_type': 'diffusion_face_port', 'variables': ['C_down', 'J_down']}],
     'general_ports': [], 'roles': mc_roles},
    glut,
    {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01, **SOLVER, 'outputs': ['J', 'J_up', 'J_down'], 'run_parameters': {},
     'invariants': [{'description': 'symmetric carrier: J = open_flag V_max (C_up/(K_m + C_up) - C_down/(K_m + C_down)), '
                                    'J_down = J, J_up = -J (mass conserved across the membrane) (exact)',
                     'expr': 'np.allclose(J, open_flag*V_max*(C_up/(K_m + C_up) - C_down/(K_m + C_down)), rtol=1e-12, '
                             'atol=0) and np.array_equal(J_down, J) and np.array_equal(J_up, -J)'},
                    {'description': 'second law: facilitated transport runs down its gradient, J ln(C_up/C_down) >= 0',
                     'expr': 'np.all(J*np.log(C_up/C_down) >= 0)'}],
     'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
    {'reviewed': False, 'notes': (
        'A membrane with facilitated-diffusion carriers between two volumes (by convention extracellular upstream, '
        'intracellular downstream): a symmetric Michaelis-Menten carrier, J = open_flag V_max (C_up/(K_m + C_up) - '
        'C_down/(K_m + C_down)), which saturates and always runs down the gradient (no energy input, no heat). Uses '
        'diffusion_face_ports, so it joins the same pools as tissue_diffusion_face. Instances: default (GLUT1, '
        'glucose, K_m 2 mM) and MCT1_lactate (K_m 4 mM; its proton coupling is not modelled).')},
    ['argus2026unpublished', 'mueckler2013slc', 'halestrap1999proton'], VERSION,
    extra_instances={'MCT1_lactate': mct})
print('written')
