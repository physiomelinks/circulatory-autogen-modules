"""
Decomposes the sympathetic neuron's SN_soma and SN_varicosity (modules/cell/neurons) into
reusable submodules, and builds the models and supermodules that use them.

Equations are not retyped: each new component copies SN_soma's or SN_varicosity's MathML
verbatim (lxml deep copies). Variable names are kept as they are in the original, so the
outputs of a decomposed model are the original's, found in another component. The only
changes to the math are where a sum of channel terms becomes one summed port variable:

  SN_membrane_soma             dV/dt:  i_Na + ... + i_leak_K   -> i_ion        (membrane_current, "sum")
  SN_Na_K_concentrations_soma  dNai/dt: j_Na + ... + j_leak_Na -> j_Na_total   (Na_flux, "sum")
                               dKi/dt:  j_Kv1_5 + ... + j_leak_K -> j_K_total  (K_flux, "sum")
  SN_Ca_handling_soma/_varicosity
                               dCai/dt: j_Ca_L|j_Ca_N + j_NaCa_Ca + j_IP3 + j_RyR -> j_Ca_in  (Ca_flux, "sum")
                               dCa_ER/dt: -(j_IP3) - j_RyR   -> -(j_ER_release)       (ER_Ca_release, "sum")

The channels' ion fluxes (j_Na = -i_Na/F, ...) move from the concentration equations into the
channel that carries the current, so a concentration pool sums whatever is connected to it.
With the neighbours listed in the original order, each sum adds the same terms in the same
order as the original, so the decomposed model reproduces SN_simple to round-off.

Writes (and rewrites on every run):
  modules/ion_channel/          the SN ion channels (<name>_SN): CellML, config, parameters,
                                units, and a spec entry for each new component (kept if present)
  modules/cell/neurons/         membrane, concentrations, Ca handling and NE release components
  modules/system/cellular/SN_simple_flat/          the decomposed SN_simple, written out
  modules/supermodules/SN_soma, SN_varicosity, sympathetic_neuron/   the supermodules
  modules/system/cellular/SN_simple_supermodules/  SN_simple as one sympathetic_neuron instance

    python tools/build_sn_modules.py
"""
import copy
import csv
import io
import json
import math
import os
import re
import shutil
import sys

import yaml
from lxml import etree

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from cam_testing import vessel_array  # noqa: E402

MODULES = os.path.join(REPO, 'modules')
NEURONS = os.path.join(MODULES, 'cell', 'neurons')
ION = os.path.join(MODULES, 'ion_channel')
SOURCE = os.path.join(NEURONS, 'neurons_modules.cellml')
SN_SIMPLE = os.path.join(MODULES, 'system', 'cellular', 'SN_simple')

C = 'http://www.cellml.org/cellml/1.1#'
M = 'http://www.w3.org/1998/Math/MathML'
NS = {'c': C, 'm': M}
GLOBALS = ('R', 'F', 'T')
LOCAL = ('epsilon_tau',)          # private constants of the original (initial_value, no interface)

# ---------------------------------------------------------------------------------------------
# ports
# ---------------------------------------------------------------------------------------------


def port(port_type, variables, multi=None):
    p = {'port_type': port_type, 'variables': list(variables)}
    if multi:
        p['multi_port'] = multi
    return p


V = port('membrane_voltage', ['V'])
I_OUT = lambda i: port('membrane_current', [i])          # noqa: E731  summed into the membrane
E_NA, E_K, E_CA = port('Na_Nernst_potential', ['E_Na']), port('E_K_potential', ['E_K']), port('Ca_potential', ['E_Ca'])
NA_FLUX = lambda j: port('Na_flux', [j])                 # noqa: E731
K_FLUX = lambda j: port('K_flux', [j])                   # noqa: E731
CA_FLUX = lambda j: port('Ca_flux', [j])                 # noqa: E731
ER_RELEASE = lambda j: port('ER_Ca_release', [j])        # noqa: E731
NAI, CAI, CA_BULK, CA_ER = (port('internal_Na', ['Nai']), port('internal_Ca', ['Cai']),
                            port('bulk_Ca', ['Ca_bulk']), port('ER_Ca', ['Ca_ER']))


def nernst(ion_out, ion_in, z=1):
    '''RT/F ln(out/in) in mV at the SN_simple constants (R 8.314, T 310.15, F 96485).'''
    return 1000 * 8.314 * 310.15 / 96485 / z * math.log(ion_out / ion_in)


# standalone values of port inputs (used when a component runs alone, e.g. in its module tests)
REST = (-70, 'definitional; clamp value when unconnected: the SN resting potential (SN_simple V_rest_approx = -70 mV)')
BC = {
    'V': REST, 'V_in': REST,
    'E_Na': (round(nernst(145, 13), 6), 'definitional; Nernst potential at SN_simple Nao 145, Nai_init 13 mM (used when unconnected)'),
    'E_K': (round(nernst(4.5, 150), 6), 'definitional; Nernst potential at SN_simple Ko 4.5, Ki_init 150 mM (used when unconnected)'),
    'E_Ca': (round(nernst(2, 1e-4), 6), 'definitional; SN_varicosity E_Ca = RT/F ln(Cao/Cai) at Cao 2 mM, Cai 1e-4 mM (used when unconnected)'),
    'Cai': (1e-4, 'definitional; SN_simple Cai_init (Hagiwara1981; used when unconnected)'),
    'Ca_bulk': (1e-4, 'definitional; SN_simple Cai_init, the initial Ca_bulk (used when unconnected)'),
    'Ca_ER': (0.6, 'definitional; SN_simple Ca_ER_init (used when unconnected)'),
    'Nai': (13, 'SN_simple: BoronBoulpaep2017_was_15 (SN_simple Nai_init / SN_varicosity Nai)'),
    'I_out': (0.0, 'definitional; no axial current when unconnected'),
    'I_in': (0.0, 'SN_simple: user_defined'),
}

# ---------------------------------------------------------------------------------------------
# components: (name, module, source component, equations by left-hand side, port inputs,
#              summed inputs {var: units}, ports, notes)
# ---------------------------------------------------------------------------------------------

CHANNEL_NOTE = ('Copied verbatim from the owner\'s SN_soma (modules/cell/neurons, the SN_simple sympathetic neuron); '
                'every gating time constant has epsilon_tau = 0.1 ms added, as in the original. ')


_N_CA = ('Vol*(c_SL*(Cai + CaB_Cai) + (1 - c_SL - c_ER - c_ER_pre)*(Ca_bulk + CaB_Ca_bulk) '
         '+ c_ER_pre*(Ca_ER_pre + CaB_Ca_ER_pre) + c_ER*Ca_ER)')
CA_TOTAL = {'description': 'alone (no channel fluxes connected) the total Ca, free + buffered in all four compartments, '
                           'is conserved (relative 1e-6)',
            'expr': f'np.max(np.abs({_N_CA} - ({_N_CA})[0])) <= 1e-6*({_N_CA})[0]'}


def gate_eqs(*names):
    return list(names)


COMPONENTS = [
    # ---- ion channels (modules/ion_channel) ---------------------------------------------------
    dict(name='i_Na_SN', module='ion_channel', source='SN_soma',
         eqs=['m_init', 'h_init', 'm_inf', 'h_inf', 'tau_m', 'tau_h', 'd/dt m', 'd/dt h', 'i_Na', 'j_Na'],
         bc=['V', 'E_Na'], general=[V, E_NA, I_OUT('i_Na'), NA_FLUX('j_Na')],
         states=['m', 'h'], current='i_Na', flux=[('j_Na', '-i_Na/F')],
         notes='Fast (NaV1.7-like) Na current of the SN soma, i_Na = g_Na g_Na_mod (V - E_Na) m^3 h, with its '
               'Na flux j_Na = -i_Na/F into the Na pool. Gates start at their steady state at V_rest_approx.'),
    dict(name='i_Na1_6_SN', module='ion_channel', source='SN_soma',
         eqs=['m_init', 'h_init', 'm_Na1_6_inf', 'h_Na1_6_inf', 'tau_m_Na1_6', 'tau_h_Na1_6', 'd/dt m_Na1_6',
              'd/dt h_Na1_6', 'i_Na1_6', 'j_Na1_6'],
         bc=['V', 'E_Na'], general=[V, E_NA, I_OUT('i_Na1_6')],
         states=['m_Na1_6', 'h_Na1_6'], current='i_Na1_6', flux=[('j_Na1_6', '-i_Na1_6/F')],
         notes='NaV1.6 current of the SN soma, i_Na1_6 = g_Na1_6 g_Na1_6_mod (V - E_Na) m_Na1_6^3 h_Na1_6 (the same '
               'law as i_Na_SN). Two features of the original are kept so that SN_simple is reproduced: the gates '
               'start at i_Na_SN\'s steady state (the original sets m_Na1_6\'s initial value to m_init, computed '
               'from V_mid_m, V_den_m, V_mid_h, V_den_h, the NaV1.7 values, which are this component\'s own '
               'parameters here), and its Na flux j_Na1_6 = -i_Na1_6/F is computed but has no port, as the original '
               'leaves it out of dNai/dt.'),
    dict(name='i_CaL_SN', module='ion_channel', source='SN_soma',
         eqs=['V_mid_hc', 'V_den_hc', 'c_init', 'hc_init', 'c_inf', 'hc_inf', 'tau_c', 'd/dt c', 'd/dt hc',
              'i_Ca_L', 'j_Ca_L'],
         bc=['V', 'Cai'], general=[V, CAI, I_OUT('i_Ca_L'), CA_FLUX('j_Ca_L')],
         states=['c', 'hc'], current='i_Ca_L', flux=[('j_Ca_L', '-i_Ca_L/(2*F)')],
         notes='L-type Ca current of the SN soma: GHK flux with permeability p_Ca (1 + cAMP), gates c^2 hc; '
               'Ca flux j_Ca_L = -i_Ca_L/(2F) into the sub-membrane Ca.'),
    dict(name='i_Kv2_SN', module='ion_channel', source='SN_soma',
         eqs=['m_Kv2_init', 'm_Kv2_inf', 'tau_m_Kv2', 'd/dt m_Kv2', 'i_Kv2', 'j_Kv2'],
         bc=['V', 'E_K'], general=[V, E_K, I_OUT('i_Kv2'), K_FLUX('j_Kv2')],
         states=['m_Kv2'], current='i_Kv2', flux=[('j_Kv2', '-i_Kv2/F')],
         notes='Kv2 delayed-rectifier K current of the SN soma, i_Kv2 = g_Kv2 (V - E_K) m_Kv2.'),
    dict(name='i_Kv1_5_SN', module='ion_channel', source='SN_soma',
         eqs=['m_Kv1_5_init', 'h_Kv1_5_init', 'm_Kv1_5_inf', 'h_Kv1_5_inf', 'tau_m_Kv1_5', 'tau_h_Kv1_5',
              'd/dt m_Kv1_5', 'd/dt h_Kv1_5', 'i_Kv1_5', 'j_Kv1_5'],
         bc=['V', 'E_K'], general=[V, E_K, I_OUT('i_Kv1_5'), K_FLUX('j_Kv1_5')],
         states=['m_Kv1_5', 'h_Kv1_5'], current='i_Kv1_5', flux=[('j_Kv1_5', '-i_Kv1_5/F')],
         notes='Kv1.5 K current of the SN soma, i_Kv1_5 = g_Kv1_5 m h (V - E_K) (Channelpedia Kv1.5 kinetics).'),
    dict(name='i_Kv4_2_SN', module='ion_channel', source='SN_soma',
         eqs=['m_Kv4_2_init', 'h_Kv4_2_init', 'm_Kv4_2_inf', 'h_Kv4_2_inf', 'tau_m_Kv4_2', 'd/dt m_Kv4_2',
              'd/dt h_Kv4_2', 'i_Kv4_2', 'j_Kv4_2'],
         bc=['V', 'E_K'], general=[V, E_K, I_OUT('i_Kv4_2'), K_FLUX('j_Kv4_2')],
         states=['m_Kv4_2', 'h_Kv4_2'], current='i_Kv4_2', flux=[('j_Kv4_2', '-i_Kv4_2/F')],
         notes='Kv4.2 A-type K current of the SN soma, i_Kv4_2 = g_Kv4_2 m^3 h (V - E_K) (Channelpedia Kv4.2 kinetics).'),
    dict(name='i_BK_SN', module='ion_channel', source='SN_soma',
         eqs=['V_mid_m_BK', 'm_BK_init', 'h_BK_init', 'z_BK_init', 'm_BK_inf', 'h_BK_inf', 'z_BK_inf', 'tau_m_BK',
              'tau_h_BK', 'd/dt m_BK', 'd/dt h_BK', 'd/dt z_BK', 'i_KCa_BK', 'j_KCa_BK'],
         bc=['V', 'E_K', 'Cai'], general=[V, E_K, CAI, I_OUT('i_KCa_BK'), K_FLUX('j_KCa_BK')],
         states=['m_BK', 'h_BK', 'z_BK'], current='i_KCa_BK', flux=[('j_KCa_BK', '-i_KCa_BK/F')],
         notes='BK (large-conductance Ca-activated K) current of the SN soma, i_KCa_BK = g_KCa_BK m^3 h z (V - E_K), '
               'z driven by the sub-membrane Cai (Khaliq 2003 kinetics). z starts at its steady state at Cai_init.'),
    dict(name='i_SK_SN', module='ion_channel', source='SN_soma',
         eqs=['m_SK_init', 'm_SK_inf', 'd/dt m_SK', 'i_KCa_SK', 'j_KCa_SK'],
         bc=['V', 'E_K', 'Ca_bulk'], general=[V, E_K, CA_BULK, I_OUT('i_KCa_SK'), K_FLUX('j_KCa_SK')],
         states=['m_SK'], current='i_KCa_SK', flux=[('j_KCa_SK', '-i_KCa_SK/F')],
         notes='SK (small-conductance Ca-activated K) current of the SN soma, i_KCa_SK = g_KCa_SK g_KCa_SK_mod m_SK '
               '(V - E_K), m_SK a Hill function of the bulk cytosolic Ca (Ca_bulk, as in the original). m_SK starts '
               'at its steady state at Cai_init.'),
    dict(name='i_M_SN', module='ion_channel', source='SN_soma',
         eqs=['w_init', 'w_inf', 'tau_w', 'd/dt w', 'i_M', 'j_KM'],
         bc=['V', 'E_K'], general=[V, E_K, I_OUT('i_M'), K_FLUX('j_KM')],
         states=['w'], current='i_M', flux=[('j_KM', '-i_M/F')],
         notes='M-type (KCNQ) K current of the SN soma, i_M = g_M_mod g_M w (V - E_K); Delta_V_mid_w shifts '
               'the activation for drug studies.'),
    dict(name='i_leak_Na_SN', module='ion_channel', source='SN_soma',
         eqs=['i_leak_Na', 'j_leak_Na'],
         bc=['V', 'E_Na'], general=[V, E_NA, I_OUT('i_leak_Na'), NA_FLUX('j_leak_Na')],
         states=[], current='i_leak_Na', flux=[('j_leak_Na', '-i_leak_Na/F')],
         notes='Na leak of the SN soma, i_leak_Na = g_leak_Na (V - E_Na).'),
    dict(name='i_leak_K_SN', module='ion_channel', source='SN_soma',
         eqs=['i_leak_K', 'j_leak_K'],
         bc=['V', 'E_K'], general=[V, E_K, I_OUT('i_leak_K'), K_FLUX('j_leak_K')],
         states=[], current='i_leak_K', flux=[('j_leak_K', '-i_leak_K/F')],
         notes='K leak of the SN soma, i_leak_K = g_leak_K (V - E_K).'),
    dict(name='i_NaCa_SN', module='ion_channel', source='SN_soma', same_in='SN_varicosity',
         eqs=['i_NaCa', 'j_NaCa_Ca', 'j_NaCa_Na'],
         bc=['V', 'Nai', 'Cai'],
         general=[V, NAI, CAI, I_OUT('i_NaCa'), NA_FLUX('j_NaCa_Na'), CA_FLUX('j_NaCa_Ca')],
         states=[], current='i_NaCa', flux=[('j_NaCa_Ca', 'i_NaCa/F'), ('j_NaCa_Na', '-3*i_NaCa/F')],
         notes='Na/Ca exchanger of the SN soma and varicosity (the same law in both; Paci 2013 form), with its Ca '
               'flux j_NaCa_Ca = i_NaCa/F and Na flux j_NaCa_Na = -3 i_NaCa/F. In the varicosity Nai is a constant '
               '(no Na pool is connected).'),
    dict(name='i_NaK_SN', module='ion_channel', source='SN_soma',
         eqs=['i_NaK', 'j_NaK_Na', 'j_NaK_K'],
         bc=['V', 'Nai'], general=[V, NAI, I_OUT('i_NaK'), NA_FLUX('j_NaK_Na'), K_FLUX('j_NaK_K')],
         states=[], current='i_NaK', flux=[('j_NaK_Na', '-3*i_NaK/F'), ('j_NaK_K', '2*i_NaK/F')],
         notes='Na/K pump of the SN soma (Lindblad 1996 form), with its fluxes j_NaK_Na = -3 i_NaK/F and '
               'j_NaK_K = 2 i_NaK/F.'),
    dict(name='IP3R_SN', module='ion_channel', source='SN_soma', same_in='SN_varicosity',
         eqs=['m_IP3_inf', 'Q_2', 'tau_hIP3', 'h_IP3_inf', 'd/dt h_IP3', 'j_IP3'],
         bc=['Cai', 'Ca_ER'], general=[CAI, CA_ER, CA_FLUX('j_IP3'), ER_RELEASE('j_IP3')],
         states=['h_IP3'], current=None, flux=[],
         notes='IP3 receptor Ca release from the ER (Li & Rinzel 1994 reduced form) of the SN soma and varicosity '
               '(the same law in both): j_IP3 = -k_IP3 m_inf^3 h^3 (Cai - Ca_ER), into the sub-membrane Ca and out '
               'of the ER.'),
    dict(name='RyR_SN', module='ion_channel', source='SN_soma', same_in='SN_varicosity',
         eqs=['c_RyR_inf', 'hc_RyR_inf', 'd/dt c_RyR', 'd/dt hc_RyR', 'j_RyR'],
         bc=['Cai', 'Ca_ER'], general=[CAI, CA_ER, CA_FLUX('j_RyR'), ER_RELEASE('j_RyR')],
         states=['c_RyR', 'hc_RyR'], current=None, flux=[],
         invariants=[{'description': 'j_RyR = k_RyR (Ca_ER - Cai) c_RyR^3 hc_RyR^3',
                      'expr': 'np.max(np.abs(j_RyR - k_RyR*(Ca_ER - Cai)*c_RyR**3*hc_RyR**3)) <= '
                              '1e-12*max(1e-30, np.max(np.abs(j_RyR)))'}],
         notes='Ryanodine receptor Ca release from the ER of the SN soma and varicosity (the same law in both): '
               'j_RyR = k_RyR (Ca_ER - Cai) c^3 hc^3, c activated by Cai and hc by Ca_ER.'),
    dict(name='i_CaN_SN', module='ion_channel', source='SN_varicosity',
         eqs=['V_mid_hc', 'V_den_hc', 'c_inf', 'hc_inf', 'tau_c', 'd/dt c', 'd/dt hc', 'i_Ca_N', 'j_Ca_N'],
         bc=['V', 'E_Ca'], general=[V, E_CA, I_OUT('i_Ca_N'), CA_FLUX('j_Ca_N')],
         states=['c', 'hc'], current='i_Ca_N', flux=[('j_Ca_N', '-i_Ca_N/(2*F)')],
         notes='N-type Ca current of the SN varicosity, i_Ca_N = g_CaN (V - E_Ca) c^2 hc (the soma\'s L-type '
               'gating with a linear driving force), with its Ca flux j_Ca_N = -i_Ca_N/(2F).'),
    # ---- neuron-specific parts (modules/cell/neurons) -----------------------------------------
    dict(name='SN_membrane_soma', module='neurons', source='SN_soma',
         eqs=['d/dt V', 'd/dt V_sensed', 'I_in_delayed'],
         bc=['I_in', 'I_out'], sums={'i_ion': 'nanoA'},
         entrance=[port('i_stim', ['I_in'])],
         exit=[port('membrane_voltage', ['V'], 'True'), port('VI_port', ['V', 'I_out'])],
         general=[port('membrane_current', ['i_ion'], 'sum')],
         collapse={'d/dt V': (['i_Na', 'i_Na1_6', 'i_Ca_L', 'i_Kv2', 'i_Kv1_5', 'i_Kv4_2', 'i_KCa_BK', 'i_KCa_SK',
                               'i_M', 'i_NaCa', 'i_NaK', 'i_leak_Na', 'i_leak_K'], 'i_ion')},
         states=['V', 'V_sensed'], run_parameters={'I_in': -0.05},
         invariants=[{'description': 'charge balance with no channels or axon connected: Cm (V - V0) + C_sensed '
                                     '(V_sensed - V_sensed0) = -1e6 I_in t (pF mV = 1e-15 C, nA s = 1e-9 C)',
                      'expr': 'np.max(np.abs(Cm*(V - V[0]) + C_sensed*(V_sensed - V_sensed[0]) + 1e6*I_in*t)) '
                              '<= 1e-6*max(1e-12, np.max(np.abs(1e6*I_in*t)))'}],
         notes='Membrane of the SN soma: Cm dV/dt = -(i_ion + I_in_delayed + I_out), i_ion the sum of the currents '
               'of the channels connected through membrane_current. The injected current I_in reaches the membrane '
               'through a sensing/electrode filter (V_sensed, C_sensed, R_diffSensed) as I_in_delayed; I_out is '
               'the axial current into the axon (VI_port). Copied verbatim from SN_soma except that the channel '
               'currents are one summed port variable.'),
    dict(name='SN_Na_K_concentrations_soma', module='neurons', source='SN_soma',
         eqs=['E_K', 'E_Na', 'd/dt Nai', 'd/dt Ki'],
         bc=[], sums={'j_Na_total': 'mol_per_s', 'j_K_total': 'mol_per_s'},
         general=[port('Na_Nernst_potential', ['E_Na'], 'True'), port('E_K_potential', ['E_K'], 'True'),
                  port('internal_Na', ['Nai'], 'True'), port('Na_flux', ['j_Na_total'], 'sum'),
                  port('K_flux', ['j_K_total'], 'sum')],
         collapse={'d/dt Nai': (['j_Na', 'j_NaCa_Na', 'j_NaK_Na', 'j_leak_Na'], 'j_Na_total'),
                   'd/dt Ki': (['j_Kv1_5', 'j_Kv4_2', 'j_Kv2', 'j_KCa_BK', 'j_KCa_SK', 'j_KM', 'j_NaK_K', 'j_leak_K'],
                               'j_K_total')},
         states=['Nai', 'Ki'],
         invariants=[{'description': 'Nernst: E_Na = 1000 RT/F ln(Nao/Nai), E_K = 1000 RT/F ln(Ko/Ki) (mV)',
                      'expr': 'np.max(np.abs(E_Na - 1000*R*T/F*np.log(Nao/Nai))) + '
                              'np.max(np.abs(E_K - 1000*R*T/F*np.log(Ko/Ki))) <= 1e-9'},
                     {'description': 'alone (no fluxes connected) Nai and Ki stay at their initial values',
                      'expr': 'np.ptp(Nai) == 0 and np.ptp(Ki) == 0'}],
         notes='Cytosolic Na and K of the SN soma: Vol dNai/dt = j_Na_total, Vol dKi/dt = j_K_total (the ion fluxes '
               'of the channels connected through Na_flux / K_flux), and their Nernst potentials E_Na, E_K.'),
    dict(name='SN_Ca_handling_soma', module='neurons', source='SN_soma',
         eqs=['c_bulk', 'j_leak_Ca', 'j_SERCA', 'd/dt B_Cai', 'CaB_Cai', 'd/dt B_Ca_bulk', 'CaB_Ca_bulk',
              'd/dt B_Ca_ER_pre', 'CaB_Ca_ER_pre', 'j_diffCa', 'j_diffCa_ER', 'd/dt Cai', 'd/dt Ca_bulk',
              'd/dt Ca_ER_pre', 'd/dt Ca_ER'],
         bc=[], sums={'j_Ca_in': 'mol_per_s', 'j_ER_release': 'mol_per_s'},
         general=[port('internal_Ca', ['Cai'], 'True'), port('bulk_Ca', ['Ca_bulk'], 'True'),
                  port('ER_Ca', ['Ca_ER'], 'True'), port('Ca_flux', ['j_Ca_in'], 'sum'),
                  port('ER_Ca_release', ['j_ER_release'], 'sum')],
         collapse={'d/dt Cai': (['j_Ca_L', 'j_NaCa_Ca', 'j_IP3', 'j_RyR'], 'j_Ca_in')}, er=True,
         invariants=[CA_TOTAL],
         states=['Cai', 'Ca_bulk', 'Ca_ER_pre', 'Ca_ER', 'B_Cai', 'B_Ca_bulk', 'B_Ca_ER_pre'],
         notes='Four-compartment Ca of the SN soma (sub-membrane Cai, bulk, pre-ER, ER) with a mobile buffer in the '
               'first three, SERCA and a leak between bulk and pre-ER, and diffusion. Cai receives j_Ca_in, the '
               'sum of the Ca fluxes of the channels connected through Ca_flux; the ER loses j_ER_release (IP3R, '
               'RyR). Copied verbatim from SN_soma except for those two summed port variables.'),
    dict(name='SN_varicosity_membrane', module='neurons', source='SN_varicosity',
         eqs=[], bc=['V_in'], units={'V_in': 'milliV'},
         invariants=[{'description': 'V is the axon\'s potential V_in', 'expr': 'np.all(V == V_in)'}], extra=[('V', 'milliV', '<eq/><ci>V</ci><ci>V_in</ci>')],
         entrance=[port('membrane_voltage', ['V_in'])], exit=[port('membrane_voltage', ['V'], 'True')],
         states=[],
         notes='Membrane of the SN varicosity: it has no dynamics of its own (the axon sets its potential), so it '
               'passes the axon\'s V_in on to the varicosity\'s channels as V (V = V_in).'),
    dict(name='SN_Ca_handling_varicosity', module='neurons', source='SN_varicosity',
         eqs=['c_bulk', 'E_Ca', 'j_leak_Ca', 'j_SERCA', 'd/dt B_Cai', 'CaB_Cai', 'd/dt B_Ca_bulk', 'CaB_Ca_bulk',
              'd/dt B_Ca_ER_pre', 'CaB_Ca_ER_pre', 'j_diffCa', 'j_diffCa_ER', 'd/dt Cai', 'd/dt Ca_bulk',
              'd/dt Ca_ER_pre', 'd/dt Ca_ER'],
         bc=[], sums={'j_Ca_in': 'mol_per_s', 'j_ER_release': 'mol_per_s'},
         general=[port('internal_Ca', ['Cai'], 'True'), port('bulk_Ca', ['Ca_bulk'], 'True'),
                  port('ER_Ca', ['Ca_ER'], 'True'), port('Ca_potential', ['E_Ca'], 'True'),
                  port('Ca_flux', ['j_Ca_in'], 'sum'), port('ER_Ca_release', ['j_ER_release'], 'sum')],
         collapse={'d/dt Cai': (['j_Ca_N', 'j_NaCa_Ca', 'j_IP3', 'j_RyR'], 'j_Ca_in')}, er=True,
         invariants=[CA_TOTAL, {'description': 'E_Ca = 1000 RT/F ln(Cao/Cai) (mV; as in the original, no 1/2)',
                                'expr': 'np.max(np.abs(E_Ca - 1000*R*T/F*np.log(Cao/Cai))) <= 1e-9'}],
         states=['Cai', 'Ca_bulk', 'Ca_ER_pre', 'Ca_ER', 'B_Cai', 'B_Ca_bulk', 'B_Ca_ER_pre'],
         notes='Four-compartment Ca of the SN varicosity, as SN_Ca_handling_soma except that the ER leak j_leak_Ca '
               'enters the sub-membrane Cai (not the bulk), as in the original SN_varicosity; also gives '
               'E_Ca = RT/F ln(Cao/Cai) (sic, no factor 1/2, as in the original).'),
    dict(name='SN_NE_release', module='neurons', source='SN_varicosity',
         eqs=['d/dt NE'], bc=['Cai'],
         general=[CAI], exit=[port('synaptic_NE', ['NE'], 'True')],
         states=['NE'], run_parameters={'Cai': 2e-4},
         invariants=[{'description': 'at a clamped Cai, NE follows the exact solution NE_init + k_NE max(Cai - '
                                     'Cai_NE_base, 0)/k_NET (1 - exp(-k_NET t)) (relative 1e-4: CVODE atol 1e-10 on a 2e-4 mM state)',
                      'expr': 'np.max(np.abs(NE - (NE_init + k_NE*max(Cai - Cai_NE_base, 0)/k_NET*(1 - np.exp(-k_NET*t))))) '
                              '<= 1e-4*max(1e-12, k_NE*max(Cai - Cai_NE_base, 0)/k_NET)'}],
         notes='Noradrenaline release of the SN varicosity: dNE/dt = k_NE max(Cai - Cai_NE_base, 0) - '
               'k_NET (NE - NE_init), release driven by sub-membrane Ca above a baseline, reuptake by NET.'),
]
BY_NAME = {c['name']: c for c in COMPONENTS}

# ---------------------------------------------------------------------------------------------
# MathML
# ---------------------------------------------------------------------------------------------


def _kids(e):
    return [k for k in e if isinstance(k.tag, str)]


def _lhs(apply):
    left = _kids(apply)[1]
    if left.tag == f'{{{M}}}ci':
        return left.text.strip()
    return 'd/dt ' + left.xpath('m:ci', namespaces=NS)[0].text.strip()


def _component(root, name):
    return root.xpath(f'//c:component[@name="{name}"]', namespaces=NS)[0]


def equations(root, name):
    comp = _component(root, name)
    return {_lhs(a): a for m in comp.xpath('m:math', namespaces=NS) for a in m.xpath('m:apply', namespaces=NS)}


def variables(root, name):
    return {v.get('name'): v for v in _component(root, name).xpath('c:variable', namespaces=NS)}


def canonical(apply):
    return re.sub(rb'\s+', b'', etree.tostring(apply, method='c14n'))


def ci(name):
    e = etree.Element(f'{{{M}}}ci')
    e.text = name
    return e


def collapse_sum(apply, terms, new):
    '''Replaces the terms of the <plus/> that holds them all by one <ci>new</ci> (at the first term).'''
    for plus in apply.iter(f'{{{M}}}apply'):
        kids = _kids(plus)
        if not kids or kids[0].tag != f'{{{M}}}plus':
            continue
        names = [k.text.strip() if k.tag == f'{{{M}}}ci' else None for k in kids[1:]]
        if not set(terms) <= set(names):
            continue
        assert names[:len(terms)] == list(terms), f'{terms} are not the leading terms of {names}'
        for k in kids[2:1 + len(terms)]:
            plus.remove(k)
        kids[1].text = new
        if len(_kids(plus)) == 2:                     # a sum of one term: the term itself
            plus.getparent().replace(plus, kids[1])
        return apply
    raise ValueError(f'no sum of {terms}')


def collapse_er_release(apply, new):
    '''(-(j_IP3) - j_RyR) -> -(new): the ER's release terms as one summed port variable.'''
    for minus in apply.iter(f'{{{M}}}apply'):
        kids = _kids(minus)
        if (len(kids) == 3 and kids[0].tag == f'{{{M}}}minus' and kids[1].tag == f'{{{M}}}apply'
                and [k.tag for k in _kids(kids[1])] == [f'{{{M}}}minus', f'{{{M}}}ci']
                and _kids(kids[1])[1].text.strip() == 'j_IP3' and kids[2].tag == f'{{{M}}}ci'
                and kids[2].text.strip() == 'j_RyR'):
            minus.remove(kids[2])
            minus.getparent().replace(minus, kids[1])
            _kids(kids[1])[1].text = new
            return apply
    raise ValueError('no -(j_IP3) - j_RyR')


def parse_math(text):
    return etree.fromstring(f'<apply xmlns="{M}" xmlns:cellml="{C}">{text}</apply>')


def cis(apply):
    return {e.text.strip() for e in apply.iter(f'{{{M}}}ci')}


# ---------------------------------------------------------------------------------------------
# building one component: CellML element, config entry, variable kinds
# ---------------------------------------------------------------------------------------------


def build(root, spec):
    src_eqs = equations(root, spec['source'])
    src_vars = variables(root, spec['source'])
    if spec.get('same_in'):
        other = equations(root, spec['same_in'])
        for key in spec['eqs']:
            if key in other and canonical(other[key]) != canonical(src_eqs[key]):
                raise ValueError(f'{spec["name"]}: {key} differs between {spec["source"]} and {spec["same_in"]}')
    eqs = []
    for key in spec['eqs']:
        a = copy.deepcopy(src_eqs[key])
        if key in (spec.get('collapse') or {}):
            terms, new = spec['collapse'][key]
            collapse_sum(a, terms, new)
        if spec.get('er') and key == 'd/dt Ca_ER':
            collapse_er_release(a, 'j_ER_release')
        eqs.append(a)
    extra_units = {}
    for name, units, text in spec.get('extra') or []:
        eqs.append(parse_math(text))
        extra_units[name] = units
    defined = set()
    for a in eqs:
        lhs = _lhs(a)
        defined.add(lhs[5:] if lhs.startswith('d/dt ') else lhs)
    used = set().union(*[cis(a) for a in eqs]) if eqs else set()
    # initial values that name a variable
    inits = {}
    for s in defined:
        v = src_vars.get(s)
        if v is not None and v.get('initial_value') is not None:
            inits[s] = v.get('initial_value')
            try:
                float(inits[s])
            except ValueError:
                used.add(inits[s])
    sums = spec.get('sums') or {}
    bcs = set(spec.get('bc') or [])
    names = sorted(used | defined)
    comp = etree.Element(f'{{{C}}}component', name=spec['name'])
    config_vars = []
    kinds = {}
    ordered = ['t'] if 't' in names else []
    ordered += [n for n in names if n != 't']
    for n in ordered:
        if n in extra_units:
            units = extra_units[n]
        elif n in (spec.get('units') or {}):
            units = spec['units'][n]
        elif n in sums:
            units = sums[n]
        elif n == 't':
            units = 'second'
        else:
            units = src_vars[n].get('units')
        attrs = {'name': n, 'units': units}
        if n in LOCAL:
            attrs['initial_value'] = src_vars[n].get('initial_value')
            etree.SubElement(comp, f'{{{C}}}variable', **attrs)
            continue
        if n == 't':
            attrs['public_interface'] = 'in'
            etree.SubElement(comp, f'{{{C}}}variable', **attrs)
            continue
        if n in inits:
            attrs['initial_value'] = inits[n]
        if n in defined:
            attrs['public_interface'], kind = 'out', 'variable'
        elif n in sums:
            attrs['public_interface'], kind = 'in', 'variable'
        elif n in bcs:
            attrs['public_interface'], kind = 'in', 'boundary_condition'
        elif n in GLOBALS:
            attrs['public_interface'], kind = 'in', 'global_constant'
        else:
            attrs['public_interface'], kind = 'in', 'constant'
        etree.SubElement(comp, f'{{{C}}}variable', **{k: attrs[k] for k in ('name', 'public_interface', 'units', 'initial_value') if k in attrs})
        config_vars.append([n, units, 'access', kind])
        kinds[n] = (kind, units)
    missing = bcs - set(kinds)
    if missing:
        raise ValueError(f'{spec["name"]}: port inputs {missing} are not used')
    math = etree.SubElement(comp, f'{{{M}}}math', nsmap={None: M})
    for a in eqs:
        math.append(a)
    entry = {
        'module_type': spec['name'], 'module_subtype': 'nn', 'module_format': 'cellml',
        'component_file': f'{"ion_channel" if spec["module"] == "ion_channel" else "neurons"}_modules.cellml',
        'component_type': spec['name'],
        'entrance_ports': spec.get('entrance') or [], 'exit_ports': spec.get('exit') or [],
        'general_ports': spec.get('general') or [],
        'variables_and_units': config_vars,
    }
    for p in entry['entrance_ports'] + entry['exit_ports'] + entry['general_ports']:
        for v in p['variables']:
            if v not in kinds:
                raise ValueError(f'{spec["name"]}: port variable {v} is not a variable')
    return comp, entry, kinds


# ---------------------------------------------------------------------------------------------
# files
# ---------------------------------------------------------------------------------------------


def write_cellml(path, comps):
    tree = etree.parse(path, etree.XMLParser(remove_blank_text=True))
    model = tree.getroot()
    new_names = {c.get('name') for c in comps}
    for old in model.xpath('c:component', namespaces=NS):
        if old.get('name') in new_names:
            model.remove(old)
    for c in comps:
        model.append(c)
    etree.indent(model, space='    ')
    old_tag = re.search(rb'<model[^>]*>', open(path, 'rb').read()).group(0)
    text = etree.tostring(tree, xml_declaration=True, encoding='UTF-8', pretty_print=True)
    text = re.sub(rb'<model[^>]*>', lambda _: old_tag, text, count=1)     # keep the file's own model tag
    with open(path, 'wb') as f:
        f.write(text)


def write_config(path, entries):
    with open(path) as f:
        config = json.load(f)
    names = {e['module_type'] for e in entries}
    config = [e for e in config if e.get('module_type') not in names] + entries
    with open(path, 'w') as f:
        json.dump(config, f, indent=2)
        f.write('\n')


def read_rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def write_parameters(path, rows):
    old = read_rows(path)
    names = {r['vessel_type'] for r in rows}
    keep = [r for r in old if r['vessel_type'] not in names]
    fields = ['vessel_type', 'BC_type', 'variable_name', 'units', 'value', 'kind', 'data_reference', 'sourced']
    eol = '\r\n' if b'\r\n' in open(path, 'rb').read() else '\n'     # keep the file's line endings
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator=eol)
        w.writeheader()
        w.writerows(keep + rows)


def units_blocks(path):
    text = open(path).read()
    return {m.group(1): m.group(0) for m in re.finditer(r'<units name="([^"]+)">.*?</units>', text, re.S)}


def add_units(target, source, needed):
    have = units_blocks(target)
    src = units_blocks(source)
    add = []
    todo = list(needed)
    while todo:                                 # units a unit is built from, too
        u = todo.pop()
        if u in have or u not in src or u in [a for a, _ in add]:
            continue
        add.append((u, src[u]))
        todo += re.findall(r'units="([^"]+)"', src[u])
    if not add:
        return []
    text = open(target).read()
    block = ''.join(f'    {b}\n' for _, b in sorted(add))
    text = text.replace('</model>', block + '</model>')
    with open(target, 'w') as f:
        f.write(text)
    return [u for u, _ in add]


def spec_entry(spec, kinds):
    outputs = [n for n in spec['states']]
    outputs += [n for n, (k, _) in kinds.items() if k == 'variable' and n not in outputs
                and n not in (spec.get('sums') or {}) and not n.endswith('_init')]
    gates = [s for s in spec['states'] if s not in ('V', 'V_sensed', 'Nai', 'Ki', 'NE') and not s.startswith(('Ca', 'B_'))]
    inv = []
    for j, expr in spec.get('flux') or []:
        inv.append({'description': f'the ion flux {j} = {expr} (mol/s; the current in nA, so the factor 1e-9)',
                    'expr': f'np.max(np.abs({j} - 1e-9*({expr}))) <= 1e-12*max(1e-30, np.max(np.abs({j})))'})
    if gates:
        inv.append({'description': 'every gating variable stays in [0, 1]',
                    'expr': ' & '.join(f'np.all(({g} >= -1e-9) & ({g} <= 1 + 1e-9))' for g in gates)})
    inv += spec.get('invariants') or []
    entry = {
        'vessel_type': spec['name'], 'BC_type': 'nn', 'module_type': spec['name'],
        'sim_time': spec.get('sim_time', 0.1), 'dt': spec.get('dt', 1.0e-4),
        'timestep': {'scheme': 'rk4', 'dts': [2.0e-06, 1.0e-06, 5.0e-07, 2.5e-07], 't_end': 0.002,
                     'min_order': 0.8, 'tol': 0.001, 'cvode_tol': 0.01},
        'outputs': outputs,
    }
    if inv:
        entry['invariants'] = inv
    if spec.get('run_parameters'):
        entry['run_parameters'] = spec['run_parameters']
    entry['notes'] = (CHANNEL_NOTE if spec['module'] == 'ion_channel' else '') + spec['notes'] + \
        ' Built by tools/build_sn_modules.py; SN_simple_flat checks that the assembled parts reproduce SN_simple.'
    reason = ('a part of the SN_simple sympathetic neuron, whose values are the owner\'s hand-tuned SN_simple values; '
              'validation belongs to the assembled neuron (SN_soma / sympathetic_neuron supermodule).')
    entry['validation'] = {'baseline': {'status': 'not_applicable', 'reason': reason},
                           'calibrate': {'status': 'not_applicable', 'reason': reason}}
    uses = sorted(g for g in ('R', 'F') if g in kinds)
    if uses:
        why = (f'generation fails alone: this component (like SN_soma) declares the global{"s" if len(uses) > 1 else ""} '
               f'{" and ".join(uses)} in the defined units J_per_mol_K / C_per_mol, but the module\'s global rows in '
               f'{spec["module"]}_parameters.csv carry the undefined J_per_MK / C_per_M of the Paci/Tao components '
               '(libcuflynx exits: units in parameters.csv do not match module_config.json). Fix: the pending review '
               'question on unifying R and F to J_per_mol_K / C_per_mol. With that fix applied to a copy of the library '
               'this component\'s run and invariant tests pass; it generates and reproduces SN_simple in SN_simple_flat.')
        entry['expected_failures'] = {t: why for t in ('run_test', 'verification_test_invariants', 'verification_test_BC',
                                                       'verification_test_timestep', 'stability_test')}
    return entry


class _NoAliases(yaml.SafeDumper):
    '''No &id anchors: an anchor name could clash with one already in the file the block is spliced into.'''
    def ignore_aliases(self, data):
        return True


def splice_spec(path, entries):
    '''Adds spec entries for new components to <module>_tests.yaml as text (leaves the rest untouched);
    a component that already has an entry keeps it.'''
    text = open(path).read()
    spec = yaml.safe_load(text)
    have = {c['vessel_type'] for c in spec.get('components') or []}
    new = [e for e in entries if e['vessel_type'] not in have]
    if not new:
        return []
    block = yaml.dump(new, Dumper=_NoAliases, sort_keys=False, width=120, default_flow_style=None, allow_unicode=True)
    lines = text.splitlines(keepends=True)
    start = next(i for i, l in enumerate(lines) if l.startswith('components:'))
    end = next((i for i in range(start + 1, len(lines)) if re.match(r'^[A-Za-z_]+:', lines[i])), len(lines))
    if not lines[end - 1].endswith('\n'):
        lines[end - 1] += '\n'
    lines[end:end] = [block]
    with open(path, 'w') as f:
        f.write(''.join(lines))
    return [e['vessel_type'] for e in new]


# ---------------------------------------------------------------------------------------------
# parameters: each new component's rows carry the original SN_soma / SN_varicosity row
# ---------------------------------------------------------------------------------------------


def module_rows(spec, kinds, neuron_params):
    by = {(r['vessel_type'], r['variable_name']): r for r in neuron_params}
    rows = []
    for n, (kind, units) in kinds.items():
        if kind == 'constant':
            r = by[(spec['source'], n)]
            rows.append({'vessel_type': spec['name'], 'BC_type': 'nn', 'variable_name': n, 'units': units,
                         'value': r['value'], 'kind': 'constant', 'data_reference': r['data_reference'],
                         'sourced': r['sourced']})
        elif kind == 'boundary_condition':
            value, ref = BC[n]
            rows.append({'vessel_type': spec['name'], 'BC_type': 'nn', 'variable_name': n, 'units': units,
                         'value': str(value), 'kind': 'boundary_condition', 'data_reference': ref,
                         'sourced': 'no' if 'definitional' in ref else by.get((spec['source'], n), {}).get('sourced', 'no')})
    return rows


# ---------------------------------------------------------------------------------------------
# the decomposed SN: submodules, internal connections, SN_simple values
# ---------------------------------------------------------------------------------------------

SOMA_CHANNELS = ['i_Na', 'i_Na1_6', 'i_CaL', 'i_Kv2', 'i_Kv1_5', 'i_Kv4_2', 'i_BK', 'i_SK', 'i_M', 'i_NaCa', 'i_NaK',
                 'i_leak_Na', 'i_leak_K']
# neighbour lists in the original's order of summation (see the module docstring)
SOMA = {   # local name: (component, [out_instances])
    'membrane': ('SN_membrane_soma', SOMA_CHANNELS),
    'Na_K': ('SN_Na_K_concentrations_soma', ['i_Na', 'i_Na1_6', 'i_NaCa', 'i_Kv1_5', 'i_Kv4_2', 'i_Kv2', 'i_BK',
                                             'i_SK', 'i_M', 'i_NaK', 'i_leak_Na', 'i_leak_K']),
    'Ca': ('SN_Ca_handling_soma', ['i_CaL', 'i_NaCa', 'IP3R', 'RyR', 'i_BK', 'i_SK']),
    **{ch: (f'{ch}_SN', []) for ch in SOMA_CHANNELS},
    'IP3R': ('IP3R_SN', []), 'RyR': ('RyR_SN', []),
}
VARICOSITY = {
    'membrane': ('SN_varicosity_membrane', ['i_CaN', 'i_NaCa']),
    'Ca': ('SN_Ca_handling_varicosity', ['i_CaN', 'i_NaCa', 'IP3R', 'RyR', 'NE']),
    'i_CaN': ('i_CaN_SN', []), 'i_NaCa': ('i_NaCa_SN', []), 'IP3R': ('IP3R_SN', []), 'RyR': ('RyR_SN', []),
    'NE': ('SN_NE_release', []),
}
PARTS = {'soma': (SOMA, 'SN_soma', 'soma_SN'), 'varicosity': (VARICOSITY, 'SN_varicosity', 'var_SN')}
# port inputs left unconnected in the assembled neuron: parameters, with the original's value
UNCONNECTED = {('soma', 'membrane'): ['I_in'], ('varicosity', 'i_NaCa'): ['Nai']}


def submodule_records(part):
    subs, _, _ = PARTS[part]
    inp = {n: [] for n in subs}
    for n, (_, outs) in subs.items():
        for o in outs:
            inp[o].append(n)
    return [{'name': n, 'module_type': comp, 'module_subtype': 'nn', 'inp_instances': inp[n], 'out_instances': outs}
            for n, (comp, outs) in subs.items()]


def part_parameters(part, kinds_of, sn_simple, prefix=None):
    '''{var}_{submodule} rows for a part (prefix None: local names; else "<prefix>_<sub>") with SN_simple's values.'''
    subs, _, old_vessel = PARTS[part]
    values = {r['variable_name'].strip(): r for r in sn_simple}
    rows = []
    for n, (comp, _) in subs.items():
        for var, (kind, units) in kinds_of[comp].items():
            if kind == 'constant' or var in UNCONNECTED.get((part, n), []):
                r = values[f'{var}_{old_vessel}']
                name = f'{var}_{prefix + "_" + n if prefix else n}'
                rows.append([name, units, r['value'].strip(), r['data_reference'].strip()])
    return rows


AXON_PARAMS = ['channel_ratio', 'C', 'R', 'V_rest_approx']


def axon_parameters(sn_simple, name):
    values = {r['variable_name'].strip(): r for r in sn_simple}
    return [[f'{v}_{name}', values[f'{v}_axon_SN']['units'].strip(), values[f'{v}_axon_SN']['value'].strip(),
             values[f'{v}_axon_SN']['data_reference'].strip()] for v in AXON_PARAMS]


def global_parameters(sn_simple):
    values = {r['variable_name'].strip(): r for r in sn_simple}
    return [[g, values[g]['units'].strip(), values[g]['value'].strip(), values[g]['data_reference'].strip()]
            for g in GLOBALS]


def write_param_csv(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['variable_name', 'units', 'value', 'data_reference'])
        w.writerows(rows)


def output_map(root, sn_simple_outputs, kinds_of):
    '''SN_simple's outputs (soma_SN/x, var_SN/x, axon_SN/x) -> the decomposed model's names.'''
    owner = {}
    for part, (subs, _, old_vessel) in PARTS.items():
        for n, (comp, _) in subs.items():
            for var, (kind, _) in kinds_of[comp].items():
                if kind == 'variable' and var not in (BY_NAME[comp].get('sums') or {}):
                    owner.setdefault((old_vessel, var), f'SN_{part}_{n}')
    out = {}
    for name in sn_simple_outputs:
        vessel, _, var = name.partition('/')
        if vessel == 'axon_SN':
            out[name] = f'SN_axon/{var}'
        elif (vessel, var) in owner:
            out[name] = f'{owner[(vessel, var)]}/{var}'
        else:
            raise KeyError(f'no component holds {name}')
    return out


def sn_simple_outputs(root):
    '''The time-varying variables of SN_simple's three components (what its model logs): the states,
    and every variable that depends on a state or on a connected port input.'''
    names = []
    for vessel, comp, inputs in (('soma_SN', 'SN_soma', {'I_out'}), ('axon_SN', 'SN_axon', {'V_in'}),
                                 ('var_SN', 'SN_varicosity', {'V'})):
        eqs = equations(root, comp)
        dynamic = {k[5:] for k in eqs if k.startswith('d/dt ')} | inputs
        changed = True
        while changed:
            changed = False
            for lhs, a in eqs.items():
                if not lhs.startswith('d/dt ') and lhs not in dynamic and cis(_kids(a)[2]) & dynamic:
                    dynamic.add(lhs)
                    changed = True
        names += [f'{vessel}/{v}' for v in sorted(dynamic - inputs)]
    return names


def copy_reference(model_dir, name):
    src, dst = os.path.join(SN_SIMPLE, 'reference'), os.path.join(model_dir, 'reference')
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(dst)
    for f in os.listdir(src):
        shutil.copy(os.path.join(src, f), os.path.join(dst, f.replace('SN_simple', name, 1)))


def write_yaml(path, data):
    with open(path, 'w') as f:
        yaml.safe_dump(data, f, sort_keys=False, width=120, default_flow_style=None, allow_unicode=True)


def main():
    root = etree.parse(SOURCE).getroot()
    neuron_params = read_rows(os.path.join(NEURONS, 'neurons_parameters.csv'))
    built = {c['name']: build(root, c) for c in COMPONENTS}
    kinds_of = {n: k for n, (_, _, k) in built.items()}

    for module, d in (('ion_channel', ION), ('neurons', NEURONS)):
        specs = [c for c in COMPONENTS if c['module'] == module]
        write_cellml(os.path.join(d, f'{module}_modules.cellml'), [built[c['name']][0] for c in specs])
        write_config(os.path.join(d, f'{module}_modules_config.json'), [built[c['name']][1] for c in specs])
        rows = [r for c in specs for r in module_rows(c, kinds_of[c['name']], neuron_params)]
        write_parameters(os.path.join(d, f'{module}_parameters.csv'), rows)
        added = splice_spec(os.path.join(d, f'{module}_tests.yaml'), [spec_entry(c, kinds_of[c['name']]) for c in specs])
        print(f'{module}: {len(specs)} components; spec entries added: {added}')
    units = {u for c in COMPONENTS if c['module'] == 'ion_channel' for _, u in kinds_of[c['name']].values()}
    units |= {cn.get(f'{{{C}}}units') for c in COMPONENTS if c['module'] == 'ion_channel'
              for cn in built[c['name']][0].iter(f'{{{M}}}cn')}                # constants' units too
    added = add_units(os.path.join(ION, 'ion_channel_units.cellml'), os.path.join(NEURONS, 'neurons_units.cellml'), units)
    print(f'ion_channel units added: {added}')

    sn_simple = [r for r in read_rows(os.path.join(SN_SIMPLE, 'SN_simple_parameters.csv'))
                 if r['variable_name'].strip() and not r['variable_name'].startswith('#')]
    omap = output_map(root, sn_simple_outputs(root), kinds_of)
    ignore = {'soma_SN/t': 'time', 'axon_SN/t': 'time', 'var_SN/t': 'time'}

    # ---- SN_simple_flat: the decomposed neuron written out (the names a sympathetic_neuron instance "SN" expands to)
    flat_dir = os.path.join(MODULES, 'system', 'cellular', 'SN_simple_flat')
    os.makedirs(flat_dir, exist_ok=True)
    records = []
    for part in ('soma', 'varicosity'):
        for rec in submodule_records(part):
            r = dict(rec, name=f'SN_{part}_{rec["name"]}',
                     inp_instances=[f'SN_{part}_{n}' for n in rec['inp_instances']],
                     out_instances=[f'SN_{part}_{n}' for n in rec['out_instances']])
            records.append(r)
        if part == 'soma':
            records[0]['out_instances'].append('SN_axon')
            records.append({'name': 'SN_axon', 'module_type': 'SN_axon', 'module_subtype': 'nn',
                            'inp_instances': ['SN_soma_membrane'], 'out_instances': ['SN_varicosity_membrane']})
        else:
            next(r for r in records if r['name'] == 'SN_varicosity_membrane')['inp_instances'].insert(0, 'SN_axon')
    stale_csv = os.path.join(flat_dir, 'SN_simple_flat_vessel_array.csv')
    if os.path.exists(stale_csv):
        os.remove(stale_csv)
    vessel_array.write_records(os.path.join(flat_dir, 'SN_simple_flat_vessel_array.json'), records)
    params = (global_parameters(sn_simple) + part_parameters('soma', kinds_of, sn_simple, 'SN_soma')
              + axon_parameters(sn_simple, 'SN_axon') + part_parameters('varicosity', kinds_of, sn_simple, 'SN_varicosity'))
    write_param_csv(os.path.join(flat_dir, 'SN_simple_flat_parameters.csv'), params)
    copy_reference(flat_dir, 'SN_simple_flat')
    write_yaml(os.path.join(flat_dir, 'SN_simple_flat_system.yaml'), {
        'model': 'SN_simple_flat', 'category': 'cellular', 'reviewed': False,
        'description': 'SN_simple (sympathetic neuron: soma, axon, varicosity) with SN_soma and SN_varicosity '
                       'decomposed into their ion channels (modules/ion_channel, *_SN) and the neuron-specific '
                       'membrane, concentration, Ca-handling and NE-release components (modules/cell/neurons), '
                       'wired explicitly. Built by tools/build_sn_modules.py; the vessel names are those a '
                       'sympathetic_neuron supermodule instance "SN" expands to (SN_soma_<part>, SN_axon, '
                       'SN_varicosity_<part>).',
        'sim_time': 2.0, 'pre_time': 0.0, 'dt': 0.01,
        'equivalence': {'reference': 'reference', 'output_map': omap, 'ignore': ignore, 'wrapped': {},
                        'reciprocal': [], 'tol': 1.0e-9, 'solver_info': {'rtol': 1.0e-10, 'atol': 1.0e-12}},
        'invariants': [],
        'notes': ['The reference is SN_simple\'s (circulatory_autogen\'s original with its bundled SN_soma, SN_axon '
                  'and SN_varicosity). Every SN_simple output is compared, under its new name (equivalence.output_map), '
                  'to 1e-9: the parts copy the original equations verbatim, and only the sums of channel currents and '
                  'fluxes become summed ports, added in the original order.'],
    })

    # ---- supermodules
    sup = os.path.join(MODULES, 'supermodules')
    for part, title in (('soma', 'SN_soma'), ('varicosity', 'SN_varicosity')):
        d = os.path.join(sup, title)
        os.makedirs(d, exist_ok=True)
        subs = submodule_records(part)
        host = ('the membrane couples to the axon: per_submodule_outputs {"membrane": ["<axon>"]} '
                '(VI_port: V out, the axial current I_out in); optional stimulus through membrane\'s i_stim'
                if part == 'soma' else
                'the membrane takes the axon\'s potential: per_submodule_inputs {"membrane": ["<axon>"]}; '
                'NE (synaptic_NE) is available from NE')
        desc = {'soma': 'Sympathetic-neuron soma (SN_simple\'s SN_soma) built from its parts: the membrane with the '
                        'sensing filter, Na/K concentrations with their Nernst potentials, four-compartment Ca handling, '
                        'and the ion channels i_Na, i_Na1_6, i_CaL, i_Kv2, i_Kv1_5, i_Kv4_2, i_BK, i_SK, i_M, i_NaCa, '
                        'i_NaK, i_leak_Na, i_leak_K, IP3R and RyR (modules/ion_channel, *_SN). ',
                'varicosity': 'Sympathetic-neuron varicosity (SN_simple\'s SN_varicosity) built from its parts: a '
                              'pass-through membrane (the axon sets V), four-compartment Ca handling with E_Ca, N-type '
                              'Ca channel, NaCa exchanger, IP3R, RyR and NE release. '}[part]
        desc += f'Each submodule becomes <instance>_<submodule> in a model; {host}.'
        entry = {'module_type': title, 'module_subtype': 'supermodule', 'module_format': 'supermodule',
                 'description': desc, 'submodules': subs, 'default_parameters': f'{title}_parameters.csv'}
        with open(os.path.join(d, f'{title}_modules_config.json'), 'w') as f:
            json.dump([entry], f, indent=2)
            f.write('\n')
        write_param_csv(os.path.join(d, f'{title}_parameters.csv'),
                        global_parameters(sn_simple) + part_parameters(part, kinds_of, sn_simple))
        write_yaml(os.path.join(d, f'{title}_supermodule.yaml'), {
            'supermodule': title,
            'description': desc + f' Default parameters are SN_simple\'s ({title}_parameters.csv). Built by '
                                  'tools/build_sn_modules.py. Tested through the sympathetic_neuron supermodule.',
            'reviewed': False, 'globals': list(GLOBALS)})

    d = os.path.join(sup, 'sympathetic_neuron')
    os.makedirs(d, exist_ok=True)
    desc = ('Sympathetic neuron (SN_simple): an SN_soma supermodule instance "soma", the SN_axon cable "axon" and an '
            'SN_varicosity supermodule instance "varicosity", coupled soma -> axon -> varicosity through the plain '
            'axon module (the soma\'s membrane gives the axon V and takes its axial current; the axon sets the '
            'varicosity membrane\'s potential). A model uses one instance, e.g. {"name": "SN", "module_type": '
            '"sympathetic_neuron", "module_subtype": "supermodule"}; nested expansion names the parts SN_soma_<part>, '
            'SN_axon and SN_varicosity_<part>.')
    entry = {'module_type': 'sympathetic_neuron', 'module_subtype': 'supermodule', 'module_format': 'supermodule',
             'description': desc,
             'submodules': [
                 {'name': 'soma', 'module_type': 'SN_soma', 'module_subtype': 'supermodule',
                  'per_submodule_outputs': {'membrane': ['axon']}},
                 {'name': 'axon', 'module_type': 'SN_axon', 'module_subtype': 'nn',
                  'inp_instances': ['soma'], 'out_instances': ['varicosity']},
                 {'name': 'varicosity', 'module_type': 'SN_varicosity', 'module_subtype': 'supermodule',
                  'per_submodule_inputs': {'membrane': ['axon']}}],
             'default_parameters': 'sympathetic_neuron_parameters.csv'}
    with open(os.path.join(d, 'sympathetic_neuron_modules_config.json'), 'w') as f:
        json.dump([entry], f, indent=2)
        f.write('\n')
    write_param_csv(os.path.join(d, 'sympathetic_neuron_parameters.csv'),
                    global_parameters(sn_simple) + axon_parameters(sn_simple, 'axon'))
    write_yaml(os.path.join(d, 'sympathetic_neuron_supermodule.yaml'), {
        'supermodule': 'sympathetic_neuron',
        'description': desc + ' The soma and varicosity parameters are the SN_soma / SN_varicosity supermodules\' '
                              'defaults; sympathetic_neuron_parameters.csv holds the axon\'s and the globals '
                              '(SN_simple\'s values). Built by tools/build_sn_modules.py.',
        'reviewed': False,
        'tests': {'equivalent': [{
            'model': 'cellular/SN_simple_supermodules', 'reproduces': 'cellular/SN_simple', 'instance': 'SN',
            'tol': 1.0e-9, 'solver_info': {'rtol': 1.0e-10, 'atol': 1.0e-12},
            'output_map': omap, 'ignore': ignore}]},
        'globals': list(GLOBALS)})

    sm_dir = os.path.join(MODULES, 'system', 'cellular', 'SN_simple_supermodules')
    os.makedirs(sm_dir, exist_ok=True)
    vessel_array.write_records(os.path.join(sm_dir, 'SN_simple_supermodules_vessel_array.json'),
                               [{'name': 'SN', 'module_type': 'sympathetic_neuron', 'module_subtype': 'supermodule'}])
    write_param_csv(os.path.join(sm_dir, 'SN_simple_supermodules_parameters.csv'), global_parameters(sn_simple))
    copy_reference(sm_dir, 'SN_simple_supermodules')
    write_yaml(os.path.join(sm_dir, 'SN_simple_supermodules_system.yaml'), {
        'model': 'SN_simple_supermodules', 'category': 'cellular', 'reviewed': False,
        'description': 'SN_simple as one sympathetic_neuron supermodule instance "SN" (modules/supermodules/'
                       'sympathetic_neuron), which nests the SN_soma and SN_varicosity supermodules. Built by '
                       'tools/build_sn_modules.py.',
        'sim_time': 2.0, 'pre_time': 0.0, 'dt': 0.01,
        'equivalence': {'reference': 'reference', 'output_map': omap, 'ignore': ignore, 'wrapped': {},
                        'reciprocal': [], 'tol': 1.0e-9, 'solver_info': {'rtol': 1.0e-10, 'atol': 1.0e-12}},
        'invariants': [],
        'notes': ['Only the globals are in the parameters file; every other value comes from the supermodules\' '
                  'defaults (SN_simple\'s values). Must reproduce circulatory_autogen\'s SN_simple (equivalence) '
                  'and the library\'s SN_simple (tests/test_supermodules.py, 1e-9).'],
    })
    print(f'wrote SN_simple_flat ({len(records)} vessels), 3 supermodules, SN_simple_supermodules; '
          f'{len(omap)} outputs mapped')


if __name__ == '__main__':
    main()
