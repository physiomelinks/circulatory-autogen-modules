"""
Library restructure proposal (2026-10-06 placement, naming and required_citations rules): a review page
for the owner to approve before anything moves.

    venv/bin/python tools/library_restructure_proposal.py      # -> reviews/library_restructure_proposal.html

Read-only: no module file is moved, renamed or edited. The page is built from the library as it is:
- the curated rename/merge table (CHANGES below) is checked against the modules/ tree;
- the equation evidence for each merge is read from the versions' CellML (MathML printed as text) and
  compared equation by equation;
- port compatibility is read from the versions' modules_config.json;
- downstream users are found by scanning system_models/ vessel arrays, supermodule configs, code,
  tests, manifests and the site;
- the required_citations proposal (CITATION_RULES) is applied to every version after the rename, with
  full references from the repo's .bib files (or NEW_BIB for papers not in any .bib yet);
- the proposed bib keys follow <surname><year><first significant title word>, all lowercase, and the
  cost of rekeying every existing key is counted.

The libcuflynx boundary-condition checks quoted in the BC-prefix section are located in the editable
install (LIBCUFLYNX below) when it exists, so their line numbers stay current.
"""
import collections
import datetime
import glob
import hashlib
import html
import json
import os
import re
import sys
import unicodedata
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(REPO, 'modules')
OUT = os.path.join(REPO, 'reviews', 'library_restructure_proposal.html')
LIBCUFLYNX = os.environ.get(
    'LIBCUFLYNX_SRC',
    '/home/farg967/Documents/git_projects/circulatory_autogen/.claude/worktrees/prediction-features/src/libcuflynx')
EXCLUDED = ('modules/system', 'modules/poiseuille_transport')   # another session's uncommitted work

# ======================================================================================== the proposal
# (old type path under modules/, old version) -> (new type path, new version, kind, evidence group, note)
# kind: 'merge' (joins a new or existing mechanism type as a version), 'rename' (type renamed, version
# renamed), 'version' (same type, version renamed), 'si' / 'old' (_2 / _OLD type folded in as a version).
CHANGES = [
    # ---- cell: one type per mechanism, versions by maths (rules 1, 2, 4, 5)
    ('cell/cardiomyocytes/Ca_dynamics_Paci_2013', 'nn', 'cell/Ca_handling', 'cardiomyocyte_Paci2013_v01', 'merge', 'Ca_handling',
     'Paci SR model (release gate g, SERCA uptake, leak, rapid buffering): different maths from both SN versions, so a version, not an instance.'),
    ('cell/neuron/soma/SN_Ca_handling_soma', 'nn', 'cell/Ca_handling', 'SN_soma_Argus2026_v01', 'merge', 'Ca_handling',
     'Khaliq sub-membrane shell (V_shell = A_mem d_shell, tau_shell exchange); leak into the bulk.'),
    ('cell/neuron/varicosity/SN_Ca_handling_varicosity', 'nn', 'cell/Ca_handling', 'SN_varicosity_Argus2026_v01', 'merge', 'Ca_handling',
     'Fixed sub-membrane fraction c_SL with D_Ca exchange; ER leak into Cai; also gives E_Ca.'),
    ('cell/cardiomyocytes/myocyte_membrane_voltage', 'nn', 'cell/membrane_potential', 'cardiomyocyte_Paci2013_v01', 'merge', 'membrane_potential',
     'dVm/dt = -(sum of the 12 named Paci currents - i_stim), currents in A/F.'),
    ('cell/neuron/sympathetic_neuron_membrane_voltage', 'nn', 'cell/membrane_potential', 'SN_Tao2011_v01', 'merge', 'membrane_potential',
     'dVm/dt = -(sum of the 5 named Tao currents - i_stim)/C_m, currents in nA.'),
    ('cell/neuron/soma/SN_membrane_soma', 'nn', 'cell/membrane_potential', 'SN_soma_Argus2026_v01', 'merge', 'membrane_potential',
     'Summed membrane_current port, electrode filter (V_sensed), axial I_out, geometry A_mem, Vol, Cm from d.'),
    ('cell/neuron/varicosity/SN_varicosity_membrane', 'nn', 'cell/membrane_potential', 'SN_varicosity_Argus2026_v01', 'merge', 'membrane_potential',
     'No ODE: V = V_in from the axon, plus geometry. Kept as a membrane_potential version (the varicosity membrane is a pass-through).'),
    ('cell/cardiomyocytes/Na_dynamics_Paci_2013', 'nn', 'cell/ion_concentrations', 'cardiomyocyte_Paci2013_v01', 'merge', 'ion_concentrations',
     'Nai only, driven by named A/F currents (i_Na, i_b_Na, 3 i_NaK, 3 i_NaCa).'),
    ('cell/neuron/soma/SN_Na_K_concentrations_soma', 'nn', 'cell/ion_concentrations', 'SN_soma_Argus2026_v01', 'merge', 'ion_concentrations',
     'Nai and Ki driven by summed mol/s fluxes, plus E_Na and E_K.'),
    ('cell/cardiomyocytes/electric_potentials_Paci_2013', 'nn', 'cell/reversal_potentials', 'cardiomyocyte_Paci2013_v01', 'merge', 'reversal_potentials',
     'Nernst/GHK reversal potentials only (E_Na, E_K, E_Ks, E_Ca), no state. New mechanism type; see open question Q6.'),
    ('cell/neuron/NE_release_Tao_2011', 'nn', 'cell/neurotransmitter_release', 'SN_varicosity_Tao2011_v01', 'merge', 'neurotransmitter_release',
     'Tao 2011 five-site Ca binding vesicle scheme (V_0Ca..V_5Ca, fusion E = gamma V_5Ca).'),
    ('cell/neuron/varicosity/SN_NE_release', 'nn', 'cell/neurotransmitter_release', 'SN_varicosity_Argus2026_v01', 'merge', 'neurotransmitter_release',
     "Owner's first-order release above a Ca baseline with NET reuptake."),
    ('cell/NKE_pump', 'nn', 'cell/ion_channels/i_NaK', 'ArgusUNPUBLISHED_v01', 'merge', 'i_NaK',
     'Four-state bond-graph Na/K-ATPase cycle (chemostats, membrane charge): different maths and ports from both current-form i_NaK versions.'),
    ('cell/ion_channels/i_piezo_Na', 'Argus_unpublished_v01', 'cell/ion_channels/i_piezo_Na', 'ArgusUNPUBLISHED_v01', 'version', None,
     'Version name to the ArgusUNPUBLISHED_vXX spelling.'),
    # ---- boundary conditions
    ('boundary_conditions/electrical_stim_input_Paci_2013', 'nn', 'boundary_conditions/i_stim_periodic', 'cardiomyocyte_Paci2013_v01', 'merge', 'i_stim_periodic',
     'Same rectangular pulse train as i_stim_periodic nn, but divided by Cm (A/F) and positive sign: different maths, so a version.'),
    # ---- control: surname-free types; all Gee 2023 modules are ArgusUNPUBLISHED for now
    *[(f'control/{n}_Gee2023', 'nn', f'control/{n}', 'ArgusUNPUBLISHED_v01', 'rename', None,
       'Gee 2023 brainstem nucleus; acronym kept (see Q9).') for n in ('DMV', 'LCN', 'NA', 'NActr', 'NTS', 'PNDMV', 'PNNA')],
    ('control/cardiopulmonary_receptor_Gee2023', 'nn', 'control/cardiopulmonary_receptor', 'ArgusUNPUBLISHED_v01', 'rename', None, ''),
    ('control/lung_stretch_receptor_Gee2023', 'nn', 'control/lung_stretch_receptor', 'ArgusUNPUBLISHED_v01', 'rename', 'lung_stretch_receptor', ''),
    ('control/lung_stretch_receptor_with_threshold_OLD', 'nn', 'control/lung_stretch_receptor', 'ArgusUNPUBLISHED_OLD', 'old', 'lung_stretch_receptor',
     'Older threshold-on-lung-volume form; _OLD type becomes an _OLD version.'),
    ('control/heart_effector_Gee2023', 'nn_lv', 'control/heart_effector', 'nn_lv_ArgusUNPUBLISHED_v01', 'rename', 'heart_effector',
     'Same maths as rv; separate versions only because the port_type names differ (Elastance_delta_port_lv / _rv); see Q8. '
     'Keeps nn_ first: "lv_"/"rv_" end in v, which libcuflynx reads as a flow-out port and rejects (BC-prefix section).'),
    ('control/heart_effector_Gee2023', 'nn_rv', 'control/heart_effector', 'nn_rv_ArgusUNPUBLISHED_v01', 'rename', 'heart_effector', ''),
    ('control/heart_effector_Gee2023_OLD', 'nn_rv', 'control/heart_effector', 'nn_rv_ArgusUNPUBLISHED_OLD', 'old', 'heart_effector',
     'Equations identical to rv; the config declares inverse (compliance) variables the CellML no longer has; see Q8.'),
    ('control/peripheral_resistance_effector_Gee2023', 'nn', 'control/peripheral_resistance_effector', 'ArgusUNPUBLISHED_v01', 'rename', None,
     'Same equations as efferent_resistance_effector_Ursino under other variable names; see Q8.'),
    ('control/sympathetic_afferent_Gee2023', 'nn', 'control/sympathetic_afferent', 'ArgusUNPUBLISHED_v01', 'rename', None, ''),
    ('control/sympathetic_efferent_Gee2023', 'nn', 'control/sympathetic_efferent', 'ArgusUNPUBLISHED_v01', 'rename', None, ''),
    ('control/unstressed_volume_effector_Gee2023', 'nn', 'control/unstressed_volume_effector', 'ArgusUNPUBLISHED_v01', 'rename', None, ''),
    ('control/efferent_resistance_effector_Ursino', 'nn', 'control/efferent_resistance_effector', 'Ursino1998_v01', 'merge', 'efferent_resistance_effector',
     'Ursino log effector; the existing nn version is the saturating (Argus) form.'),
    ('control/efferent_respiratory_effector_Gee_Argus_2026', 'nn', 'control/efferent_respiratory_effector', 'GeeArgus2026_v01', 'merge', 'efferent_respiratory_effector',
     'Multiplicative muscle-pressure drive (alpha_c alpha_p, relaxing to 1) vs the additive nn form.'),
    # ---- respiratory
    ('respiratory/simple_lung_bg', 'nn', 'respiratory/lung_mechanics', 'ArgusUNPUBLISHED_v01', 'merge', 'lung_mechanics',
     'Albanese 2016 airway tree with owner changes (pleural compliance through R_plwall, alpha_LR).'),
    ('respiratory/lung_bg_Argus_Chase', 'nn', 'respiratory/lung_mechanics', 'ventilated_ArgusUNPUBLISHED_v01', 'merge', 'lung_mechanics',
     "Owner's mechanically ventilated lung (ventilator PI control, tube, mouth leak, B-spline effort) around the same airway tree."),
    # ---- _2 (SI duplicates) and _OLD
    ('coupling/sum_blood_vol_2', 'nn', 'coupling/sum_blood_vol', 'nn_SI', 'si', 'sum_blood_vol', 'Identical equations; SI units.'),
    *[(f'vessels/compartments/simple_2', v, 'vessels/compartments/simple', f'{v}_SI', 'si', 'simple', 'Identical equations; SI units.')
      for v in ('pp_noI', 'pv_noI', 'vp_noI', 'vv_noI')],
    *[(f'vessels/junctions/Min_junction_2', v, 'vessels/junctions/Min_junction', f'{v}_SI', 'si', 'Min_junction', 'Identical equations; SI units.')
      for v in ('vp_noI', 'vv_noI')],
    ('vessels/junctions/MinNout_junction_2', 'vv_noI', 'vessels/junctions/MinNout_junction', 'vv_noI_SI', 'si', 'MinNout_junction', 'Identical equations; SI units.'),
    *[(f'vessels/junctions/Nout_junction_2', v, 'vessels/junctions/Nout_junction', f'{v}_SI', 'si', 'Nout_junction', 'Identical equations; SI units.')
      for v in ('pv_noI', 'vv_noI')],
    ('heart', 'vp_simple_2', 'heart', 'vp_simple_SI', 'version', None,
     'An SI duplicate kept as a version with a _2 suffix; renamed for the _SI convention (the rule names _2 types; see Q10).'),
]

# Curated verdicts per merge group, shown beside the computed equation comparison.
GROUP_NOTES = {
    'Ca_handling': (
        'Three versions. The two SN versions share 15 of their 20 equations but differ in the sub-membrane compartment '
        '(soma: a Khaliq 2003 shell of volume A_mem d_shell exchanging with time constant tau_shell; varicosity: a fixed '
        'volume fraction c_SL exchanging through D_Ca), in where the ER leak lands (bulk vs Cai) and in the varicosity\'s E_Ca. '
        'Paci 2013 is a different model altogether: one cytosolic and one SR pool with rapid-equilibrium buffering, '
        'a release current gated by d and the Ca-inactivation gate g, a Hill-2 SERCA and a linear leak. '
        'The owner expected "probably a cardiomyocyte instance for now", i.e. an instance if only the parameters differed; '
        'the equations differ, so by rule 4 it is a version, cardiomyocyte_Paci2013_v01.'),
    'membrane_potential': (
        'Four versions, all different. Paci and Tao both integrate a fixed list of named currents (12 in A/F; 5 in nA divided by C_m); '
        'the SN soma sums one multi-port and adds the electrode filter and the geometry; the SN varicosity has no ODE. '
        'Ports differ throughout: the named-current versions take one port per current, the SN versions one summed membrane_current port.'),
    'ion_concentrations': (
        'Two versions. Paci integrates Nai only, from named A/F currents with fixed stoichiometry; the SN soma integrates Nai and Ki '
        'from summed mol/s fluxes and gives E_Na, E_K. Different maths and ports.'),
    'reversal_potentials': (
        'One version. Placement alternatives: (a) a new mechanism type cell/reversal_potentials (proposed); '
        '(b) a version of ion_concentrations (the SN soma computes E_Na and E_K there), which would put two ion_concentrations '
        'versions in one Paci cell. See Q6.'),
    'neurotransmitter_release': (
        'Two versions, different maths: Tao 2011 tracks vesicles through five Ca-binding states; the SN varicosity release is a '
        'first-order rate above a Ca baseline with NET reuptake.'),
    'i_NaK': (
        'Three versions. Paci (Luo-Rudy-type voltage factor) and Argus2026 (Nygren/Lindblad form) are current laws on the same kind of '
        'ports; NKE_pump is a bond-graph enzyme cycle on chemostats and a membrane charge, with no current port. It joins i_NaK by '
        'mechanism (rule 5), but cannot be swapped for the other two in a cell without a bond-graph interface; see the port table.'),
    'i_stim_periodic': (
        'Same pulse logic (start, end, period, pulse duration); the Paci form divides by Cm (A/F) and has the opposite sign, '
        'and its frequency is an input (period = 1/frequency).'),
    'efferent_resistance_effector': (
        'Different maths: the existing nn version is a saturating first-order effector (G sigma/(f_mid + sigma)); '
        'Ursino 1998 is the logarithmic static gain followed by a first-order lag.'),
    'efferent_respiratory_effector': (
        'Same structure; the GeeArgus form relaxes the muscle-pressure multipliers towards 1 and multiplies them '
        '(alpha_c alpha_p) where nn relaxes to 0 and adds them.'),
    'heart_effector': (
        'lv and rv have identical equations and differ only in parameter values and in their port_type names; the OLD rv has '
        'identical equations too. Under rule 4 lv/rv would be instances of one version if the port_types were made chamber-neutral (Q8).'),
    'lung_stretch_receptor': 'Different input: the current form takes q_str; the OLD form thresholds the lung volume q_A + q_b.',
    'lung_mechanics': (
        'Both are bond-graph airway trees (larynx/trachea-bronchi-alveoli); the ventilated version adds the ventilator, tube, '
        'mouth compartment and spline effort and splits R_tot/C_tot. Different maths, so two versions.'),
    'sum_blood_vol': 'Equations identical; only the units differ (SI).',
    'simple': 'Equations identical; only the units differ (SI).',
    'Min_junction': 'Equations identical; only the units differ (SI).',
    'MinNout_junction': 'Equations identical; only the units differ (SI).',
    'Nout_junction': 'Equations identical; only the units differ (SI).',
}

# Same-meaning or identical-maths candidates outside the proposal (open question Q11): shown with the
# computed comparison so the owner can decide; nothing here is in the rename table.
AMBIGUOUS = [
    ('Ursino log effector under five names',
     [('control/efferent_heart_elastance_effector', 'nn_lv'), ('control/efferent_resistance_effector_Ursino', 'nn'),
      ('control/peripheral_resistance_effector_Gee2023', 'nn'), ('control/unstressed_volume_effector_Gee2023', 'nn')],
     'sigma = G ln(f - f_min + 1) for f >= f_min, then a first-order lag. Identical for the first two; the Gee two use other '
     'variable names (Delta_Rsp, f_esp; Delta_Vusv, f_esv). One generic efferent_effector type with instances would need the '
     'variable and port names unified, which changes the system models that connect them.'),
    ('Sympathetic heart effectors',
     [('control/heart_effector_Gee2023', 'nn_lv'), ('control/efferent_heart_period_effector', 'nn')],
     'Both add a logarithmic sympathetic and a linear vagal first-order term; Gee also scales the sum.'),
    ('Kidney left / right',
     [('organs/BVC_Kidney_left', 'nn'), ('organs/BVC_Kidney_right', 'nn')],
     'Identical equations: by rule 4 one type (e.g. kidney_volume_control) with left/right instances, if the ports agree.'),
    ('Capillary volume control intake / uptake',
     [('vessels/compartments/capillary_simple_vol_ctl_intake', 'vv'), ('vessels/compartments/capillary_simple_vol_ctl_uptake', 'vv')],
     'Identical equations: instances of one version, if the ports agree.'),
    ('Blood-volume sums',
     [('coupling/volume_sum', 'nn'), ('coupling/sum_blood_vol', 'nn')],
     'Same mechanism (sum the volume ports of a network). volume_sum also exposes the total on an exit port.'),
    ('Background Na current',
     [('cell/ion_channels/i_b_Na', 'Paci2013_v01'), ('cell/ion_channels/i_leak_Na', 'Argus2026_v01')],
     'Both are ohmic background Na currents. If one mechanism, i_b_Na Paci2013_v01 becomes a version of one of them.'),
    ('Tao 2011 channels beside their molecular counterparts',
     [('cell/ion_channels/i_KV', 'Tao2011_v01'), ('cell/ion_channels/i_Kv2', 'Argus2026_v01')],
     'i_KV (Tao delayed rectifier) and i_Kv2 are both Belluzzi & Sacchi 1988 delayed rectifiers; similarly i_KCa vs i_BK/i_SK and '
     'i_Ca vs i_CaL/i_CaN. Kept separate unless you count them as one mechanism.'),
    ('Local flow control brain / muscle',
     [('control/local_flow_control_brain', 'nn'), ('control/local_flow_control_muscle', 'nn')],
     'Named after the organ (rule 1). Different maths; could become versions of one local_flow_control type (brain_..., muscle_...).'),
]

# Version names that are kept but questioned (shown in the open questions).
# ---------------------------------------------------------------------------- required_citations
# Ordered rules: (regex on the post-rename "<type path>/<version>", [keys], uncertain, basis). First match wins.
# Keys are in the proposed convention (see make_key); ARGUS_SN etc. are placeholders the owner decides on (Q1).
ARGUS_SN = 'argus2026sympathetic'          # placeholder: the planned sympathetic-neuron paper
ARGUS_UNP = 'argus2026unpublished'          # placeholder: owner-modified, no publication
GEE_ARGUS = 'gee2026unpublished'            # placeholder: the Gee & Argus respiratory-control work
CA_SOFT = 'argus2026circulatory'            # circulatory_autogen (libcuflynx) software citation (existing key circulatory_autogen)
CITATION_RULES = [
    # ---------------------------------------------------------------- cell: ion channels
    (r'cell/ion_channels/[^/]+/Paci2013_v01$', ['paci2013computational'], False, 'Paci 2013 hiPSC-CM model, channel as published.'),
    (r'cell/ion_channels/i_Na/Tao2011_v01$', ['tao2011model', 'belluzzi1986quantitative'], False,
     'Tao 2011 uses the Belluzzi & Sacchi 1986 Na current.'),
    (r'cell/ion_channels/i_Ca/Tao2011_v01$', ['tao2011model', 'belluzzi1989calcium'], False, 'Tao 2011 Eqs. 91-94 = Belluzzi & Sacchi 1989.'),
    (r'cell/ion_channels/i_KV/Tao2011_v01$', ['tao2011model', 'belluzzi1988interactions'], False, 'Tao 2011 Eqs. 96-97 after Belluzzi & Sacchi 1988.'),
    (r'cell/ion_channels/i_A/Tao2011_v01$', ['tao2011model', 'belluzzi1985fast'], True,
     'Tao 2011 A-current; Belluzzi, Sacchi & Wanke 1985 is the likely primary (not checked equation by equation).'),
    (r'cell/ion_channels/i_KCa/Tao2011_v01$', ['tao2011model'], True, 'Tao 2011 KCa; primary source not traced.'),
    (r'cell/ion_channels/i_Na/Argus2026_v01$', ['belluzzi1986quantitative', 'tao2011model', ARGUS_SN], False,
     'Belluzzi & Sacchi 1986 kinetics via Tao 2011 (owner example), modified by the owner.'),
    (r'cell/ion_channels/i_Na1_6/Argus2026_v01$', ['belluzzi1986quantitative', 'tao2011model', ARGUS_SN], True,
     'Same law as i_Na (BS1986/Tao); Herzog 2003 supplies NaV1.6 values only.'),
    (r'cell/ion_channels/i_CaL/Argus2026_v01$', ['belluzzi1989calcium', 'tao2011model', ARGUS_SN], False,
     'Belluzzi & Sacchi 1989 gating via Tao 2011 (GHK driving force added by the owner).'),
    (r'cell/ion_channels/i_CaN/Argus2026_v01$', ['belluzzi1989calcium', ARGUS_SN], True,
     'Gating shared with i_CaL (BS1989). The original CellML says "fit to data from Huang1998", which is untraced.'),
    (r'cell/ion_channels/i_Kv2/Argus2026_v01$', ['belluzzi1988interactions', 'tao2011model', ARGUS_SN], False,
     'Belluzzi & Sacchi 1988 delayed rectifier via Tao 2011.'),
    (r'cell/ion_channels/i_Kv1_5/Argus2026_v01$', ['philipson1991sequence', 'ranjanchannelpedia', ARGUS_SN], True,
     'Channelpedia Kv1.5 model 21 (R. Ranjan), after Philipson 1991; whether to cite Channelpedia (no year) is open.'),
    (r'cell/ion_channels/i_Kv4_2/Argus2026_v01$', ['bekkers2000properties', ARGUS_SN], False,
     'Bekkers 2000 (owner example), via Channelpedia Kv4.2 model 40.'),
    (r'cell/ion_channels/i_BK/Argus2026_v01$', ['khaliq2003contribution', ARGUS_SN], False, 'Khaliq 2003 BK kinetics.'),
    (r'cell/ion_channels/i_SK/Argus2026_v01$', ['hirschberg1998gating', ARGUS_SN], True,
     'Hirschberg 1998 tau and Hill coefficient; the original CA token cited Gupta & Manchanda 2016 (unread).'),
    (r'cell/ion_channels/i_M/Argus2026_v01$', ['martin2023modelling', ARGUS_SN], True,
     'Martin & Pedersen (owner example "2023"): the bib holds both the 2023 bioRxiv and the 2024 PLOS Comput Biol article; '
     'their w equations cite Zhou et al. 2018.'),
    (r'cell/ion_channels/i_NaCa/Argus2026_v01$', ['paci2013computational', ARGUS_SN], False, 'Paci 2013 exchanger form.'),
    (r'cell/ion_channels/i_NaK/Argus2026_v01$', ['nygren1998mathematical', ARGUS_SN], True,
     'Nygren 1998 (owner example); the module comment says Lindblad 1996, the form Nygren inherited - add lindblad1996model?'),
    (r'cell/ion_channels/i_NaK/ArgusUNPUBLISHED_v01$', ['smith2004development', 'terkildsen2007balance', 'pan2020cardiac', ARGUS_UNP], True,
     'Four-state lumping as in Smith & Crampin 2004 / Terkildsen 2007, constants from Pan 2020; who built it ("Alfred_database") is untraced.'),
    (r'cell/ion_channels/i_leak_(Na|K)/Argus2026_v01$', [ARGUS_SN], False, "Owner's ohmic leak; no model paper."),
    (r'cell/ion_channels/IP3R/Argus2026_v01$', ['li1994equations', 'deyoung1992single', ARGUS_SN], True,
     'Li & Rinzel 1994 reduction of De Young & Keizer 1992; whether to cite DYK too is open.'),
    (r'cell/ion_channels/RyR/Argus2026_v01$', [ARGUS_SN], False, '"Ryanodine receptor model, developed by Finbar Argus" (CellML).'),
    (r'cell/ion_channels/RyR/HernandezCruz1997_v01$', ['hernandezcruz1997ca'], False, 'Hernandez-Cruz 1997 Appendix (owner example).'),
    (r'cell/ion_channels/i_piezo_Na/ArgusUNPUBLISHED_v01$', [ARGUS_UNP], False, "Owner's unpublished stretch-activated current."),
    # ---------------------------------------------------------------- cell: mechanisms and neuron parts
    (r'cell/[A-Za-z_]+/cardiomyocyte_Paci2013_v01$', ['paci2013computational'], False, 'Paci 2013.'),
    (r'cell/Ca_handling/SN_soma_Argus2026_v01$', [ARGUS_SN, 'khaliq2003contribution', 'shannon2004mathematical'], True,
     "Owner's four-compartment model; the shell is Khaliq 2003; the reversible SERCA form looks like Shannon 2004 (to confirm)."),
    (r'cell/Ca_handling/SN_varicosity_Argus2026_v01$', [ARGUS_SN, 'shannon2004mathematical'], True, 'As the soma, without the shell.'),
    (r'cell/membrane_potential/SN_Tao2011_v01$', ['tao2011model'], False, 'Tao 2011 membrane equation.'),
    (r'cell/neurotransmitter_release/SN_varicosity_Tao2011_v01$', ['tao2011model', 'schneggenburger2000intracellular'], True,
     'Tao 2011 Eqs. 107-113; the five-site Ca-binding scheme is Schneggenburger & Neher 2000 (add?).'),
    (r'cell/.*SN_(soma|varicosity)_Argus2026_v01$', [ARGUS_SN], False, "Owner's SN model part."),
    (r'cell/neuron/.*/(sympathetic|sympathetic_monolithic_v01)$|cell/neuron/sympathetic$', [ARGUS_SN], False,
     "Owner's SN model; a supermodule's report and exports take the union of its submodules' citations."),
    # ---------------------------------------------------------------- boundary conditions
    (r'boundary_conditions/i_stim_periodic/cardiomyocyte_Paci2013_v01$', ['paci2013computational'], False, 'Paci 2013 stimulus.'),
    (r'boundary_conditions/inlet_flow/nn_(adan|adan_2|aorticbif)$', ['boileau2015benchmark', CA_SOFT], True,
     'Prescribed waveform; the model is the Boileau 2015 benchmark inflow (a data source rather than a model paper).'),
    (r'boundary_conditions/', [CA_SOFT], True, 'Generic prescribed input; no model paper (Q1: software citation or none).'),
    # ---------------------------------------------------------------- control
    (r'control/(DMV|LCN|NA|NActr|PNDMV|PNNA|cardiopulmonary_receptor|lung_stretch_receptor)/ArgusUNPUBLISHED_(v01|OLD)$', ['gee2023closed', ARGUS_UNP], False,
     'Gee 2023 (owner example), modified by the owner.'),
    (r'control/NTS/ArgusUNPUBLISHED_v01$', ['gee2023closed', 'park2020investigating', ARGUS_UNP], True,
     'Gee 2023; the NTS notes cite Park 2020, which Gee builds on (bib entry to be added and checked).'),
    (r'control/(sympathetic_afferent|sympathetic_efferent|peripheral_resistance_effector|unstressed_volume_effector)/ArgusUNPUBLISHED_v01$',
     ['gee2023closed', 'ursino1998interaction', ARGUS_UNP], False, 'Gee 2023 built on the Ursino 1998 afferent/efferent/effector laws.'),
    (r'control/heart_effector/nn_(lv|rv)_ArgusUNPUBLISHED_(v01|OLD)$', ['gee2023closed', 'ursino1998interaction', ARGUS_UNP], False,
     'Gee 2023 heart effector (Ursino 1998 form).'),
    (r'control/(baroreceptor|efferent_heart_elastance_effector|efferent_heart_period_effector)/', ['ursino1998interaction'], False, 'Ursino 1998.'),
    (r'control/efferent_resistance_effector/Ursino1998_v01$', ['ursino1998interaction'], False, 'Ursino 1998.'),
    (r'control/afferent_to_(sympathetic|vagal)_efferent/', ['ursino1998interaction', 'ursino2000acute'], True,
     'Ursino 1998 efferent law with the chemoreceptor input of Ursino & Magosso 2000.'),
    (r'control/chemoreceptor/', ['ursino2000acute', 'spencer1979computational'], True, 'Ursino & Magosso 2000 chemoreceptor; Spencer 1979 blood-gas contents.'),
    (r'control/(efferent_resistance_effector|efferent_venous_us_volume_and_compliance_effector)/nn$', ['ursino1998interaction', ARGUS_UNP], True,
     "Owner's saturating effector (component *_Argus_type) replacing Ursino's log gain."),
    (r'control/efferent_respiratory_effector/GeeArgus2026_v01$', [GEE_ARGUS, 'albanese2016integrated'], True, 'Gee & Argus work in progress.'),
    (r'control/(efferent_respiratory_effector/nn|afferent_to_respiratory_efferent/)', ['albanese2016integrated', 'ursino2001integrated'], True,
     'Respiratory chemoreflex after Albanese 2016 / Ursino 2001 (to confirm).'),
    (r'control/local_flow_control_(brain|muscle)/', ['ursino2000acute', 'magosso2001mathematical'], True, 'Local metabolic control after Ursino & Magosso 2000.'),
    (r'control/', [CA_SOFT], True, 'Generic (observer, PID); no model paper.'),
    # ---------------------------------------------------------------- heart
    (r'heart/chamber/', ['liang2009multi'], False, 'Liang 2009 time-varying elastance chamber and timing.'),
    (r'heart/valve/(pp|pp_rmod)$', ['mynard2012simple'], False, 'Mynard 2012 valve.'),
    (r'heart/valve/pp_linear$', [CA_SOFT], True, 'Linear valve; no model paper.'),
    (r'heart/cardiac_clock/', ['liang2009multi', CA_SOFT], True, 'Phase timing as in Liang 2009.'),
    (r'heart/Argus2026_v01$', [ARGUS_UNP], True, 'Supermodule: its own entry plus the union of chamber/valve/clock citations.'),
    (r'heart/vp_simple(_SI|_OLD)?$', [CA_SOFT], True, 'Time-varying elastance with diode valves; generic textbook model (Suga-type), no single source.'),
    (r'heart/', ['liang2009multi', 'mynard2012simple'], False, 'Liang 2009 chambers with Mynard 2012 valves (notes).'),
    # ---------------------------------------------------------------- vessels, coupling, organs, respiratory, transport
    (r'vessels/properties/material_prop_visco', ['olufsen2000numerical', 'toro2022cerebrospinal', 'alastruey2011pulse'], True,
     'Empirical wave speed (Olufsen), tube law (Toro 2022 Eq. 84), Kelvin-Voigt viscoelasticity (Alastruey 2011?).'),
    (r'vessels/properties/', ['olufsen2000numerical', 'toro2022cerebrospinal'], True, 'Empirical wave speed and tube law.'),
    (r'vessels/terminals/', ['safaei2016roadmap', 'westerhof2009arterial'], True, 'Bond-graph Windkessel terminals (RCR after Westerhof).'),
    (r'vessels/microvasculature/', ['safaei2016roadmap', CA_SOFT], True, 'Lumped bond-graph micro-vessels; R, C given directly.'),
    (r'vessels/.*_nonlinear', ['safaei2016roadmap', 'toro2022cerebrospinal'], True, 'Bond-graph vessel with a nonlinear tube law.'),
    (r'vessels/', ['safaei2016roadmap', CA_SOFT], True, 'Lumped bond-graph vessels (Safaei 2016 roadmap) as implemented in circulatory_autogen.'),
    (r'coupling/FV1D_vessel/', ['muller2014global', CA_SOFT], True, '1D finite-volume solver coupling.'),
    (r'coupling/', [CA_SOFT], True, 'Generic coupling/sum; no model paper.'),
    (r'organs/', [ARGUS_UNP], True, "Owner's blood-volume-control bond graphs; no publication found."),
    (r'respiratory/lung_mechanics/ventilated_ArgusUNPUBLISHED_v01$', [ARGUS_UNP, 'albanese2016integrated'], True, "Owner's (Argus & Chase) ventilated lung."),
    (r'respiratory/lung_mechanics/ArgusUNPUBLISHED_v01$', ['albanese2016integrated', ARGUS_UNP], False, 'Albanese 2016 with owner changes.'),
    (r'respiratory/tissue_GE_simple/', ['ursino2001integrated', 'spencer1979computational'], True, 'Tissue gas exchange after Ursino 2001.'),
    (r'respiratory/', ['albanese2016integrated', 'spencer1979computational'], True, 'Albanese 2016 gas exchange; Spencer 1979 blood-gas contents.'),
    (r'transport/', ['fang2008oxygen', 'secomb2020mass'], True, 'Oxygen advection-diffusion (Fang 2008) with the Secomb 2020 mass-transfer coefficient.'),
    # ---------------------------------------------------------------- benchmarks
    (r'benchmarks/FitzHugh_Nagumo/', ['fitzhugh1961impulses', 'nagumo1962active'], False, 'FitzHugh 1961, Nagumo 1962.'),
    (r'benchmarks/Lotka_Volterra/', ['lotka1925elements', 'volterra1926fluctuations'], False, 'Lotka 1925, Volterra 1926.'),
    (r'benchmarks/VanDerPol/', ['vanderpol1926relaxation', 'hairer1996solving'], False, 'Van der Pol 1926; the stiff VDPOL form of Hairer & Wanner 1996.'),
    (r'benchmarks/', [CA_SOFT], False, 'circulatory_autogen test fixture.'),
]

# Papers needed by the proposal that are in no .bib yet (or placeholders). 'verify' marks details to check.
NEW_BIB = {
    ARGUS_SN: dict(author='Argus, Finbar and others', year='2026', title='Sympathetic neuron model (title to be decided)',
                   journal='in preparation', placeholder=True),
    ARGUS_UNP: dict(author='Argus, Finbar', year='2026', title='Unpublished modules of the circulatory-autogen module library',
                    journal='unpublished (placeholder until a publication exists)', placeholder=True),
    GEE_ARGUS: dict(author='Gee, Michelle M. and Argus, Finbar', year='2026', title='Unpublished respiratory control work (title to be decided)',
                    journal='in preparation', placeholder=True),
    CA_SOFT: dict(author='Argus, Finbar and Maso Talou, Gonzalo D. and others', year='2026', title='circulatory_autogen (libcuflynx)',
                  journal='GitHub, https://github.com/physiomelinks/circulatory_autogen (existing key circulatory_autogen; year = version cited)',
                  placeholder=True),
    'ranjanchannelpedia': dict(author='Ranjan, R.', year='', title='Channelpedia Kv1.5 model 21 (Kv1.5_5), after Philipson et al. 1991',
                               journal='Channelpedia, EPFL (no year: key has none)', verify=True),
    'park2020investigating': dict(author='Park, James H. and Gorky, Jonathan and Ogunnaike, Babatunde and Vadigepalli, Rajanikanth and Schwaber, James S.',
                                  year='2020', title='Investigating the effects of brainstem neuronal adaptation on cardiovascular homeostasis',
                                  journal='Frontiers in Neuroscience 14:470', verify=True),
    'westerhof2009arterial': dict(author='Westerhof, Nico and Lankhaar, Jan-Willem and Westerhof, Berend E.', year='2009',
                                  title='The arterial Windkessel', journal='Medical & Biological Engineering & Computing 47(2):131-141', verify=True),
    'fitzhugh1961impulses': dict(author='FitzHugh, Richard', year='1961', title='Impulses and physiological states in theoretical models of nerve membrane',
                                 journal='Biophysical Journal 1(6):445-466'),
    'nagumo1962active': dict(author='Nagumo, J. and Arimoto, S. and Yoshizawa, S.', year='1962', title='An active pulse transmission line simulating nerve axon',
                             journal='Proceedings of the IRE 50(10):2061-2070'),
    'lotka1925elements': dict(author='Lotka, Alfred J.', year='1925', title='Elements of Physical Biology', journal='Williams & Wilkins, Baltimore'),
    'volterra1926fluctuations': dict(author='Volterra, Vito', year='1926', title='Fluctuations in the abundance of a species considered mathematically',
                                     journal='Nature 118:558-560'),
    'vanderpol1926relaxation': dict(author='van der Pol, Balthasar', year='1926', title='On "relaxation-oscillations"',
                                    journal='The London, Edinburgh, and Dublin Philosophical Magazine 2(11):978-992'),
}

STOPWORDS = {'a', 'an', 'the', 'on', 'of', 'in', 'for', 'to', 'and', 'at', 'by', 'with', 'from', 'is', 'are'}
PARTICLES_JOINED = True   # "ten Tusscher" -> tentusscher, "De Young" -> deyoung, "van der Pol" -> vanderpol

OPEN_QUESTIONS = [
    ('Q1', 'Citing the owner\'s original models',
     'Versions with no publication behind them (the SN parts, RyR Argus2026_v01, the leaks, i_piezo_Na, the organs, the owner-modified '
     'control modules) need something in required_citations. Proposed: placeholder @unpublished entries, '
     f'<code>{ARGUS_SN}</code> (the planned SN paper; replaced by the real key on publication), <code>{ARGUS_UNP}</code> (owner-modified, '
     f'no paper) and <code>{GEE_ARGUS}</code>; generic mechanisms (constant BCs, sums, observers, PID, vessels without a model paper) cite '
     f'the software, <code>{CA_SOFT}</code> (today\'s key <code>circulatory_autogen</code>). Alternatives: (b) cite the software for every '
     'owner-original version too; (c) allow an empty list when <code>creator</code> is set. Which?'),
    ('Q2', 'Key casing',
     'Proposed: all lowercase, <code>belluzzi1986quantitative</code> (Google Scholar style; no case questions for "ten Tusscher" or '
     '"De Young"; accents and hyphens dropped: <code>hernandezcruz1997ca</code>). Alternative: <code>Belluzzi1986Quantitative</code>. '
     'Collisions add the next significant title word. Which casing?'),
    ('Q3', 'Rekey every existing key?',
     'See the rekey section for the count. Recommended: yes, in one scripted commit (bib files, CSV data_reference prefixes, '
     'reference_proposals and also_cites, source_figures.json, obs_data reference_key), with the structure test extended to require '
     'the convention, because new required_citations keys and old parameter keys would otherwise mix two conventions in one .bib.'),
    ('Q4', 'Boundary-condition prefix',
     'See the BC-prefix section. Recommended: (A) a small libcuflynx fix (treat a neighbour as a vessel only when its subtype starts '
     'with pp/pv/vp/vv, and read port kinds only from such names) so mechanism-style version names are safe everywhere; until it is '
     'merged, (B) control, organ, respiratory and coupling versions that can sit next to a junction keep an <code>nn_</code> prefix '
     '(e.g. <code>nn_ArgusUNPUBLISHED_v01</code>). Either way, no non-vessel version may start with two letters ending in p or v '
     '(hence <code>nn_lv_</code>, not <code>lv_</code>). Cell versions are safe without a prefix. A or B (or both)?'),
    ('Q5', 'Paci Ca handling: version, not instance',
     'You expected "probably a cardiomyocyte instance for now". The equations differ from both SN versions (see the evidence), so '
     'rule 4 makes it the version <code>cardiomyocyte_Paci2013_v01</code>. Agree?'),
    ('Q6', 'Where the Paci reversal potentials go',
     'Proposed: a new mechanism type <code>cell/reversal_potentials</code>. Alternative: a version of <code>ion_concentrations</code> '
     '(the SN soma computes E_Na, E_K there).'),
    ('Q7', 'Tao 2011 NE release version name',
     'Rule 6 keeps the source name; rule 4 adds system and compartment. Proposed <code>SN_varicosity_Tao2011_v01</code> (beside '
     '<code>SN_varicosity_Argus2026_v01</code>); alternatives <code>SN_Tao2011_v01</code> or <code>Tao2011_v01</code>. Likewise the Tao '
     'membrane: <code>SN_Tao2011_v01</code>.'),
    ('Q8', 'Same maths under different port names (control)',
     'heart_effector lv/rv (and the OLD rv) have identical equations; they are separate versions only because their port_type names '
     'carry the chamber. Proposed now: versions <code>nn_lv_/nn_rv_ArgusUNPUBLISHED_v01</code> and <code>nn_rv_ArgusUNPUBLISHED_OLD</code> '
     '(the nn_ is required: see Q4). '
     'Later: chamber-neutral port_types would make lv/rv two instances of one version. Also: peripheral_resistance_effector, '
     'unstressed_volume_effector and the Ursino effectors share one law (Q11).'),
    ('Q9', 'Brainstem nuclei names',
     'DMV, LCN, NA, NActr, NTS, PNDMV, PNNA name anatomical nuclei, not mechanisms. Proposed: keep the acronyms (as approved on '
     '2026-10-02). Alternative: nest them as parts of a <code>control/brainstem</code> supermodule.'),
    ('Q10', 'heart vp_simple_2 and the other _2 version names',
     'The rule names _2 <em>types</em>. heart <code>vp_simple_2</code> is a version that is an SI duplicate, so it is proposed as '
     '<code>vp_simple_SI</code>. inlet_flow <code>nn_adan_2</code> (interpolated waveform) and <code>nn_constant_2</code> (two inflows) '
     'are not SI duplicates and are left as they are.'),
    ('Q11', 'Further merges (not in the table)',
     'The ambiguous-merges section lists same-meaning or identical-maths types found outside the proposal (Ursino effectors, '
     'kidney left/right, the identical vessel compartments arterial_simple/capillary_simple/venous/arteriole/venule/vein, junction '
     'twins, background Na current, Tao channels). Proposed: a second phase after this one, because most of them change the CVS system '
     'models\' port wiring. Which, if any, now?'),
    ('Q12', 'Version names that stay',
     'Kept unchanged: the ion-channel <code>Argus2026_v01</code> versions (named after the planned SN paper), the existing '
     '<code>nn</code>/<code>nn_*</code>/BC-prefixed vessel versions, the neuron parts\' <code>sympathetic</code> and '
     '<code>sympathetic_monolithic_v01</code>, and <code>heart/Argus2026_v01</code> (a supermodule: libcuflynx expands it before any '
     'BC_type check, so it needs no prefix). The existing <code>efferent_resistance_effector/nn</code> (owner\'s saturating form) could '
     'become <code>ArgusUNPUBLISHED_v01</code> now that the type has two versions; proposed: leave it, to avoid touching five system models.'),
]


# ============================================================================================ helpers
M = '{http://www.w3.org/1998/Math/MathML}'
_OPS = {'plus': '+', 'minus': '-', 'times': '*', 'divide': '/', 'power': '^', 'eq': '=', 'lt': '<', 'gt': '>', 'leq': '<=',
        'geq': '>=', 'and': ' and ', 'or': ' or ', 'neq': '!='}


def _tag(e):
    return e.tag.replace(M, '')


def _wrap(a):
    return a if re.fullmatch(r'[\w.\-]+|\(.*\)|\w+\(.*\)', a) else '(' + a + ')'


def _tx(e):
    t = _tag(e)
    if t == 'ci':
        return (e.text or '').strip()
    if t == 'cn':
        txt = (e.text or '').strip()
        if len(e) and _tag(e[0]) == 'sep':
            txt += 'e' + (e[0].tail or '').strip()
        return txt
    if t in ('true', 'false', 'pi', 'exponentiale', 'infinity'):
        return t
    if t == 'apply':
        op = _tag(e[0])
        args = [_tx(a) for a in e[1:] if _tag(a) not in ('bvar', 'degree', 'logbase')]
        if op == 'diff':
            return f'd({args[0]})/dt'
        if op in _OPS:
            if op == 'minus' and len(args) == 1:
                return f'-({args[0]})'
            if op == 'eq':
                return f'{args[0]} = {args[1]}'
            if op in ('plus', 'minus'):
                return '(' + _OPS[op].join(args) + ')'
            return _OPS[op].join(_wrap(a) for a in args)
        if op == 'root':
            return f'sqrt({args[0]})'
        return f'{op}({", ".join(args)})'
    if t == 'piecewise':
        parts = []
        for c in e:
            if _tag(c) == 'piece':
                parts.append(f'{_tx(c[0])} if {_tx(c[1])}')
            else:
                parts.append(f'{_tx(c[0])} otherwise')
        return '{' + '; '.join(parts) + '}'
    return f'<{t}>'


def cellml_equations(path):
    '''Every equation of every component of a CellML file, as text.'''
    if not path or not os.path.isfile(path):
        return []
    root = ET.parse(path).getroot()
    out = []
    for comp in root.iter():
        if comp.tag.endswith('}component'):
            for m in comp:
                if m.tag == M + 'math':
                    out += [_tx(a) for a in m]
    return out


def strip_latex(s):
    s = re.sub(r'\\[\'"`^~=.uvHcdbtk]\s*\{?([A-Za-z])\}?', r'\1', s)      # \'{a} -> a, \v{z} -> z
    s = s.replace('{', '').replace('}', '').replace('\\', '')
    return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode()


def parse_bib(path):
    '''{key: {type, fields}} for a .bib file (brace-balanced fields).'''
    txt = open(path, encoding='utf-8', errors='replace').read()
    out = {}
    for m in re.finditer(r'@(\w+)\s*\{\s*([^,\s]+)\s*,', txt):
        start, depth, i = m.end(), 1, m.end()
        while i < len(txt) and depth:
            depth += {'{': 1, '}': -1}.get(txt[i], 0)
            i += 1
        body = txt[start:i - 1]
        fields = {}
        for f in re.finditer(r'(\w+)\s*=\s*', body):
            j = f.end()
            if j < len(body) and body[j] == '{':
                d, k = 1, j + 1
                while k < len(body) and d:
                    d += {'{': 1, '}': -1}.get(body[k], 0)
                    k += 1
                val = body[j + 1:k - 1]
            elif j < len(body) and body[j] == '"':
                k = body.index('"', j + 1)
                val = body[j + 1:k]
            else:
                val = re.match(r'[^,\n]*', body[j:])[0]
            fields.setdefault(f[1].lower(), re.sub(r'\s+', ' ', val).strip())
        out[m[2]] = {'type': m[1].lower(), **fields}
    return out


def first_surname(author):
    first = re.split(r'\s+and\s+', author.strip())[0]
    if ',' in first:
        sur = first.split(',')[0]
    else:
        parts = first.split()
        sur = ' '.join(parts[-1:]) if parts else first
        if first.startswith('{') and first.endswith('}'):
            sur = first
    sur = strip_latex(sur)
    return re.sub(r'[^a-z]', '', sur.lower())


def first_title_word(title):
    words = strip_latex(title).split()
    for w in words:
        m = re.match(r'[A-Za-z]+', w)
        if m and m[0].lower() not in STOPWORDS:
            return m[0].lower()
    return ''


def make_key(author, year, title):
    return f'{first_surname(author)}{re.sub(r"[^0-9]", "", year or "")[:4]}{first_title_word(title)}'


def esc(s):
    return html.escape(str(s), quote=True)


def rel(p):
    return os.path.relpath(p, REPO)


def load_json(path):
    with open(path) as f:
        return json.load(f)


# ========================================================================================== library
class Version:
    def __init__(self, cfg_path):
        self.dir = os.path.dirname(cfg_path)
        self.name = os.path.basename(self.dir)
        self.type_dir = os.path.dirname(os.path.dirname(self.dir))
        self.type_path = os.path.relpath(self.type_dir, MODULES)
        self.type_name = os.path.basename(self.type_dir)
        self.cfg = load_json(cfg_path)[0]
        self.format = self.cfg.get('module_format', 'cellml')
        cellml = os.path.join(self.dir, f'{self.type_name}_{self.name}_modules.cellml')
        self.cellml = cellml if os.path.isfile(cellml) else None
        idir = os.path.join(self.dir, 'instances')
        self.instances = sorted(os.listdir(idir)) if os.path.isdir(idir) else []
        self.has_results = os.path.isdir(os.path.join(self.dir, 'results'))
        self.has_plots = os.path.isdir(os.path.join(self.dir, 'plots'))
        self.bib = {}
        for b in glob.glob(os.path.join(self.dir, '*_references*.bib')):
            for k, v in parse_bib(b).items():
                self.bib.setdefault(k, v)

    @property
    def key(self):
        return f'{self.type_path}/{self.name}'

    def ports(self):
        out = {}
        for d in ('entrance_ports', 'exit_ports', 'general_ports'):
            for p in self.cfg.get(d) or []:
                out[(d.split('_')[0], p['port_type'])] = tuple(p.get('variables') or [])
        return out

    def equations(self):
        return cellml_equations(self.cellml)


def load_library():
    versions = []
    for cfg in sorted(glob.glob(os.path.join(MODULES, '**', 'versions', '*', '*_modules_config.json'), recursive=True)):
        if any(rel(cfg).startswith(e) for e in EXCLUDED):
            continue
        versions.append(Version(cfg))
    return versions


def norm_equations(eqs, names=()):
    '''Equations with component-name prefixes removed, for comparing versions.'''
    return [e for e in eqs]


# ======================================================================================== analysis
def check_changes(lib):
    by_key = {v.key: v for v in lib}
    problems, rows = [], []
    for old_t, old_v, new_t, new_v, kind, group, note in CHANGES:
        v = by_key.get(f'{old_t}/{old_v}')
        if v is None:
            problems.append(f'{old_t}/{old_v} not found in modules/')
        rows.append(dict(old_t=old_t, old_v=old_v, new_t=new_t, new_v=new_v, kind=kind, group=group, note=note, version=v))
    return rows, problems


def post_rename_versions(lib, rows):
    '''[(new type path, new version, Version)] for the whole library after the proposal.'''
    m = {(r['old_t'], r['old_v']): (r['new_t'], r['new_v']) for r in rows}
    return [(m.get((v.type_path, v.name), (v.type_path, v.name)) + (v,)) for v in lib]


def group_members(post, group_type_path):
    return [(t, n, v) for t, n, v in post if t == group_type_path]


def compare_equations(members):
    '''Equation-set comparison across versions: shared by all, and per version the ones only it has.'''
    sets = {n: set(v.equations()) for _, n, v in members if v.cellml}
    if not sets:
        return None
    shared = set.intersection(*sets.values()) if len(sets) > 1 else set()
    only = {n: sorted(s - shared) for n, s in sets.items()}
    identical = len(sets) > 1 and all(not o for o in only.values())
    return dict(sets=sets, shared=sorted(shared), only=only, identical=identical)


def port_table(members):
    allp = collections.OrderedDict()
    per = {}
    for _, n, v in members:
        p = v.ports()
        per[n] = p
        for k in p:
            allp.setdefault(k, None)
    names = [n for _, n, _ in members]
    sets = [set(per[n]) for n in names]
    same = all(s == sets[0] for s in sets)
    return names, list(allp), per, same


def scan_users(lib, rows):
    '''Users of every changed (module_type, version) pair: vessel arrays, supermodule configs.'''
    pairs = {(os.path.basename(r['old_t']), r['old_v']): r for r in rows}
    users = collections.defaultdict(list)
    for f in sorted(glob.glob(os.path.join(REPO, 'system_models', '**', '*_vessel_array.json'), recursive=True)):
        if '/reference/' in f:      # circulatory_autogen originals, kept as they are
            continue
        try:
            recs = load_json(f)
        except Exception:  # noqa: BLE001
            continue
        for r in recs if isinstance(recs, list) else []:
            k = (r.get('module_type'), r.get('module_subtype'))
            if k in pairs:
                users[k].append(('vessel array', rel(f), r.get('name')))
    for v in lib:
        for s in v.cfg.get('submodules') or []:
            k = (s.get('module_type'), s.get('module_subtype'))
            if k in pairs:
                users[k].append(('supermodule', v.key, s.get('name')))
    return users


def grep_counts(names, roots, exts, skip_dirs=('results', 'node_modules', '__pycache__', 'venv', '.git')):
    rx = re.compile(r'(?<![A-Za-z0-9_])(' + '|'.join(re.escape(n) for n in sorted(names, key=len, reverse=True)) + r')(?![A-Za-z0-9_])')
    out = collections.Counter()
    files = collections.defaultdict(set)
    for root in roots:
        base = os.path.join(REPO, root)
        if os.path.isfile(base):
            walk = [(os.path.dirname(base), [], [os.path.basename(base)])]
        else:
            walk = os.walk(base)
        for d, dirs, fs in walk:
            dirs[:] = [x for x in dirs if x not in skip_dirs]
            if any(rel(os.path.join(d, '')).startswith(e) for e in EXCLUDED):
                continue
            for f in fs:
                if not f.endswith(exts):
                    continue
                p = os.path.join(d, f)
                try:
                    txt = open(p, encoding='utf-8', errors='ignore').read()
                except OSError:
                    continue
                for m in rx.findall(txt):
                    out[m] += 1
                    files[m].add(rel(p))
    return out, files


def libcuflynx_bc_lines():
    '''(file, line, text) of the BC_type checks that matter for version names.'''
    hits = []
    pats = [("generators/CVSCellMLGenerator.py", r"_BC!='nn'|_BC != 'nn'"),
            ("generators/CVSCellMLGenerator.py", r"BC_type\.startswith\('v'\)"),
            ("generators/CVSCellMLGenerator.py", r"startswith\(\('vp','pp'\)\)"),
            ("generators/CVSCellMLGenerator.py", r"if out_vessel_BC_type.startswith\('nn'\)|temp_main_vessel_BC_type.endswith|self.__check_input_output_modules\("),
            ("generators/CVSCppGenerator.py", r"BC_type\[1\] == 'v'"),
            ("parsers/ModelParsers.py", r'BC_type"\]=="nn"|startswith\("nn"\)'),
            ("utilities/supermodules.py", r"SUPERMODULE_FORMAT")]
    for f, p in pats:
        path = os.path.join(LIBCUFLYNX, f)
        if not os.path.isfile(path):
            continue
        for i, line in enumerate(open(path, encoding='utf-8', errors='ignore'), 1):
            if re.search(p, line):
                hits.append((f, i, line.strip()))
    return hits


def is_vessel_type(t, post):
    names = [n for tt, n, _ in post if tt == t]
    return bool(names) and all(n[:2] in ('pp', 'pv', 'vp', 'vv') for n in names)


def bc_risk(t, n, post):
    '''A non-vessel version name that libcuflynx would read as a port convention.'''
    if is_vessel_type(t, post) or n.startswith('nn') or n[:2] in ('pp', 'pv', 'vp', 'vv'):
        return False
    return n[:2].endswith(('p', 'v'))


def junction_neighbours():
    '''Non-vessel modules wired next to junctions/terminals in the system models today.'''
    cnt = collections.Counter()
    for f in glob.glob(os.path.join(REPO, 'system_models', '**', '*_vessel_array.json'), recursive=True):
        if '/reference/' in f:
            continue
        try:
            recs = load_json(f)
        except Exception:  # noqa: BLE001
            continue
        by = {r['name']: r for r in recs}
        for r in recs:
            mt = r['module_type']
            if 'junction' in mt or 'FV1D' in mt:
                for n in r.get('inp_instances', []) + r.get('out_instances', []):
                    o = by.get(n)
                    if o and o['module_subtype'][:2] not in ('vv', 'vp', 'pv', 'pp'):
                        cnt[(mt, o['module_type'], o['module_subtype'])] += 1
    return cnt


def identical_groups(lib):
    '''Groups of versions of different module types with identical equations (outside cell/).'''
    g = collections.defaultdict(list)
    for v in lib:
        if v.type_path.startswith('cell/') or not v.cellml:
            continue
        eqs = sorted(v.equations())
        if eqs:
            g[hashlib.md5('\n'.join(eqs).encode()).hexdigest()].append(v)
    return [vs for vs in g.values() if len({v.type_path for v in vs}) > 1]


# ========================================================================================== citations
def all_bib_entries(lib):
    out = {}
    for v in lib:
        for k, e in v.bib.items():
            out.setdefault(k, e)
    return out


def citations_for(key):
    for rx, keys, unc, basis in CITATION_RULES:
        if re.search(rx, key):
            return keys, unc, basis
    return [], True, 'no rule'


def citation_text(newkey, by_newkey):
    if newkey in NEW_BIB:
        e = NEW_BIB[newkey]
        tag = ' [placeholder]' if e.get('placeholder') else (' [details to verify]' if e.get('verify') else '')
        return f"{e['author']} ({e['year'] or 'n.d.'}). {e['title']}. {e['journal']}.{tag}", None
    if newkey in by_newkey:
        old, e = by_newkey[newkey][0]
        venue = e.get('journal') or e.get('booktitle') or e.get('publisher') or e.get('howpublished') or ''
        vol = e.get('volume', '')
        pages = e.get('pages', '')
        bits = ', '.join(x for x in (venue, vol and f'vol. {vol}', pages and f'pp. {pages}') if x)
        return f"{strip_latex(e.get('author', ''))} ({e.get('year', 'n.d.')}). {strip_latex(e.get('title', ''))}. {strip_latex(bits)}.", old
    return '(not in any .bib: entry to be added)', None


# ============================================================================================== HTML
CSS = '''
:root { --bg:#ffffff; --fg:#111827; --muted:#4b5563; --card:#f9fafb; --line:#e5e7eb; --accent:#b4440b; --ok:#e7f5ec; --warn:#fff4e0; --bad:#fdecec; --code:#f3f4f6; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; --ok:#16301f; --warn:#3a2c12; --bad:#3a1a1a; --code:#1f252d; } }
:root[data-theme="dark"] { --bg:#111418; --fg:#e5e7eb; --muted:#9ca3af; --card:#1a1f26; --line:#2d3540; --accent:#ff8a4c; --ok:#16301f; --warn:#3a2c12; --bad:#3a1a1a; --code:#1f252d; }
* { box-sizing: border-box; }
body { background: var(--bg); color: var(--fg); font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; margin: 0; padding: 0 16px 64px; }
main { max-width: 1180px; margin: 0 auto; }
h1 { font-size: 1.7rem; margin: 28px 0 4px; } h2 { margin-top: 40px; border-bottom: 1px solid var(--line); padding-bottom: 4px; }
h3 { margin-top: 26px; } .muted { color: var(--muted); } a { color: var(--accent); }
code, pre { background: var(--code); border-radius: 4px; font: 12.5px/1.45 ui-monospace, SFMono-Regular, Menlo, monospace; }
code { padding: 1px 4px; } pre { padding: 8px 10px; overflow-x: auto; white-space: pre-wrap; word-break: break-word; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13.5px; }
th, td { border: 1px solid var(--line); padding: 5px 7px; text-align: left; vertical-align: top; }
th { background: var(--card); position: sticky; top: 0; }
.scroll { overflow-x: auto; } .card { background: var(--card); border: 1px solid var(--line); border-radius: 8px; padding: 12px 16px; margin: 12px 0; }
.ok { background: var(--ok); } .warn { background: var(--warn); } .bad { background: var(--bad); }
.tag { display: inline-block; font-size: 11.5px; padding: 0 6px; border-radius: 9px; border: 1px solid var(--line); margin-left: 4px; }
details { margin: 6px 0; } summary { cursor: pointer; font-weight: 600; }
.toc a { margin-right: 14px; white-space: nowrap; } .stats span { display: inline-block; margin-right: 18px; }
.stats b { font-size: 1.25rem; } ol.rules > li { margin-bottom: 8px; }
'''


def html_table(head, rows, cls=None):
    out = ['<div class="scroll"><table>', '<tr>' + ''.join(f'<th>{h}</th>' for h in head) + '</tr>']
    for i, r in enumerate(rows):
        c = f' class="{cls[i]}"' if cls and cls[i] else ''
        out.append(f'<tr{c}>' + ''.join(f'<td>{x}</td>' for x in r) + '</tr>')
    out.append('</table></div>')
    return '\n'.join(out)


def eq_block(eqs, limit=45):
    if not eqs:
        return '<span class="muted">(no CellML: supermodule or C++)</span>'
    shown = eqs[:limit]
    more = f'\n… {len(eqs) - limit} more' if len(eqs) > limit else ''
    return '<pre>' + esc('\n'.join(shown) + more) + '</pre>'


RULE_TEXT = '''
<ol class="rules">
<li><b>Name by mechanism.</b> A module type is named after the mechanism it models (<code>Ca_handling</code>,
<code>membrane_potential</code>, <code>ion_concentrations</code>, <code>neurotransmitter_release</code>, <code>i_Na</code>, ...),
never after the cell, compartment, system or author it was built for. Type names never carry surnames or years.</li>
<li><b>Placement.</b> A type sits in the most general category where its mechanism is meaningful: Ca handling exists in every
cell, so <code>cell/Ca_handling</code>; ion channels, pumps and exchangers are <code>cell/ion_channels/&lt;channel&gt;</code>.</li>
<li><b>Nesting.</b> Nest a type only when it is an anatomical part of its parent: <code>cell/neuron/{soma,axon,varicosity}</code>,
<code>heart/{chamber,valve,cardiac_clock}</code>. Parts are supermodules whose versions are cell or organ types
(<code>soma/sympathetic</code>). No directory name appears twice in the tree, and no two sibling categories or types mean the same
thing (one <code>heart/</code>, no <code>cardiac/</code>).</li>
<li><b>Version or instance: compare the equations.</b>
<ul>
<li><b>Different maths is a version</b>, named <code>&lt;system&gt;_&lt;compartment&gt;_&lt;Source&gt;&lt;Year&gt;_vXX</code>
(<code>Ca_handling/versions/SN_soma_Argus2026_v01</code>), or <code>&lt;system&gt;_&lt;Source&gt;&lt;Year&gt;_vXX</code> without a
compartment (<code>membrane_potential/versions/cardiomyocyte_Paci2013_v01</code>), or <code>&lt;Source&gt;&lt;Year&gt;_vXX</code> where
the system adds nothing (ion channels: <code>i_Na/versions/Tao2011_v01</code>). <code>Source</code> is the first author's surname
in ASCII CamelCase (<code>HernandezCruz1997</code>); the owner's unpublished or modified models use <code>ArgusUNPUBLISHED</code>;
an SI-unit duplicate adds <code>_SI</code>; a superseded formulation ends <code>_OLD</code>.</li>
<li><b>Same maths with different parameter values is an instance</b> of the existing version
(<code>versions/&lt;v&gt;/instances/&lt;instance&gt;</code>).</li>
<li>Vessel-network modules keep their port prefix first (<code>vp_wCont_ASD</code>), because libcuflynx reads a module subtype's
first two letters as its inlet/outlet port kinds (see "Boundary-condition prefix" for the exact rule and for non-vessel modules).</li>
</ul></li>
<li><b>One mechanism across cell types.</b> Modules of the same mechanism built for different cells or compartments join one type
as versions or instances (rule 4). Versions of one type expose compatible ports where the mechanism allows; where they cannot
(e.g. a bond-graph pump beside current-law pumps), the type's report says which ports differ.</li>
<li><b>required_citations.</b> Every version's <code>modules_config.json</code> entry has <code>"required_citations": [bib keys]</code>:
the papers whose <em>model</em> the version's equations come from (not the per-parameter value sources, which stay in the CSV
<code>data_reference</code>). Each key is in the version's <code>references.bib</code> and follows
<code>&lt;surname&gt;&lt;year&gt;&lt;first significant title word&gt;</code>, all lowercase (<code>belluzzi1986quantitative</code>).
A supermodule lists its own and its report and exports add its submodules'.</li>
</ol>
'''


def build(lib):
    now = datetime.date.today().isoformat()
    rows, problems = check_changes(lib)
    post = post_rename_versions(lib, rows)
    users = scan_users(lib, rows)
    changed_types = {r['old_t'] for r in rows}
    all_types = sorted({v.type_path for v in lib})
    # a type changes when its path changes or any of its versions are renamed
    renamed_types = sorted({r['old_t'] for r in rows if r['old_t'] != r['new_t']})
    version_only = sorted({r['old_t'] for r in rows if r['old_t'] == r['new_t']})
    new_types = sorted({r['new_t'] for r in rows} - set(all_types))
    merged_into_existing = sorted({r['new_t'] for r in rows if r['new_t'] in all_types and r['old_t'] != r['new_t']})
    unchanged = [t for t in all_types if t not in changed_types and t not in merged_into_existing]
    remaining_dirs = sorted({'cell/cardiomyocytes'} if all(t.startswith('cell/cardiomyocytes/') is False or t in changed_types
                                                            for t in all_types) else set())
    # post-rename tree: duplicate directory names
    new_paths = sorted({t for t, _, _ in post})
    seg_count = collections.Counter()
    seg_where = collections.defaultdict(set)
    for t in new_paths:
        for i, s in enumerate(t.split('/')):
            seg_where[s].add('/'.join(t.split('/')[:i + 1]))
    dup_names = {s: sorted(p) for s, p in seg_where.items() if len(p) > 1}

    out = []
    A = out.append
    A('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">')
    A('<title>Library restructure proposal</title>')
    A(f'<style>{CSS}</style></head><body><main>')
    A('<h1>Library restructure proposal</h1>')
    A(f'<p class="muted">Placement, naming and required_citations under the rules of 2026-10-06. Proposal only: no module file has moved. '
      f'Generated {now} by <code>tools/library_restructure_proposal.py</code> from the library at this commit; regenerate after any change.</p>')

    n_cites = len(post)
    cite_rows = []
    by_newkey = collections.defaultdict(list)
    entries = all_bib_entries(lib)
    for k, e in entries.items():
        by_newkey[make_key(e.get('author', ''), e.get('year', ''), e.get('title', ''))].append((k, e))
    uncertain = 0
    for t, n, v in sorted(post, key=lambda x: (x[0], x[1])):
        keys, unc, basis = citations_for(f'{t}/{n}')
        uncertain += bool(unc)
        cite_rows.append((t, n, keys, unc, basis))

    kinds = collections.defaultdict(set)
    for r in rows:
        if r['old_t'] != r['new_t']:
            kinds[r['kind']].add(r['old_t'])
    A('<div class="card stats">'
      f'<span><b>{len(renamed_types)}</b> types change path</span>'
      f'<span><b>{len(kinds["merge"])}</b> merged into a mechanism type as versions</span>'
      f'<span><b>{len(kinds["rename"])}</b> renamed (surname dropped)</span>'
      f'<span><b>{len(kinds["si"]) + len(kinds["old"])}</b> _2/_OLD folded in</span>'
      f'<span><b>{len(new_types)}</b> new types</span>'
      f'<span><b>{len(version_only)}</b> types with a version rename only</span>'
      f'<span><b>{len(unchanged)}</b> types unchanged</span>'
      f'<span><b>{n_cites}</b> versions need required_citations</span>'
      f'<span><b>{uncertain}</b> flagged uncertain</span></div>')
    if problems:
        A('<div class="card bad"><b>Table problems:</b><ul>' + ''.join(f'<li>{esc(p)}</li>' for p in problems) + '</ul></div>')
    A('<p class="toc"><a href="#rules">1 Rules</a><a href="#table">2 Rename/merge table</a><a href="#evidence">2a Evidence</a>'
      '<a href="#ports">2b Ports</a><a href="#moves">2c What moves</a><a href="#unchanged">2d Unchanged</a>'
      '<a href="#bc">BC prefix</a><a href="#impact">3 Downstream</a><a href="#citations">4 required_citations</a>'
      '<a href="#keys">Keys and rekey</a><a href="#ambiguous">Ambiguous merges</a><a href="#questions">5 Open questions</a></p>')

    # ------------------------------------------------------------------------- 1 rules
    A('<h2 id="rules">1. Rule text for modules/README.md "Placement and naming"</h2>')
    A('<div class="card">' + RULE_TEXT + '</div>')
    A('<p class="muted">Superseded: the "most specific place covering its uses" rule, the SN_ part names, and the 2026-10-02 rename '
      'table where it conflicts (the SN parts were to become <code>soma_…</code>/<code>varicosity_…</code> types with version '
      '<code>Argus2026_v01</code>; they now join mechanism types as <code>SN_soma_/SN_varicosity_Argus2026_v01</code>; '
      '<code>lung_bg_Argus_Chase → lung/ArgusChase_v01</code> becomes <code>lung_mechanics/ventilated_ArgusUNPUBLISHED_v01</code>; '
      'Gee 2023 versions become <code>ArgusUNPUBLISHED</code> instead of <code>Gee2023_v01</code>; <code>Ca_dynamics</code>, '
      '<code>Na_dynamics</code>, <code>electric_potentials</code> become <code>Ca_handling</code>, <code>ion_concentrations</code>, '
      '<code>reversal_potentials</code>). Kept from it: surname-free control types, <code>Ursino1998</code> joining '
      '<code>efferent_resistance_effector</code>, <code>GeeArgus2026_v01</code>, _2 → _SI, _OLD → version.</p>')

    # ------------------------------------------------------------------------- 2 table
    A('<h2 id="table">2. Rename / merge table</h2>')
    A('<p>Every module type that changes. Kind: <b>merge</b> joins a mechanism type as a version; <b>rename</b> renames the type '
      '(surname removed) and its version; <b>si</b>/<b>old</b> folds a _2/_OLD type into its base type; <b>version</b> renames a '
      'version in place. Users: system models (vessel arrays) and supermodules that name the pair.</p>')
    trs, cls = [], []
    for r in rows:
        v = r['version']
        u = users.get((os.path.basename(r['old_t']), r['old_v']), [])
        ustr = '<br>'.join(esc(f'{k}: {p}' + (f' ({n})' if n else '')) for k, p, n in u[:6]) + (f'<br>… {len(u) - 6} more' if len(u) > 6 else '')
        trs.append([f'<code>{esc(r["old_t"])}</code> / <code>{esc(r["old_v"])}</code>',
                    f'<code>{esc(r["new_t"])}</code>', f'<code>{esc(r["new_v"])}</code>', esc(r['kind']),
                    (f'<a href="#ev-{esc(r["group"])}">evidence</a> ' if r['group'] in GROUP_NOTES else '') + esc(r['note']),
                    ustr or '<span class="muted">none</span>'])
        cls.append('bad' if v is None else None)
    A(html_table(['Current type / version', 'New type', 'New version', 'Kind', 'Why', 'Users to update'], trs, cls))
    A('<p>Category removed when empty: <code>cell/cardiomyocytes</code> (all four of its types move). New mechanism types: '
      + ', '.join(f'<code>{esc(t)}</code>' for t in new_types) + '. Existing types gaining versions: '
      + ', '.join(f'<code>{esc(t)}</code>' for t in merged_into_existing) + '.</p>')
    if dup_names:
        gen = {'versions', 'instances', 'risk', 'plots', 'results', 'reference'}
        real = {s: p for s, p in dup_names.items() if s not in gen}
        if real:
            A('<div class="card warn"><b>Directory names that appear twice after the rename:</b> '
              + '; '.join(f'<code>{esc(s)}</code> ({", ".join(esc(x) for x in p)})' for s, p in real.items()) + '</div>')
        else:
            A('<p class="muted">Checked: no directory name appears twice in the post-rename tree.</p>')

    # ------------------------------------------------------------------------- 2a evidence
    A('<h2 id="evidence">2a. Equation evidence for each merge</h2>')
    A('<p>Equations are read from each version\'s CellML (MathML printed as text). "Shared" are equations that appear verbatim '
      'in every version of the group; each version\'s column lists the rest. The verdict is the curated reading of the comparison.</p>')
    groups = []
    for r in rows:
        if r['group'] and r['group'] not in [g for g, _ in groups]:
            groups.append((r['group'], r['new_t']))
    si_groups = {r['group'] for r in rows if r['kind'] == 'si'}
    by_key_all = {v.key: v for v in lib}
    si_rows = []
    for r in rows:
        if r['kind'] == 'si' and r['version']:
            base = by_key_all.get(f"{r['new_t']}/{r['new_v'][:-3]}")
            eq_same = base is not None and sorted(base.equations()) == sorted(r['version'].equations())
            port_same = base is not None and base.ports() == r['version'].ports()
            si_rows.append([f"<code>{esc(r['old_t'])}/{esc(r['old_v'])}</code>",
                            f"<code>{esc(r['new_t'])}/{esc(r['new_v'][:-3])}</code>" if base else '—',
                            'identical' if eq_same else 'DIFFER', 'identical' if port_same else 'differ (units/variables)'])
    A('<h3 id="ev-si">_2 types folded in as _SI versions</h3>')
    A('<p>Each _2 version against the base version it duplicates (computed from the CellML and the configs):</p>')
    A(html_table(['_2 version', 'Base version', 'Equations', 'Ports'], si_rows,
                 ['ok' if x[2] == 'identical' else 'bad' for x in si_rows]))
    for g in si_groups:
        A(f'<span id="ev-{esc(g)}"></span>')
    for g, t in groups:
        if g in si_groups:
            continue
        members = sorted(group_members(post, t), key=lambda x: x[1])
        cmp = compare_equations(members)
        A(f'<h3 id="ev-{esc(g)}"><code>{esc(t)}</code></h3>')
        verdict = GROUP_NOTES.get(g, '')
        ident = cmp and cmp['identical']
        A(f'<div class="card {"ok" if ident else ""}"><b>Verdict:</b> {esc(verdict)}'
          + (' <span class="tag">computed: identical equations</span>' if ident else
             (' <span class="tag">computed: equations differ</span>' if cmp and len(cmp['sets']) > 1 else '')) + '</div>')
        cols = [n for _, n, _ in members]
        src = {n: v for _, n, v in members}
        if cmp and len(cmp['sets']) > 2:
            ns = list(cmp['sets'])
            A('<p class="muted">Equations in common, pairwise (of each row version\'s total): '
              + '; '.join(f'{esc(a)} &amp; {esc(b)}: {len(cmp["sets"][a] & cmp["sets"][b])} of {len(cmp["sets"][a])}/{len(cmp["sets"][b])}'
                          for i, a in enumerate(ns) for b in ns[i + 1:]) + '</p>')
        A('<details><summary>Equations side by side (' + ', '.join(esc(c) for c in cols) + ')</summary>')
        if cmp and cmp['shared']:
            A(f'<p><b>Shared by all ({len(cmp["shared"])}):</b></p>' + eq_block(cmp['shared'], 30))
        A(html_table([f'<code>{esc(c)}</code><br><span class="muted">{esc(src[c].key)}</span>' for c in cols],
                     [[eq_block(cmp['only'].get(c, []) if cmp else src[c].equations()) if src[c].cellml else eq_block([]) for c in cols]]))
        A('</details>')

    # ------------------------------------------------------------------------- 2b ports
    A('<h2 id="ports">2b. Port compatibility of the merged types</h2>')
    A('<p>Ports by (direction, port_type), with their variables, for every version of each type after the merge. '
      'A type is <b>compatible</b> when all its versions expose the same port types.</p>')
    for g, t in groups:
        members = sorted(group_members(post, t), key=lambda x: x[1])
        if len(members) < 2 or g in si_groups:
            continue
        names, ports, per, same = port_table(members)
        A(f'<details><summary><code>{esc(t)}</code>: {"compatible" if same else "ports differ"}</summary>')
        trs = []
        for p in ports:
            cells = []
            for n in names:
                vv = per[n].get(p)
                cells.append('—' if vv is None else '<code>' + esc(', '.join(vv)) + '</code>')
            trs.append([f'{esc(p[0])}: <code>{esc(p[1])}</code>'] + cells)
        A(html_table(['Port'] + [f'<code>{esc(n)}</code>' for n in names], trs))
        A('</details>')
    A('<div class="card warn"><b>Where ports cannot be made compatible:</b><ul>'
      '<li><code>i_NaK</code>: the bond-graph NKE cycle has chemostat and membrane-charge ports, not a current port; it cannot replace the '
      'current-law versions inside a cell without a bond-graph membrane.</li>'
      '<li><code>membrane_potential</code>, <code>ion_concentrations</code>: the Paci and Tao versions take one port per named current '
      '(fixed cell composition); the SN versions take one summed multi-port. Making them compatible means rewriting Paci/Tao to the '
      'summed form (a maths change, not part of this proposal).</li>'
      '<li><code>Ca_handling</code>: Paci uses A/F currents and the d gate; SN uses mol/s fluxes (Ca_flux) and ER release ports.</li>'
      '<li><code>i_stim_periodic</code>: A/F (Paci) vs nA current.</li>'
      '<li><code>heart_effector</code>: lv/rv differ only by port_type name (by design).</li>'
      '<li>The _SI versions are compared with their base versions in 2a (same ports, SI units).</li></ul></div>')

    # ------------------------------------------------------------------------- 2c moves
    A('<h2 id="moves">2c. What moves where</h2>')
    A('<p>Each version directory moves whole to <code>modules/&lt;new type&gt;/versions/&lt;new version&gt;/</code>; its files are '
      'renamed from <code>&lt;old type&gt;_&lt;old version&gt;_*</code> to <code>&lt;new type&gt;_&lt;new version&gt;_*</code> '
      '(config, CellML, units, tests, verification config, references.bib and references_proposed.bib, html); inside the config '
      '<code>module_type</code>, <code>module_subtype</code> and <code>component_file</code> change (the CellML component name may stay; '
      'renaming components that carry surnames, e.g. <code>NA_Gee2023</code>, is optional). Instances keep their names; '
      '<code>results/</code> and <code>plots/</code> are regenerated rather than moved.</p>')
    trs = []
    for r in rows:
        v = r['version']
        if not v:
            continue
        trs.append([f'<code>{esc(rel(v.dir))}</code>', f'<code>modules/{esc(r["new_t"])}/versions/{esc(r["new_v"])}</code>',
                    esc(', '.join(v.instances) or '—'), 'yes' if v.has_results else '—', 'yes' if v.has_plots else '—',
                    str(len(v.bib))])
    A('<details><summary>Per version: directories, instances, results, references</summary>'
      + html_table(['From', 'To', 'Instances (kept)', 'results/', 'plots/', 'bib entries'], trs) + '</details>')

    # ------------------------------------------------------------------------- 2d unchanged
    A(f'<h2 id="unchanged">2d. Unchanged types ({len(unchanged)})</h2>')
    bycat = collections.defaultdict(list)
    for t in unchanged:
        bycat[t.rsplit('/', 1)[0] if '/' in t else '(top level)'].append(t.rsplit('/', 1)[-1])
    A('<details><summary>Show all unchanged types by category</summary>'
      + html_table(['Category', 'Types'], [[f'<code>{esc(c)}</code>', ', '.join(f'<code>{esc(x)}</code>' for x in ts)]
                                           for c, ts in sorted(bycat.items())]) + '</details>')

    # ------------------------------------------------------------------------- BC prefix
    A('<h2 id="bc">Boundary-condition prefix: does it matter for the new version names?</h2>')
    hits = libcuflynx_bc_lines()
    A('<p>libcuflynx puts the module_subtype in the vessel array\'s <code>BC_type</code> column and reads its first letters in '
      'these places (from the editable install, <code>' + esc(LIBCUFLYNX) + '</code>):</p>')
    if hits:
        A(html_table(['File', 'Line', 'Code'], [[f'<code>{esc(f)}</code>', str(i), f'<code>{esc(s[:150])}</code>'] for f, i, s in hits]))
    else:
        A('<p class="muted">(libcuflynx source not found at generation time.)</p>')
    neigh = junction_neighbours()
    risky = sorted({(t, n) for t, n, v in post if bc_risk(t, n, post)})
    A('<div class="card ' + ('bad' if risky else 'ok') + '"><b>Computed check of every post-rename version name</b> '
      '(non-vessel types only; a name is at risk when it does not start with <code>nn</code> and its first two letters end in '
      '<code>p</code> or <code>v</code>): '
      + (', '.join(f'<code>{esc(t)}/{esc(n)}</code>' for t, n in risky) if risky else 'none at risk.') + '</div>')
    A('<div class="card"><b>What this means.</b><ul>'
      '<li><b>Every connection</b> (CVSCellMLGenerator <code>__check_input_output_modules</code>, called for each module pair): if '
      'neither side starts with <code>nn</code>, the first two letters are read as port kinds: a module whose <code>BC_type[:2]</code> '
      'ends in <code>v</code> must feed one starting with <code>p</code>, one ending in <code>p</code> must feed one starting with '
      '<code>v</code>, otherwise generation stops (<code>exit()</code>). So a version named <code>lv_…</code>, <code>rv_…</code>, '
      '<code>Up…</code> or <code>Sp…</code> breaks its connections; <code>SN_</code>, <code>ca</code>rdiomyocyte, <code>Ar</code>gus, '
      '<code>Pa</code>ci, <code>Ta</code>o, <code>He</code>rnandezCruz, <code>Ur</code>sino, <code>Ge</code>eArgus, <code>ve</code>ntilated '
      'pass silently. This is why the heart_effector versions keep <code>nn_</code>.</li>'
      '<li><b>Min_junction / Nout_junction / MinNout_junction</b> (CVSCellMLGenerator): every neighbour whose <code>BC_type[:2]</code> is '
      'not <code>nn</code> is treated as a vessel joined at the junction ("otherwise K_tube modules or other non-vessel modules will '
      'also be included"). A non-vessel module named e.g. <code>ArgusUNPUBLISHED_v01</code>, <code>SN_soma_…</code>, '
      '<code>cardiomyocyte_…</code>, <code>lv_…</code> wired to such a junction would be taken for a vessel. The schema note "exempts '
      'names starting nn" is therefore the whole rule: anything else counts.</li>'
      '<li><b>Terminal inputs</b>: a vessel with a terminal among its inputs and <code>BC_type.startswith(\'v\')</code> gets a venous '
      'connection. Only vessels have terminal inputs.</li>'
      '<li><b>1D coupling</b> (ModelParsers, CVSCppGenerator): inlet BCs next to a 1D vessel must start <code>nn</code>; heart vessel '
      'types are read by <code>BC_type[1]</code>.</li>'
      '<li><b>Supermodules</b> are expanded before any of this (utilities/supermodules.py), so <code>soma/sympathetic</code>, '
      '<code>heart/Argus2026_v01</code> etc. are never read.</li></ul>'
      '<p><b>Cell modules:</b> never wired to junctions, terminals or 1D vessels, so the proposed <code>SN_…</code>, '
      '<code>cardiomyocyte_…</code>, <code>ArgusUNPUBLISHED_v01</code> cell versions are safe (as the existing <code>Argus2026_v01</code> '
      'ion channels already are). <b>Control, organ, respiratory and coupling modules</b> can sit next to junctions: today only '
      '<code>nn</code>-prefixed ones do (below), so nothing breaks now, but a renamed control version wired to a Min/Nout junction '
      'later would.</p>'
      '<p><b>Safe naming, two options (Q4):</b> (A) fix libcuflynx: one helper <code>is_vessel_bc(bc) = bc[:2] in {"pp","pv","vp","vv"}</code> '
      'used by the junction neighbour tests, <code>__check_input_output_modules</code> (return unless both sides are vessel BCs) and the '
      '<code>startswith("nn")</code>/<code>=="nn"</code> checks in ModelParsers; then every version name is safe. (B) Without the fix, '
      'non-cell, non-vessel versions keep <code>nn_</code> first: <code>nn_ArgusUNPUBLISHED_v01</code>, <code>nn_lv_ArgusUNPUBLISHED_v01</code>, '
      '<code>nn_Ursino1998_v01</code>, <code>nn_GeeArgus2026_v01</code>, <code>nn_ventilated_ArgusUNPUBLISHED_v01</code> (the existing '
      '<code>nn_SI</code>, <code>nn_lv</code> style). Recommended: A, with B as the interim spelling if the rename lands first.</p></div>')
    A('<details><summary>Non-vessel modules next to junctions/1D vessels in the system models today</summary>'
      + html_table(['Junction type', 'Neighbour type', 'Neighbour subtype', 'Count'],
                   [[f'<code>{esc(a)}</code>', f'<code>{esc(b)}</code>', f'<code>{esc(c)}</code>', str(n)] for (a, b, c), n in sorted(neigh.items())])
      + '</details>')

    # ------------------------------------------------------------------------- 3 downstream
    A('<h2 id="impact">3. Downstream impact</h2>')
    names = sorted({os.path.basename(r['old_t']) for r in rows if r['old_t'] != r['new_t']} | {r['old_v'] for r in rows if r['kind'] == 'version'})
    code_counts, code_files = grep_counts(names, ['cam_testing', 'tests', 'tools', 'conftest.py', 'Makefile', '.github'], ('.py', '.yml', '.yaml', 'Makefile'))
    self_counts, _ = grep_counts(names, ['tools/library_restructure_proposal.py'], ('.py',))
    code_counts.subtract(self_counts)
    code_files = {k: {f for f in v if not f.endswith('library_restructure_proposal.py')} for k, v in code_files.items()}
    code_counts = collections.Counter({k: c for k, c in code_counts.items() if code_files.get(k) and c > 0})
    man_counts, man_files = grep_counts(names, ['manifests'], ('.json',))
    site_counts, site_files = grep_counts(names, ['site'], ('.html',))
    sys_counts, sys_files = grep_counts(names, ['system_models'], ('.json', '.yaml', '.csv'))
    n_va = sum(1 for u in users.values() for k, _, _ in u if k == 'vessel array')
    n_sm = sum(1 for u in users.values() for k, _, _ in u if k == 'supermodule')
    va_files = sorted({p for u in users.values() for k, p, _ in u if k == 'vessel array'})
    sm_files = sorted({p for u in users.values() for k, p, _ in u if k == 'supermodule'})
    impact = [
        ('System model vessel arrays', f'{n_va} records in {len(va_files)} files name a changed pair',
         '<br>'.join(f'<code>{esc(f)}</code>' for f in va_files) + '<br>Rewrite <code>module_type</code>/<code>module_subtype</code>; '
         'the <code>reference/</code> copies stay as they are (they are the CA originals). Other system_model files mentioning old names: '
         f'{len(set().union(*sys_files.values())) if sys_files else 0}.'),
        ('Supermodule configs', f'{n_sm} submodule entries in {len(sm_files)} versions',
         '<br>'.join(f'<code>{esc(f)}</code>' for f in sm_files) + '<br>Each submodule\'s <code>module_type</code>/<code>module_subtype</code> '
         'changes; submodule <em>names</em> (membrane, Ca, Na_K, NE) and so the generated variable prefixes stay the same, so '
         'instance parameter files, obs_data and calibration rows are unaffected. Descriptions that name the old types are updated.'),
        ('Old-name lookup', '<code>cam_testing/library.py legacy_renames()</code>',
         'Extend it with this table (old pair → new pair) so circulatory_autogen\'s own models and <code>tools/import_systems.py</code> / '
         '<code>cam_testing/vessel_array.py</code> keep resolving the old names; it reads <code>tools/restructure_map.yaml</code> + '
         '<code>HEART_VERSION_RENAMES</code> today.'),
        ('Tests and code', f'{sum(code_counts.values())} mentions in {len(set().union(*code_files.values())) if code_files else 0} files',
         '<br>'.join(f'<code>{esc(f)}</code>' for f in sorted(set().union(*code_files.values()))) if code_files else ''),
        ('Structure test and schema', '<code>tests/test_structure.py</code>, <code>modules/directory_schema.json</code>',
         'Add: no surname/year in type names; version-name pattern; required_citations present and in the version bib; key convention; '
         'no duplicate directory names. The schema\'s version_names note is corrected (only <code>nn</code> is exempt).'),
        ('Manifests (PhLynx)', f'{sum(man_counts.values())} mentions in {len(set().union(*man_files.values())) if man_files else 0} files',
         'Rebuild with <code>tools/build_manifests.py</code>. PhLynx loads <code>manifests/vitalworkshop.json</code> from <code>main</code> '
         'through jsDelivr, so it changes only when this reaches main.'),
        ('Report site', f'{sum(site_counts.values())} mentions in {len(set().union(*site_files.values())) if site_files else 0} pages',
         'Regenerate <code>site/</code> (module pages move with their types; old pages are deleted, not redirected).'),
        ('module_library_dirs users', 'cam_testing (system.py, checks.py, harness.py), README, the calibration handover, CUFLynx '
         'workflow_manager/settings (calibration-workflow worktree), PhLynx export bridge',
         'They pass the <code>modules/</code> root; libcuflynx collects every config below it and matches by (module_type, module_subtype), '
         'not by path. Directory moves need no change; renamed pairs do (in the workflows/system models that name them).'),
        ('CUFLynx per-instance OMEX', '<code>tools/build_instance_omex.py</code>, <code>cam_testing/omex.py</code>',
         'Archive names derive from the version key, so every moved version gets new OMEX names; rebuild all and drop the old ones. '
         'Known-failure lists keyed by version key (test_instance_omex) need the new keys. required_citations should be written into '
         'the OMEX metadata (the rule says exports carry them).'),
        ('Sister repo sympathetic_neuron', 'no functional reference',
         'Its models use circulatory_autogen\'s built-in names (<code>SN_soma</code>, <code>SN_axon</code>, <code>SN_varicosity</code>, '
         '<code>SN_soma_Argus2026</code>) through its own module_config_user, not this library\'s paths. Only prose in '
         '<code>SN_full/model_changes/steps.json</code> (review/sn-model-changes) names library types; later steps should use the new names.'),
        ('circulatory_autogen built-ins', 'none required',
         'libcuflynx\'s bundled modules keep their names; the library\'s renamed pairs do not collide with them. When '
         '<code>use_builtin_modules</code> is on, the legacy table maps old names.'),
    ]
    A(html_table(['Where', 'Scale', 'What changes'], [[esc(a), b, c] for a, b, c in impact]))

    # ------------------------------------------------------------------------- 4 citations
    A('<h2 id="citations">4. required_citations for every version</h2>')
    A(f'<p>{n_cites} versions after the rename; {uncertain} rows are flagged <b>uncertain</b> (amber). Keys are in the proposed convention; '
      'the full reference comes from the repo\'s .bib files (old key shown) or is new. Placeholders are the open question Q1. '
      'Rows of the same family are grouped; open a family to see every version.</p>')
    fam = collections.OrderedDict()
    for t, n, keys, unc, basis in cite_rows:
        f = t.split('/')[0] + ('/' + t.split('/')[1] if t.split('/')[0] in ('cell', 'vessels') and '/' in t else '')
        fam.setdefault(f, []).append((t, n, keys, unc, basis))
    used_keys = collections.OrderedDict()
    for f, items in fam.items():
        nu = sum(1 for x in items if x[3])
        A(f'<details{" open" if f.startswith("cell") else ""}><summary>{esc(f)}: {len(items)} versions'
          + (f' <span class="tag">{nu} uncertain</span>' if nu else '') + '</summary>')
        trs, cls = [], []
        for t, n, keys, unc, basis in items:
            for k in keys:
                used_keys.setdefault(k, None)
            trs.append([f'<code>{esc(t)}</code>', f'<code>{esc(n)}</code>', ', '.join(f'<code>{esc(k)}</code>' for k in keys), esc(basis)])
            cls.append('warn' if unc else None)
        A(html_table(['Type', 'Version', 'required_citations', 'Basis'], trs, cls))
        A('</details>')
    A('<h3>Full references of the proposed keys</h3>')
    trs, cls = [], []
    for k in used_keys:
        text, old = citation_text(k, by_newkey)
        trs.append([f'<code>{esc(k)}</code>', f'<code>{esc(old)}</code>' if old else ('<span class="muted">new</span>'), esc(text)])
        cls.append('warn' if (k in NEW_BIB and (NEW_BIB[k].get('verify') or NEW_BIB[k].get('placeholder'))) or text.startswith('(not') else None)
    A(html_table(['Key', 'Existing key', 'Reference'], trs, cls))
    missing = [k for k in used_keys if k not in NEW_BIB and k not in by_newkey]
    if missing:
        A('<div class="card bad"><b>Proposed keys with no entry anywhere:</b> ' + ', '.join(f'<code>{esc(k)}</code>' for k in missing) + '</div>')
    A('<div class="card"><b>Notes on specific citations.</b><ul>'
      '<li><code>i_M</code>: the owner\'s example says Martin &amp; Pedersen 2023. The repo has both the bioRxiv preprint '
      '(<code>MartinPedersen2023</code>, 2023) and the PLOS Comput Biol article (<code>Martin2024</code>, 2024, same title). '
      'Cite the published 2024 article (<code>martin2024modelling</code>) or keep 2023?</li>'
      '<li><code>i_NaK Argus2026_v01</code>: reviewed as the Nygren 1998 form; the CellML comment and the provenance trace say Lindblad 1996 '
      '(which Nygren adopted). Add <code>lindblad1996model</code>?</li>'
      '<li><code>i_CaN</code>: "fit to data from Huang1998" is untraced; closest found is Huang 1997 (JPET).</li>'
      '<li>Duplicated entries under two keys today (merge in the rekey): <code>Bennett1998</code>/<code>BennettCheungBrain1998</code>, '
      '<code>Hilgemann1991</code>/<code>HilgemannNicollPhilipson1991</code>, <code>Landowne1970</code>/<code>LandowneRitchie1970</code>, '
      '<code>Mynard2015</code>/<code>Mynard2015ABME</code>.</li>'
      '<li>Vessel, terminal and property families: the model papers are not stated anywhere in the modules (their bib entries are '
      'parameter sources). Safaei 2016 / circulatory_autogen is a proposal to confirm.</li></ul></div>')

    # ------------------------------------------------------------------------- keys and rekey
    A('<h2 id="keys">Key convention and rekeying all existing keys</h2>')
    rk = []
    newkeys = collections.Counter()
    for k, e in sorted(entries.items(), key=lambda x: x[0].lower()):
        nk = make_key(e.get('author', ''), e.get('year', ''), e.get('title', ''))
        newkeys[nk] += 1
        rk.append((k, nk, e))
    counts, files = grep_counts([k for k, _, _ in rk], ['modules', 'system_models', 'tests', 'cam_testing', 'reviews'],
                                ('.csv', '.yaml', '.json', '.bib', '.py'))
    total_refs = 0
    for k, _, _ in rk:
        total_refs += counts.get(k, 0)
    nfiles = len(set().union(*files.values())) if files else 0
    collisions = {nk for nk, c in newkeys.items() if c > 1}
    same = sum(1 for k, nk, _ in rk if k == nk)
    A(f'<p><b>{len(rk)}</b> distinct keys in the repo\'s .bib files today (references.bib and references_proposed.bib); '
      f'<b>{same}</b> already match the convention. Rekeying touches about <b>{total_refs}</b> occurrences in <b>{nfiles}</b> files '
      '(whole-word count over .bib, .csv data_reference, tests.yaml reference_proposals/also_cites, source_figures.json, obs_data '
      'reference_key, reviews/*.yaml, and code). What has to change:</p><ul>'
      '<li>every <code>*_references.bib</code> / <code>*_references_proposed.bib</code> entry key;</li>'
      '<li>CSV <code>data_reference</code> strings that start with <code>&lt;key&gt;;</code> (the report and the structure test parse the '
      'key from there: <code>cam_testing/bib.py reference_key</code>);</li>'
      '<li><code>reference_proposals</code> and <code>also_cites</code> in tests.yaml, obs_data <code>reference_key</code>, '
      '<code>source_figures.json</code> <code>source</code>, review yaml files;</li>'
      '<li>tests that assert on keys: none hard-code a key; <code>test_structure</code> checks that cited keys exist, so it guards '
      'the rewrite. Add a key-format check.</li></ul>'
      '<p>Recommended: yes, rekey everything in one scripted commit (old → new map written to <code>tools/bib_rekey_map.json</code> '
      'for traceability), merging the duplicate entries. Collisions under the convention get the next title word: '
      + (', '.join(f'<code>{esc(c)}</code>' for c in sorted(collisions)) if collisions else 'none') + '.</p>')
    A('<details><summary>Old key → proposed key (all entries)</summary>'
      + html_table(['Old key', 'Proposed', 'Uses', 'First author / year / title'],
                   [[f'<code>{esc(k)}</code>', f'<code>{esc(nk)}</code>' + (' <span class="tag">collision</span>' if nk in collisions else ''),
                     str(counts.get(k, 0)), esc(f"{strip_latex(e.get('author', ''))[:50]} / {e.get('year', '')} / {strip_latex(e.get('title', ''))[:80]}")]
                    for k, nk, e in rk]) + '</details>')

    # ------------------------------------------------------------------------- ambiguous
    A('<h2 id="ambiguous">Ambiguous merges (not in the table)</h2>')
    A('<p>Same-meaning or identical-maths module types found while checking the rules. None is proposed for this round (Q11).</p>')
    by_key = {v.key: v for v in lib}
    for title, keys, note in AMBIGUOUS:
        vs = [by_key[f'{a}/{b}'] for a, b in keys if f'{a}/{b}' in by_key]
        members = [(v.type_path, v.key, v) for v in vs]
        cmp = compare_equations([(t, k, v) for t, k, v in members])
        A(f'<details><summary>{esc(title)}' + (' <span class="tag">identical equations</span>' if cmp and cmp['identical'] else '')
          + '</summary><p>' + esc(note) + '</p>')
        if cmp:
            if cmp['shared']:
                A(f'<p><b>Shared ({len(cmp["shared"])}):</b></p>' + eq_block(cmp['shared'], 20))
            A(html_table([f'<code>{esc(v.key)}</code>' for v in vs], [[eq_block(cmp['only'].get(v.key, []), 25) for v in vs]]))
        A('</details>')
    ig = identical_groups(lib)
    A('<details><summary>Every group of versions with identical equations across different types (computed, outside cell/)</summary>'
      + html_table(['Versions (type/version)'], [[', '.join(f'<code>{esc(v.type_path.split("/")[-1])}/{esc(v.name)}</code>' for v in g)] for g in ig])
      + '<p>Identical equations means identical variable names too; whether the ports also agree is the next check before any merge. '
        'The vessel compartment groups (arterial_simple, capillary_simple, venous, arteriole, venule, vein per port prefix) would '
        'become one type with instances under rule 4 and touch every CVS system model.</p></details>')

    # ------------------------------------------------------------------------- questions
    A('<h2 id="questions">5. Open questions</h2>')
    A(html_table(['', 'Question', 'Proposal'], [[q, esc(t), b] for q, t, b in OPEN_QUESTIONS]))
    A('<p class="muted">After approval: one commit per category (cell, BCs, control, respiratory, _SI/_OLD), each regenerating the '
      'affected reports, manifests and OMEX, running the structure and system tests, and confirming that every moved version generates '
      'the same model as before.</p>')
    A('</main></body></html>')
    return '\n'.join(out), dict(renamed=len(renamed_types), merged_types=len(kinds['merge']), renamed_only=len(kinds['rename']),
                                si_old=len(kinds['si']) + len(kinds['old']), version_only=len(version_only), new_types=len(new_types),
                                merges=len([r for r in rows if r['kind'] == 'merge']), unchanged=len(unchanged),
                                versions=n_cites, uncertain=uncertain, problems=problems, rekey_keys=len(rk),
                                rekey_refs=total_refs, rekey_files=nfiles)


def main():
    lib = load_library()
    page, stats = build(lib)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as f:
        f.write(page)
    print(rel(OUT))
    print(json.dumps(stats, indent=1))


if __name__ == '__main__':
    main()
