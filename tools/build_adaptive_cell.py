"""Builds the generic adaptive cell of AdaptOrDeity (module type cell/generic_adaptive_cell) and the modules it is made
of, from modules/cell/generic_adaptive_cell/parameters.json (every physiological value and its source; nothing is
hardcoded here: values in this script are either read from that file or derived from it by a stated calibration).

    python tools/build_adaptive_cell.py

Supermodules (one flat supermodule per level; a level extends another by adding submodules and links):

  cell/generic_adaptive_cell/empty_membrane   every general mechanism, no transport proteins
  cell/generic_adaptive_cell/fermenter        AdaptOrDeity's starting organism: empty_membrane plus its proteins

Module versions (equations as infix text, tools/cellml_writer.py):

  cell/energy_metabolism/adaptive_ArgusUNPUBLISHED_v01   bond-graph cytosol: glucose, lactate, ATP, ADP, Pi, glycogen;
                                                         glycolysis (scaled by its enzymes' expression), maintenance
                                                         ATPase; heat; H+ from ATP use; the ATP budget by reaction
  cell/membrane_potential/adaptive_ArgusUNPUBLISHED_v01  the membrane as a capacitor (disc geometry)
  cell/Ca_handling/adaptive_ArgusUNPUBLISHED_v01         cytosolic Ca with a fixed buffer (no ER)
  cell/H_handling/adaptive_ArgusUNPUBLISHED_v01          buffered cytosolic H+, pH and E_H
  cell/biosynthesis/adaptive_ArgusUNPUBLISHED_v01        RNA, ribosomes, lipids, damage and repair, growth to a daughter
  cell/protein_pool/ArgusUNPUBLISHED_v01                 one protein: expression, upkeep, build cost, footprint
  cell/protein_expression/adaptive_ArgusUNPUBLISHED_v01  protein synthesis serving the pools (capacity, ATP throttle)
  cell/ion_channels/i_leak_Ca/ArgusUNPUBLISHED_v01       Ca leak with its Nernst potential
  cell/ion_channels/i_leak_H/ArgusUNPUBLISHED_v01        H+ leak
  cell/ion_channels/i_H_ATPase/FoF1_ArgusUNPUBLISHED_v01 F1Fo ATPase pumping H+ out (10/3 H+ per ATP)
  cell/ion_channels/i_NaH/NhaA_ArgusUNPUBLISHED_v01      NhaA Na+/H+ antiporter (2 H+ in per Na+ out)
  cell/ion_channels/i_KH/KefC_ArgusUNPUBLISHED_v01       KefC K+/H+ exchanger (electroneutral)
  cell/ion_channels/i_leak_K/adaptive_ArgusUNPUBLISHED_v01  i_leak_K Argus2026_v01 times expression (the Trk channel)
  cell/ion_channels/i_NaCa/adaptive_ArgusUNPUBLISHED_v01    i_NaCa Argus2026_v01 times expression
  cell/ion_channels/i_NaK/ATP_ArgusUNPUBLISHED_v01       the Argus2026_v01 pump plus ATP use (eukaryotic; not in the
                                                         fermenter, kept for later bases)
  transport/cell_membrane_face/adaptive_ArgusUNPUBLISHED_v01  Glf glucose facilitator and lactate-/H+ symporter
  transport/glucose_carrier/ATP_ArgusUNPUBLISHED_v01     glucose uptake driven by n ATP equivalents (instances PTS: the
                                                         phosphotransferase system, 1 PEP; ABC: binding-protein importer,
                                                         2 ATP)
  transport/glucose_carrier/H_symport_ArgusUNPUBLISHED_v01  glucose/H+ symport (GalP-type MFS carrier, 1 H+)

and reuses, unchanged, the sympathetic-neuron maths for the Na/K concentrations (ion_concentrations
SN_soma_Argus2026_v01) and the Na leak (i_leak_Na Argus2026_v01), through the supermodules' parameters.

Calibration (all at T = 20 C, the temperature Lynch & Marinov normalise to): the fermenter's default instance is an
exact steady state at V = -105 mV and pH 7.6 with the target concentrations. Total ATP turnover equals Lynch & Marinov's
maintenance cost for the cell's volume; the H+-ATPase spends H_ATPase_share of it; of the protons it pumps, NhaA_share
return through NhaA (removing Na+), KefC_share through KefC (removing the K+ that Trk brings in) and the rest through
the H+ leak; the Na/Ca exchanger runs at NCX_share of NhaA's cycles. Conductances, densities and capacities follow from
the balances. The medium pH (7.3) matters: above about 7.8 NhaA can no longer export Na+ at -105 mV.
"""
import csv
import json
import math
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from module_writer import BIB, REPO, register_bib  # noqa: E402
from module_writer import write_version as _write_version  # noqa: E402


def write_version(category, module_type, component, variables, equations, config_extra, params, verification, tests,
                  bib_keys, version, extra_instances=None):
    params = sourced_rows(f'{module_type}/{version}', params)
    extra = {k: sourced_rows(f'{module_type}/{version}:{k}', v) for k, v in (extra_instances or {}).items()}
    return _write_version(category, module_type, component, variables, equations, config_extra, params, verification,
                          tests, bib_keys, version, extra or None)

V_NEW = 'adaptive_ArgusUNPUBLISHED_v01'   # modules of the generic adaptive cell
PH = 'unreferenced_placeholder'

# ------------------------------------------------------------------ parameters (modules/cell/generic_adaptive_cell/parameters.json)
PARAM_FILE = os.path.join(REPO, 'modules', 'cell', 'generic_adaptive_cell', 'parameters.json')
PARAMS = json.load(open(PARAM_FILE))['parameters']


class PVal(float):
    """A value read from parameters.json; it carries its key, so a module CSV row holding it takes the JSON's source.
    Arithmetic gives plain floats: anything computed from it is derived (calibrated) and documented as such."""
    key = None


def par(key):
    v = PVal(PARAMS[key]['value'])
    v.key = key
    return v


def par_ref(key):
    e = PARAMS[key]
    src = e['source'] if e['source'] != 'uniprot' else 'definitional'
    return f"{src}; {e['note']}" if src not in ('definitional',) else f"definitional; {e['note']}"


class DVal(float):
    """A value derived from parameters.json entries by a stated rule: its CSV row names the rule and every input's source."""
    key = None


def der(value, how, *keys):
    v = DVal(value)
    v.ref = 'derived: ' + how + ' (' + '; '.join(f"{k}: {PARAMS[k]['source']}, {PARAMS[k]['status']}" for k in keys) + ')'
    v.sourced = all(par_sourced(k) for k in keys)
    return v


def par_sourced(key):
    return PARAMS[key]['status'] in ('sourced',)


UNSOURCED = []   # module CSV rows whose value is neither from parameters.json nor derived (reported at the end)
DERIVED_PREFIXES = ('calibrated', 'definitional', 'no ', 'geometry', 'operating point', 'harness', 'Nernst potential',
                    'the fermenter cell\'s calibrated', 'steady state', 'the calibrated rest')


def sourced_rows(where, rows):
    """Rows (name, value, ref, sourced): a parameters.json value takes its source from the JSON; anything else must be
    derived (its ref says how) or a definition."""
    out = []
    for name, value, ref, sourced in rows:
        key = getattr(value, 'key', None)
        if key is not None:
            out.append((name, float(value), par_ref(key), par_sourced(key)))
            continue
        if isinstance(value, DVal):
            out.append((name, float(value), value.ref, value.sourced))
            continue
        if not str(ref).startswith(DERIVED_PREFIXES):
            UNSOURCED.append(f'{where}: {name} = {value!r} ({ref})')
        out.append((name, float(value) if isinstance(value, float) else value, ref, sourced))
    return out


R, F, NA = par('R'), par('F'), par('N_A')
E_CHARGE_NC = par('e_charge')
T = par('T')
R_CELL, H_CELL = par('R_cell'), par('h_cell')
A_MEM = 2 * math.pi * R_CELL * H_CELL          # um2: in-plane membrane only (top and bottom are insulated in the game)
VOL = math.pi * R_CELL ** 2 * H_CELL * 1e-18   # m3
C_M = par('c_m')
V_REST = par('V_rest')
NAO, KO, CAO = par('Nao'), par('Ko'), par('Cao')
NAI, KI, CAI = par('Nai'), par('Ki'), par('Cai')
PH_I, PH_O = par('pH_i'), par('pH_o')
H_I = der(1000.0 * 10 ** -PH_I, 'free H+ = 10^-pH_i M', 'pH_i')
H_O = der(1000.0 * 10 ** -PH_O, 'free H+ = 10^-pH_o M', 'pH_o')
BETA_H = par('beta_H')
N_H_ATPASE = par('n_H')
ATP, ADP, PI = par('ATP'), par('ADP'), par('Pi')
GLC_E, LAC_E = par('glc_e'), par('lac_e')
H_ATPASE_SHARE = par('H_ATPase_share')
NHA_SHARE, KEF_SHARE = par('NhaA_share'), par('KefC_share')
NCX_SHARE = par('NCX_share')
K_H_ACT = der(1000.0 * 10 ** -par('pH_act_H_ATPase'), 'H+ at the half-activation pH', 'pH_act_H_ATPase')
V_NAK = par('V_NaK')
PUMP_SHARE = par('NaK_share')
DG0_FERM, DG0_HYD = par('dG0_ferm'), par('dG0_hyd')
H_LAC_GLY, H_ATP = par('h_lac_gly'), par('h_ATP')
C_STD = par('c_std')
RT = R * T


def lynch_maintenance_mol_s(vol_m3):
    """Lynch & Marinov (2015): C_M = 0.39 V^0.88 x 1e9 ATP/cell/h at 20 C, V in um3 -> mol ATP/s."""
    return par('maintenance_coefficient') * (vol_m3 * 1e18) ** par('maintenance_exponent') / 3600.0 / NA


# Nernst potentials (mV)
E_NA = 1000 * RT / F * math.log(NAO / NAI)
E_K = 1000 * RT / F * math.log(KO / KI)
E_CA = 1000 * RT / (2 * F) * math.log(CAO / CAI)

E_H = 1000 * RT / F * math.log(H_O / H_I)

# ---- ATP budget and ion cycles at rest (mol/s per cell)
R_TOT = lynch_maintenance_mol_s(VOL)
nA = 1e9  # A -> nA
R_F = H_ATPASE_SHARE * R_TOT                 # ATP hydrolysed by the F1Fo H+-ATPase
H_PUMPED = N_H_ATPASE * R_F                  # H+ pumped out
R_NHA = NHA_SHARE * H_PUMPED / 2             # NhaA cycles (2 H+ in, 1 Na+ out)
R_KEF = KEF_SHARE * H_PUMPED                 # KefC cycles (1 H+ in, 1 K+ out)
J_HLEAK = (1 - NHA_SHARE - KEF_SHARE) * H_PUMPED   # H+ leak into the cell
N_CYC = NCX_SHARE * R_NHA                    # NCX cycles (3 Na+ in, 1 Ca2+ out)

# leak currents (nA, positive outward) and conductances (uS) from the ion balances; Trk is the K+ uptake channel
I_K = -R_KEF * F * nA                        # K+ enters by Trk what KefC exports
I_NA = -(R_NHA - 3 * N_CYC) * F * nA         # Na+ leaks in what NhaA removes, less what NCX brings in
I_CA = -2 * N_CYC * F * nA                   # Ca2+ leaks in what NCX removes
I_H = -J_HLEAK * F * nA
G_K, G_NA, G_CA, G_H = (I_K / (V_REST - E_K), I_NA / (V_REST - E_NA), I_CA / (V_REST - E_CA), I_H / (V_REST - E_H))
GAMMA_K, GAMMA_NA, GAMMA_CA = par('gamma_K'), par('gamma_Na'), par('gamma_Ca')
RHO_K, RHO_NA, RHO_CA = (G / (A_MEM * GAM * 1e-6) for G, GAM in ((G_K, GAMMA_K), (G_NA, GAMMA_NA), (G_CA, GAMMA_CA)))
GS_H = G_H / A_MEM                           # uS/um2: proton leak per membrane area
# electroneutral and electrogenic exchangers: linear flux-force laws, r = g dE / F (g in uS, dE in mV -> nA)
DRIVE_NHA = -(V_REST + E_NA - 2 * E_H)        # mV, > 0 drives Na+ out
DRIVE_KEF = E_H - E_K                         # mV, > 0 drives K+ out
assert DRIVE_NHA > 0 and DRIVE_KEF > 0, 'NhaA and KefC must run forward at the resting point'
GS_NHA = R_NHA * F * nA / DRIVE_NHA / A_MEM   # uS/um2
GS_KEF = R_KEF * F * nA / DRIVE_KEF / A_MEM

A_ATP_REST = -DG0_HYD - RT * math.log((ADP / C_STD) * (PI / C_STD) / (ATP / C_STD))
# F1Fo H+-ATPase in hydrolysis mode: k ATP/(K_m + ATP) H_i/(K_act + H_i) f_thermo
FATP = dict(r_H_ATPase=par('r_H_ATPase'), K_m_ATP_H=par('K_m_ATP_H'))
W_H_REST = N_H_ATPASE * F * (E_H - V_REST) * 1e-3
F_THERMO_H = 1 - math.exp(-(A_ATP_REST - W_H_REST) / RT)
assert F_THERMO_H > 0, 'the H+-ATPase would violate the second law at the resting point'
I_HATP = H_PUMPED * F * nA
h_factor = ATP / (FATP['K_m_ATP_H'] + ATP) * H_I / (K_H_ACT + H_I) * F_THERMO_H
RHO_HATP = I_HATP / h_factor / (A_MEM * FATP['r_H_ATPase'] * N_H_ATPASE * E_CHARGE_NC)

# Na/K pump (eukaryotic; not in the fermenter cell): Argus2026_v01 law times ATP/(K_m_ATP + ATP), calibrated at its
# own operating point (V_NAK, 25% of the turnover) for its default instance
NAK = {k: par(k) for k in ('Ko_mid_NaK', 'Nai_mid_NaK', 'Nai_power_NaK', 'r_NaK', 'K_m_ATP_NaK')}
P_CYC = PUMP_SHARE * R_TOT
W_ION_NAK = F * (3 * (E_NA - V_NAK) + 2 * (V_NAK - E_K)) * 1e-3
F_THERMO_NAK = 1 - math.exp(-(A_ATP_REST - W_ION_NAK) / RT)
pump_factor = ((V_NAK + 150) / (V_NAK + 200) / (1 + (NAK['Nai_mid_NaK'] / NAI) ** NAK['Nai_power_NaK'])
               * KO / (KO + NAK['Ko_mid_NaK']) * ATP / (NAK['K_m_ATP_NaK'] + ATP) * F_THERMO_NAK)
RHO_NAK = P_CYC * F * nA / pump_factor / (A_MEM * NAK['r_NaK'] * E_CHARGE_NC)

# Na/Ca exchanger: Argus2026_v01 law
NCX = {k: par(k) for k in ('KmCa', 'KmNai', 'Ksat', 'alpha_NaCa', 'gamma_NaCa', 'r_NaCa')}
phi = V_REST * 1e-3 * F / RT
ncx_factor = ((math.exp(NCX['gamma_NaCa'] * phi) * NAI ** 3 * CAO
               - math.exp((NCX['gamma_NaCa'] - 1) * phi) * NAO ** 3 * CAI * NCX['alpha_NaCa'])
              / ((NCX['KmNai'] ** 3 + NAO ** 3) * (NCX['KmCa'] + CAO)
                 * (1 + NCX['Ksat'] * math.exp((NCX['gamma_NaCa'] - 1) * phi))))
I_NCX = -N_CYC * F * nA                     # inward: 3 Na in, 1 Ca out per cycle
assert ncx_factor < 0, 'the exchanger must run forward (Ca out) at the resting point'
RHO_NCX = I_NCX / ncx_factor / (A_MEM * NCX['r_NaCa'] * E_CHARGE_NC)

# ---- metabolism (rates in mM/s = mol/m3/s)
Q10, T_REF, DT_Q10 = par('Q10'), T, par('dT_Q10')
V_GLY_SS = R_TOT / (2 * VOL)                # glycolysis makes all the ATP
V_HYD_SS = None                              # set below, after the protein costs
GLY = {k: par(k) for k in ('K_glc', 'K_ADP_gly', 'K_Pi_gly', 'K_lac', 'K_ATP_gly')}
HYD = {k: par(k) for k in ('K_ATP_hyd', 'K_ADP_hyd', 'K_Pi_hyd')}
GLYC = {k: par(k) for k in ('K_glc_store', 'K_ATP_store', 'C_glyc_max', 'K_glyc_mob', 'K_ADP_mob', 'K_Pi_mob', 'K_i_ATP_mob',
                            'dG0_store', 'h_glyc_poly')}
GLYC['V_store_max'] = der(par('glycogen_store_rate') * V_GLY_SS, 'glycogen_store_rate x resting glycolytic flux', 'glycogen_store_rate', 'maintenance_coefficient')
GLYC['V_mob'] = der(par('glycogen_mobilise_rate') * V_GLY_SS, 'glycogen_mobilise_rate x resting glycolytic flux', 'glycogen_mobilise_rate', 'maintenance_coefficient')
# Membrane carriers and glycolytic capacity. A growing fermenter runs glycolysis 10-30x faster than its maintenance flux
# (glucose uptake ~10-30 mmol/gDW/h), so the enzymes have spare capacity and the cell is uptake-limited: internal glucose
# sits low (GLC_I) and glycolysis has about 20x headroom at rest. Glf then works near its full driving force from 1 mM
# outside; growth needs richer food (packets) or more carriers. The lactate exporter runs at 4% of its capacity at rest.
K_GLC_T, K_LAC_T = par('K_glc_t'), par('K_lac_t')
GLC_I = par('glc_i')
J_GLC = V_GLY_SS * VOL
VMAX_GLC = J_GLC / (GLC_E / (K_GLC_T + GLC_E) - GLC_I / (K_GLC_T + GLC_I))
J_LAC = 2 * V_GLY_SS * VOL
VMAX_LAC = par('lactate_exporter_headroom') * J_LAC
# lactate leaves with a proton (lactate-/H+ symport): the carrier sees X = C_lac H / H_ref (lactic acid activity)
H_REF = par('H_ref')
_x_lac = 1.0 / par('lactate_exporter_headroom')                  # exporter saturation at rest
LAC_I = K_LAC_T * _x_lac / (1 - _x_lac) * H_REF / H_I

# ---- protein expression: ATP to keep (turn over) and to build each transport protein and the rest of the proteome
C_SYN, C_POLY, C_DEG, C_RECYCLE = par('c_syn'), par('c_poly'), par('c_deg'), par('c_recycle')
K_DEG = par('k_deg')
PROTEOME_PER_UM3, L_PROTEOME = par('proteome_per_um3'), par('L_proteome')
GLYCOLYTIC_SHARE = par('glycolytic_share')
# key: (label, amino acids per functional unit, copies at rest, how the copies were estimated)
PROTEINS = {
    'H_ATPase': ('F1Fo H+-ATPase', par('L_H_ATPase'), RHO_HATP * A_MEM, 'rho_H_ATPase x A_mem (calibrated density)'),
    'NhaA': ('NhaA', par('L_NhaA'), R_NHA * NA / par('turnover_NhaA'), 'cycles at rest / turnover_NhaA'),
    'KefC': ('KefC', par('L_KefC'), R_KEF * NA / par('turnover_KefC'), 'cycles at rest / turnover_KefC'),
    'Trk': ('TrkH2-TrkA4', par('L_Trk'), R_KEF * NA / par('turnover_Trk'),
            'K+ entry at rest / turnover_Trk (the channel law lumps conductance)'),
    'NCX': ('NCX_Mj-type exchanger', par('L_NCX'), RHO_NCX * A_MEM, 'rho_NaCa x A_mem (calibrated density)'),
    'Glf': ('Glf glucose facilitator', par('L_Glf'), VMAX_GLC * NA / par('turnover_Glf'), 'capacity / turnover_Glf'),
    'LacT': ('lactate-/H+ symporter (LldP-sized)', par('L_LacT'), VMAX_LAC * NA / par('turnover_LacT'),
             'capacity / turnover_LacT'),
    'Flagellum': ('flagellum (filament, hook, basal body, MotAB stators)', par('L_Flagellum'),
                  der(par('n_flagella_init'), 'the starting number of flagella', 'n_flagella_init'),
                  'one flagellum per unit of expression: m is the number of flagella'),
    'Glyc': ('glycolytic enzymes', par('L_glycolytic'),
             der(GLYCOLYTIC_SHARE * PROTEOME_PER_UM3 * VOL * 1e18, 'glycolytic_share x proteome_per_um3 x cell volume',
                 'glycolytic_share', 'proteome_per_um3'), 'glycolytic_share x proteome_per_um3 x cell volume'),
    'proteome': ('rest of the proteome (ribosomal, nucleoid and other proteins)', L_PROTEOME,
                 der((1 - GLYCOLYTIC_SHARE) * PROTEOME_PER_UM3 * VOL * 1e18,
                     '(1 - glycolytic_share) x proteome_per_um3 x cell volume', 'glycolytic_share', 'proteome_per_um3'),
                 'proteome_per_um3 x cell volume, less the glycolytic enzymes'),
}
C_MAINT, C_BUILD = C_POLY + C_DEG + C_RECYCLE, C_SYN + C_POLY
R_PROT = sum(n * L for _, L, n, _ in PROTEINS.values()) * C_MAINT * K_DEG / NA     # mol ATP/s at rest
# ---- biosynthesis: RNA, ribosomes, lipids, DNA, growth, damage and repair (ATP in mol/s, costs in ~P)
VOL_UM3 = VOL * 1e18
A_SURFACE = 2 * math.pi * R_CELL ** 2 + 2 * math.pi * R_CELL * H_CELL     # um2: the whole surface (lipids cover it all)
BIO = {k: par(k) for k in ('L_mRNA', 't_half_mRNA', 'L_rRNA', 'k_ribo', 'c_nt_recycle', 'c_nt_denovo', 'area_per_lipid',
                            'c_lipid', 'k_lipid', 'genome_nt', 'K_D', 'repair_cost_fraction', 'Q10_dmg')}
BIO['n_mRNA'] = der(par('mRNA_per_um3') * VOL_UM3, 'mRNA_per_um3 x cell volume', 'mRNA_per_um3')
BIO['n_ribo'] = der(par('ribosomes_per_um3') * VOL_UM3, 'ribosomes_per_um3 x cell volume', 'ribosomes_per_um3')
BIO['r_dmg'] = der(1.0 / par('damage_lifetime'), '1 / damage_lifetime', 'damage_lifetime')
BIO['r_rep_ratio'] = par('repair_capacity_ratio')
BIO['n_lipid'] = 2 * A_SURFACE / BIO['area_per_lipid']
J_RNA = BIO['n_mRNA'] * BIO['L_mRNA'] * BIO['c_nt_recycle'] * math.log(2) / BIO['t_half_mRNA'] / NA
J_RIBO = BIO['n_ribo'] * BIO['L_rRNA'] * BIO['c_nt_recycle'] * BIO['k_ribo'] / NA
J_LIPID = BIO['n_lipid'] * BIO['c_lipid'] * BIO['k_lipid'] / NA
M_BUILD_PROT = sum(n * L for _, L, n, _ in PROTEINS.values()) * C_BUILD / NA
M_BUILD_RNA = (BIO['n_mRNA'] * BIO['L_mRNA'] + BIO['n_ribo'] * BIO['L_rRNA']) * BIO['c_nt_denovo'] / NA
M_BUILD_LIPID = BIO['n_lipid'] * BIO['c_lipid'] / NA
M_BUILD_DNA = BIO['genome_nt'] * BIO['c_nt_denovo'] / NA
M_DAUGHTER = M_BUILD_PROT + M_BUILD_RNA + M_BUILD_LIPID + M_BUILD_DNA
C_REPAIR = der(BIO['repair_cost_fraction'] * M_DAUGHTER, 'repair_cost_fraction x the daughter\'s itemised cost', 'repair_cost_fraction')                   # mol ATP to repair a whole unit of damage
R_REP_MAX = der(BIO['r_rep_ratio'] * BIO['r_dmg'], 'repair_capacity_ratio / damage_lifetime', 'repair_capacity_ratio', 'damage_lifetime')
# steady damage at full repair: R_rep_max D/(K_D + D) = r_dmg
D_REST = BIO['K_D'] * BIO['r_dmg'] / (R_REP_MAX - BIO['r_dmg'])
J_REPAIR = BIO['r_dmg'] * C_REPAIR
R_BIO = J_RNA + J_RIBO + J_LIPID + J_REPAIR
V_HYD_SS = (R_TOT - R_F - R_PROT - R_BIO) / VOL   # the maintenance ATPase: what the H+-ATPase, proteins and biosynthesis don't use
assert V_HYD_SS > 0, 'the itemised costs exceed the Lynch & Marinov total'


def keq(dg0):
    return math.exp(-dg0 / RT)


KEQ_GLY = keq(DG0_FERM - 2 * DG0_HYD)
KEQ_HYD = keq(DG0_HYD)


def cm_gly(glc, lac, atp, adp, pi):
    g = lambda c, k: c / k
    s = g(glc, GLY['K_glc']) * g(adp, GLY['K_ADP_gly']) ** 2 * g(pi, GLY['K_Pi_gly']) ** 2
    gamma = (lac / C_STD) ** 2 * (atp / C_STD) ** 2 / ((glc / C_STD) * (adp / C_STD) ** 2 * (pi / C_STD) ** 2)
    d = ((1 + g(glc, GLY['K_glc'])) * (1 + g(adp, GLY['K_ADP_gly'])) ** 2 * (1 + g(pi, GLY['K_Pi_gly'])) ** 2
         + (1 + g(lac, GLY['K_lac'])) ** 2 * (1 + g(atp, GLY['K_ATP_gly'])) ** 2 - 1)
    return s * (1 - gamma / KEQ_GLY) / d


def cm_hyd(atp, adp, pi):
    s = atp / HYD['K_ATP_hyd']
    gamma = (adp / C_STD) * (pi / C_STD) / (atp / C_STD)
    d = (1 + s) + (1 + adp / HYD['K_ADP_hyd']) * (1 + pi / HYD['K_Pi_hyd']) - 1
    return s * (1 - gamma / KEQ_HYD) / d


V_GLY = V_GLY_SS / cm_gly(GLC_I, LAC_I, ATP, ADP, PI)
GLY_HEADROOM = V_GLY * cm_gly(1e3, LAC_I, ATP, ADP, PI) / V_GLY_SS   # glycolysis if glucose inside were saturating
V_HYD = V_HYD_SS / cm_hyd(ATP, ADP, PI)
H_GLC = 2 * H_LAC_GLY                        # J/mol relative to lactate (lactate is this cell's end product)
A_HYD = RT * math.log(KEQ_HYD / ((ADP / C_STD) * (PI / C_STD) / (ATP / C_STD)))
Q_SS = VOL * V_GLY_SS * 2 * H_LAC_GLY       # W: at steady ATP, 63 kJ per mol lactate

print(f'cell: V {VOL * 1e18:.4g} um3, A_mem {A_MEM:.4g} um2, Cm {C_M * A_MEM * 0.01:.4g} pF')
print(f'Nernst: E_Na {E_NA:.2f}, E_K {E_K:.2f}, E_Ca {E_CA:.2f}, E_H {E_H:.2f} mV; rest {V_REST} mV, pH {PH_I} in, '
      f'{PH_O} out')
print(f'ATP turnover {R_TOT * NA:.4g} ATP/s ({R_TOT / VOL:.4g} mM/s); H+-ATPase {R_F * NA:.4g} ATP/s pumping '
      f'{H_PUMPED * NA:.4g} H+/s; NhaA {R_NHA * NA:.4g}, KefC {R_KEF * NA:.4g}, H+ leak {J_HLEAK * NA:.4g}, '
      f'NCX {N_CYC * NA:.4g} /s')
print(f'densities: Trk {RHO_K:.4g}, Na leak {RHO_NA:.4g}, Ca leak {RHO_CA:.4g}, H+-ATPase {RHO_HATP:.4g}, '
      f'NCX {RHO_NCX:.4g} /um2; g_H {GS_H:.4g}, g_NhaA {GS_NHA:.4g}, g_KefC {GS_KEF:.4g} uS/um2')
print(f'protein turnover {R_PROT * NA:.4g} ATP/s ({R_PROT / R_TOT:.1%} of the turnover): ' + ', '.join(f'{k} {n:.3g} x {L:.0f} aa' for k, (_, L, n, _) in PROTEINS.items()))
print(f'biosynthesis: mRNA {J_RNA * NA:.4g}, ribosomes {J_RIBO * NA:.4g}, lipids {J_LIPID * NA:.4g}, repair {J_REPAIR * NA:.4g} '
      f'ATP/s ({R_BIO / R_TOT:.1%}); daughter {M_DAUGHTER * NA:.4g} ATP (protein {M_BUILD_PROT / M_DAUGHTER:.1%}, RNA '
      f'{M_BUILD_RNA / M_DAUGHTER:.1%}, lipid {M_BUILD_LIPID / M_DAUGHTER:.1%}, DNA {M_BUILD_DNA / M_DAUGHTER:.2%}); Lynch & '
      f'Marinov growth cost {par("growth_cost_coefficient") * VOL_UM3 ** par("growth_cost_exponent"):.4g}; other maintenance '
      f'{V_HYD_SS * VOL / R_TOT:.1%}')
print(f'glycolytic headroom {GLY_HEADROOM:.3g}x the resting flux; Glf capacity from the background ({GLC_E} mM): '
      f'{VMAX_GLC * GLC_E / (K_GLC_T + GLC_E) / J_GLC:.3g}x, from a 10 mM packet: {VMAX_GLC * 10 / (K_GLC_T + 10) / J_GLC:.3g}x')
print(f'glycolysis {V_GLY_SS:.4g} mM/s (V_gly {V_GLY:.4g}), ATPase {V_HYD_SS:.4g} mM/s (V_hyd {V_HYD:.4g}); '
      f'glc_i {GLC_I:.4g} mM, lac_i {LAC_I:.4g} mM; heat {Q_SS:.4g} W ({Q_SS / VOL:.4g} W/m3)')
print(f'H+-ATPase thermodynamics: ATP affinity {A_HYD / 1e3:.2f} kJ/mol vs H+ work {W_H_REST / 1e3:.2f} kJ/mol per ATP; '
      f'thermodynamic factor at rest {F_THERMO_H:.4f}; NhaA drive {DRIVE_NHA:.1f} mV, KefC drive {DRIVE_KEF:.1f} mV')

# ------------------------------------------------------------------ references
register_bib({
    'meister1987proton': dict(type='article', author='Meister, Markus and Lowe, Gordon and Berg, Howard C.',
                              title='The proton flux through the bacterial flagellar motor', journal='Cell', year='1987',
                              volume='49', number='5', pages='643--650', doi='10.1016/0092-8674(87)90540-X',
                              note='About 1,200 protons per revolution of the motor'),
    'argus2026unpublished': dict(
        type='unpublished', author='Argus, Finbar', title='Unpublished modules of the circulatory-autogen module library',
        year='2026', note='Placeholder citation (required_citations, 2026-10-06): replace with the real reference when it '
                          'exists: the owner\'s unpublished or modified models, no publication'),
    'lynch2015bioenergetic': dict(type='article', author='Lynch, Michael and Marinov, Georgi K.',
                                  title='The bioenergetic costs of a gene',
                                  journal='Proceedings of the National Academy of Sciences', year='2015', volume='112',
                                  number='51', pages='15690--15695', doi='10.1073/pnas.1514974112',
                                  note='Basal maintenance cost C_M = 0.39 V^0.88 x 1e9 ATP/cell/h at 20 C (V in um3)'),
    'rolfe1997cellular': dict(type='article', author='Rolfe, D. F. S. and Brown, Guy C.',
                              title='Cellular energy utilization and molecular origin of standard metabolic rate in mammals',
                              journal='Physiological Reviews', year='1997', volume='77', number='3', pages='731--758',
                              note='Na+/K+-ATPase uses 19-28% of the ATP-coupled oxygen consumption in mammals'),
    'millero2008composition': dict(type='article',
                                   author='Millero, Frank J. and Feistel, Rainer and Wright, Daniel G. and McDougall, '
                                          'Trevor J.',
                                   title='The composition of Standard Seawater and the definition of the '
                                         'Reference-Composition Salinity Scale',
                                   journal='Deep Sea Research Part I: Oceanographic Research Papers', year='2008',
                                   volume='55', number='1', pages='50--72', doi='10.1016/j.dsr.2007.10.001',
                                   note='Reference composition at S = 35: Na 469.0, K 10.21, Ca 10.28 mmol/kg'),
    'gangola1987maintenance': dict(type='article', author='Gangola, Prem and Rosen, Barry P.',
                                   title='Maintenance of intracellular calcium in Escherichia coli',
                                   journal='Journal of Biological Chemistry', year='1987', volume='262', number='26',
                                   pages='12570--12574',
                                   note='Free cytosolic Ca2+ 90 +- 10 nM (fura-2), independent of external Ca'),
    'bennett2009absolute': dict(type='article',
                                author='Bennett, Bryson D. and Kimball, Elizabeth H. and Gao, Melissa and Osterhout, '
                                       'Robin and Van Dien, Stephen J. and Rabinowitz, Joshua D.',
                                title='Absolute metabolite concentrations and implied enzyme active site occupancy in '
                                      'Escherichia coli',
                                journal='Nature Chemical Biology', year='2009', volume='5', number='8', pages='593--599',
                                doi='10.1038/nchembio.186',
                                note='Glucose-fed, exponentially growing E. coli: ATP 9.6 mM'),
    'weisser1995functional': dict(type='article',
                                  author='Weisser, Peter and Kr{\\"a}mer, Reinhard and Sahm, Hermann and Sprenger, '
                                         'Georg A.',
                                  title='Functional expression of the glucose transporter of Zymomonas mobilis leads to '
                                        'restoration of glucose and fructose uptake in Escherichia coli mutants and '
                                        'provides evidence for its facilitator action',
                                  journal='Journal of Bacteriology', year='1995', volume='177', number='11',
                                  pages='3351--3354',
                                  note='Glf facilitated diffusion, Km 4.1 mM for glucose'),
    'liebermeister2010modular': dict(type='article',
                                     author='Liebermeister, Wolfram and Uhlendorf, Jannis and Klipp, Edda',
                                     title='Modular rate laws for enzymatic reactions: thermodynamics, elasticities and '
                                           'implementation',
                                     journal='Bioinformatics', year='2010', volume='26', number='12', pages='1528--1534',
                                     doi='10.1093/bioinformatics/btq141',
                                     note='Common modular rate law: thermodynamically consistent (Haldane) rate law for '
                                          'enzymes of unknown mechanism'),
    'gnaiger1990anaerobic': dict(type='article', author='Gnaiger, Erich and Kemp, Richard B.',
                                 title='Anaerobic metabolism in aerobic mammalian cells: information from the ratio of '
                                       'calorimetric heat flux and respirometric oxygen flux',
                                 journal='Biochimica et Biophysica Acta (BBA) - Bioenergetics', year='1990',
                                 volume='1016', number='3', pages='328--332',
                                 note='Net glucose to lactate in bicarbonate buffer -63 kJ/mol lactate'),
    'gajewski1986thermodynamics': dict(type='article',
                                       author='Gajewski, Ellen and Steckler, Donald K. and Goldberg, Robert N.',
                                       title='Thermodynamics of the hydrolysis of adenosine 5\'-triphosphate to adenosine '
                                             '5\'-diphosphate',
                                       journal='Journal of Biological Chemistry', year='1986', volume='261',
                                       number='27', pages='12733--12737',
                                       note='Delta H = -20.5 +- 0.4 kJ/mol at 25 C'),
    'milo2010bionumbers': dict(type='article',
                               author='Milo, Ron and Jorgensen, Paul and Moran, Uri and Weber, Griffin and Springer, '
                                      'Michael',
                               title='BioNumbers--the database of key numbers in molecular and cell biology',
                               journal='Nucleic Acids Research', year='2010', volume='38', number='suppl 1',
                               pages='D750--D753', doi='10.1093/nar/gkp889',
                               note='BNID 107114: E. coli cytosolic K+ 180-200 mM; BNID 103386: E. coli membrane '
                                    'potential -140 to -105 mV (fermentative growth at the less negative end)'),
    'slonczewski1981ph': dict(type='article', author='Slonczewski, J. L. and Rosen, B. P. and Alger, J. R. and Macnab, R. M.',
                              title='pH homeostasis in Escherichia coli: measurement by 31P nuclear magnetic resonance of '
                                    'methylphosphonate and phosphate',
                              journal='Proceedings of the National Academy of Sciences', year='1981', volume='78',
                              number='10', pages='6271--6275', doi='10.1073/pnas.78.10.6271',
                              note='Cytoplasmic pH = 7.6 + 0.1 (external pH - 7.6) over external pH 5.5-9'),
    'slonczewski1982effects': dict(type='article',
                                   author='Slonczewski, Joan L. and Macnab, Robert M. and Alger, Jeffry R. and Castle, '
                                          'Anna M.',
                                   title='Effects of pH and repellent tactic stimuli on protein methylation levels in '
                                         'Escherichia coli',
                                   journal='Journal of Bacteriology', year='1982', volume='152', number='1',
                                   pages='384--399', doi='10.1128/jb.152.1.384-399.1982',
                                   note='Cytoplasmic buffering capacity of E. coli about 50 mM per pH unit'),
    'jiang2001preferred': dict(type='article', author='Jiang, Weiping and Hermolin, Joe and Fillingame, Robert H.',
                               title='The preferred stoichiometry of c subunits in the rotary motor sector of '
                                     'Escherichia coli ATP synthase is 10',
                               journal='Proceedings of the National Academy of Sciences', year='2001', volume='98',
                               number='9', pages='4966--4971', doi='10.1073/pnas.081424898',
                               note='c10 ring with 3 catalytic beta subunits: 10/3 H+ per ATP'),
    'taglicht1993proton': dict(type='article', author='Taglicht, D. and Padan, E. and Schuldiner, S.',
                               title='Proton-sodium stoichiometry of NhaA, an electrogenic antiporter from Escherichia coli',
                               journal='Journal of Biological Chemistry', year='1993', volume='268', number='8',
                               pages='5382--5387', doi='10.1016/s0021-9258(18)53333-0',
                               note='2 H+ per Na+, independent of pH between 7.2 and 8.1'),
    'halevy2017geologic': dict(type='article', author='Halevy, I. and Bachan, A.',
                               title='The geologic history of seawater pH', journal='Science', year='2017', volume='355',
                               number='6329', pages='1069--1071', doi='10.1126/science.aal4151',
                               note='Seawater pH about 6.5-7.0 in the early Archean rising to 7.5-9.0 in the Phanerozoic'),
    'kobayashi1985proton': dict(type='article', author='Kobayashi, H.',
                                title='A proton-translocating ATPase regulates pH of the bacterial cytoplasm',
                                journal='Journal of Biological Chemistry', year='1985', volume='260', number='1',
                                pages='72--76', doi='10.1016/s0021-9258(18)89694-6',
                                note='Streptococcus faecalis (Enterococcus hirae): the F1Fo ATPase hydrolyses ATP made by '
                                     'fermentation to pump H+ out and regulate cytoplasmic pH'),
    'gulati2024structure': dict(type='article',
                                author='Gulati, Ashutosh and Kokane, Surabhi and Perez-Boerema, Annemarie and Alleva, '
                                       'Claudia and Meier, Pascal F. and Matsuoka, Rei and Drew, David',
                                title='Structure and mechanism of the K+/H+ exchanger KefC',
                                journal='Nature Communications', year='2024', volume='15', pages='4751',
                                doi='10.1038/s41467-024-49082-7',
                                note='KefC, the glutathione-gated K+/H+ exchanger of E. coli'),
    'bernstein2002global': dict(type='article',
                                author='Bernstein, Jonathan A. and Khodursky, Arkady B. and Lin, Pei-Hsun and Lin-Chao, Sue '
                                       'and Cohen, Stanley N.',
                                title='Global analysis of mRNA decay and abundance in Escherichia coli at single-gene '
                                      'resolution using two-color fluorescent DNA microarrays',
                                journal='Proceedings of the National Academy of Sciences', year='2002', volume='99',
                                number='15', pages='9697--9702', doi='10.1073/pnas.112318199',
                                note='About 80% of E. coli mRNA half-lives between 3 and 8 min'),
    'akashi2002metabolic': dict(type='article', author='Akashi, Hiroshi and Gojobori, Takashi',
                                title='Metabolic efficiency and amino acid composition in the proteomes of Escherichia coli '
                                      'and Bacillus subtilis',
                                journal='Proceedings of the National Academy of Sciences', year='2002', volume='99',
                                number='6', pages='3695--3700', doi='10.1073/pnas.062526999',
                                note='Biosynthetic cost 12-74 high-energy phosphate bonds per amino acid in E. coli'),
    'cao2011crystal': dict(type='article',
                           author='Cao, Yu and Jin, Xiangshu and Huang, Hua and Derebe, Mehabaw Getahun and Levin, Elena J. '
                                  'and Kabaleeswaran, Venkataraman and Pan, Yaping and Punta, Marco and Love, James and '
                                  'Weng, Jun and Quick, Matthias and Ye, Sheng and Kloss, Brian and Bruni, Renato and '
                                  'Martinez-Hackert, Erik and Hendrickson, Wayne A. and Rost, Burkhard and Javitch, '
                                  'Jonathan A. and Rajashankar, Kanagalaghatta R. and Jiang, Youxing and Zhou, Ming',
                           title='Crystal structure of a potassium ion transporter, TrkH', journal='Nature', year='2011',
                           volume='471', number='7338', pages='336--340', doi='10.1038/nature09731',
                           note='TrkH, the bacterial K+ uptake protein, has a channel-like selectivity filter'),
    'deutscher2006how': dict(type='article', author='Deutscher, Josef and Francke, Christof and Postma, Pieter W.',
                             title='How phosphotransferase system-related protein phosphorylation regulates carbohydrate '
                                   'metabolism in bacteria',
                             journal='Microbiology and Molecular Biology Reviews', year='2006', volume='70', number='4',
                             pages='939--1031', doi='10.1128/MMBR.00024-06',
                             note='The PTS: one PEP per sugar taken up and phosphorylated; glucose PTS Km in the '
                                  'micromolar range'),
    'davidson2008structure': dict(type='article',
                                  author='Davidson, Amy L. and Dassa, Elie and Orelle, Cedric and Chen, Jue',
                                  title='Structure, function, and evolution of bacterial ATP-binding cassette systems',
                                  journal='Microbiology and Molecular Biology Reviews', year='2008', volume='72',
                                  number='2', pages='317--364', doi='10.1128/MMBR.00031-07',
                                  note='Binding-protein-dependent ABC importers hydrolyse two ATP per transport cycle'),
    'pao1998major': dict(type='article', author='Pao, Stephanie S. and Paulsen, Ian T. and Saier, Milton H., Jr.',
                         title='Major facilitator superfamily', journal='Microbiology and Molecular Biology Reviews',
                         year='1998', volume='62', number='1', pages='1--34', doi='10.1128/MMBR.62.1.1-34.1998',
                         note='MFS uniporters, symporters and antiporters; sugar/H+ symporters (GalP, LacY)'),
    'ferenci1996adaptation': dict(type='article', author='Ferenci, Thomas',
                                  title='Adaptation to life at micromolar nutrient levels: the regulation of Escherichia '
                                        'coli glucose transport by endoinduction and cAMP',
                                  journal='FEMS Microbiology Reviews', year='1996', volume='18', number='4',
                                  pages='301--317', doi='10.1111/j.1574-6976.1996.tb00226.x',
                                  note='E. coli glucose uptake at micromolar glucose: the Mgl binding-protein system has '
                                       'sub-micromolar affinity'),
})



TEMPERATURE_SWEEP = {
    'ranges': {'T': [273.15, 313.15], 'T_ref': [283.15, 303.15]},
    'rationale': 'temperatures from 0 to 40 C (the default factors of 0.5 and 2 would put the cell at 147 K and 586 K, '
                 'where Q10 = 2 changes every rate by 2^15 to 2^29 and nothing physical remains to check)',
}

# ------------------------------------------------------------------ 1. bond-graph cytosol (energy_metabolism)
cyto_species = ['glc', 'lac', 'ATP', 'ADP', 'Pi']
cyto_vars = [('t', 'second', 'in', None, None)]
cyto_vars += [(f'C_{s}_init', 'millimolar', 'in', None, None) for s in cyto_species]
cyto_vars += [(f'C_{s}', 'millimolar', 'out', f'C_{s}_init',
               'species stores (bond-graph C elements): mu = mu0 + RT ln(C/c_std)' if s == 'glc' else None)
              for s in cyto_species]
cyto_vars += [(f'J_{s}_in', 'mol_per_s', 'in', None, 'summed fluxes of the faces on the entrance and exit cytosol_face_ports'
               if s == 'glc' else None) for s in cyto_species]
cyto_vars += [(f'J_{s}_out', 'mol_per_s', 'in', None, None) for s in cyto_species]
cyto_vars += [(f'j_{s}_ext', 'mol_per_s', 'in', None, 'summed ATP-consuming transporters (ATP_hydrolysis_port)'
               if s == 'ATP' else None) for s in ('ATP', 'ADP', 'Pi')]
cyto_vars += [
    ('V_c', 'm3', 'in', None, 'volume (compartment_volume port)'),
    ('T', 'kelvin', 'in', None, None), ('R', 'J_per_mol_K', 'in', None, None),
    ('T_ref', 'kelvin', 'in', None, None), ('dT_Q10', 'kelvin', 'in', None, None), ('Q10', 'dimensionless', 'in', None, None),
    ('c_std', 'millimolar', 'in', None, None),
    ('dG0_ferm', 'J_per_mol', 'in', None, 'standard Gibbs energies (glucose -> 2 lactate; ATP hydrolysis)'),
    ('dG0_hyd', 'J_per_mol', 'in', None, None),
    ('V_gly', 'millimolar_per_s', 'in', None, None), ('K_glc', 'millimolar', 'in', None, None),
    ('K_ADP_gly', 'millimolar', 'in', None, None), ('K_Pi_gly', 'millimolar', 'in', None, None),
    ('K_lac', 'millimolar', 'in', None, None), ('K_ATP_gly', 'millimolar', 'in', None, None),
    ('V_hyd', 'millimolar_per_s', 'in', None, None), ('K_ATP_hyd', 'millimolar', 'in', None, None),
    ('K_ADP_hyd', 'millimolar', 'in', None, None), ('K_Pi_hyd', 'millimolar', 'in', None, None),
    ('h_lac_gly', 'J_per_mol', 'in', None, 'enthalpies released (positive = heat out)'), ('h_ATP', 'J_per_mol', 'in', None, None),
    ('f_T', 'dimensionless', 'out', None, None),
    ('K_eq_gly', 'dimensionless', 'out', None, None), ('K_eq_hyd', 'dimensionless', 'out', None, None),
    ('Gamma_gly', 'dimensionless', 'out', None, 'mass-action ratios'), ('Gamma_hyd', 'dimensionless', 'out', None, None),
    ('v_gly', 'millimolar_per_s', 'out', None, 'glucose + 2 ADP + 2 Pi <=> 2 lactate + 2 ATP (common modular rate law)'),
    ('v_hyd', 'millimolar_per_s', 'out', None, 'ATP <=> ADP + Pi, the cell\'s maintenance work'),
    ('A_gly', 'J_per_mol', 'out', None, 'affinities A = RT ln(K_eq/Gamma)'), ('A_hyd', 'J_per_mol', 'out', None, None),
    ('P_diss', 'watt', 'out', None, 'free energy dissipated by the reactions (>= 0, second law)'),
    ('Q_met', 'watt', 'out', None, 'heat released: reaction enthalpies, including ATP hydrolysed by transporters'),
    ('J_ATP_turnover', 'mol_per_s', 'out', None, 'ATP hydrolysed by the ATPase and transporters'),
    ('j_H_met', 'mol_per_s', 'out', None, 'H+ released into the cytosol by all ATP hydrolysed, the maintenance ATPase and '
                                          'the ATP_hydrolysis_port users (H_production port): '
                                          'ATP4- + H2O -> ADP3- + HPO4 2- + H+'),
    # glycogen: a polymer store of glucose units (activity 1, so it does not enter the mass-action ratios)
    ('C_glyc_init', 'millimolar', 'in', None, None),
    ('C_glyc', 'millimolar', 'out', 'C_glyc_init', 'glycogen, as mM of glucose units'),
    ('store_alloc', 'dimensionless', 'in', None, 'fraction of the storage capacity in use (set by the game, 0-1)'),
    ('V_store_max', 'millimolar_per_s', 'in', None, None), ('K_glc_store', 'millimolar', 'in', None, None),
    ('K_ATP_store', 'millimolar', 'in', None, None), ('C_glyc_max', 'millimolar', 'in', None, 'storage ceiling'),
    ('V_mob', 'millimolar_per_s', 'in', None, None), ('K_glyc_mob', 'millimolar', 'in', None, None),
    ('K_ADP_mob', 'millimolar', 'in', None, None), ('K_Pi_mob', 'millimolar', 'in', None, None),
    ('K_i_ATP_mob', 'millimolar', 'in', None, 'ATP inhibition of mobilisation (phosphorylase)'),
    ('dG0_store', 'J_per_mol', 'in', None, 'glucose + 2 ATP -> glycogen unit + 2 ADP + 2 Pi'),
    ('h_glyc_poly', 'J_per_mol', 'in', None, 'enthalpy released by adding a glucose to glycogen'),
    ('K_eq_store', 'dimensionless', 'out', None, None), ('K_eq_mob', 'dimensionless', 'out', None, None),
    ('Gamma_store', 'dimensionless', 'out', None, None), ('Gamma_mob', 'dimensionless', 'out', None, None),
    ('v_store', 'millimolar_per_s', 'out', None, 'glycogen synthesis: glucose + 2 ATP -> glycogen + 2 ADP + 2 Pi'),
    ('v_mob', 'millimolar_per_s', 'out', None, 'glycogenolysis to lactate: glycogen + 3 ADP + 3 Pi -> 2 lactate + 3 ATP'),
    ('A_store', 'J_per_mol', 'out', None, None), ('A_mob', 'J_per_mol', 'out', None, None),
    ('m_gly', 'dimensionless', 'in', None, 'expression of the glycolytic enzymes, 1 = the calibrated cell (enzyme_expression '
                                          'port, from their protein_pool)'),
    ('J_ATP_gly', 'mol_per_s', 'out', None, 'ATP made by glycolysis (2 per glucose)'),
    ('J_ATP_mob', 'mol_per_s', 'out', None, 'ATP made from glycogen (3 per glucose unit)'),
    ('J_ATP_store', 'mol_per_s', 'out', None, 'ATP spent storing glucose as glycogen (2 per glucose)'),
    ('J_ATP_maint', 'mol_per_s', 'out', None, 'ATP hydrolysed by the maintenance ATPase (everything not itemised)'),
    ('M_glyc_ATP', 'mole', 'out', None, 'ATP the glycogen store would give if spent (3 per glucose unit; glycogen_store port)'),
]
cyto_eqs = [
    ('f_T', 'Q10 ** ((T - T_ref) / dT_Q10)'),
    ('K_eq_gly', 'exp(-(dG0_ferm - 2.0 * dG0_hyd) / (R * T))'),
    ('K_eq_hyd', 'exp(-dG0_hyd / (R * T))'),
    ('Gamma_gly', '((C_lac / c_std) ** 2.0 * (C_ATP / c_std) ** 2.0) / '
                  '((C_glc / c_std) * (C_ADP / c_std) ** 2.0 * (C_Pi / c_std) ** 2.0)'),
    ('Gamma_hyd', '((C_ADP / c_std) * (C_Pi / c_std)) / (C_ATP / c_std)'),
    ('v_gly', 'V_gly * m_gly * f_T * (C_glc / K_glc) * (C_ADP / K_ADP_gly) ** 2.0 * (C_Pi / K_Pi_gly) ** 2.0 '
              '* (1.0 - Gamma_gly / K_eq_gly) / ((1.0 + C_glc / K_glc) * (1.0 + C_ADP / K_ADP_gly) ** 2.0 '
              '* (1.0 + C_Pi / K_Pi_gly) ** 2.0 + (1.0 + C_lac / K_lac) ** 2.0 * (1.0 + C_ATP / K_ATP_gly) ** 2.0 - 1.0)'),
    ('v_hyd', 'V_hyd * f_T * (C_ATP / K_ATP_hyd) * (1.0 - Gamma_hyd / K_eq_hyd) / ((1.0 + C_ATP / K_ATP_hyd) '
              '+ (1.0 + C_ADP / K_ADP_hyd) * (1.0 + C_Pi / K_Pi_hyd) - 1.0)'),
    # ln of |Gamma| + 1e-300: a substrate that runs out sits at round-off (about 1e-14 mM, either sign) and its ratio
    # would make ln undefined; the rate law itself self-corrects (its sign follows 1 - Gamma/K_eq)
    ('A_gly', 'R * T * (ln(K_eq_gly) - ln(abs(Gamma_gly) + 1e-300))'),
    ('A_hyd', 'R * T * (ln(K_eq_hyd) - ln(abs(Gamma_hyd) + 1e-300))'),
    ('K_eq_store', 'exp(-dG0_store / (R * T))'),
    ('K_eq_mob', 'exp(-(dG0_ferm - dG0_store - dG0_hyd) / (R * T))'),
    ('Gamma_store', '((C_ADP / c_std) ** 2.0 * (C_Pi / c_std) ** 2.0) / ((C_glc / c_std) * (C_ATP / c_std) ** 2.0)'),
    ('Gamma_mob', '((C_lac / c_std) ** 2.0 * (C_ATP / c_std) ** 3.0) / ((C_ADP / c_std) ** 3.0 * (C_Pi / c_std) ** 3.0)'),
    # storage saturates near the ceiling (smooth, always >= 0 so its sign follows the affinity)
    ('v_store', 'V_store_max * store_alloc * f_T * C_glc / (K_glc_store + C_glc) * C_ATP / (K_ATP_store + C_ATP) '
                '/ (1.0 + (C_glyc / C_glyc_max) ** 8.0) * (1.0 - Gamma_store / K_eq_store)'),
    ('v_mob', 'V_mob * f_T * C_glyc / (K_glyc_mob + C_glyc) * C_ADP / (K_ADP_mob + C_ADP) * C_Pi / (K_Pi_mob + C_Pi) '
              '* K_i_ATP_mob / (K_i_ATP_mob + C_ATP) * (1.0 - Gamma_mob / K_eq_mob)'),
    ('A_store', 'R * T * (ln(K_eq_store) - ln(abs(Gamma_store) + 1e-300))'),
    ('A_mob', 'R * T * (ln(K_eq_mob) - ln(abs(Gamma_mob) + 1e-300))'),
    ('P_diss', 'V_c * (v_gly * A_gly + v_hyd * A_hyd + v_store * A_store + v_mob * A_mob)'),
    ('J_ATP_turnover', 'V_c * v_hyd - j_ATP_ext'),
    # one H+ per ATP hydrolysed (2 per glucose stored); glycogenolysis to lactate takes up one H+ (charge balance)
    ('j_H_met', 'V_c * (v_hyd + 2.0 * v_store - v_mob) - j_ATP_ext'),
    ('Q_met', 'V_c * (v_gly * (2.0 * h_lac_gly - 2.0 * h_ATP) + v_hyd * h_ATP + v_store * (2.0 * h_ATP + h_glyc_poly) '
              '+ v_mob * (2.0 * h_lac_gly - h_glyc_poly - 3.0 * h_ATP)) - j_ATP_ext * h_ATP'),
    ('d/dt C_glc', '(J_glc_in + J_glc_out) / V_c - v_gly - v_store'),
    ('d/dt C_lac', '(J_lac_in + J_lac_out) / V_c + 2.0 * v_gly + 2.0 * v_mob'),
    ('d/dt C_ATP', '(J_ATP_in + J_ATP_out + j_ATP_ext) / V_c + 2.0 * v_gly - v_hyd - 2.0 * v_store + 3.0 * v_mob'),
    ('d/dt C_ADP', '(J_ADP_in + J_ADP_out + j_ADP_ext) / V_c - 2.0 * v_gly + v_hyd + 2.0 * v_store - 3.0 * v_mob'),
    ('d/dt C_Pi', '(J_Pi_in + J_Pi_out + j_Pi_ext) / V_c - 2.0 * v_gly + v_hyd + 2.0 * v_store - 3.0 * v_mob'),
    ('d/dt C_glyc', 'v_store - v_mob'),
    # the ATP budget by reaction, from the stoichiometry above
    ('J_ATP_gly', '2.0 * V_c * v_gly'), ('J_ATP_mob', '3.0 * V_c * v_mob'), ('J_ATP_store', '2.0 * V_c * v_store'),
    ('J_ATP_maint', 'V_c * v_hyd'), ('M_glyc_ATP', '3.0 * V_c * C_glyc'),
]
cyto_roles = {v[0]: 'constant' for v in cyto_vars if v[0] != 't'}
cyto_roles.update({f'C_{s}': 'variable' for s in cyto_species})
cyto_roles.update({f'J_{s}_{d}': 'variable' for s in cyto_species for d in ('in', 'out')})
cyto_roles.update({f'j_{s}_ext': 'variable' for s in ('ATP', 'ADP', 'Pi')})
cyto_roles.update({n: 'variable' for n in ('f_T', 'K_eq_gly', 'K_eq_hyd', 'Gamma_gly', 'Gamma_hyd', 'v_gly', 'v_hyd',
                                            'A_gly', 'A_hyd', 'P_diss', 'Q_met', 'J_ATP_turnover', 'V_c', 'j_H_met',
                                            'C_glyc', 'K_eq_store', 'K_eq_mob', 'Gamma_store', 'Gamma_mob', 'v_store',
                                            'v_mob', 'A_store', 'A_mob', 'J_ATP_gly', 'J_ATP_mob', 'J_ATP_store',
                                            'J_ATP_maint', 'M_glyc_ATP')})
cyto_roles['store_alloc'] = 'boundary_condition'   # set by the game
cyto_roles['m_gly'] = 'boundary_condition'   # from the glycolytic enzymes' protein_pool
cyto_roles.update({n: 'global_constant' for n in ('T', 'R')})
cyto_roles['V_c'] = 'boundary_condition'   # a single port input, as the library's SN modules
face_vars_in = [f'C_{s}' for s in cyto_species]
cyto_ports = {
    'entrance_ports': [{'port_type': 'cytosol_face_port', 'variables': face_vars_in + [f'J_{s}_in' for s in cyto_species],
                        'multi_port': ['True'] * 5 + ['sum'] * 5}],
    'exit_ports': [{'port_type': 'cytosol_face_port', 'variables': face_vars_in + [f'J_{s}_out' for s in cyto_species],
                    'multi_port': ['True'] * 5 + ['sum'] * 5}],
    'general_ports': [
        {'port_type': 'ATP_hydrolysis_port', 'variables': ['C_ATP', 'C_ADP', 'C_Pi', 'j_ATP_ext', 'j_ADP_ext', 'j_Pi_ext'],
         'multi_port': ['True', 'True', 'True', 'sum', 'sum', 'sum']},
        {'port_type': 'compartment_volume', 'variables': ['V_c']},
        {'port_type': 'heat_source_port', 'variables': ['T', 'Q_met']},
        {'port_type': 'H_production', 'variables': ['j_H_met']},
        {'port_type': 'glycogen_store', 'variables': ['M_glyc_ATP'], 'multi_port': 'True'},
        # the glycolytic enzymes' pool (a port type of its own: the membrane face, linked to the cytosol, has
        # protein_expression ports)
        {'port_type': 'enzyme_expression', 'variables': ['m_gly']},
    ],
}
cyto_params = [
    ('C_glc_init', GLC_I, 'calibrated: steady state of the lumped cell (uptake by the membrane carrier = glycolysis)', False),
    ('C_lac_init', LAC_I, 'calibrated: steady state of the lumped cell (export by the membrane carrier = production)', False),
    ('C_ATP_init', ATP, 'bennett2009absolute; ATP 9.6 mM in glucose-fed, exponentially growing E. coli', True),
    ('C_ADP_init', ADP, f'{PH}; ADP of order 0.5 mM in E. coli (reported about 0.56 mM by Bennett et al. 2009; to confirm)',
     False),
    ('C_Pi_init', PI, f'{PH}; cytoplasmic Pi kept near 10 mM by the Pho regulon in E. coli (1-10 mM reported)', False),
    ('J_glc_in', 0.0, 'no face connected: closed volume', False), ('J_glc_out', 0.0, 'no face connected', False),
    ('J_lac_in', 0.0, 'no face connected', False), ('J_lac_out', 0.0, 'no face connected', False),
    ('J_ATP_in', 0.0, 'no face connected', False), ('J_ATP_out', 0.0, 'no face connected', False),
    ('J_ADP_in', 0.0, 'no face connected', False), ('J_ADP_out', 0.0, 'no face connected', False),
    ('J_Pi_in', 0.0, 'no face connected', False), ('J_Pi_out', 0.0, 'no face connected', False),
    ('j_ATP_ext', 0.0, 'no transporter connected', False), ('j_ADP_ext', 0.0, 'no transporter connected', False),
    ('j_Pi_ext', 0.0, 'no transporter connected', False),
    ('V_c', VOL, 'geometry: disc of radius 1 um, 2 um deep (the starting cell)', False),
    ('T', T, 'operating point: 20 C, the temperature Lynch & Marinov normalise maintenance costs to', False),
    ('R', R, 'definitional; gas constant, as the library\'s SN modules', True),
    ('T_ref', T_REF, 'reference temperature of the rate capacities (20 C)', False),
    ('dT_Q10', DT_Q10, 'definition of Q10 (rate ratio per 10 K)', True),
    ('Q10', Q10, f'{PH}; metabolic Q10 is typically 2-3; 2 chosen pending a source', False),
    ('c_std', C_STD, 'definitional; standard concentration 1 M for the thermodynamic factors', True),
    ('dG0_ferm', DG0_FERM, f'{PH}; textbook standard Gibbs energy of glucose -> 2 lactate (-196 kJ/mol; BioNumbers BNID '
                           '106539 lists -183.6)', False),
    ('dG0_hyd', DG0_HYD, f'{PH}; textbook standard Gibbs energy of ATP hydrolysis (-30.5 kJ/mol)', False),
    ('V_gly', V_GLY, 'calibrated: glycolysis supplies the total ATP turnover (lynch2015bioenergetic maintenance at 20 C) '
                     'at the default steady state', False),
    ('K_glc', GLY['K_glc'], f'{PH}; lumped glycolysis affinity for glucose (glucokinase order 0.1 mM)', False),
    ('K_ADP_gly', GLY['K_ADP_gly'], f'{PH}; lumped glycolysis affinity for ADP', False),
    ('K_Pi_gly', GLY['K_Pi_gly'], f'{PH}; lumped glycolysis affinity for Pi', False),
    ('K_lac', GLY['K_lac'], f'{PH}; lumped glycolysis product constant for lactate', False),
    ('K_ATP_gly', GLY['K_ATP_gly'], f'{PH}; lumped glycolysis product constant for ATP', False),
    ('V_hyd', V_HYD, 'calibrated: the maintenance ATPase uses the ATP turnover that neither the H+-ATPase (25%, '
                     'placeholder) nor protein turnover (protein_expression) uses', False),
    ('K_ATP_hyd', HYD['K_ATP_hyd'], f'{PH}; ATPase affinity for ATP, order 0.1-0.5 mM', False),
    ('K_ADP_hyd', HYD['K_ADP_hyd'], f'{PH}; ATPase product constant for ADP', False),
    ('K_Pi_hyd', HYD['K_Pi_hyd'], f'{PH}; ATPase product constant for Pi', False),
    ('h_lac_gly', H_LAC_GLY, 'gnaiger1990anaerobic; glucose to lactate in bicarbonate buffer -63 kJ/mol lactate', True),
    ('h_ATP', H_ATP, 'gajewski1986thermodynamics; Delta H of ATP hydrolysis -20.5 kJ/mol at 25 C', True),
    ('C_glyc_init', 0.0, 'definitional; the starting cell stores nothing', False),
    ('store_alloc', 0.0, 'definitional; storage off until the player allocates energy to it', False),
    ('V_store_max', GLYC['V_store_max'], f'{PH}; at full allocation the cell stores glucose at about half its resting '
                                         'glycolytic flux', False),
    ('K_glc_store', GLYC['K_glc_store'], f'{PH}; glucokinase-order affinity', False),
    ('K_ATP_store', GLYC['K_ATP_store'], f'{PH}; ADP-glucose pyrophosphorylase ATP affinity, order 1 mM', False),
    ('C_glyc_max', GLYC['C_glyc_max'], f'{PH}; about 10% of dry weight (cell 1.1 g/mL, 30% dry): 200 mM glucose units',
     False),
    ('V_mob', GLYC['V_mob'], f'{PH}; mobilisation capacity twice the resting glycolytic flux', False),
    ('K_glyc_mob', GLYC['K_glyc_mob'], f'{PH}; phosphorylase affinity for glycogen', False),
    ('K_ADP_mob', GLYC['K_ADP_mob'], f'{PH}; lumped affinity for ADP', False),
    ('K_Pi_mob', GLYC['K_Pi_mob'], f'{PH}; phosphorylase affinity for Pi, order 5 mM', False),
    ('K_i_ATP_mob', GLYC['K_i_ATP_mob'], f'{PH}; ATP inhibits mobilisation, so glycogen is spent when ATP falls', False),
    ('dG0_store', GLYC['dG0_store'], f'{PH}; textbook steps glucose -> G6P -> G1P -> UDP-glucose -> glycogen with PPi '
                                     'hydrolysis, about -42 kJ/mol', False),
    ('h_glyc_poly', GLYC['h_glyc_poly'], f'{PH}; polymerisation enthalpy taken as zero', False),
    ('m_gly', 1.0, 'definitional; baseline expression of the glycolytic enzymes (the calibrated cell)', False),
]
h_glc = 2 * H_LAC_GLY
cyto_invariants = [
    {'description': 'carbon balance with the membrane carriers: d/dt V_c (6 C_glc + 3 C_lac + 6 C_glyc) = 6 J_glc + '
                    '3 J_lac, integrated (relative 1e-6)',
     'expr': 'abs(V_c*((6*C_glc[-1] + 3*C_lac[-1] + 6*C_glyc[-1]) - (6*C_glc[0] + 3*C_lac[0] + 6*C_glyc[0])) - '
             'np.trapz(6*mf__J_glc + 3*mf__J_lac, t)) <= 1e-6*np.trapz(np.abs(6*mf__J_glc) + np.abs(3*mf__J_lac), t)'},
    {'description': 'adenine is conserved (bond-graph structure; ATP, ADP do not cross): C_ATP + C_ADP constant '
                    '(relative 1e-9)',
     'expr': 'np.ptp(C_ATP + C_ADP) <= 1e-9*np.max(C_ATP + C_ADP)'},
    {'description': 'phosphoryl groups are conserved: C_ATP + C_Pi constant (relative 1e-9)',
     'expr': 'np.ptp(C_ATP + C_Pi) <= 1e-9*np.max(C_ATP + C_Pi)'},
    {'description': f'first law: the enthalpy brought in with glucose equals the change in stored enthalpy plus the heat, '
                    f'h_glc int J_glc dt = V_c (h_glc dC_glc + h_ATP dC_ATP) + int Q_met dt, h_glc = 2 h_lac_gly = '
                    f'{h_glc:g} J/mol relative to lactate, ADP + Pi (relative 2e-5 of the heat: trapezoid integration of the output)',
     'expr': 'abs(2*h_lac_gly*np.trapz(mf__J_glc, t) - V_c*(2*h_lac_gly*(C_glc[-1] - C_glc[0]) + h_ATP*(C_ATP[-1] - '
             'C_ATP[0]) + (2*h_lac_gly - h_glyc_poly)*(C_glyc[-1] - C_glyc[0])) - np.trapz(Q_met, t)) <= '
             '2e-5*np.trapz(np.abs(Q_met), t)'},   # trapezoid integration of the logged heat over fast transients
    {'description': 'second law: each reaction runs down its affinity, v A >= 0, so P_diss >= 0',
     'expr': 'np.all(v_gly*A_gly >= 0) and np.all(v_hyd*A_hyd >= 0) and np.all(v_store*A_store >= 0) and '
             'np.all(v_mob*A_mob >= 0) and np.all(P_diss >= 0)'},
    {'description': 'the ATP budget by reaction follows the stoichiometry (exact): 2 per glycolysis, 3 per glycogen '
                    'unit mobilised, 2 per glucose stored, 1 per maintenance turnover; M_glyc_ATP = 3 V_c C_glyc',
     'expr': 'np.allclose(J_ATP_gly, 2*V_c*v_gly, rtol=1e-12, atol=0) and np.allclose(J_ATP_mob, 3*V_c*v_mob, rtol=1e-12, '
             'atol=0) and np.allclose(J_ATP_store, 2*V_c*v_store, rtol=1e-12, atol=0) and np.allclose(J_ATP_maint, '
             'V_c*v_hyd, rtol=1e-12, atol=0) and np.allclose(M_glyc_ATP, 3*V_c*C_glyc, rtol=1e-12, atol=0)'},
    {'description': 'glycogen stays between 0 and its ceiling (plus 5%)',
     'expr': 'np.all(C_glyc >= -1e-9) and np.all(C_glyc <= 1.05*C_glyc_max)'},
]
cyto_notes = (
    'The cytosol of a fermenting prokaryote-grade cell (AdaptOrDeity\'s starting organism) as a bond graph in concentration '
    'form: species stores (C elements, mu = mu0 + RT ln(C/c_std)) for glucose, lactate, ATP, ADP and Pi, joined by two '
    'reactions (R elements): homolactic glycolysis, glucose + 2 ADP + 2 Pi <=> 2 lactate + 2 ATP, and the maintenance '
    'ATPase, ATP <=> ADP + Pi; plus a glycogen store (a polymer of glucose units, activity 1): synthesis, glucose + 2 ATP '
    '-> glycogen + 2 ADP + 2 Pi (scaled by the game\'s store_alloc and slowed near the ceiling C_glyc_max), and glycogenolysis '
    'to lactate, glycogen + 3 ADP + 3 Pi -> 2 lactate + 3 ATP, inhibited by ATP so the store is spent when ATP falls. '
    'Storing a glucose costs 2 ATP and gives 3 later (instead of 2 now): about 75% round trip. The junction structure is the stoichiometry, so carbon (6 glc + 3 lac), adenine (ATP + ADP) '
    'and phosphoryl groups (ATP + Pi) are conserved by construction. Each rate is the common modular rate law '
    '(Liebermeister et al. 2010) written with the thermodynamic factor (1 - Gamma/K_eq): it saturates in each reactant and '
    'its sign is always that of the affinity A = RT ln(K_eq/Gamma), so the dissipation v A is never negative. Heat is the '
    'reaction enthalpy released: 2 h_lac_gly - 2 h_ATP per glycolysis turnover, h_ATP per ATP hydrolysed (including ATP '
    'hydrolysed by transporters on the ATP_hydrolysis_port, e.g. the H+-ATPase), so at steady ATP the heat is 63 kJ per '
    'mol lactate (Gnaiger and Kemp 1990). Protons: glycolysis to lactate written with the ionised species (glucose + 2 '
    'ADP3- + 2 HPO4 2- -> 2 lactate- + 2 ATP4-) releases none, and each ATP hydrolysed releases one H+ (ATP4- + H2O -> '
    'ADP3- + HPO4 2- + H+), so the cytosol exports j_H_met = V_c v_hyd - j_ATP_ext on the H_production port for all ATP hydrolysed, '
    'its own and that of the ATP_hydrolysis_port users (transporters, protein turnover); at steady ATP the cell makes 2 H+ per glucose, which leave with lactate. Ports: cytosol_face_port (entrance and exit; the five concentrations shared, the '
    'five fluxes summed) for diffusion faces and the cell membrane; ATP_hydrolysis_port for ATP-driven transporters; '
    'compartment_volume (from the membrane); heat_source_port to a thermal_volume. The default instance is the steady state '
    'of the whole fermenter cell (generic_adaptive_cell fermenter) with total ATP turnover at the Lynch & Marinov (2015) '
    'maintenance cost for its volume at 20 C. Its tests feed it through cell_membrane_face from a 1 mM glucose medium. '
    'When a substrate runs out its concentration sits at solver round-off, so the affinities use ln(|Gamma| + 1e-300).\n\n'
    'Ports of the versions of energy_metabolism (they differ; allowed and documented, see modules/README.md, Versions):\n'
    'general O2/glucose/lactate/CO2 exchange ports, heat_source_port: mitochondrial_ArgusUNPUBLISHED_v01\n'
    'entrance/exit cytosol_face_port, general ATP_hydrolysis_port, compartment_volume, heat_source_port: '
    'adaptive_ArgusUNPUBLISHED_v01')
cyto_verif = {
    'sim_time': 120.0, 'pre_time': 0.0, 'dt': 0.01, 'solver': 'CVODE_myokit',
    'solver_info': {'rtol': 1e-10, 'atol': 1e-14, 'MaximumStep': 0.5},
    'outputs': [f'C_{s}' for s in cyto_species] + ['C_glyc', 'v_store', 'v_mob', 'A_store', 'A_mob'] +
               ['v_gly', 'v_hyd', 'A_gly', 'A_hyd', 'P_diss', 'Q_met', 'f_T', 'J_ATP_gly', 'J_ATP_mob', 'J_ATP_store',
                'J_ATP_maint', 'M_glyc_ATP',
                                                     'mf/J_glc', 'mf/J_lac'],
    'run_parameters': {'store_alloc': 0.5, 'C_glyc_init': 50.0},   # exercise storage and mobilisation
    # fed through the cell membrane from the background glucose of the medium (the medium side is open: C_glc_e, C_lac_e parameters);
    # without the H+-ATPase the maintenance ATPase is the only ATP user, so ATP settles a little above its whole-cell value
    'harness': {
        'module_array': [['mf', V_NEW, 'cell_membrane_face', '', 'mod', 'default'],
                         ['mod', V_NEW, 'energy_metabolism', 'mf', '', 'default']],
        'parameters': [['C_glc_e_mf', 'millimolar', GLC_E, 'harness: the background glucose of the medium'],
                       ['C_lac_e_mf', 'millimolar', LAC_E, 'harness: no lactate in the medium'],
                       ['A_mem_mf', 'um2', A_MEM, 'harness: in-plane membrane of the 2 um cell']],
    },
    'invariants': cyto_invariants,
    'bc_sweep': TEMPERATURE_SWEEP,
    'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]},
}
write_version('cell', 'energy_metabolism', 'energy_metabolism_fermenter', cyto_vars, cyto_eqs,
              {'required_citations': ['argus2026unpublished', 'liebermeister2010modular'], **cyto_ports,
               'roles': cyto_roles},
              cyto_params, cyto_verif, {'reviewed': False, 'notes': cyto_notes},
              ['argus2026unpublished', 'liebermeister2010modular', 'lynch2015bioenergetic', 'rolfe1997cellular',
               'bennett2009absolute', 'gnaiger1990anaerobic', 'gajewski1986thermodynamics'], V_NEW)

# ------------------------------------------------------------------ 2. membrane potential (capacitor, disc geometry)
mem_vars = [
    ('t', 'second', 'in', None, None),
    ('R_cell', 'um', 'in', None, 'disc radius'), ('h_cell', 'um', 'in', None, 'disc depth'),
    ('c_m', 'uF_per_cm2', 'in', None, None),
    ('A_mem', 'um2', 'out', None, 'in-plane membrane area 2 pi R h (top and bottom exchange nothing in AdaptOrDeity)'),
    ('Vol', 'm3', 'out', None, None), ('Cm', 'picoF', 'out', None, None),
    ('V_rest_approx', 'milliV', 'in', None, None),
    ('V', 'milliV', 'out', 'V_rest_approx', None),
    ('i_ion', 'nanoA', 'in', None, 'summed membrane currents (positive outward)'),
    ('I_in', 'nanoA', 'in', None, 'injected current (positive depolarises)'),
]
mem_eqs = [
    ('A_mem', '2.0 * pi * R_cell * h_cell'),
    ('Vol', 'pi * R_cell ** 2.0 * h_cell * u(1e-18, m3_per_um3)'),
    ('Cm', 'c_m * A_mem * u(0.01, picoF_cm2_per_uF_um2)'),
    ('d/dt V', 'u(1e6, milliV_per_kiloV) * (I_in - i_ion) / Cm'),
]
mem_roles = {'R_cell': 'constant', 'h_cell': 'constant', 'c_m': 'constant', 'A_mem': 'variable', 'Vol': 'variable',
             'Cm': 'variable', 'V_rest_approx': 'constant', 'V': 'variable', 'i_ion': 'variable',
             'I_in': 'boundary_condition'}
write_version('cell', 'membrane_potential', 'membrane_fermenter', mem_vars, mem_eqs,
              {'required_citations': ['argus2026unpublished'],
               'entrance_ports': [{'port_type': 'i_stim', 'variables': ['I_in']}],
               'exit_ports': [{'port_type': 'membrane_voltage', 'variables': ['V'], 'multi_port': 'True'}],
               'general_ports': [{'port_type': 'membrane_current', 'variables': ['i_ion'], 'multi_port': 'sum'},
                                 {'port_type': 'membrane_area', 'variables': ['A_mem'], 'multi_port': 'True'},
                                 {'port_type': 'compartment_volume', 'variables': ['Vol'], 'multi_port': 'True'},
                                 {'port_type': 'H_compartment_volume', 'variables': ['Vol'], 'multi_port': 'True'}],
               'roles': mem_roles},
              [('R_cell', R_CELL, 'geometry: 2 um cell (typical prokaryote)', False),
               ('h_cell', H_CELL, 'geometry: depth of the 2D sheet equal to the cell diameter', False),
               ('c_m', C_M, 'gentet2000direct; specific membrane capacitance about 0.9 uF/cm2; 1 uF/cm2 is the '
                            'conventional value (as membrane_potential SN_soma_Argus2026_v01)', True),
               ('V_rest_approx', V_REST, f'{PH}; starting voltage, the calibrated rest of the fermenter cell (E. coli '
                                         'rests at -140 to -105 mV through proton pumping, milo2010bionumbers BNID 103386)',
                False),
               ('i_ion', 0.0, 'no channel connected', False), ('I_in', 0.0, 'no stimulus', False)],
              {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.001, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-14}, 'outputs': ['V', 'A_mem', 'Vol', 'Cm'],
               'run_parameters': {'I_in': 1e-6},
               'invariants': [{'description': 'a capacitor: Cm (V - V0) = I_in t with no channels (relative 1e-9); '
                                              'A_mem = 2 pi R h, Vol = pi R^2 h (exact)',
                               'expr': 'np.allclose(Cm*(V - V[0]), 1e6*I_in*t, rtol=1e-9, atol=1e-12) and '
                                       'np.allclose(A_mem, 2*np.pi*R_cell*h_cell, rtol=1e-12, atol=0) and '
                                       'np.allclose(Vol, np.pi*R_cell**2*h_cell*1e-18, rtol=1e-12, atol=0)',
                               'applies': 'run'}],
               'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
              {'reviewed': False, 'notes': (
                  'The membrane of a cell without dendrites, axon or electrode: a capacitor Cm = c_m A_mem charged by the '
                  'summed membrane currents, dV/dt = (I_in - i_ion)/Cm (the maths of SN_soma_Argus2026_v01 without its '
                  'sensed-voltage filter and axial current). Geometry is a disc of radius R_cell and depth h_cell, as a cell '
                  'sits in AdaptOrDeity\'s 2D sheet: A_mem = 2 pi R h is the in-plane membrane only (the top and bottom are '
                  'insulated), Vol = pi R^2 h. Same port types as SN_soma_Argus2026_v01 (membrane_current, membrane_voltage, '
                  'membrane_area, compartment_volume, i_stim), so the library\'s ion channels connect unchanged.\n\n'
                  'The volume is also exported on H_compartment_volume for the proton pool (H_handling), which is connected '
                  'to the cytosol and so cannot take compartment_volume a second time.\n\n'
                  'Ports of the versions of membrane_potential (they differ; allowed and documented): '
                  'adaptive_ArgusUNPUBLISHED_v01 has no VI_port and no shell_membrane_area, and adds H_compartment_volume.')},
              ['argus2026unpublished', 'milo2010bionumbers'], V_NEW)

# ------------------------------------------------------------------ 3. Ca handling (cytosol + fixed buffer, no ER)
KD_BUF = 8e3 / 0.1e6
B_INIT = 8.0 * KD_BUF / (KD_BUF + CAI)
ca_vars = [
    ('t', 'second', 'in', None, None),
    ('Cai_init', 'millimolar', 'in', None, None), ('Cai', 'millimolar', 'out', 'Cai_init', 'free cytosolic Ca'),
    ('B_Ca_init', 'millimolar', 'in', None, None), ('B_Cai', 'millimolar', 'out', 'B_Ca_init', 'free buffer'),
    ('CaB_Cai', 'millimolar', 'out', None, 'Ca-bound buffer'),
    ('B_Ca_total', 'millimolar', 'in', None, None),
    ('k_bind_Ca', 'per_millimolar_s', 'in', None, None), ('k_unbind_Ca', 'per_s', 'in', None, None),
    ('j_Ca_in', 'mol_per_s', 'in', None, 'summed Ca fluxes into the cell (Ca_flux port)'),
    ('Vol', 'm3', 'in', None, None),
    ('Ca_total', 'millimolar', 'out', None, 'free + bound'),
]
ca_eqs = [
    ('CaB_Cai', 'B_Ca_total - B_Cai'),
    ('d/dt B_Cai', 'k_unbind_Ca * CaB_Cai - k_bind_Ca * Cai * B_Cai'),
    ('d/dt Cai', 'j_Ca_in / Vol + k_unbind_Ca * CaB_Cai - k_bind_Ca * Cai * B_Cai'),
    ('Ca_total', 'Cai + CaB_Cai'),
]
write_version('cell', 'Ca_handling', 'Ca_fermenter', ca_vars, ca_eqs,
              {'required_citations': ['argus2026unpublished'],
               'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'internal_Ca', 'variables': ['Cai'], 'multi_port': 'True'},
                                 {'port_type': 'Ca_flux', 'variables': ['j_Ca_in'], 'multi_port': 'sum'},
                                 {'port_type': 'compartment_volume', 'variables': ['Vol']}],
               'roles': {'Cai_init': 'constant', 'Cai': 'variable', 'B_Ca_init': 'constant', 'B_Cai': 'variable',
                         'CaB_Cai': 'variable', 'B_Ca_total': 'constant', 'k_bind_Ca': 'constant',
                         'k_unbind_Ca': 'constant', 'j_Ca_in': 'variable', 'Vol': 'boundary_condition',
                         'Ca_total': 'variable'}},
              [('Cai_init', CAI, 'gangola1987maintenance; free cytosolic Ca2+ in E. coli 90 +- 10 nM', True),
               ('B_Ca_init', B_INIT, 'definitional; free buffer in equilibrium with Cai_init: B_tot Kd/(Kd + Cai), '
                                     'Kd = k_unbind/k_bind = 80 uM', False),
               ('B_Ca_total', par('B_Ca_total'), f'{PH}; the SN soma\'s fixed buffer (8 mM, matveev2002new); bacterial Ca buffering is '
                                   'not characterised here', False),
               ('k_bind_Ca', par('k_bind_Ca'), f'{PH}; the SN soma\'s buffer on-rate (matveev2002new, quoted set)', False),
               ('k_unbind_Ca', par('k_unbind_Ca'), f'{PH}; the SN soma\'s buffer off-rate (matveev2002new, actual set)', False),
               ('j_Ca_in', 0.0, 'no Ca pathway connected', False),
               ('Vol', VOL, 'geometry: disc of radius 1 um, 2 um deep', False)],
              {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 1e-4, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-16}, 'outputs': ['Cai', 'B_Cai', 'CaB_Cai', 'Ca_total'],
               'run_parameters': {'Cai_init': 2e-4},
               'invariants': [{'description': 'closed cell: total Ca (free + bound) is conserved while the buffer '
                                              're-equilibrates (relative 1e-9); buffer conserved (exact)',
                               'expr': 'np.ptp(Ca_total) <= 1e-9*np.max(np.abs(Ca_total)) and '
                                       'np.allclose(B_Cai + CaB_Cai, B_Ca_total, rtol=1e-12, atol=0)'},
                              {'description': 'the buffer reaches equilibrium: k_bind Cai B = k_unbind CaB at the end '
                                              '(relative 1e-6)',
                               'expr': 'abs(k_bind_Ca*Cai[-1]*B_Cai[-1] - k_unbind_Ca*CaB_Cai[-1]) <= '
                                       '1e-6*k_unbind_Ca*CaB_Cai[-1]', 'applies': 'run'}],
               'bc_sweep': {'constraints': ['B_Ca_init <= B_Ca_total'],
                            'rationale': 'the free buffer cannot exceed the total buffer'},
               'timestep': {'scheme': 'rk4', 'dts': [2e-7, 1e-7, 5e-8], 't_end': 1e-4, 'min_order': 0.8, 'tol': 1e-3,
                            'cvode_tol': 1e-2},
               'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-10}], 't_end': 0.01,
                             'time_budget': 10}},
              {'reviewed': False, 'notes': (
                  'Cytosolic Ca of a cell without an ER (a prokaryote): one free Ca pool with a fixed buffer, the bulk '
                  'equations of SN_soma_Argus2026_v01 without its sub-membrane shell, pre-ER and ER compartments. dCai/dt '
                  '= j_Ca_in/Vol + k_unbind CaB - k_bind Cai B, dB/dt = k_unbind CaB - k_bind Cai B. The buffer equilibrates '
                  'in about 1/(k_bind B) ~ 1 us, the stiffest part of the fermenter cell (fixed-step solvers need very small '
                  'steps; a rapid-buffer approximation is a candidate for the game).\n\n'
                  'Ports of the versions of Ca_handling (they differ): adaptive_ArgusUNPUBLISHED_v01 has internal_Ca, '
                  'Ca_flux and compartment_volume only.')},
              ['argus2026unpublished'], V_NEW)

# ------------------------------------------------------------------ 4. i_leak_Ca (new type)
lca_vars = [
    ('t', 'second', 'in', None, None),
    ('V', 'milliV', 'in', None, None), ('Cai', 'millimolar', 'in', None, None), ('Cao', 'millimolar', 'in', None, None),
    ('R', 'J_per_mol_K', 'in', None, None), ('T', 'kelvin', 'in', None, None), ('F', 'C_per_mol', 'in', None, None),
    ('A_mem', 'um2', 'in', None, None), ('rho_leak_Ca', 'per_um2', 'in', None, None),
    ('gamma_leak_Ca', 'picoS', 'in', None, None),
    ('g_leak_Ca', 'microS', 'out', None, None), ('E_Ca', 'milliV', 'out', None, 'Ca Nernst potential'),
    ('i_leak_Ca', 'nanoA', 'out', None, None), ('j_leak_Ca', 'mol_per_s', 'out', None, 'Ca into the cell'),
]
lca_eqs = [
    ('g_leak_Ca', 'A_mem * rho_leak_Ca * gamma_leak_Ca * u(1e-6, microS_per_picoS)'),
    ('E_Ca', 'R * T / (2.0 * F) * ln(Cao / Cai) * u(1000.0, milliV_per_V)'),
    ('i_leak_Ca', 'g_leak_Ca * (V - E_Ca)'),
    ('j_leak_Ca', '-i_leak_Ca / (2.0 * F) * u(1e-9, A_per_nanoA)'),
]
write_version('cell/ion_channels', 'i_leak_Ca', 'i_leak_Ca', lca_vars, lca_eqs,
              {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'membrane_voltage', 'variables': ['V']},
                                 {'port_type': 'internal_Ca', 'variables': ['Cai']},
                                 {'port_type': 'membrane_current', 'variables': ['i_leak_Ca']},
                                 {'port_type': 'Ca_flux', 'variables': ['j_leak_Ca']},
                                 {'port_type': 'membrane_area', 'variables': ['A_mem']}],
               'roles': {'V': 'boundary_condition', 'Cai': 'boundary_condition', 'Cao': 'constant', 'R': 'global_constant',
                         'T': 'global_constant', 'F': 'global_constant', 'A_mem': 'boundary_condition',
                         'rho_leak_Ca': 'constant', 'gamma_leak_Ca': 'constant', 'g_leak_Ca': 'variable',
                         'E_Ca': 'variable', 'i_leak_Ca': 'variable', 'j_leak_Ca': 'variable'}},
              [('V', V_REST, 'operating point', False), ('Cai', CAI, 'gangola1987maintenance; 90 nM', True),
               ('Cao', CAO, 'millero2008composition; standard seawater Ca 10.28 mmol/kg', True),
               ('R', R, 'definitional; as the library\'s SN modules', True),
               ('T', T, 'operating point: 20 C', False),
               ('F', F, 'definitional; Faraday constant, as the library\'s SN modules', True),
               ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
               ('rho_leak_Ca', RHO_CA, 'calibrated: Ca entry balances Na/Ca exchange (1% of NhaA\'s cycles) at rest',
                False),
               ('gamma_leak_Ca', GAMMA_CA, f'{PH}; single-channel conductance of a Ca leak pathway', False)],
              {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-14}, 'outputs': ['g_leak_Ca', 'E_Ca', 'i_leak_Ca', 'j_leak_Ca'],
               'run_parameters': {},
               'invariants': [{'description': 'Ohmic leak at the Ca Nernst potential; j = -i/(2F) (exact)',
                               'expr': 'np.allclose(E_Ca, 1000*R*T/(2*F)*np.log(Cao/Cai), rtol=1e-12, atol=0) and '
                                       'np.allclose(i_leak_Ca, g_leak_Ca*(V - E_Ca), rtol=1e-12, atol=0) and '
                                       'np.allclose(j_leak_Ca, -i_leak_Ca/(2*F)*1e-9, rtol=1e-12, atol=0)'},
                              {'description': 'second law: the leak runs down its electrochemical gradient, '
                                              'i (V - E_Ca) >= 0', 'expr': 'np.all(i_leak_Ca*(V - E_Ca) >= 0)'}],
               'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
              {'reviewed': False, 'notes': (
                  'A Ca leak, the maths of i_leak_Na / i_leak_K Argus2026_v01 for a divalent ion: i = g (V - E_Ca), '
                  'E_Ca = (RT/2F) ln(Cao/Cai), j = -i/(2F), g = A_mem rho gamma. It computes its own Nernst potential from '
                  'the internal_Ca port, since Ca_handling versions do not export one.')},
              ['argus2026unpublished', 'millero2008composition', 'gangola1987maintenance'], 'ArgusUNPUBLISHED_v01')

# ------------------------------------------------------------------ 5. i_NaK with ATP use
nak_vars = [
    ('t', 'second', 'in', None, None),
    ('V', 'milliV', 'in', None, None), ('Nai', 'millimolar', 'in', None, None), ('Ko', 'millimolar', 'in', None, None),
    ('Ko_mid_NaK', 'millimolar', 'in', None, None), ('Nai_mid_NaK', 'millimolar', 'in', None, None),
    ('Nai_power_NaK', 'dimensionless', 'in', None, None),
    ('A_mem', 'um2', 'in', None, None), ('rho_NaK', 'per_um2', 'in', None, None), ('r_NaK', 'per_s', 'in', None, None),
    ('F', 'C_per_mol', 'in', None, None),
    ('ATP', 'millimolar', 'in', None, 'cytosolic ATP (ATP_hydrolysis_port)'), ('ADP', 'millimolar', 'in', None, None),
    ('Pi', 'millimolar', 'in', None, None), ('K_m_ATP_NaK', 'millimolar', 'in', None, None),
    ('k_NaK', 'nanoA', 'out', None, None), ('i_NaK', 'nanoA', 'out', None, None),
    ('j_NaK_Na', 'mol_per_s', 'out', None, None), ('j_NaK_K', 'mol_per_s', 'out', None, None),
    ('j_NaK_ATP', 'mol_per_s', 'out', None, '1 ATP hydrolysed per cycle'),
    ('j_NaK_ADP', 'mol_per_s', 'out', None, None), ('j_NaK_Pi', 'mol_per_s', 'out', None, None),
    ('E_Na', 'milliV', 'in', None, 'Nernst potentials (from ion_concentrations), for the ion work per cycle'),
    ('E_K', 'milliV', 'in', None, None),
    ('R', 'J_per_mol_K', 'in', None, None), ('T', 'kelvin', 'in', None, None),
    ('dG0_hyd', 'J_per_mol', 'in', None, None), ('c_std', 'millimolar', 'in', None, None),
    ('A_ATP', 'J_per_mol', 'out', None, 'free energy of ATP hydrolysis available per cycle'),
    ('W_ion', 'J_per_mol', 'out', None, 'electrochemical work per cycle: 3 Na out, 2 K in'),
    ('f_thermo', 'dimensionless', 'out', None, 'thermodynamic factor: 0 at equilibrium, negative (reverse) below it'),
]
nak_eqs = [
    ('k_NaK', 'A_mem * rho_NaK * r_NaK * u(1.602176634e-10, nanoC)'),
    ('A_ATP', '-dG0_hyd - R * T * ln(abs((ADP / c_std) * (Pi / c_std) / (ATP / c_std)) + 1e-300)'),
    ('W_ion', 'F * (3.0 * (E_Na - V) + 2.0 * (V - E_K)) * u(1e-3, V_per_milliV)'),
    ('f_thermo', '1.0 - exp(-(A_ATP - W_ion) / (R * T))'),
    ('i_NaK', 'k_NaK * ((V + u(150.0, milliV)) / (V + u(200.0, milliV))) / (1.0 + (Nai_mid_NaK / Nai) ** Nai_power_NaK) '
              '* Ko / (Ko + Ko_mid_NaK) * ATP / (K_m_ATP_NaK + ATP) * f_thermo'),
    ('j_NaK_Na', '-3.0 * i_NaK / F * u(1e-9, A_per_nanoA)'),
    ('j_NaK_K', '2.0 * i_NaK / F * u(1e-9, A_per_nanoA)'),
    ('j_NaK_ATP', '-i_NaK / F * u(1e-9, A_per_nanoA)'),
    ('j_NaK_ADP', 'i_NaK / F * u(1e-9, A_per_nanoA)'),
    ('j_NaK_Pi', 'i_NaK / F * u(1e-9, A_per_nanoA)'),
]
V_NAK_VERSION = 'ATP_ArgusUNPUBLISHED_v01'
write_version('cell/ion_channels', 'i_NaK', 'i_NaK_ATP', nak_vars, nak_eqs,
              {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'membrane_voltage', 'variables': ['V']},
                                 {'port_type': 'internal_Na', 'variables': ['Nai']},
                                 {'port_type': 'membrane_current', 'variables': ['i_NaK']},
                                 {'port_type': 'Na_flux', 'variables': ['j_NaK_Na']},
                                 {'port_type': 'K_flux', 'variables': ['j_NaK_K']},
                                 {'port_type': 'membrane_area', 'variables': ['A_mem']},
                                 {'port_type': 'ATP_hydrolysis_port',
                                  'variables': ['ATP', 'ADP', 'Pi', 'j_NaK_ATP', 'j_NaK_ADP', 'j_NaK_Pi']},
                                 {'port_type': 'Na_Nernst_potential', 'variables': ['E_Na']},
                                 {'port_type': 'E_K_potential', 'variables': ['E_K']}],
               'roles': {'V': 'boundary_condition', 'Nai': 'boundary_condition', 'Ko': 'constant',
                         'Ko_mid_NaK': 'constant', 'Nai_mid_NaK': 'constant', 'Nai_power_NaK': 'constant',
                         'A_mem': 'boundary_condition', 'rho_NaK': 'constant', 'r_NaK': 'constant',
                         'F': 'global_constant', 'ATP': 'boundary_condition', 'ADP': 'boundary_condition',
                         'Pi': 'boundary_condition', 'K_m_ATP_NaK': 'constant', 'k_NaK': 'variable',
                         'i_NaK': 'variable', 'j_NaK_Na': 'variable', 'j_NaK_K': 'variable', 'j_NaK_ATP': 'variable',
                         'j_NaK_ADP': 'variable', 'j_NaK_Pi': 'variable', 'E_Na': 'boundary_condition',
                         'E_K': 'boundary_condition', 'R': 'global_constant', 'T': 'global_constant',
                         'dG0_hyd': 'constant', 'c_std': 'constant', 'A_ATP': 'variable', 'W_ion': 'variable',
                         'f_thermo': 'variable'}},
              [('V', V_NAK, 'operating point of an animal-like cell at rest', False),
               ('Nai', NAI, f'{PH}; cytosolic Na, 20 mM assumed', False),
               ('Ko', KO, 'millero2008composition; standard seawater K 10.21 mmol/kg', True),
               ('Ko_mid_NaK', NAK['Ko_mid_NaK'], 'as i_NaK Argus2026_v01 (default instance)', False),
               ('Nai_mid_NaK', NAK['Nai_mid_NaK'], 'as i_NaK Argus2026_v01 (default instance)', False),
               ('Nai_power_NaK', NAK['Nai_power_NaK'], 'as i_NaK Argus2026_v01 (default instance)', False),
               ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
               ('rho_NaK', RHO_NAK, 'calibrated: the pump uses 25% of a 2 um cell\'s maintenance ATP at -70 mV (rolfe1997cellular, '
                                    'mammals 19-28%)', False),
               ('r_NaK', NAK['r_NaK'], 'as i_NaK Argus2026_v01 (default instance): cycles per pump per second', False),
               ('F', F, 'definitional; as the library\'s SN modules', True),
               ('ATP', ATP, 'bennett2009absolute; ATP 9.6 mM in E. coli', True),
               ('ADP', ADP, f'{PH}; about 0.56 mM in E. coli (to confirm)', False),
               ('Pi', PI, f'{PH}; near 10 mM in E. coli', False),
               ('K_m_ATP_NaK', NAK['K_m_ATP_NaK'], f'{PH}; pump affinity for ATP, order 0.1-0.5 mM', False),
               ('E_Na', E_NA, 'Nernst potential at the default concentrations', False),
               ('E_K', E_K, 'Nernst potential at the default concentrations', False),
               ('R', R, 'definitional; as the library\'s SN modules', True), ('T', T, 'operating point: 20 C', False),
               ('dG0_hyd', DG0_HYD, f'{PH}; textbook standard Gibbs energy of ATP hydrolysis (-30.5 kJ/mol)', False),
               ('c_std', C_STD, 'definitional; standard concentration 1 M', True)],
              {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-14},
               'outputs': ['i_NaK', 'j_NaK_Na', 'j_NaK_K', 'j_NaK_ATP', 'j_NaK_ADP', 'j_NaK_Pi', 'A_ATP', 'W_ion',
                           'f_thermo'], 'run_parameters': {},
               'invariants': [{'description': 'stoichiometry: 3 Na out, 2 K in, 1 ATP -> ADP + Pi per cycle; the current '
                                              'is the net charge, j_Na + j_K = -i/F (exact)',
                               'expr': 'np.allclose(j_NaK_Na, 3*j_NaK_ATP, rtol=1e-12, atol=0) and '
                                       'np.allclose(j_NaK_K, -2*j_NaK_ATP, rtol=1e-12, atol=0) and '
                                       'np.allclose(j_NaK_ADP, -j_NaK_ATP, rtol=1e-12, atol=0) and '
                                       'np.allclose(j_NaK_Pi, -j_NaK_ATP, rtol=1e-12, atol=0) and '
                                       'np.allclose(j_NaK_Na + j_NaK_K, -i_NaK/F*1e-9, rtol=1e-12, atol=0)'},
                              {'description': 'second law: the pump runs forward only while the ATP free energy exceeds the '
                                              'ion work per cycle, i_NaK (A_ATP - W_ion) >= 0',
                               'expr': 'np.all(i_NaK*(A_ATP - W_ion) >= 0)'}],
               'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
              {'reviewed': False, 'notes': (
                  'The Na/K pump of i_NaK Argus2026_v01 (same voltage, Na and K dependence) made to spend ATP and obey the '
                  'second law: the rate is multiplied by ATP/(K_m_ATP + ATP) and by the thermodynamic factor f_thermo = 1 - '
                  'exp(-(A_ATP - W_ion)/RT), where A_ATP is the free energy of ATP hydrolysis and W_ion = F(3(E_Na - V) + '
                  '2(V - E_K)) the electrochemical work of one cycle, so the pump stops at equilibrium and can never run '
                  'forward against it (the Argus2026_v01 law alone keeps pumping when ATP runs low, which creates energy). '
                  'Every cycle hydrolyses one ATP to ADP + Pi, taken from the cytosol through the ATP_hydrolysis_port. So keeping the membrane potential has an energy cost, and '
                  'the pump slows as ATP runs out. The Na/K-ATPase is a eukaryotic protein (AdaptOrDeity unlocks it with '
                  'eukaryotes; the prokaryotic fermenter cell uses an F1Fo H+-ATPase, NhaA and Trk/KefC instead). The '
                  'default instance is calibrated at -70 mV with the pump using 25% of a 2 um cell\'s maintenance ATP '
                  '(rolfe1997cellular).\n\n'
                  'Ports of the versions of i_NaK (they differ): ATP_ArgusUNPUBLISHED_v01 adds ATP_hydrolysis_port, '
                  'Na_Nernst_potential and E_K_potential to Argus2026_v01\'s ports.')},
              ['argus2026unpublished', 'millero2008composition', 'bennett2009absolute', 'rolfe1997cellular'], V_NAK_VERSION)

# ------------------------------------------------------------------ 6. cell membrane face (glucose and lactate carriers)
mf_vars = [
    ('t', 'second', 'in', None, None),
    ('C_glc_e', 'millimolar', 'in', None, 'medium side (entrance medium_face_port)'), ('C_lac_e', 'millimolar', 'in', None, None),
    ('C_glc', 'millimolar', 'in', None, 'cytosol side (exit cytosol_face_port)'), ('C_lac', 'millimolar', 'in', None, None),
    ('C_ATP', 'millimolar', 'in', None, 'impermeant, carried for the port'), ('C_ADP', 'millimolar', 'in', None, None),
    ('C_Pi', 'millimolar', 'in', None, None),
    ('A_mem', 'um2', 'in', None, None),
    ('jmax_glc', 'mol_per_s_per_um2', 'in', None, 'carrier capacity per membrane area'),
    ('K_glc_t', 'millimolar', 'in', None, None),
    ('jmax_lac', 'mol_per_s_per_um2', 'in', None, None), ('K_lac_t', 'millimolar', 'in', None, None),
    ('J_glc', 'mol_per_s', 'out', None, 'into the cell'), ('J_lac', 'mol_per_s', 'out', None, 'into the cell'),
    ('J_glc_e', 'mol_per_s', 'out', None, None), ('J_lac_e', 'mol_per_s', 'out', None, None),
    ('J_glc_c', 'mol_per_s', 'out', None, None), ('J_lac_c', 'mol_per_s', 'out', None, None),
    ('J_ATP_c', 'mol_per_s', 'out', None, None), ('J_ADP_c', 'mol_per_s', 'out', None, None),
    ('J_Pi_c', 'mol_per_s', 'out', None, None),
    ('H_i', 'millimolar', 'in', None, 'free cytosolic H+ (internal_H port)'),
    ('H_e', 'millimolar', 'in', None, 'free H+ in the medium'),
    ('H_ref', 'millimolar', 'in', None, 'reference H+ for the lactic-acid activity X = C_lac H / H_ref'),
    ('X_lac_e', 'millimolar', 'out', None, None), ('X_lac_i', 'millimolar', 'out', None, None),
    ('J_H_lac', 'mol_per_s', 'out', None, 'H+ into the cell with lactate (H_flux port)'),
    ('m_Glf', 'dimensionless', 'in', None, 'glucose carrier expression, 1 = the calibrated cell (first protein_expression '
                                          'port)'),
    ('m_LacT', 'dimensionless', 'in', None, 'lactate carrier expression (second protein_expression port)'),
]
mf_eqs = [
    ('J_glc', 'jmax_glc * m_Glf * A_mem * (C_glc_e / (K_glc_t + C_glc_e) - C_glc / (K_glc_t + C_glc))'),
    ('X_lac_e', 'C_lac_e * H_e / H_ref'), ('X_lac_i', 'C_lac * H_i / H_ref'),
    ('J_lac', 'jmax_lac * m_LacT * A_mem * (X_lac_e / (K_lac_t + X_lac_e) - X_lac_i / (K_lac_t + X_lac_i))'),
    ('J_H_lac', 'J_lac'),
    ('J_glc_e', '-J_glc'), ('J_lac_e', '-J_lac'), ('J_glc_c', 'J_glc'), ('J_lac_c', 'J_lac'),
    ('J_ATP_c', 'u(0.0, mol_per_s)'), ('J_ADP_c', 'u(0.0, mol_per_s)'), ('J_Pi_c', 'u(0.0, mol_per_s)'),
]
mf_roles = {n: 'boundary_condition' for n in ('C_glc_e', 'C_lac_e', 'C_glc', 'C_lac', 'C_ATP', 'C_ADP', 'C_Pi', 'A_mem',
                                               'H_i', 'm_Glf', 'm_LacT')}
mf_roles.update({n: 'constant' for n in ('jmax_glc', 'K_glc_t', 'jmax_lac', 'K_lac_t', 'H_e', 'H_ref')})
mf_roles.update({n: 'variable' for n in ('X_lac_e', 'X_lac_i', 'J_H_lac')})
mf_roles.update({n: 'variable' for n in ('J_glc', 'J_lac', 'J_glc_e', 'J_lac_e', 'J_glc_c', 'J_lac_c', 'J_ATP_c', 'J_ADP_c',
                                         'J_Pi_c')})
write_version('transport', 'cell_membrane_face', 'cell_membrane_face_fermenter', mf_vars, mf_eqs,
              {'required_citations': ['argus2026unpublished'],
               'entrance_ports': [{'port_type': 'medium_face_port', 'variables': ['C_glc_e', 'C_lac_e', 'J_glc_e', 'J_lac_e']}],
               'exit_ports': [{'port_type': 'cytosol_face_port',
                               'variables': ['C_glc', 'C_lac', 'C_ATP', 'C_ADP', 'C_Pi', 'J_glc_c', 'J_lac_c', 'J_ATP_c',
                                             'J_ADP_c', 'J_Pi_c']}],
               'general_ports': [{'port_type': 'membrane_area', 'variables': ['A_mem']},
                                 {'port_type': 'internal_H', 'variables': ['H_i']},
                                 {'port_type': 'H_flux', 'variables': ['J_H_lac']},
                                 # two pools connect here in submodule order: the Glf pool first, then the LacT pool
                                 {'port_type': 'protein_expression', 'variables': ['m_Glf']},
                                 {'port_type': 'protein_expression', 'variables': ['m_LacT']}], 'roles': mf_roles},
              [('C_glc_e', GLC_E, '', False),
               ('C_lac_e', LAC_E, 'no lactate in the medium', False),
               ('C_glc', GLC_I, 'the fermenter cell\'s calibrated steady state', False),
               ('C_lac', LAC_I, 'the fermenter cell\'s calibrated steady state', False),
               ('C_ATP', ATP, 'bennett2009absolute; ATP 9.6 mM in E. coli (impermeant here)', True),
               ('C_ADP', ADP, f'{PH}; impermeant here', False), ('C_Pi', PI, f'{PH}; impermeant here', False),
               ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
               ('jmax_glc', VMAX_GLC / A_MEM, 'calibrated: at the default steady state the uptake-limited cell (5 uM glucose '
                                             'inside) takes up its resting need; from a richer medium it takes up more; '
                                             'carrier density is a game trait', False),
               ('K_glc_t', K_GLC_T, 'weisser1995functional; Glf glucose facilitator of Zymomonas mobilis, Km 4.1 mM', True),
               ('jmax_lac', VMAX_LAC / A_MEM, 'calibrated: lactate export at 4% of capacity at the default steady state, so '
                                             'it keeps up with fast growth',
                False),
               ('K_lac_t', K_LAC_T, f'{PH}; lactate carrier affinity for lactic-acid activity (no prokaryotic '
                                    'lactate exporter characterised here)', False),
               ('H_i', H_I, 'slonczewski1981ph; cytoplasmic pH 7.6 in E. coli', True),
               ('H_e', H_O, f'{PH}; medium pH 7.3, a Proterozoic ocean between the Archean 6.5-7.0 and Phanerozoic '
                            '7.5-9.0 of halevy2017geologic', False),
               ('H_ref', H_REF, 'definitional; pH 7, sets the scale of X only', False),
               ('m_Glf', 1.0, 'definitional; baseline expression (the calibrated cell)', False),
               ('m_LacT', 1.0, 'definitional; baseline expression (the calibrated cell)', False)],
              {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-14},
               'outputs': ['J_glc', 'J_lac', 'J_glc_e', 'J_lac_e', 'J_glc_c', 'J_lac_c', 'J_ATP_c', 'J_H_lac', 'X_lac_e',
                           'X_lac_i'], 'run_parameters': {},
               'invariants': [{'description': 'symmetric carriers; mass conserved across the membrane; ATP, ADP and Pi do '
                                              'not cross (exact)',
                               'expr': 'np.allclose(J_glc, jmax_glc*m_Glf*A_mem*(C_glc_e/(K_glc_t + C_glc_e) - C_glc/(K_glc_t + '
                                       'C_glc)), rtol=1e-12, atol=0) and np.array_equal(J_glc_e, -J_glc_c) and '
                                       'np.array_equal(J_lac_e, -J_lac_c) and np.all(J_ATP_c == 0) and np.array_equal(J_H_lac, J_lac)'},
                              {'description': 'second law: each carrier runs down its gradient, J ln(C_out/C_in) >= 0',
                               'expr': 'np.all(J_glc*np.log(C_glc_e/C_glc) >= 0) and np.all(J_lac*(X_lac_e - X_lac_i) >= 0)'}],
               'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}},
              {'reviewed': False, 'notes': (
                  'The plasma membrane of the fermenter cell between the medium (entrance medium_face_port: glucose, '
                  'lactate) and the cytosol (exit cytosol_face_port, the energy_metabolism fermenter species): a glucose '
                  'facilitator (Glf-type, Km 4.1 mM) and a lactate-/H+ symporter, each a symmetric carrier J = jmax A_mem '
                  '(S_out/(K + S_out) - S_in/(K + S_in)) (the maths of membrane_carrier_face ArgusUNPUBLISHED_v01, per '
                  'membrane area). For glucose S is the concentration; lactate leaves with one proton (electroneutral, as '
                  'lactic acid), so its carrier sees the lactic-acid activity S = X = C_lac H/H_ref and stops when '
                  'C_lac_e H_e = C_lac H_i, its electrochemical equilibrium. The protons go to the cytosolic proton pool '
                  '(J_H_lac on the H_flux port, H_i from the internal_H port). This export removes the 2 H+ per glucose '
                  'that ATP hydrolysis releases. ATP, ADP and Pi cannot cross. No heat (facilitated diffusion is entropic). '
                  'Ion currents cross the membrane in the ion-channel modules, not here.')},
              ['argus2026unpublished', 'weisser1995functional', 'bennett2009absolute', 'slonczewski1981ph',
               'halevy2017geologic'], V_NEW)

# ------------------------------------------------------------------ 6b. protons: the cytosolic pool and the prokaryotic ion set
V_PROK = 'ArgusUNPUBLISHED_v01'
OP = lambda name: (name, V_REST, 'operating point of the fermenter cell', False)
GEN = [('R', R, 'definitional; gas constant, as the library\'s SN modules', True),
       ('T', T, 'operating point: 20 C', False),
       ('F', F, 'definitional; Faraday constant, as the library\'s SN modules', True)]
STAB = {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-6, 'atol': 1e-8}]}
RUN1 = {'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01, 'solver': 'CVODE_myokit',
        'solver_info': {'rtol': 1e-10, 'atol': 1e-14}}

# --- H_handling: free H+ with a fixed (rapid, linear) buffer
h_vars = [
    ('t', 'second', 'in', None, None),
    ('H_i_init', 'millimolar', 'in', None, None),
    ('H_i', 'millimolar', 'out', 'H_i_init', 'free cytosolic H+'),
    ('H_o', 'millimolar', 'in', None, 'free H+ in the medium'),
    ('beta_H', 'millimolar', 'in', None, 'buffering capacity: mM of H+ bound per pH unit'),
    ('c_std', 'millimolar', 'in', None, None),
    ('j_H_in', 'mol_per_s', 'in', None, 'summed H+ fluxes into the cytosol across the membrane (H_flux port)'),
    ('j_H_prod', 'mol_per_s', 'in', None, 'summed H+ released by reactions in the cytosol (H_production port)'),
    ('Vol', 'm3', 'in', None, None),
    ('R', 'J_per_mol_K', 'in', None, None), ('T', 'kelvin', 'in', None, None), ('F', 'C_per_mol', 'in', None, None),
    ('pH_i', 'dimensionless', 'out', None, None), ('pH_o', 'dimensionless', 'out', None, None),
    ('E_H', 'milliV', 'out', None, 'H+ Nernst potential'),
    ('H_total_rel', 'millimolar', 'out', None, 'free + bound H+ relative to a reference: H_i - beta_H pH_i'),
]
h_eqs = [
    ('d/dt H_i', '(j_H_in + j_H_prod) / Vol / (1.0 + beta_H / (ln(10.0) * H_i))'),
    ('pH_i', '-ln(H_i / c_std) / ln(10.0)'), ('pH_o', '-ln(H_o / c_std) / ln(10.0)'),
    ('E_H', 'R * T / F * ln(H_o / H_i) * u(1000.0, milliV_per_V)'),
    ('H_total_rel', 'H_i - beta_H * pH_i'),
]
write_version('cell', 'H_handling', 'H_fermenter', h_vars, h_eqs,
              {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'internal_H', 'variables': ['H_i'], 'multi_port': 'True'},
                                 {'port_type': 'H_Nernst_potential', 'variables': ['E_H'], 'multi_port': 'True'},
                                 {'port_type': 'H_flux', 'variables': ['j_H_in'], 'multi_port': 'sum'},
                                 {'port_type': 'H_production', 'variables': ['j_H_prod'], 'multi_port': 'sum'},
                                 {'port_type': 'H_compartment_volume', 'variables': ['Vol']}],
               'roles': {'H_i_init': 'constant', 'H_i': 'variable', 'H_o': 'constant', 'beta_H': 'constant',
                         'c_std': 'constant', 'j_H_in': 'variable', 'j_H_prod': 'variable', 'Vol': 'boundary_condition',
                         'R': 'global_constant', 'T': 'global_constant', 'F': 'global_constant', 'pH_i': 'variable',
                         'pH_o': 'variable', 'E_H': 'variable', 'H_total_rel': 'variable'}},
              [('H_i_init', H_I, 'slonczewski1981ph; cytoplasmic pH 7.6 in E. coli', True),
               ('H_o', H_O, f'{PH}; pH 7.3, a Proterozoic ocean between the Archean 6.5-7.0 and Phanerozoic 7.5-9.0 of '
                            'halevy2017geologic (an environment parameter of the game)', False),
               ('beta_H', BETA_H, 'slonczewski1982effects; E. coli cytoplasmic buffering about 50 mM per pH unit', True),
               ('c_std', C_STD, 'definitional; standard concentration 1 M (pH = -log10(H/1 M))', True),
               ('j_H_in', 0.0, 'no H+ pathway connected', False),
               ('j_H_prod', 0.0, 'no reaction connected', False),
               ('Vol', VOL, 'geometry: disc of radius 1 um, 2 um deep', False)] + GEN,
              {'sim_time': 10.0, 'pre_time': 0.0, 'dt': 0.01, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-18},
               'outputs': ['H_i', 'pH_i', 'pH_o', 'E_H', 'H_total_rel', 'leak/j_leak_H'],
               'run_parameters': {},
               # fed by the H+ leak at -105 mV: protons flow in and the cytosol acidifies
               'harness': {
                   'module_array': [['leak', V_PROK, 'i_leak_H', '', 'mod', 'default'],
                                    ['mod', V_NEW, 'H_handling', 'leak', '', 'default']],
                   'parameters': [['V_leak', 'milliV', V_REST, 'harness: the fermenter cell\'s resting potential'],
                                  ['A_mem_leak', 'um2', A_MEM, 'harness: in-plane membrane of the 2 um cell']]},
               'invariants': [{'description': 'proton conservation with the buffer: Vol d(H_i - beta_H pH_i) = int j_H dt '
                                              '(relative 1e-6)',
                               'expr': 'abs(Vol*(H_total_rel[-1] - H_total_rel[0]) - np.trapz(leak__j_leak_H, t)) <= '
                                       '1e-6*np.trapz(np.abs(leak__j_leak_H), t) + 1e-30'},
                              {'description': 'pH and Nernst potential definitions (exact)',
                               'expr': 'np.allclose(pH_i, -np.log10(H_i/c_std), rtol=1e-12, atol=0) and '
                                       'np.allclose(E_H, 1000*R*T/F*np.log(H_o/H_i), rtol=1e-10, atol=1e-12)'}],
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'Cytosolic protons of a prokaryote: free H+ with a fixed linear buffer of capacity beta_H (mM of H+ '
                  'bound per pH unit), in rapid equilibrium, so the total H_i + beta_H (pH_ref - pH_i) changes only by the '
                  'summed fluxes: dH_i/dt = ((j_H_in + j_H_prod)/Vol)/(1 + beta_H/(ln10 H_i)), with membrane fluxes on H_flux and '
                  'protons released by cytosolic reactions on H_production (separate port types, because a membrane carrier '
                  'and the cytosol are themselves connected). Outputs pH_i, the medium pH_o and the H+ '
                  'Nernst potential E_H = (RT/F) ln(H_o/H_i) for the proton pathways. The buffer is algebraic, so unlike a '
                  'kinetic buffer it adds no stiffness. Its fluxes come from the H+-ATPase, NhaA, KefC, the H+ leak, the '
                  'lactate-/H+ symporter (H_flux) and ATP hydrolysis in the cytosol (H_production).')},
              ['argus2026unpublished', 'slonczewski1981ph', 'slonczewski1982effects', 'halevy2017geologic'], V_NEW)

# --- i_leak_H: proton leak
lh_vars = [
    ('t', 'second', 'in', None, None),
    ('V', 'milliV', 'in', None, None), ('E_H', 'milliV', 'in', None, None),
    ('A_mem', 'um2', 'in', None, None), ('g_spec_leak_H', 'microS_per_um2', 'in', None, 'conductance per membrane area'),
    ('F', 'C_per_mol', 'in', None, None),
    ('g_leak_H', 'microS', 'out', None, None), ('i_leak_H', 'nanoA', 'out', None, None),
    ('j_leak_H', 'mol_per_s', 'out', None, 'H+ into the cell'),
]
write_version('cell/ion_channels', 'i_leak_H', 'i_leak_H', lh_vars,
              [('g_leak_H', 'A_mem * g_spec_leak_H'), ('i_leak_H', 'g_leak_H * (V - E_H)'),
               ('j_leak_H', '-i_leak_H / F * u(1e-9, A_per_nanoA)')],
              {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'membrane_voltage', 'variables': ['V']},
                                 {'port_type': 'H_Nernst_potential', 'variables': ['E_H']},
                                 {'port_type': 'membrane_current', 'variables': ['i_leak_H']},
                                 {'port_type': 'H_flux', 'variables': ['j_leak_H']},
                                 {'port_type': 'membrane_area', 'variables': ['A_mem']}],
               'roles': {'V': 'boundary_condition', 'E_H': 'boundary_condition', 'A_mem': 'boundary_condition',
                         'g_spec_leak_H': 'constant', 'F': 'global_constant', 'g_leak_H': 'variable',
                         'i_leak_H': 'variable', 'j_leak_H': 'variable'}},
              [OP('V'), ('E_H', E_H, 'Nernst potential at pH 7.6 in, 7.3 out', False),
               ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
               ('g_spec_leak_H', GS_H, 'calibrated: H+ leak takes the share of pumped H+ not re-entering by NhaA and '
                                       'KefC (50%, placeholder)', False),
               ('F', F, 'definitional; as the library\'s SN modules', True)],
              {**RUN1, 'outputs': ['g_leak_H', 'i_leak_H', 'j_leak_H'], 'run_parameters': {},
               'invariants': [{'description': 'Ohmic leak at the H+ Nernst potential; j = -i/F (exact)',
                               'expr': 'np.allclose(i_leak_H, g_leak_H*(V - E_H), rtol=1e-12, atol=0) and '
                                       'np.allclose(j_leak_H, -i_leak_H/F*1e-9, rtol=1e-12, atol=0)'},
                              {'description': 'second law: i (V - E_H) >= 0', 'expr': 'np.all(i_leak_H*(V - E_H) >= 0)'}],
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'Proton leak across the membrane (lipid permeability and unspecific pathways), the maths of i_leak_Na '
                  'Argus2026_v01 per membrane area: i = A_mem g_spec (V - E_H), j = -i/F, with E_H from the proton pool '
                  '(H_Nernst_potential port).')},
              ['argus2026unpublished'], V_PROK)

# --- i_H_ATPase: F1Fo ATPase in hydrolysis mode, an H+ pump
ha_vars = [
    ('t', 'second', 'in', None, None),
    ('V', 'milliV', 'in', None, None), ('E_H', 'milliV', 'in', None, None),
    ('H_i', 'millimolar', 'in', None, 'free cytosolic H+ (internal_H port)'),
    ('K_act_H', 'millimolar', 'in', None, 'H+ for half activation (pH regulation)'),
    ('A_mem', 'um2', 'in', None, None), ('rho_H_ATPase', 'per_um2', 'in', None, 'density at baseline expression'),
    ('m_expr', 'dimensionless', 'in', None, 'expression, 1 = the calibrated cell (protein_expression port)'),
    ('r_H_ATPase', 'per_s', 'in', None, 'ATP hydrolysed per enzyme per second at full activity'),
    ('n_H', 'dimensionless', 'in', None, 'H+ pumped per ATP'),
    ('ATP', 'millimolar', 'in', None, None), ('ADP', 'millimolar', 'in', None, None), ('Pi', 'millimolar', 'in', None, None),
    ('K_m_ATP_H', 'millimolar', 'in', None, None),
    ('R', 'J_per_mol_K', 'in', None, None), ('T', 'kelvin', 'in', None, None), ('F', 'C_per_mol', 'in', None, None),
    ('dG0_hyd', 'J_per_mol', 'in', None, None), ('c_std', 'millimolar', 'in', None, None),
    ('k_H_ATPase', 'nanoA', 'out', None, None),
    ('A_ATP', 'J_per_mol', 'out', None, 'free energy of ATP hydrolysis'),
    ('W_H', 'J_per_mol', 'out', None, 'electrochemical work of pumping n_H H+ out'),
    ('f_thermo', 'dimensionless', 'out', None, 'thermodynamic factor: 0 at equilibrium, negative (synthesis) below it'),
    ('i_H_ATPase', 'nanoA', 'out', None, 'outward current of the pumped H+'),
    ('v_ATP_H', 'mol_per_s', 'out', None, 'ATP hydrolysed'),
    ('j_H_pumped', 'mol_per_s', 'out', None, 'H+ into the cell through the rotor (negative: pumped out)'),
    ('j_H_ATPase', 'mol_per_s', 'out', None, 'net H+ into the cytosol (H_flux port)'),
    ('j_ATP_H', 'mol_per_s', 'out', None, None), ('j_ADP_H', 'mol_per_s', 'out', None, None),
    ('j_Pi_H', 'mol_per_s', 'out', None, None),
]
ha_eqs = [
    ('k_H_ATPase', 'A_mem * rho_H_ATPase * m_expr * r_H_ATPase * n_H * u(1.602176634e-10, nanoC)'),
    ('A_ATP', '-dG0_hyd - R * T * ln(abs((ADP / c_std) * (Pi / c_std) / (ATP / c_std)) + 1e-300)'),
    ('W_H', 'n_H * F * (E_H - V) * u(1e-3, V_per_milliV)'),
    ('f_thermo', '1.0 - exp(-(A_ATP - W_H) / (R * T))'),
    ('i_H_ATPase', 'k_H_ATPase * ATP / (K_m_ATP_H + ATP) * H_i / (K_act_H + H_i) * f_thermo'),
    ('v_ATP_H', 'i_H_ATPase / (n_H * F) * u(1e-9, A_per_nanoA)'),
    ('j_H_pumped', '-n_H * v_ATP_H'), ('j_H_ATPase', 'j_H_pumped'),
    ('j_ATP_H', '-v_ATP_H'), ('j_ADP_H', 'v_ATP_H'), ('j_Pi_H', 'v_ATP_H'),
]
ha_roles = {v[0]: 'constant' for v in ha_vars if v[0] != 't'}
ha_roles.update({n: 'boundary_condition' for n in ('V', 'E_H', 'H_i', 'A_mem', 'ATP', 'ADP', 'Pi', 'm_expr')})
ha_roles.update({n: 'global_constant' for n in ('R', 'T', 'F')})
ha_roles.update({v[0]: 'variable' for v in ha_vars if v[2] == 'out'})
write_version('cell/ion_channels', 'i_H_ATPase', 'i_H_ATPase_FoF1', ha_vars, ha_eqs,
              {'required_citations': ['argus2026unpublished'],
               'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'protein_expression', 'variables': ['m_expr']},
                                 {'port_type': 'ATP_hydrolysis_port',
                                  'variables': ['ATP', 'ADP', 'Pi', 'j_ATP_H', 'j_ADP_H', 'j_Pi_H']},
                                 {'port_type': 'membrane_voltage', 'variables': ['V']},
                                 {'port_type': 'H_Nernst_potential', 'variables': ['E_H']},
                                 {'port_type': 'internal_H', 'variables': ['H_i']},
                                 {'port_type': 'membrane_current', 'variables': ['i_H_ATPase']},
                                 {'port_type': 'H_flux', 'variables': ['j_H_ATPase']},
                                 {'port_type': 'membrane_area', 'variables': ['A_mem']}],
               'roles': ha_roles},
              [OP('V'), ('E_H', E_H, 'Nernst potential at pH 7.6 in, 7.3 out', False),
               ('H_i', H_I, 'slonczewski1981ph; cytoplasmic pH 7.6', True),
               ('K_act_H', K_H_ACT, f'{PH}; half activation at pH 7.6: Kobayashi (kobayashi1985proton) found the enzyme '
                                    'regulates cytoplasmic pH at 7.6-7.8, more active as the cytoplasm acidifies', False),
               ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
               ('m_expr', 1.0, 'definitional; baseline expression (the calibrated cell)', False),
               ('rho_H_ATPase', RHO_HATP, 'calibrated: the H+-ATPase uses 25% of the cell\'s ATP turnover at rest '
                                          '(placeholder share)', False),
               ('r_H_ATPase', FATP['r_H_ATPase'], f'{PH}; ATP per enzyme per second (order 100/s)', False),
               ('n_H', N_H_ATPASE, 'jiang2001preferred; E. coli c10 ring and 3 catalytic sites: 10/3 H+ per ATP', True),
               ('ATP', ATP, 'bennett2009absolute; ATP 9.6 mM in E. coli', True),
               ('ADP', ADP, f'{PH}; about 0.56 mM in E. coli (to confirm)', False), ('Pi', PI, f'{PH}; near 10 mM', False),
               ('K_m_ATP_H', FATP['K_m_ATP_H'], f'{PH}; affinity for ATP, order 0.1-0.5 mM', False)] + GEN + [
               ('dG0_hyd', DG0_HYD, f'{PH}; textbook standard Gibbs energy of ATP hydrolysis (-30.5 kJ/mol)', False),
               ('c_std', C_STD, 'definitional; standard concentration 1 M', True)],
              {**RUN1, 'outputs': ['i_H_ATPase', 'v_ATP_H', 'j_H_pumped', 'j_H_ATPase', 'j_ATP_H', 'j_ADP_H',
                                   'j_Pi_H', 'A_ATP', 'W_H', 'f_thermo'], 'run_parameters': {},
               'invariants': [{'description': 'stoichiometry: n_H H+ out and 1 ATP -> ADP + Pi per turnover; the current '
                                              'is the pumped charge, i = -F j_H_pumped (exact)',
                               'expr': 'np.allclose(j_H_pumped, n_H*j_ATP_H, rtol=1e-12, atol=0) and '
                                       'np.allclose(j_ADP_H, -j_ATP_H, rtol=1e-12, atol=0) and '
                                       'np.allclose(j_Pi_H, -j_ATP_H, rtol=1e-12, atol=0) and '
                                       'np.allclose(i_H_ATPase, -F*j_H_pumped*1e9, rtol=1e-12, atol=0)'},
                              {'description': 'second law: it hydrolyses ATP only while the ATP free energy exceeds the '
                                              'proton work, v (A_ATP - W_H) >= 0 (and synthesises ATP otherwise)',
                               'expr': 'np.all(v_ATP_H*(A_ATP - W_H) >= 0)'}],
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'The F1Fo ATPase of a fermenting bacterium, which hydrolyses the ATP made by glycolysis to pump protons '
                  'out and so sets the proton motive force and regulates cytoplasmic pH (Enterococcus hirae / Streptococcus '
                  'faecalis, kobayashi1985proton). n_H = 10/3 H+ per ATP (the E. coli c10 rotor over 3 catalytic sites, '
                  'jiang2001preferred). Rate i = k ATP/(K_m + ATP) H_i/(K_act + H_i) f_thermo with f_thermo = 1 - exp(-(A_ATP '
                  '- W_H)/RT), W_H = n_H F (E_H - V): the same thermodynamic brake as the Na/K pump (i_NaK '
                  'ATP_ArgusUNPUBLISHED_v01), so below equilibrium it runs backwards and synthesises ATP, as the enzyme is '
                  'reversible. The H+ activation term is a placeholder for the pH regulation Kobayashi reports. Each ATP '
                  'hydrolysed also releases one H+, which the cytosol (energy_metabolism) accounts for.')},
              ['argus2026unpublished', 'kobayashi1985proton', 'jiang2001preferred', 'slonczewski1981ph',
               'bennett2009absolute'], 'FoF1_ArgusUNPUBLISHED_v01')

# --- exchangers: linear flux-force laws in the Nernst potentials
def exchanger(module_type, component, version, desc_vars, drive, drive_desc, j_vars, current, params, notes, bib,
              ports, invariants):
    vars_ = [('t', 'second', 'in', None, None)] + desc_vars + [
        ('A_mem', 'um2', 'in', None, None), ('g_spec', 'microS_per_um2', 'in', None, 'flux-force coefficient per area'),
        ('m_expr', 'dimensionless', 'in', None, 'expression, 1 = the calibrated cell (protein_expression port)'),
        ('F', 'C_per_mol', 'in', None, None),
        ('drive', 'milliV', 'out', None, drive_desc),
        ('r_cycle', 'mol_per_s', 'out', None, 'exchange cycles per second (forward positive)')] + j_vars
    eqs = [('drive', drive), ('r_cycle', 'A_mem * g_spec * m_expr * drive / F * u(1e-9, A_per_nanoA)')] + current
    roles = {v[0]: 'boundary_condition' for v in desc_vars}
    roles.update({'A_mem': 'boundary_condition', 'm_expr': 'boundary_condition', 'g_spec': 'constant',
                  'F': 'global_constant', 'drive': 'variable',
                  'r_cycle': 'variable'})
    roles.update({v[0]: 'variable' for v in j_vars})
    write_version('cell/ion_channels', module_type, component, vars_, eqs,
                  {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
                   'general_ports': ports + [{'port_type': 'membrane_area', 'variables': ['A_mem']}, {'port_type': 'protein_expression', 'variables': ['m_expr']}],
                   'roles': roles},
                  params + [('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
                            ('m_expr', 1.0, 'definitional; baseline expression (the calibrated cell)', False),
                            ('F', F, 'definitional; as the library\'s SN modules', True)],
                  {**RUN1, 'outputs': ['drive', 'r_cycle'] + [v[0] for v in j_vars], 'run_parameters': {},
                   'invariants': invariants + [{'description': 'second law: the exchange runs down its free energy, '
                                                               'r drive >= 0', 'expr': 'np.all(r_cycle*drive >= 0)'}],
                   'stability': STAB},
                  {'reviewed': False, 'notes': notes}, bib, version)


exchanger('i_NaH', 'i_NaH_NhaA', 'NhaA_ArgusUNPUBLISHED_v01',
          [('V', 'milliV', 'in', None, None), ('E_Na', 'milliV', 'in', None, None), ('E_H', 'milliV', 'in', None, None)],
          '-(V + E_Na - 2.0 * E_H)', 'free energy released per cycle over F: -(V + E_Na - 2 E_H)',
          [('i_NaH', 'nanoA', 'out', None, 'net +1 charge in per cycle (outward positive)'),
           ('j_NaH_Na', 'mol_per_s', 'out', None, 'Na+ into the cell'), ('j_NaH_H', 'mol_per_s', 'out', None, 'H+ into the cell')],
          [('i_NaH', '-r_cycle * F * u(1e9, nanoA_per_A)'), ('j_NaH_Na', '-r_cycle'), ('j_NaH_H', '2.0 * r_cycle')],
          [OP('V'), ('E_Na', E_NA, 'Nernst potential at the default concentrations', False),
           ('E_H', E_H, 'Nernst potential at pH 7.6 in, 7.3 out', False),
           ('g_spec', GS_NHA, 'calibrated: NhaA returns 30% of the pumped H+ (placeholder share), removing the Na+ that '
                              'leaks in', False)],
          'NhaA, the main Na+/H+ antiporter of E. coli: 1 Na+ out for 2 H+ in (taglicht1993proton), electrogenic (+1 '
          'charge in per cycle), driven by the proton motive force. The cycle free energy over F is V + E_Na - 2 E_H (Na+ '
          'out costs E_Na - V, each H+ in gives V - E_H), and the rate is the linear flux-force law r = A_mem g_spec '
          'drive/F, the exchanger analogue of an Ohmic leak (placeholder kinetics; thermodynamics exact). It cannot hold '
          'Na+ down when V + E_Na > 2 E_H: in modern seawater (pH 8.1) at -105 mV with 20 mM Na+ inside it would run '
          'backwards, which is why marine bacteria add primary Na+ pumps. NhaA is strongly activated at alkaline pH; '
          'that regulation is not modelled.',
          ['argus2026unpublished', 'taglicht1993proton'],
          [{'port_type': 'membrane_voltage', 'variables': ['V']}, {'port_type': 'Na_Nernst_potential', 'variables': ['E_Na']},
           {'port_type': 'H_Nernst_potential', 'variables': ['E_H']},
           {'port_type': 'membrane_current', 'variables': ['i_NaH']}, {'port_type': 'Na_flux', 'variables': ['j_NaH_Na']},
           {'port_type': 'H_flux', 'variables': ['j_NaH_H']}],
          [{'description': 'stoichiometry: 1 Na+ out, 2 H+ in; current = net charge in, i = -F (j_Na + j_H) (exact)',
            'expr': 'np.allclose(j_NaH_H, -2*j_NaH_Na, rtol=1e-12, atol=0) and '
                    'np.allclose(i_NaH, -F*(j_NaH_Na + j_NaH_H)*1e9, rtol=1e-12, atol=0)'}])

exchanger('i_KH', 'i_KH_KefC', 'KefC_ArgusUNPUBLISHED_v01',
          [('E_K', 'milliV', 'in', None, None), ('E_H', 'milliV', 'in', None, None)],
          'E_H - E_K', 'free energy released per cycle over F: E_H - E_K (electroneutral)',
          [('j_KH_K', 'mol_per_s', 'out', None, 'K+ into the cell'), ('j_KH_H', 'mol_per_s', 'out', None, 'H+ into the cell')],
          [('j_KH_K', '-r_cycle'), ('j_KH_H', 'r_cycle')],
          [('E_K', E_K, 'Nernst potential at the default concentrations', False),
           ('E_H', E_H, 'Nernst potential at pH 7.6 in, 7.3 out', False),
           ('g_spec', GS_KEF, 'calibrated: KefC returns 20% of the pumped H+ (placeholder share), exporting the K+ that '
                              'Trk brings in', False)],
          'KefC-type K+/H+ exchanger (gulati2024structure): 1 K+ out for 1 H+ in, electroneutral, so it is driven by the '
          'concentration gradients alone: free energy per cycle over F = E_H - E_K. With the membrane more negative than '
          'E_K, K+ keeps entering through the Trk channel; this exchanger returns it, so V can sit below E_K. Linear '
          'flux-force law r = A_mem g_spec (E_H - E_K)/F (placeholder kinetics; thermodynamics exact). In E. coli KefC is '
          'gated by glutathione (mostly closed); its basal activity stands in for all K+ efflux here.',
          ['argus2026unpublished', 'gulati2024structure'],
          [{'port_type': 'E_K_potential', 'variables': ['E_K']}, {'port_type': 'H_Nernst_potential', 'variables': ['E_H']},
           {'port_type': 'K_flux', 'variables': ['j_KH_K']}, {'port_type': 'H_flux', 'variables': ['j_KH_H']}],
          [{'description': 'stoichiometry: 1 K+ out, 1 H+ in, no net charge (exact)',
            'expr': 'np.allclose(j_KH_H, -j_KH_K, rtol=1e-12, atol=0)'}])

# --- i_H_flagellum: the proton-driven flagellar motor and the swimming it gives
DRIVE_H_REST = E_H - V_REST          # mV: proton motive force at rest (inward), what the motor runs on
fl_vars = [
    ('t', 'second', 'in', None, None),
    ('V', 'milliV', 'in', None, None), ('E_H', 'milliV', 'in', None, None),
    ('m_expr', 'dimensionless', 'in', None, 'number of flagella (protein_expression port, from the Flagellum pool)'),
    ('a_motor', 'dimensionless', 'in', None, 'motor switch: on (1) or off (0), set by the game'),
    ('tau_bundle', 'second', 'in', None, 'time constant of spinning up into a bundle or unwinding'),
    ('a_init', 'dimensionless', 'in', None, None),
    ('a_eff', 'dimensionless', 'out', 'a_init', 'motor and bundle activity: follows the switch with tau_bundle'),
    ('f_ref', 'per_s', 'in', None, 'rotation rate at the reference proton motive force'),
    ('pmf_ref', 'milliV', 'in', None, 'reference proton motive force'),
    ('n_H_rev', 'dimensionless', 'in', None, 'H+ through the stators per revolution'),
    ('v_one', 'm_per_s', 'in', None, 'swimming speed with one flagellum at f_ref'),
    ('n_half', 'dimensionless', 'in', None, 'flagella at which extra flagella add half as much again'),
    ('eta', 'pascal_second', 'in', None, 'viscosity of the water'),
    ('R_body', 'um', 'in', None, 'radius of the cell body (Stokes drag)'),
    ('N_A', 'per_mol', 'in', None, None), ('F', 'C_per_mol', 'in', None, None),
    ('xi', 'newton_s_per_m', 'out', None, 'Stokes drag coefficient of the cell body, 6 pi eta R'),
    ('F_one', 'newton', 'out', None, 'thrust of each flagellum (its share falls as flagella crowd: n_half)'),
    ('F_thrust', 'newton', 'out', None, 'thrust of all the flagella'),
    ('P_mech', 'watt', 'out', None, 'mechanical power of swimming alone, F_thrust v_swim (all of it ends as heat in the water)'),
    ('pmf', 'milliV', 'out', None, 'proton motive force driving H+ in: E_H - V'),
    ('f_rot', 'per_s', 'out', None, 'rotation rate: proportional to the PMF (E. coli motors)'),
    ('j_H_flag', 'mol_per_s', 'out', None, 'H+ into the cell through the motors (H_flux port)'),
    ('i_flag', 'nanoA', 'out', None, 'membrane current of those H+ (outward positive)'),
    ('P_flag', 'watt', 'out', None, 'proton-motive power the motors use'),
    ('v_swim', 'm_per_s', 'out', None, 'swimming speed: one flagellum gives v_one at f_ref; more saturate'),
]
fl_eqs = [
    ('pmf', 'E_H - V'),
    ('d/dt a_eff', '(a_motor - a_eff) / tau_bundle'),
    ('f_rot', 'a_eff * f_ref * pmf / pmf_ref'),
    ('j_H_flag', 'm_expr * f_rot * n_H_rev / N_A'),
    ('i_flag', '-j_H_flag * F * u(1e9, nanoA_per_A)'),
    ('P_flag', 'j_H_flag * F * pmf * u(1e-3, V_per_milliV)'),
    ('xi', '6.0 * pi * eta * R_body * u(1e-6, m_per_um)'),
    ('F_one', 'xi * v_one * f_rot / f_ref * (1.0 + n_half) / (m_expr + n_half)'),
    ('F_thrust', 'm_expr * F_one'),
    ('v_swim', 'F_thrust / xi'),
    ('P_mech', 'F_thrust * v_swim'),
]
fl_roles = {v[0]: ('variable' if v[2] == 'out' else 'constant') for v in fl_vars if v[0] != 't'}
fl_roles.update({n: 'boundary_condition' for n in ('V', 'E_H', 'm_expr', 'a_motor')})
fl_roles['F'] = 'global_constant'
write_version('cell/ion_channels', 'i_H_flagellum', 'i_H_flagellum_MotAB', fl_vars, fl_eqs,
              {'required_citations': ['argus2026unpublished'], 'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'membrane_voltage', 'variables': ['V']},
                                 {'port_type': 'H_Nernst_potential', 'variables': ['E_H']},
                                 {'port_type': 'membrane_current', 'variables': ['i_flag']},
                                 {'port_type': 'H_flux', 'variables': ['j_H_flag']},
                                 {'port_type': 'protein_expression', 'variables': ['m_expr']}],
               'roles': fl_roles},
              [OP('V'), ('E_H', E_H, 'Nernst potential at pH 7.6 in, 7.3 out', False),
               ('m_expr', par('n_flagella_init'), '', False),
               ('a_motor', 0.0, 'definitional; motor off until the game runs it', False),
               ('tau_bundle', par('tau_bundle'), '', False),
               ('a_init', 0.0, 'definitional; at rest the motors are still', False),
               ('f_ref', par('f_flag_ref'), '', False),
               ('pmf_ref', DRIVE_H_REST, 'calibrated: the proton motive force of the fermenter at rest (E_H - V)', False),
               ('n_H_rev', par('n_H_rev'), '', False), ('v_one', par('v_swim_one'), '', False),
               ('n_half', par('flagella_half'), '', False), ('eta', par('eta_water'), '', False),
               ('R_body', R_CELL, '', False), ('N_A', NA, '', False), ('F', F, '', False)],
              {**RUN1, 'sim_time': 3.0, 'outputs': ['pmf', 'a_eff', 'f_rot', 'j_H_flag', 'i_flag', 'P_flag', 'v_swim', 'xi',
                                                    'F_one', 'F_thrust', 'P_mech'],
               'run_parameters': {'a_motor': 1.0},
               'invariants': [
                   {'description': 'motor law (exact): f = a f_ref pmf/pmf_ref, j = m f n_H_rev / N_A, i = -F j, '
                                   'v = v_one (f/f_ref) m (1 + n_half)/(m + n_half)',
                    'expr': 'np.allclose(f_rot, a_eff*f_ref*(E_H - V)/pmf_ref, rtol=1e-12, atol=0) and '
                            'np.allclose(j_H_flag, m_expr*f_rot*n_H_rev/N_A, rtol=1e-12, atol=0) and '
                            'np.allclose(i_flag, -j_H_flag*F*1e9, rtol=1e-12, atol=0) and '
                            'np.allclose(v_swim, v_one*f_rot/f_ref*m_expr*(1 + n_half)/(m_expr + n_half), rtol=1e-12, atol=0)'},
                   {'description': 'the motors spin up smoothly: from rest a_eff = a_motor (1 - exp(-t/tau_bundle)) '
                                   '(relative 1e-6), and stays within 0-1',
                    'expr': 'np.allclose(a_eff, a_motor*(1 - np.exp(-t/tau_bundle)), rtol=1e-6, atol=1e-9) and '
                            'np.all(a_eff >= -1e-12) and np.all(a_eff <= 1 + 1e-12)', 'applies': 'run'},
                   {'description': 'second law: protons only run down the proton motive force, j pmf >= 0; the power '
                                   'used is never negative',
                    'expr': 'np.all(j_H_flag*pmf >= 0) and np.all(P_flag >= 0)'},
                   {'description': 'one flagellum at the reference PMF swims at v_one; speed rises with the number but '
                                   'saturates below v_one (1 + n_half)',
                    'expr': 'np.all(v_swim <= v_one*(1 + n_half)*np.abs(f_rot)/f_ref*(1 + 1e-12))'},
                   {'description': 'thrust and drag balance for a free cell: F_thrust = xi v_swim with xi = 6 pi eta R '
                                   '(Stokes), F_thrust = m F_one',
                    'expr': 'np.allclose(xi, 6*np.pi*eta*R_body*1e-6, rtol=1e-12, atol=0) and np.allclose(F_thrust, '
                            'xi*v_swim, rtol=1e-12, atol=0) and np.allclose(F_thrust, m_expr*F_one, rtol=1e-12, atol=0)'},
                   {'description': 'energy: swimming never does more work than the motors take from the proton motive '
                                   'force, P_mech <= P_flag (flagellar propulsion is ~0.1-2% efficient; the rest is heat)',
                    'expr': 'np.all(P_mech <= P_flag*(1 + 1e-12) + 1e-30)'}],
               'timestep': {'scheme': 'rk4', 'dts': [0.04, 0.02, 0.01], 't_end': 1.0, 'min_order': 0.8, 'tol': 1e-3,
                            'cvode_tol': 1e-2},
               'bc_sweep': {'ranges': {'m_expr': [0.0, 8.0], 'a_motor': [0.0, 1.0], 'V': [-150.0, -60.0]},
                            'rationale': 'zero to eight flagella, motor off to on, membrane potentials from hyperpolarised '
                                         'to a cell running down'},
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'The bacterial flagellar motor as a proton-driven rotor. The game switches the motors on or off '
                  '(a_motor); the motors and the bundle follow with time constant tau_bundle (a_eff), so thrust rises '
                  'and falls smoothly. Each flagellum turns at f = a_eff f_ref pmf/pmf_ref (rotation proportional to the proton motive force, '
                  'as measured in E. coli) and passes n_H_rev H+ per revolution into the cell (Meister et al. 1987: about '
                  '1,200), on the H_flux port to the proton pool and as an inward membrane current. The power is that '
                  'proton flux times the PMF. Swimming speed: one flagellum at f_ref gives v_one; more flagella add with '
                  'diminishing returns, v = v_one (f/f_ref) m (1 + n_half)/(m + n_half), as a thrust: each flagellum '
                  'pushes F_one, all of them F_thrust = xi v_swim against the Stokes drag xi = 6 pi eta R of the body (the '
                  'game moves the whole creature as one rigid body with these forces). Its mechanical power F v is far '
                  'below the proton-motive power P_flag (~0.2%, as flagellar propulsion is ~0.1-2% efficient): the rest '
                  'is heat in the motor, filament and water. The number of flagella m comes '
                  'from the Flagellum protein_pool (expression = count), so building one costs its ~1e7 amino acids. '
                  'Placeholders: f_ref, v_one, n_half. Tumbling, chemotaxis and the bundle are not modelled; the game '
                  'decides when the motor runs and in which direction the cell swims.')},
              ['argus2026unpublished', 'meister1987proton'], 'MotAB_ArgusUNPUBLISHED_v01')

# ------------------------------------------------------------------ 6b2. glucose carriers the cell can evolve (after level 1)
# Each sits between the medium (medium_face_port) and the cytosol (cytosol_face_port), like the Glf face, with its
# capacity per membrane area times the expression of its protein_pool. m = 1 is as many copies as the calibrated cell has
# Glf (n_0 below), so the carriers compare at equal copy number: capacity n_0 turnover / N_A. In the fermenter their pools
# start at 0, so the level-1 cell is unchanged.
GC = 'glucose_carrier'
V_GC_ATP, V_GC_H = 'ATP_ArgusUNPUBLISHED_v01', 'H_symport_ArgusUNPUBLISHED_v01'
N0_CARRIER = der(PROTEINS['Glf'][2], 'the calibrated cell\'s Glf copies (Glf capacity / turnover_Glf)', 'turnover_Glf')
# key: (label, module version, instance, turnover key, length key, footprint key)
CARRIERS = {
    'PTS': ('glucose phosphotransferase system (PtsG, EIIA-Glc, EI, HPr)', V_GC_ATP, 'PTS'),
    'ABC': ('binding-protein ABC glucose importer (MglABC)', V_GC_ATP, 'ABC'),
    'GalP': ('GalP-type glucose/H+ symporter (MFS)', V_GC_H, 'default'),
}


def carrier_jmax(key):
    return der(N0_CARRIER * par(f'turnover_{key}') / NA / A_MEM,
               f'n_0 x turnover_{key} / N_A / A_mem: m = 1 has as many copies as the calibrated cell\'s Glf',
               f'turnover_{key}', 'turnover_Glf')


face_side_vars = [
    ('t', 'second', 'in', None, None),
    ('C_glc_e', 'millimolar', 'in', None, 'medium side (entrance medium_face_port)'),
    ('C_lac_e', 'millimolar', 'in', None, 'carried for the port'),
    ('C_glc', 'millimolar', 'in', None, 'cytosol side (exit cytosol_face_port)'),
    ('C_lac', 'millimolar', 'in', None, 'carried for the port'), ('C_ATP', 'millimolar', 'in', None, None),
    ('C_ADP', 'millimolar', 'in', None, None), ('C_Pi', 'millimolar', 'in', None, None),
    ('A_mem', 'um2', 'in', None, None),
    ('m_expr', 'dimensionless', 'in', None, 'expression (protein_expression port): n_0 copies at 1'),
    ('jmax_glc', 'mol_per_s_per_um2', 'in', None, 'capacity per membrane area at m_expr = 1'),
    ('K_glc_t', 'millimolar', 'in', None, 'affinity for glucose outside'),
    ('R', 'J_per_mol_K', 'in', None, None), ('T', 'kelvin', 'in', None, None),
    ('J_glc', 'mol_per_s', 'out', None, 'glucose into the cell'),
    ('A_glc', 'J_per_mol', 'out', None, 'affinity of one transport cycle (> 0 drives glucose in)'),
    ('J_glc_e', 'mol_per_s', 'out', None, None), ('J_lac_e', 'mol_per_s', 'out', None, None),
    ('J_glc_c', 'mol_per_s', 'out', None, None), ('J_lac_c', 'mol_per_s', 'out', None, None),
    ('J_ATP_c', 'mol_per_s', 'out', None, 'ATP, ADP and Pi do not cross'), ('J_ADP_c', 'mol_per_s', 'out', None, None),
    ('J_Pi_c', 'mol_per_s', 'out', None, None),
]
face_side_eqs = [
    ('J_glc_e', '-J_glc'), ('J_lac_e', 'u(0.0, mol_per_s)'), ('J_glc_c', 'J_glc'), ('J_lac_c', 'u(0.0, mol_per_s)'),
    ('J_ATP_c', 'u(0.0, mol_per_s)'), ('J_ADP_c', 'u(0.0, mol_per_s)'), ('J_Pi_c', 'u(0.0, mol_per_s)'),
]
face_side_ports = {
    'entrance_ports': [{'port_type': 'medium_face_port', 'variables': ['C_glc_e', 'C_lac_e', 'J_glc_e', 'J_lac_e']}],
    'exit_ports': [{'port_type': 'cytosol_face_port',
                    'variables': ['C_glc', 'C_lac', 'C_ATP', 'C_ADP', 'C_Pi', 'J_glc_c', 'J_lac_c', 'J_ATP_c', 'J_ADP_c',
                                  'J_Pi_c']}],
}


def face_side_roles(extra_bcs, extra_constants):
    roles = {n: 'boundary_condition' for n in ('C_glc_e', 'C_lac_e', 'C_glc', 'C_lac', 'C_ATP', 'C_ADP', 'C_Pi', 'A_mem',
                                               'm_expr', *extra_bcs)}
    roles.update({n: 'constant' for n in ('jmax_glc', 'K_glc_t', *extra_constants)})
    roles.update({'R': 'global_constant', 'T': 'global_constant'})
    roles.update({n: 'variable' for n in ('J_glc', 'A_glc', 'J_glc_e', 'J_lac_e', 'J_glc_c', 'J_lac_c', 'J_ATP_c',
                                          'J_ADP_c', 'J_Pi_c')})
    return roles


def face_side_rows(key):
    return [('C_glc_e', GLC_E, '', False), ('C_lac_e', LAC_E, 'no lactate in the medium', False),
            ('C_glc', GLC_I, 'the fermenter cell\'s calibrated steady state', False),
            ('C_lac', LAC_I, 'the fermenter cell\'s calibrated steady state', False),
            ('C_ATP', ATP, '', False), ('C_ADP', ADP, '', False), ('C_Pi', PI, '', False),
            ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
            ('m_expr', 1.0, 'definitional; n_0 copies, as many as the calibrated cell has Glf', False),
            ('jmax_glc', carrier_jmax(key), '', False), ('K_glc_t', par(f'K_glc_{key}'), '', False),
            ('R', R, 'definitional; gas constant, as the library\'s SN modules', True),
            ('T', T, 'operating point: 20 C', False)]


FACE_SIDE_INV = ('np.array_equal(J_glc_e, -J_glc) and np.array_equal(J_glc_c, J_glc) and np.all(J_lac_e == 0) and '
                 'np.all(J_lac_c == 0) and np.all(J_ATP_c == 0) and np.all(J_ADP_c == 0) and np.all(J_Pi_c == 0)')
# ln of |C/c_std| + 1e-300, as the cytosol's affinities: an empty side sits at round-off
A_GLC_EXPR = 'R * T * (ln(abs(C_glc_e / c_std) + 1e-300) - ln(abs(C_glc / c_std) + 1e-300))'
CARRIER_SWEEP = {'ranges': {'C_glc_e': {'min': 1e-4, 'max': 50.0, 'points': 7, 'scale': 'log'},
                            'C_glc': [0.0, 5.0], 'm_expr': [0.0, 4.0]},
                 'rationale': 'glucose outside from 0.1 uM to 50 mM (the game\'s range, both sides of every Km), inside '
                              'from empty to above the outside at the low end (so the carrier runs backwards there), '
                              'expression from none to the game\'s maximum'}

# --- glucose_carrier ATP: uptake driven by n_P ATP equivalents (PTS: 1 PEP; ABC: 2 ATP)
gca_vars = face_side_vars + [
    ('n_P', 'dimensionless', 'in', None, 'ATP equivalents hydrolysed per glucose taken up'),
    ('K_m_ATP', 'millimolar', 'in', None, 'ATP for half-rate'),
    ('ATP', 'millimolar', 'in', None, 'cytosol ATP (ATP_hydrolysis_port)'), ('ADP', 'millimolar', 'in', None, None),
    ('Pi', 'millimolar', 'in', None, None),
    ('dG0_hyd', 'J_per_mol', 'in', None, None), ('c_std', 'millimolar', 'in', None, None),
    ('ATP_eq', 'millimolar', 'out', None, 'ATP at which hydrolysis would be at equilibrium with the present ADP and Pi'),
    ('A_ATP', 'J_per_mol', 'out', None, 'free energy of ATP hydrolysis'),
    ('J_ATP_cost', 'mol_per_s', 'out', None, 'ATP hydrolysed, n_P per glucose'),
    ('j_ATP', 'mol_per_s', 'out', None, None), ('j_ADP', 'mol_per_s', 'out', None, None),
    ('j_Pi', 'mol_per_s', 'out', None, None),
]
gca_eqs = [
    ('ATP_eq', 'ADP * Pi / c_std * exp(dG0_hyd / (R * T))'),
    ('A_ATP', '-dG0_hyd - R * T * ln(abs((ADP / c_std) * (Pi / c_std) / (ATP / c_std)) + 1e-300)'),
    # glc_e + n_P ATP <=> glc_i + n_P (ADP + Pi): saturating in glucose outside and in ATP, zero at equilibrium
    ('J_glc', 'jmax_glc * m_expr * A_mem * (C_glc_e * (ATP / K_m_ATP) ** n_P - C_glc * (ATP_eq / K_m_ATP) ** n_P) '
              '/ ((K_glc_t + C_glc_e) * (1.0 + ATP / K_m_ATP) ** n_P)'),
    ('A_glc', f'{A_GLC_EXPR} + n_P * A_ATP'),
    ('J_ATP_cost', 'n_P * J_glc'), ('j_ATP', '-J_ATP_cost'), ('j_ADP', 'J_ATP_cost'), ('j_Pi', 'J_ATP_cost'),
] + face_side_eqs
gca_roles = face_side_roles(('ATP', 'ADP', 'Pi'), ('n_P', 'K_m_ATP', 'dG0_hyd', 'c_std'))
gca_roles.update({n: 'variable' for n in ('ATP_eq', 'A_ATP', 'J_ATP_cost', 'j_ATP', 'j_ADP', 'j_Pi')})


def gca_rows(key, n_key):
    return face_side_rows(key) + [
        ('n_P', par(n_key), '', True), ('K_m_ATP', par('K_m_ATP_carrier'), '', False),
        ('ATP', ATP, '', False), ('ADP', ADP, '', False), ('Pi', PI, '', False),
        ('dG0_hyd', DG0_HYD, f'{PH}; textbook standard Gibbs energy of ATP hydrolysis (-30.5 kJ/mol)', False),
        ('c_std', C_STD, 'definitional; standard concentration 1 M', True)]


write_version('transport', GC, 'glucose_carrier_ATP', gca_vars, gca_eqs,
              {'required_citations': ['argus2026unpublished'], **face_side_ports,
               'general_ports': [{'port_type': 'membrane_area', 'variables': ['A_mem']},
                                 {'port_type': 'protein_expression', 'variables': ['m_expr']},
                                 {'port_type': 'ATP_hydrolysis_port',
                                  'variables': ['ATP', 'ADP', 'Pi', 'j_ATP', 'j_ADP', 'j_Pi']}],
               'roles': gca_roles},
              gca_rows('PTS', 'n_P_PTS'),
              {**RUN1, 'outputs': ['J_glc', 'A_glc', 'A_ATP', 'ATP_eq', 'J_ATP_cost', 'j_ATP', 'j_ADP', 'j_Pi', 'J_glc_e',
                                   'J_lac_e', 'J_glc_c', 'J_lac_c', 'J_ATP_c', 'J_ADP_c', 'J_Pi_c'],
               'run_parameters': {},
               'invariants': [
                   {'description': 'rate law (relative 1e-9): J = jmax m A_mem (C_e a^n - C_i a_eq^n)/((K + C_e)(1 + '
                                   'a)^n), a = ATP/K_m_ATP, a_eq = ATP_eq/K_m_ATP, ATP_eq = ADP Pi/(c_std K_eq_hyd)',
                    'expr': 'np.allclose(J_glc, jmax_glc*m_expr*A_mem*(C_glc_e*(ATP/K_m_ATP)**n_P - C_glc*(ADP*Pi/c_std*'
                            'np.exp(dG0_hyd/(R*T))/K_m_ATP)**n_P)/((K_glc_t + C_glc_e)*(1 + ATP/K_m_ATP)**n_P), '
                            'rtol=1e-9, atol=1e-12*abs(jmax_glc*A_mem*max(m_expr, 1e-30)))'},
                   {'description': 'stoichiometry (exact): n_P ATP -> ADP + Pi per glucose; glucose is conserved across '
                                   'the membrane; lactate, ATP, ADP and Pi do not cross',
                    'expr': 'np.array_equal(J_ATP_cost, n_P*J_glc) and np.array_equal(j_ATP, -J_ATP_cost) and '
                            'np.array_equal(j_ADP, J_ATP_cost) and np.array_equal(j_Pi, J_ATP_cost) and ' + FACE_SIDE_INV},
                   {'description': 'second law: glucose moves down the affinity of the coupled reaction, J A >= 0, A = '
                                   'RT ln(C_e/C_i) + n_P A_ATP',
                    'expr': 'np.all(J_glc*A_glc >= 0)'}],
               'bc_sweep': CARRIER_SWEEP, 'stability': STAB},
              {'reviewed': False, 'notes': (
                  'A glucose carrier driven by n_P ATP equivalents, between the medium (entrance medium_face_port) and '
                  'the cytosol (exit cytosol_face_port, as the Glf face cell_membrane_face adaptive), the ATP taken from '
                  'and ADP and Pi returned to the cytosol on the ATP_hydrolysis_port. Reaction glc_e + n_P ATP <=> glc_i '
                  '+ n_P (ADP + Pi), affinity A = RT ln(C_e/C_i) + n_P A_ATP, rate J = jmax m A_mem (C_e a^n - C_i '
                  'a_eq^n)/((K + C_e)(1 + a)^n) with a = ATP/K_m_ATP and a_eq = ATP_eq/K_m_ATP, ATP_eq = ADP Pi/(c_std '
                  'K_eq_hyd): Michaelis-Menten in glucose outside (J = jmax m A_mem C_e/(K + C_e) (a/(1 + a))^n with '
                  'nothing inside), saturating in ATP, zero at equilibrium and always of the sign of A, so J A >= 0. At '
                  'the cell\'s ATP (A_ATP ~49 kJ/mol) one ATP holds a 10^8-fold gradient: the carrier is irreversible in '
                  'practice. J_ATP_cost = n_P J is the ATP it spends; the cytosol counts the H+ and heat of that '
                  'hydrolysis. Instances: PTS (phosphotransferase system, Km 5 uM, n_P = 1: one PEP per glucose, '
                  'deutscher2006how) and ABC (binding-protein importer, Km 0.5 uM, ferenci1996adaptation; n_P = 2 ATP, '
                  'davidson2008structure). The lumped cytosol has no PEP or glucose 6-phosphate: the PTS PEP is '
                  'charged as one ATP, and its glucose enters as free glucose, which the lumped glycolysis phosphorylates '
                  'again (hexokinase), so a PTS glucose nets 1 ATP here where a real PTS cell keeps 2 (the PEP replaces '
                  'the hexokinase ATP). The driving forces agree: PEP to pyruvate with glucose phosphorylation releases '
                  'about 48 kJ/mol, one ATP at the cell\'s ratios about 49. m = 1 is as many copies as the calibrated '
                  'cell\'s Glf; turnover numbers and K_m_ATP are placeholders.')},
              ['argus2026unpublished', 'deutscher2006how', 'davidson2008structure', 'ferenci1996adaptation'], V_GC_ATP,
              extra_instances={'PTS': gca_rows('PTS', 'n_P_PTS'), 'ABC': gca_rows('ABC', 'n_ATP_ABC')})

# --- glucose_carrier H_symport: glucose with n_H H+ down the proton motive force (GalP-type MFS)
gch_vars = face_side_vars + [
    ('n_H', 'dimensionless', 'in', None, 'H+ per glucose'),
    ('V', 'milliV', 'in', None, None), ('E_H', 'milliV', 'in', None, 'H+ Nernst potential (H_Nernst_potential port)'),
    ('F', 'C_per_mol', 'in', None, None), ('c_std', 'millimolar', 'in', None, None),
    ('W_H', 'J_per_mol', 'out', None, 'free energy released by n_H H+ entering: n_H F (E_H - V)'),
    ('j_H', 'mol_per_s', 'out', None, 'H+ into the cell (H_flux port)'),
    ('i_glc_H', 'nanoA', 'out', None, 'current of the H+ carried in (positive outward)'),
]
gch_eqs = [
    ('W_H', 'n_H * F * (E_H - V) * u(1e-3, V_per_milliV)'),
    # glc_e + n_H H_e <=> glc_i + n_H H_i: Michaelis-Menten in glucose outside, zero at equilibrium
    ('J_glc', 'jmax_glc * m_expr * A_mem * (C_glc_e - C_glc * exp(-W_H / (R * T))) / (K_glc_t + C_glc_e)'),
    ('A_glc', f'{A_GLC_EXPR} + W_H'),
    ('j_H', 'n_H * J_glc'), ('i_glc_H', '-F * j_H * u(1e9, nanoA_per_A)'),
] + face_side_eqs
gch_roles = face_side_roles(('V', 'E_H'), ('n_H', 'c_std'))
gch_roles.update({'F': 'global_constant', 'W_H': 'variable', 'j_H': 'variable', 'i_glc_H': 'variable'})
write_version('transport', GC, 'glucose_carrier_H_symport', gch_vars, gch_eqs,
              {'required_citations': ['argus2026unpublished'], **face_side_ports,
               'general_ports': [{'port_type': 'membrane_area', 'variables': ['A_mem']},
                                 {'port_type': 'protein_expression', 'variables': ['m_expr']},
                                 {'port_type': 'membrane_voltage', 'variables': ['V']},
                                 {'port_type': 'H_Nernst_potential', 'variables': ['E_H']},
                                 {'port_type': 'membrane_current', 'variables': ['i_glc_H']},
                                 {'port_type': 'H_flux', 'variables': ['j_H']}],
               'roles': gch_roles},
              face_side_rows('GalP') + [
                  ('n_H', par('n_H_GalP'), '', True), OP('V'),
                  ('E_H', E_H, 'Nernst potential at pH 7.6 in, 7.3 out', False),
                  ('F', F, 'definitional; as the library\'s SN modules', True),
                  ('c_std', C_STD, 'definitional; standard concentration 1 M', True)],
              {**RUN1, 'outputs': ['J_glc', 'A_glc', 'W_H', 'j_H', 'i_glc_H', 'J_glc_e', 'J_lac_e', 'J_glc_c', 'J_lac_c',
                                   'J_ATP_c', 'J_ADP_c', 'J_Pi_c'],
               'run_parameters': {},
               'invariants': [
                   {'description': 'rate law (relative 1e-9): J = jmax m A_mem (C_e - C_i exp(-W_H/RT))/(K + C_e), W_H = '
                                   'n_H F (E_H - V)',
                    'expr': 'np.allclose(J_glc, jmax_glc*m_expr*A_mem*(C_glc_e - C_glc*np.exp(-n_H*F*(E_H - V)*1e-3/'
                            '(R*T)))/(K_glc_t + C_glc_e), rtol=1e-9, atol=1e-12*abs(jmax_glc*A_mem*max(m_expr, 1e-30)))'},
                   {'description': 'stoichiometry (exact): n_H H+ per glucose, carrying the current i = -F j_H; glucose is '
                                   'conserved across the membrane; lactate, ATP, ADP and Pi do not cross',
                    'expr': 'np.array_equal(j_H, n_H*J_glc) and np.allclose(i_glc_H, -F*j_H*1e9, rtol=1e-12, atol=0) and '
                            + FACE_SIDE_INV},
                   {'description': 'second law: J A >= 0, A = RT ln(C_e/C_i) + n_H F (E_H - V)',
                    'expr': 'np.all(J_glc*A_glc >= 0)'}],
               'bc_sweep': {**CARRIER_SWEEP, 'ranges': {**CARRIER_SWEEP['ranges'], 'V': [-200.0, 50.0]},
                            'rationale': CARRIER_SWEEP['rationale'] + '; membrane potential from hyperpolarised to '
                                                                      'depolarised past E_H, where the carrier reverses'},
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'A glucose/H+ symporter of the major facilitator superfamily (GalP, LacY type; pao1998major), between '
                  'the medium (entrance medium_face_port) and the cytosol (exit cytosol_face_port, as the Glf face '
                  'cell_membrane_face adaptive). Reaction glc_e + n_H H_e <=> glc_i + n_H H_i, affinity A = RT '
                  'ln(C_e/C_i) + W_H with W_H = n_H F (E_H - V) the proton motive force it uses (about 12 kJ/mol at '
                  '-105 mV and pH 7.6 in, 7.3 out: a 130-fold glucose gradient). Rate J = jmax m A_mem (C_e - C_i '
                  'exp(-W_H/RT))/(K + C_e): Michaelis-Menten in glucose outside with nothing inside, zero at '
                  'equilibrium and of the sign of A, so J A >= 0. The protons go to the cytosolic pool (j_H on H_flux, '
                  'E_H from H_Nernst_potential) and carry a membrane current i = -F j_H (membrane_current): the '
                  'H+-ATPase pumps them out again, about 0.3 ATP per glucose (3 ATP per 10 H+). Km 0.1 mM, the middle of '
                  'the MFS sugar/H+ symporters\' 10 uM-1 mM; m = 1 is as many copies as the calibrated cell\'s Glf; the '
                  'turnover number is a placeholder.')},
              ['argus2026unpublished', 'pao1998major'], V_GC_H)

# ------------------------------------------------------------------ 6c. proteins: one protein_pool per protein, and their synthesis
# Each protein of the cell is a protein_pool: its expression m (copies n = n_0 m) approaches the game's target m_target
# with time constant tau. Keeping the copies costs their turnover; raising m costs building the new copies, at most what the
# cell's protein synthesis allows (protein_expression: a capacity, throttled as ATP falls). The pool passes m to the
# protein it scales (a transporter's density or capacity, or the glycolytic enzymes' V_gly) on its protein_expression
# port, and its costs, build demand, daughter share and membrane area to protein_expression on protein_synthesis.
POOL_V = 'ArgusUNPUBLISHED_v01'
C_BUILD_EXPR = '(c_syn + c_poly)'
pool_vars = [
    ('t', 'second', 'in', None, None),
    ('m_init', 'dimensionless', 'in', None, None),
    ('m', 'dimensionless', 'out', 'm_init', 'expression: copies relative to n_0 (protein_expression port)'),
    ('m_target', 'dimensionless', 'in', None, 'expression the cell is moving to (set by the game)'),
    ('tau', 'second', 'in', None, 'time constant of the approach'),
    ('n_0', 'dimensionless', 'in', None, 'copies per cell at m = 1'), ('L', 'dimensionless', 'in', None,
                                                                       'amino acids per functional unit'),
    ('footprint', 'nm2', 'in', None, 'membrane area one copy takes (0 for a cytosolic protein)'),
    ('c_syn', 'dimensionless', 'in', None, '~P per amino acid made from glucose'),
    ('c_poly', 'dimensionless', 'in', None, '~P per peptide bond'),
    ('c_deg', 'dimensionless', 'in', None, '~P per residue for ATP-dependent proteolysis'),
    ('c_recycle', 'dimensionless', 'in', None, '~P per amino acid returned to the pool'),
    ('k_deg', 'per_s', 'in', None, 'fraction of the protein degraded and remade per second'),
    ('N_A', 'per_mol', 'in', None, None),
    ('K_m_ATP', 'millimolar', 'in', None, 'ATP for half-rate synthesis'),
    ('ATP_ref', 'millimolar', 'in', None, 'ATP at which the costs are stated (rest)'),
    ('ATP', 'millimolar', 'in', None, 'cytosol ATP (protein_synthesis port)'),
    ('s_build', 'dimensionless', 'in', None, 'share of the build demand the cell can serve (protein_synthesis port)'),
    ('n', 'dimensionless', 'out', None, 'copies per cell'),
    ('r_rel', 'per_s', 'out', None, '(m_target - m)/tau'),
    ('r_up', 'per_s', 'out', None, 'rising part of r_rel (copies to build)'),
    ('r_down', 'per_s', 'out', None, 'falling part of r_rel (copies left to degrade, no ATP)'),
    ('f_ATP', 'dimensionless', 'out', None, 'synthesis slows as ATP runs out: 1 at ATP_ref'),
    ('f_build', 'dimensionless', 'out', None, 'f_ATP capped at 1: building never runs faster than asked'),
    ('J_keep', 'mol_per_s', 'out', None, 'ATP to keep the copies (turnover)'),
    ('J_demand', 'mol_per_s', 'out', None, 'ATP the rise asks for: building r_up n_0 copies per second'),
    ('J_build', 'mol_per_s', 'out', None, 'ATP spent building'),
    ('J_ATP', 'mol_per_s', 'out', None, 'all ATP this protein uses'),
    ('M_build', 'mole', 'out', None, 'ATP to build the present copies from scratch (its share of a daughter)'),
    ('a_unit', 'um2', 'out', None, 'membrane area per unit of m'),
    ('A_protein', 'um2', 'out', None, 'membrane area the copies take'),
]
pool_eqs = [
    ('n', 'n_0 * m'),
    ('r_rel', '(m_target - m) / tau'),
    ('r_up', '(r_rel + abs(r_rel)) / 2.0'), ('r_down', '(r_rel - abs(r_rel)) / 2.0'),
    ('f_ATP', 'ATP / (K_m_ATP + ATP) * (K_m_ATP + ATP_ref) / ATP_ref'),
    ('f_build', '(f_ATP + 1.0 - abs(f_ATP - 1.0)) / 2.0'),
    ('J_keep', 'n * L * (c_poly + c_deg + c_recycle) * k_deg / N_A * f_ATP'),
    ('J_demand', f'r_up * n_0 * L * {C_BUILD_EXPR} / N_A'),
    ('J_build', 'J_demand * s_build * f_build'),
    ('J_ATP', 'J_keep + J_build'),
    ('M_build', f'n * L * {C_BUILD_EXPR} / N_A'),
    ('a_unit', 'n_0 * footprint * u(1e-6, um2_per_nm2)'),
    ('A_protein', 'a_unit * m'),
    # copies rise as fast as the ATP actually spent builds them; falling needs no ATP (degradation, not refunded)
    ('d/dt m', 'r_up * s_build * f_build + r_down'),
]
pool_roles = {v[0]: ('variable' if v[2] == 'out' else 'constant') for v in pool_vars if v[0] != 't'}
pool_roles.update({n: 'boundary_condition' for n in ('m_target', 'ATP', 's_build')})


def pool_rows(key):
    if key in CARRIERS:   # a carrier the level-1 cell lacks: none expressed until the game sets a target
        start = [('m_init', 0.0, 'definitional; the level-1 cell has none (the game unlocks it)', False),
                 ('m_target', 0.0, 'definitional; none until the game sets a target', False),
                 ('tau', par('tau_expression'), '', False), ('n_0', N0_CARRIER, '', False),
                 ('L', par(f'L_{key}'), '', False)]
    else:
        label, L, n, how = PROTEINS[key]
        start = [('m_init', 1.0, 'definitional; the calibrated cell', False),
                 ('m_target', 1.0, 'definitional; at rest the cell keeps its baseline', False),
                 ('tau', par('tau_expression'), '', False),
                 ('n_0', n, f'calibrated: {how}', False), ('L', L, '', False)]
    return start + [
        ('footprint', par(f'footprint_{key}'), '', False),
        ('c_syn', C_SYN, '', False), ('c_poly', C_POLY, '', False), ('c_deg', C_DEG, '', False),
        ('c_recycle', C_RECYCLE, '', False), ('k_deg', K_DEG, '', False), ('N_A', NA, '', False),
        ('K_m_ATP', par('K_m_ATP_prot'), '', False),
        ('ATP_ref', ATP, '', False), ('ATP', ATP, '', False),
        ('s_build', 1.0, 'definitional; synthesis serves the whole demand when unconnected', False),
    ]


pool_inv = [
    {'description': 'costs (exact): J_keep = n L (c_poly + c_deg + c_recycle) k_deg f_ATP / N_A, J_demand = r_up n_0 L '
                    '(c_syn + c_poly) / N_A, J_build = J_demand s_build min(f_ATP, 1), M_build = n L (c_syn + c_poly) / '
                    'N_A, A_protein = n_0 footprint m',
     'expr': 'np.allclose(J_keep, n_0*m*L*(c_poly + c_deg + c_recycle)*k_deg/N_A*f_ATP, rtol=1e-12, atol=0) and '
             'np.allclose(J_demand, np.maximum(m_target - m, 0)/tau*n_0*L*(c_syn + c_poly)/N_A, rtol=1e-9, atol=1e-30) and '
             'np.allclose(J_build, J_demand*s_build*np.minimum(f_ATP, 1), rtol=1e-9, atol=1e-30) and '
             'np.allclose(M_build, n_0*m*L*(c_syn + c_poly)/N_A, rtol=1e-12, atol=0) and '
             'np.allclose(A_protein, n_0*footprint*1e-6*m, rtol=1e-12, atol=0)'},
    {'description': 'the copies built are paid for: n_0 L (c_syn + c_poly) / N_A (m(t) - m(0)) = int J_build dt while m '
                    'rises (relative 1e-4: trapezoid integration of the logged flux)',
     'expr': 'abs(n_0*L*(c_syn + c_poly)/N_A*(m[-1] - m[0]) - np.trapz(J_build, t)) <= 1e-4*abs(np.trapz(J_build, t)) + '
             '1e-30', 'applies': 'run'},
    {'description': 'expression approaches its target and never overshoots: m moves monotonically towards m_target',
     'expr': 'np.all(np.diff(m)*np.sign(m_target - m[0]) >= -1e-12) and np.all((m - m_target)*(m[0] - m_target) >= -1e-12)'},
    {'description': 'costs are never negative', 'expr': 'np.all(J_keep >= 0) and np.all(J_build >= 0) and np.all(M_build >= 0)'},
]
pool_extra = {k: pool_rows(k) for k in [*PROTEINS, *CARRIERS]}
write_version('cell', 'protein_pool', 'protein_pool', pool_vars, pool_eqs,
              {'required_citations': ['argus2026unpublished'],
               'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'protein_expression', 'variables': ['m']},
                                 {'port_type': 'enzyme_expression', 'variables': ['m']},
                                 {'port_type': 'daughter_share', 'variables': ['M_build']},
                                 {'port_type': 'protein_synthesis',
                                  'variables': ['ATP', 's_build', 'J_demand', 'J_ATP', 'J_keep', 'J_build', 'M_build',
                                                'A_protein']}],
               'roles': pool_roles},
              pool_rows('NhaA'),
              {'sim_time': 60.0, 'pre_time': 0.0, 'dt': 0.1, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-14},
               'outputs': ['m', 'n', 'r_up', 'r_down', 'f_ATP', 'J_keep', 'J_demand', 'J_build', 'J_ATP', 'M_build',
                           'A_protein'],
               'run_parameters': {'m_target': 2.0},
               'invariants': pool_inv,
               'bc_sweep': {'ranges': {'m_target': [0.25, 4.0], 'ATP': [0.5, 12.0], 's_build': [0.0, 1.0]},
                            'rationale': 'expression over the range the game offers (x0.25-4); ATP from near exhaustion '
                                         'to above rest; synthesis from stalled to unconstrained'},
               'timestep': {'scheme': 'rk4', 'dts': [2.0, 1.0, 0.5], 't_end': 20.0, 'min_order': 0.8, 'tol': 1e-3,
                            'cvode_tol': 1e-2},
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'One protein of an adaptive cell as a pool of copies: n = n_0 m, with expression m approaching the '
                  'target m_target (set by the game) at rate (m_target - m)/tau. Keeping the copies costs their turnover, '
                  'J_keep = n L (c_poly + c_deg + c_recycle) k_deg / N_A (as protein_expression did for the whole cell), '
                  'slowed by f_ATP = ATP/(K_m + ATP) (K_m + ATP_ref)/ATP_ref as ATP runs out. Raising m asks for the '
                  'ATP to build r_up n_0 copies per second from glucose, J_demand = r_up n_0 L (c_syn + c_poly) / N_A; the '
                  'cell\'s protein synthesis (protein_expression adaptive) serves the share s_build of it, and the copies '
                  'rise only as fast as the ATP actually spent: dm/dt = r_up s_build min(f_ATP, 1) + r_down. Falling '
                  'needs no ATP (degradation) and is not refunded. M_build = n L (c_syn + c_poly) / N_A is the ATP to '
                  'build the present copies: its share of a daughter. A_protein = n footprint is the membrane area the '
                  'copies take (0 for cytosolic proteins). m goes to the protein it scales on the protein_expression '
                  'port (its density or capacity times m), or on enzyme_expression to the cytosol\'s enzymes; M_build '
                  'goes to biosynthesis (daughter_share, summed there); the costs, demand, daughter share and area go to '
                  'protein_expression on protein_synthesis, which returns ATP and s_build. Instances: one per protein of '
                  'the fermenter cell (copy numbers, lengths and footprints in each; see parameters.json).'),
               'instances': {**{k: f'{PROTEINS[k][0]}: {PROTEINS[k][3]}' for k in PROTEINS},
                             **{k: f'{CARRIERS[k][0]}: none at start; n_0 = the calibrated cell\'s Glf copies' for k in CARRIERS}}},
              ['argus2026unpublished', 'akashi2002metabolic', 'lynch2015bioenergetic', 'milo2010bionumbers'], POOL_V,
              extra_instances=pool_extra)

# --- protein_expression adaptive: the cell's protein synthesis, serving its pools
pe_vars = [
    ('t', 'second', 'in', None, None),
    ('ATP', 'millimolar', 'in', None, 'cytosol ATP (ATP_hydrolysis_port), passed to the pools'),
    ('ADP', 'millimolar', 'in', None, None), ('Pi', 'millimolar', 'in', None, None),
    # the pools' sums arrive on protein_synthesis (*_in, summed by libcuflynx) and are named here
    ('J_demand_in', 'mol_per_s', 'in', None, 'summed build demand of the pools (protein_synthesis port)'),
    ('J_ATP_in', 'mol_per_s', 'in', None, None), ('J_keep_in', 'mol_per_s', 'in', None, None),
    ('J_build_in', 'mol_per_s', 'in', None, None), ('M_build_in', 'mole', 'in', None, None),
    ('A_protein_in', 'um2', 'in', None, None),
    ('J_demand_total', 'mol_per_s', 'out', None, 'build demand of all pools'),
    ('J_ATP_total', 'mol_per_s', 'out', None, 'ATP use of all pools'),
    ('J_keep_total', 'mol_per_s', 'out', None, 'turnover cost of all pools'),
    ('J_build_total', 'mol_per_s', 'out', None, 'build spending of all pools'),
    ('M_build_protein', 'mole', 'out', None, 'ATP to build every pool\'s copies: the protein share of a daughter'),
    ('A_protein_total', 'um2', 'out', None, 'membrane area of all pools'),
    ('J_ref', 'mol_per_s', 'in', None, 'resting ATP turnover of the cell: the scale of synthesis and growth'),
    ('build_capacity', 'dimensionless', 'in', None, 'most ATP synthesis can take, in multiples of J_ref'),
    ('ATP_floor', 'millimolar', 'in', None, 'ATP at which synthesis and growth stop'),
    ('ATP_ref', 'millimolar', 'in', None, 'ATP at which synthesis runs at full speed (rest)'),
    ('J_eps', 'mol_per_s', 'in', None, 'numerical floor on the demand (so s_build is defined with nothing to build)'),
    ('A_mem', 'um2', 'in', None, 'in-plane membrane area (membrane_area port)'),
    ('surface_factor', 'dimensionless', 'in', None, 'whole cell surface / in-plane membrane area'),
    ('area_max_fraction', 'dimensionless', 'in', None, 'most of the cell surface the proteins may take'),
    ('x_ATP', 'dimensionless', 'out', None, '(ATP - ATP_floor)/(ATP_ref - ATP_floor)'),
    ('throttle', 'dimensionless', 'out', None, 'x_ATP clamped to 0-1 (the stringent response)'),
    ('J_cap', 'mol_per_s', 'out', None, 'what synthesis can spend now: build_capacity J_ref throttle'),
    ('s_build', 'dimensionless', 'out', None, 'share of the demand served: min(1, J_cap/J_demand_total)'),
    ('A_surface', 'um2', 'out', None, 'whole cell surface'),
    ('A_protein_max', 'um2', 'out', None, 'most membrane area the proteins may take'),
    ('f_area', 'dimensionless', 'out', None, 'share of A_protein_max in use'),
    ('j_ATP_prot', 'mol_per_s', 'out', None, None), ('j_ADP_prot', 'mol_per_s', 'out', None, None),
    ('j_Pi_prot', 'mol_per_s', 'out', None, None),
]
pe_eqs = [
    ('J_demand_total', 'J_demand_in'), ('J_ATP_total', 'J_ATP_in'), ('J_keep_total', 'J_keep_in'),
    ('J_build_total', 'J_build_in'), ('M_build_protein', 'M_build_in'), ('A_protein_total', 'A_protein_in'),
    ('x_ATP', '(ATP - ATP_floor) / (ATP_ref - ATP_floor)'),
    ('throttle', '((x_ATP + abs(x_ATP)) / 2.0 + 1.0 - abs((x_ATP + abs(x_ATP)) / 2.0 - 1.0)) / 2.0'),
    ('J_cap', 'build_capacity * J_ref * throttle'),
    # min(1, J_cap/J_demand) as (J_cap + eps)/(max(J_cap, J_demand) + eps): 1 with nothing to build, 0 with no capacity
    ('s_build', '(J_cap + J_eps) / ((J_cap + J_demand_total + abs(J_cap - J_demand_total)) / 2.0 + J_eps)'),
    ('A_surface', 'A_mem * surface_factor'),
    ('A_protein_max', 'area_max_fraction * A_surface'),
    ('f_area', 'A_protein_total / A_protein_max'),
    ('j_ATP_prot', '-J_ATP_total'), ('j_ADP_prot', 'J_ATP_total'), ('j_Pi_prot', 'J_ATP_total'),
]
pe_roles = {v[0]: ('variable' if v[2] == 'out' else 'constant') for v in pe_vars if v[0] != 't'}
pe_roles.update({n: 'boundary_condition' for n in ('ATP', 'ADP', 'Pi', 'A_mem')})
pe_roles.update({n: 'variable' for n in ('J_demand_in', 'J_ATP_in', 'J_keep_in', 'J_build_in', 'M_build_in',
                                         'A_protein_in')})
pe_params = [
    ('ATP', ATP, '', False), ('ADP', ADP, '', False), ('Pi', PI, '', False),
    ('J_demand_in', 0.0, 'no pool connected (an open sum port is 0)', False),
    ('J_ATP_in', 0.0, 'no pool connected', False), ('J_keep_in', 0.0, 'no pool connected', False),
    ('J_build_in', 0.0, 'no pool connected', False), ('M_build_in', 0.0, 'no pool connected', False),
    ('A_protein_in', 0.0, 'no pool connected', False),
    ('J_ref', der(R_TOT, 'maintenance_coefficient x cell volume ^ maintenance_exponent (Lynch & Marinov at 20 C)',
                  'maintenance_coefficient', 'maintenance_exponent'), '', False),
    ('build_capacity', par('build_capacity'), '', False),
    ('ATP_floor', par('build_ATP_floor'), '', False),
    ('ATP_ref', ATP, '', False),
    ('J_eps', 1e-30, 'definitional; numerical floor, 1e-30 mol/s (12 orders below the resting turnover)', False),
    ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
    ('surface_factor', der(A_SURFACE / A_MEM, '(2 pi R^2 + 2 pi R h)/(2 pi R h): a disc\'s whole surface over its rim',
                           'R_cell', 'h_cell'), '', False),
    ('area_max_fraction', par('protein_area_max_fraction'), '', False),
]
pe_inv = [
    {'description': 'throttle = clamp((ATP - ATP_floor)/(ATP_ref - ATP_floor), 0, 1); J_cap = build_capacity J_ref '
                    'throttle; s_build = min(1, J_cap/J_demand_total) (exact up to the 1e-30 mol/s floor)',
     'expr': 'np.allclose(throttle, np.clip((cyto__C_ATP - ATP_floor)/(ATP_ref - ATP_floor), 0, 1), rtol=1e-12, '
             'atol=1e-15) and np.allclose(J_cap, build_capacity*J_ref*throttle, rtol=1e-12, atol=0) and '
             'np.allclose(s_build, np.where(J_demand_total > J_cap, (J_cap + J_eps)/(J_demand_total + J_eps), 1.0), '
             'rtol=1e-9, atol=1e-12)'},
    {'description': 'the pools never spend more than synthesis can: J_build_total <= J_cap (relative 1e-9)',
     'expr': 'np.all(J_build_total <= J_cap*(1 + 1e-9) + 1e-30)', 'applies': 'run'},
    {'description': 'the ATP port takes exactly what the pools use: j_ATP = -J_ATP_total = -(J_keep_total + '
                    'J_build_total) (exact; relative 1e-9 with pools)',
     'expr': 'np.allclose(j_ATP_prot, -J_ATP_total, rtol=1e-12, atol=0) and np.allclose(J_ATP_total, J_keep_total + '
             'J_build_total, rtol=1e-9, atol=1e-30)', 'applies': 'run'},
]
write_version('cell', 'protein_expression', 'protein_expression_adaptive', pe_vars, pe_eqs,
              {'required_citations': ['argus2026unpublished'],
               'entrance_ports': [], 'exit_ports': [],
               'general_ports': [{'port_type': 'ATP_hydrolysis_port',
                                  'variables': ['ATP', 'ADP', 'Pi', 'j_ATP_prot', 'j_ADP_prot', 'j_Pi_prot']},
                                 {'port_type': 'membrane_area', 'variables': ['A_mem']},
                                 {'port_type': 'protein_synthesis',
                                  'variables': ['ATP', 's_build', 'J_demand_in', 'J_ATP_in', 'J_keep_in', 'J_build_in',
                                                'M_build_in', 'A_protein_in'],
                                  'multi_port': ['True', 'True', 'sum', 'sum', 'sum', 'sum', 'sum', 'sum']}],
               'roles': pe_roles},
              pe_params,
              {'sim_time': 60.0, 'pre_time': 0.0, 'dt': 0.1, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-14},
               'outputs': ['throttle', 'J_cap', 's_build', 'J_demand_total', 'J_build_total', 'J_keep_total',
                           'J_ATP_total', 'j_ATP_prot', 'M_build_protein', 'A_protein_total', 'f_area', 'x_ATP',
                           'pool_a/m', 'pool_b/m', 'cyto/C_ATP'],
               # two pools rising at once ask for more than synthesis can give: both are served the same share. A
               # cytosol supplies the ATP (the pools read it through this module); with nothing feeding it, its ATP
               # runs down over the minute, so the throttle is exercised from full speed to stopped
               'harness': {
                   'module_array': [['pool_a', POOL_V, 'protein_pool', '', 'mod', 'H_ATPase'],
                                    ['pool_b', POOL_V, 'protein_pool', '', 'mod', 'proteome'],
                                    ['cyto', V_NEW, 'energy_metabolism', '', 'mod', 'default'],
                                    ['mod', V_NEW, 'protein_expression', 'pool_a pool_b cyto', '', 'default']],
                   'parameters': [['m_target_pool_a', 'dimensionless', 3.0, 'harness: triple the H+-ATPase'],
                                  ['m_target_pool_b', 'dimensionless', 1.2, 'harness: a fifth more of the proteome']]},
               'run_parameters': {},
               'invariants': pe_inv,
               'bc_sweep': {'ranges': {'build_capacity': [0.5, 10.0], 'ATP_floor': [1.0, 8.0]},
                            'exclude': ['ATP', 'ADP', 'Pi'],   # supplied by the harness cytosol
                            'rationale': 'synthesis capacity from below to well above the two pools\' demand; the ATP '
                                         'floor over the range where the cytosol\'s falling ATP crosses it'},
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'The protein synthesis of an adaptive cell, serving its protein_pools (one per protein, on the '
                  'protein_synthesis port: their ATP use, turnover, build spending, daughter share and membrane area are '
                  'summed here, and each receives the cytosol ATP and s_build). Building is limited to J_cap = '
                  'build_capacity J_ref throttle: a growing cell builds protein at several times its maintenance rate, '
                  'and building backs off linearly as ATP falls from rest (ATP_ref) to ATP_floor (the stringent '
                  'response), so it never starves the cell. When the pools ask for more, each gets the same share s_build '
                  '= J_cap/J_demand_total. The summed ATP goes to the cytosol on the ATP_hydrolysis_port. The membrane '
                  'area the pools take is compared with the most the cell surface allows (area_max_fraction of it); the '
                  'game keeps targets within it. Biosynthesis throttles growth the same way (its own ATP_floor and J_ref, '
                  'the same values of parameters.json).')},
              ['argus2026unpublished', 'lynch2015bioenergetic'], V_NEW)


# --- adaptive versions of library channels: their density times the expression m_expr (protein_expression port)
def adaptive_library_version(module_type, src_version, density, note, defaults):
    """Copies a library version with its density scaled by m_expr; defaults: {variable: (value, ref, sourced)}."""
    import shutil
    from supermodule_writer import module_dir
    src = module_dir(module_type, src_version)
    dst = os.path.join(os.path.dirname(src), V_NEW)
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(os.path.join(dst, 'instances', 'default'))
    s_stem, d_stem = f'{module_type}_{src_version}', f'{module_type}_{V_NEW}'
    cellml = open(os.path.join(src, f'{s_stem}_modules.cellml')).read().replace(s_stem, d_stem)
    cfg_src = json.load(open(os.path.join(src, f'{s_stem}_modules_config.json')))[0]
    # a component name of its own: libcuflynx finds a module's component by name across the library
    comp_src, comp_dst = cfg_src['component_type'], f'{module_type}_adaptive'
    assert cellml.count(f'<component name="{comp_src}">') == 1
    cellml = cellml.replace(f'<component name="{comp_src}">', f'<component name="{comp_dst}">')
    decl = f'<variable name="{density}" public_interface="in" units="per_um2"/>'
    assert cellml.count(decl) == 1
    cellml = cellml.replace(decl, decl + '\n        <variable name="m_expr" public_interface="in" units="dimensionless"/>')
    ci = f'<ci>{density}</ci>'
    n_ci = cellml.count(ci)
    assert n_ci == 1, f'{module_type}: {density} appears {n_ci} times in the maths'
    cellml = cellml.replace(ci, f'<apply><times/>{ci}<ci>m_expr</ci></apply>')
    open(os.path.join(dst, f'{d_stem}_modules.cellml'), 'w').write(cellml)
    shutil.copy(os.path.join(src, f'{s_stem}_units.cellml'), os.path.join(dst, f'{d_stem}_units.cellml'))
    cfg = json.load(open(os.path.join(src, f'{s_stem}_modules_config.json')))
    cfg[0].update({'module_subtype': V_NEW, 'component_file': f'{d_stem}_modules.cellml', 'component_type': comp_dst})
    cfg[0]['general_ports'] = (cfg[0].get('general_ports') or []) + [{'port_type': 'protein_expression',
                                                                       'variables': ['m_expr']}]
    cfg[0]['variables_and_units'].append(['m_expr', 'dimensionless', 'access', 'boundary_condition'])
    cfg[0]['required_citations'] = sorted(set(cfg[0].get('required_citations', [])) | {'argus2026unpublished'})
    json.dump(cfg, open(os.path.join(dst, f'{d_stem}_modules_config.json'), 'w'), indent=2)
    rows = []
    for r in csv.DictReader(open(os.path.join(src, 'instances', 'default', 'default_parameters.csv'))):
        if r['variable_name'] in defaults:
            v, ref, src_ok = defaults[r['variable_name']]
            r.update(value=repr(float(v)), data_reference=ref, sourced='yes' if src_ok else 'no')
        rows.append(r)
    rows.append({'variable_name': 'm_expr', 'units': 'dimensionless', 'value': '1.0',
                 'data_reference': 'definitional; baseline expression (the calibrated cell)', 'sourced': 'no'})
    with open(os.path.join(dst, 'instances', 'default', 'default_parameters.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['variable_name', 'units', 'value', 'data_reference', 'sourced'],
                           lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    ver = json.load(open(os.path.join(src, f'{s_stem}_verification_config.json')))
    ver['version'] = V_NEW
    for inv in ver.get('invariants', []):
        inv['expr'] = inv['expr'].replace(density, f'({density}*m_expr)')
    json.dump(ver, open(os.path.join(dst, f'{d_stem}_verification_config.json'), 'w'), indent=2)
    tests = yaml.safe_load(open(os.path.join(src, f'{s_stem}_tests.yaml')))
    tests = {'module_type': module_type, 'version': V_NEW, 'reviewed': False,
             'description': tests.get('description', ''),
             'notes': (f'{src_version} with its density {density} times the expression m_expr (protein_expression '
                       f'port; 1 = the calibrated cell), for the generic adaptive cell: {note} Written by '
                       'tools/build_adaptive_cell.py from the library version; the maths is otherwise unchanged.')}
    yaml.safe_dump(tests, open(os.path.join(dst, f'{d_stem}_tests.yaml'), 'w'), sort_keys=False, width=120)
    shutil.copy(os.path.join(src, f'{s_stem}_references.bib'), os.path.join(dst, f'{d_stem}_references.bib'))
    bib = open(os.path.join(dst, f'{d_stem}_references.bib')).read()
    if '{argus2026unpublished,' not in bib:
        open(os.path.join(dst, f'{d_stem}_references.bib'), 'a').write('\n' + bib_text(['argus2026unpublished']))
    return dst


from module_writer import bib_text  # noqa: E402
adaptive_library_version('i_leak_K', 'Argus2026_v01', 'rho_leak_K', 'here the Trk K+ uptake channel of the fermenter.', {})
adaptive_library_version('i_NaCa', 'Argus2026_v01', 'rho_NaCa', 'here the archaeal NCX_Mj-type exchanger.', {})

# ------------------------------------------------------------------ 6d. biosynthesis: RNA, ribosomes, lipids, damage, growth
bio_vars = [
    ('t', 'second', 'in', None, None),
    ('ATP', 'millimolar', 'in', None, 'cytosol ATP (ATP_hydrolysis_port)'), ('ADP', 'millimolar', 'in', None, None),
    ('Pi', 'millimolar', 'in', None, None),
    ('A_mem', 'um2', 'in', None, 'in-plane membrane area (membrane_area port)'),
    ('surface_factor', 'dimensionless', 'in', None, 'whole cell surface / in-plane membrane area'),
    ('n_mRNA', 'dimensionless', 'in', None, 'mRNA molecules per cell'), ('L_mRNA', 'dimensionless', 'in', None, 'nt'),
    ('t_half_mRNA', 'second', 'in', None, 'mRNA half-life'),
    ('n_ribo', 'dimensionless', 'in', None, 'ribosomes per cell'), ('L_rRNA', 'dimensionless', 'in', None, 'nt per ribosome'),
    ('k_ribo', 'per_s', 'in', None, 'rRNA turnover'),
    ('c_nt_recycle', 'dimensionless', 'in', None, '~P to remake a nucleotide from its NMP'),
    ('c_nt_denovo', 'dimensionless', 'in', None, '~P to make a nucleotide from glucose'),
    ('area_per_lipid', 'um2', 'in', None, None), ('c_lipid', 'dimensionless', 'in', None, '~P per phospholipid made'),
    ('k_lipid', 'per_s', 'in', None, 'lipid turnover'),
    ('genome_nt', 'dimensionless', 'in', None, 'nucleotides in one genome (both strands)'),
    ('N_A', 'per_mol', 'in', None, None),
    ('K_m_ATP_bio', 'millimolar', 'in', None, None), ('ATP_ref', 'millimolar', 'in', None, None),
    ('T', 'kelvin', 'in', None, None), ('T_ref_dmg', 'kelvin', 'in', None, None),
    ('dT_Q10', 'kelvin', 'in', None, None), ('Q10_dmg', 'dimensionless', 'in', None, None),
    ('r_dmg', 'per_s', 'in', None, 'damage accrual: fraction of a lethal load per second'),
    ('r_rep_max', 'per_s', 'in', None, 'repair capacity at full repair level'),
    ('K_D', 'dimensionless', 'in', None, 'damage at half repair rate'),
    ('C_repair', 'mole', 'in', None, 'ATP to repair one whole lethal load of damage'),
    ('repair_level', 'dimensionless', 'in', None, 'how much the cell invests in repair, 0-1 (set by the game)'),
    ('growth_level', 'dimensionless', 'in', None, 'ATP the cell puts into a daughter, in multiples of J_ref (set by the '
                                                 'game)'),
    ('J_ref', 'mol_per_s', 'in', None, 'resting ATP turnover of the cell: growth is asked in multiples of it'),
    ('ATP_floor', 'millimolar', 'in', None, 'ATP at which growth stops (as protein synthesis)'),
    ('x_ATP', 'dimensionless', 'out', None, '(ATP - ATP_floor)/(ATP_ref - ATP_floor)'),
    ('throttle', 'dimensionless', 'out', None, 'x_ATP clamped to 0-1: growth backs off as ATP falls (stringent response)'),
    ('M_build_in', 'mole', 'in', None, 'the pools\' build costs, summed on the daughter_share port'),
    ('M_build_protein', 'mole', 'out', None, 'ATP to build the cell\'s proteins'),
    ('M_glyc_ATP', 'mole', 'in', None, 'ATP the glycogen store would give (glycogen_store port)'),
    ('J_growth', 'mol_per_s', 'out', None, 'ATP asked for the daughter: growth_level J_ref throttle'),
    ('M_daughter', 'mole', 'out', None, 'ATP to build a daughter: proteins, RNA, lipids and DNA'),
    ('growth_progress', 'dimensionless', 'out', None, 'share of a daughter built'),
    ('glycogen_progress', 'dimensionless', 'out', None, 'share of a daughter the glycogen store could finish'),
    ('D_init', 'dimensionless', 'in', None, None), ('P_init', 'mole', 'in', None, None),
    ('D', 'dimensionless', 'out', 'D_init', 'damage (1 = lethal)'),
    ('P_growth', 'mole', 'out', 'P_init', 'ATP invested in the daughter so far'),
    ('n_lipid', 'dimensionless', 'out', None, 'phospholipids in the membrane (both leaflets of the whole surface)'),
    ('f_ATP_bio', 'dimensionless', 'out', None, 'synthesis slows as ATP runs out: 1 at ATP_ref'),
    ('f_T_dmg', 'dimensionless', 'out', None, None),
    ('J_RNA', 'mol_per_s', 'out', None, 'mRNA turnover'), ('J_ribo', 'mol_per_s', 'out', None, 'ribosome (rRNA) upkeep'),
    ('J_lipid', 'mol_per_s', 'out', None, 'membrane lipid turnover, scales with membrane area'),
    ('R_repair', 'per_s', 'out', None, 'damage repaired per second'),
    ('J_repair', 'mol_per_s', 'out', None, 'ATP spent on repair'),
    ('J_growth_used', 'mol_per_s', 'out', None, 'ATP actually put into the daughter'),
    ('J_ATP_bio', 'mol_per_s', 'out', None, 'all ATP used here'),
    ('M_build_bio', 'mole', 'out', None, 'ATP to build a daughter\'s RNA, lipids and DNA (proteins: protein_expression)'),
    ('j_ATP_bio', 'mol_per_s', 'out', None, None), ('j_ADP_bio', 'mol_per_s', 'out', None, None),
    ('j_Pi_bio', 'mol_per_s', 'out', None, None),
]
bio_eqs = [
    ('n_lipid', '2.0 * A_mem * surface_factor / area_per_lipid'),
    ('f_ATP_bio', 'ATP / (K_m_ATP_bio + ATP) * (K_m_ATP_bio + ATP_ref) / ATP_ref'),
    ('f_T_dmg', 'Q10_dmg ** ((T - T_ref_dmg) / dT_Q10)'),
    ('J_RNA', 'n_mRNA * L_mRNA * c_nt_recycle * ln(2.0) / t_half_mRNA / N_A * f_ATP_bio'),
    ('J_ribo', 'n_ribo * L_rRNA * c_nt_recycle * k_ribo / N_A * f_ATP_bio'),
    ('J_lipid', 'n_lipid * c_lipid * k_lipid / N_A * f_ATP_bio'),
    ('R_repair', 'repair_level * r_rep_max * D / (K_D + D) * f_ATP_bio'),
    ('J_repair', 'R_repair * C_repair'),
    ('x_ATP', '(ATP - ATP_floor) / (ATP_ref - ATP_floor)'),
    ('throttle', '((x_ATP + abs(x_ATP)) / 2.0 + 1.0 - abs((x_ATP + abs(x_ATP)) / 2.0 - 1.0)) / 2.0'),
    ('M_build_protein', 'M_build_in'),
    ('J_growth', 'growth_level * J_ref * throttle'),
    ('J_growth_used', 'J_growth * f_ATP_bio'),
    ('M_daughter', 'M_build_bio + M_build_protein'),
    ('growth_progress', 'P_growth / M_daughter'), ('glycogen_progress', 'M_glyc_ATP / M_daughter'),
    ('J_ATP_bio', 'J_RNA + J_ribo + J_lipid + J_repair + J_growth_used'),
    ('M_build_bio', '((n_mRNA * L_mRNA + n_ribo * L_rRNA + genome_nt) * c_nt_denovo + n_lipid * c_lipid) / N_A'),
    ('j_ATP_bio', '-J_ATP_bio'), ('j_ADP_bio', 'J_ATP_bio'), ('j_Pi_bio', 'J_ATP_bio'),
    ('d/dt D', 'r_dmg * f_T_dmg - R_repair'),
    ('d/dt P_growth', 'J_growth_used'),
]
bio_roles = {v[0]: ('variable' if v[2] == 'out' else 'constant') for v in bio_vars if v[0] != 't'}
bio_roles.update({n: 'boundary_condition' for n in ('ATP', 'ADP', 'Pi', 'A_mem', 'repair_level', 'growth_level',
                                                    'M_glyc_ATP')})
bio_roles['M_build_in'] = 'variable'   # summed from the pools
bio_roles['T'] = 'global_constant'
bio_params = [
    ('ATP', ATP, 'bennett2009absolute; 9.6 mM', True), ('ADP', ADP, f'{PH}; port only', False),
    ('Pi', PI, f'{PH}; port only', False),
    ('A_mem', A_MEM, 'geometry: in-plane membrane of the 2 um cell', False),
    ('surface_factor', A_SURFACE / A_MEM, 'geometry: a disc of radius 1 um and depth 2 um has 1.5x its in-plane area '
                                          'as surface', False),
    ('n_mRNA', BIO['n_mRNA'], 'milo2010bionumbers; about 1,380 mRNAs in an E. coli cell of about 1 um3, scaled by volume',
     True),
    ('L_mRNA', BIO['L_mRNA'], f'{PH}; mean mRNA length about 1,000 nt', False),
    ('t_half_mRNA', BIO['t_half_mRNA'], 'bernstein2002global; 80% of E. coli mRNA half-lives 3-8 min (median ~5 min)', True),
    ('n_ribo', BIO['n_ribo'], f'{PH}; 7,000-70,000 ribosomes per E. coli cell with growth rate; the slow end per um3', False),
    ('L_rRNA', BIO['L_rRNA'], 'definitional; E. coli 16S + 23S + 5S rRNA: 1542 + 2904 + 120 nt', True),
    ('k_ribo', BIO['k_ribo'], f'{PH}; rRNA is stable in growing cells; 1%/h', False),
    ('c_nt_recycle', BIO['c_nt_recycle'], 'definitional; NMP -> NDP -> NTP costs 2 ~P', False),
    ('c_nt_denovo', BIO['c_nt_denovo'], f'{PH}; of order 50 ~P per nucleotide from glucose (lynch2015bioenergetic)', False),
    ('area_per_lipid', BIO['area_per_lipid'], 'milo2010bionumbers; BNID 106993, about 0.5 nm2 per lipid', True),
    ('c_lipid', BIO['c_lipid'], f'{PH}; two C16 acyl chains (ATP and NADPH counted as ~P) plus head group: ~100 ~P', False),
    ('k_lipid', BIO['k_lipid'], f'{PH}; phospholipid turnover 1%/h', False),
    ('genome_nt', BIO['genome_nt'], f'{PH}; an E. coli-sized genome, 4.6 Mb', False),
    ('N_A', NA, 'definitional; Avogadro constant', True),
    ('K_m_ATP_bio', par('K_m_ATP_bio'), f'{PH}; polymerase and synthetase ATP affinity, order 0.1 mM', False),
    ('ATP_ref', ATP, 'definitional; costs stated at the resting ATP', False),
    ('T', T, 'operating point: 20 C', False), ('T_ref_dmg', T, 'reference temperature of the damage rate (20 C)', False),
    ('dT_Q10', DT_Q10, 'definition of Q10', True), ('Q10_dmg', BIO['Q10_dmg'], f'{PH}; damage doubles per 10 K', False),
    ('r_dmg', BIO['r_dmg'], f'{PH}; an unrepaired cell accumulates a lethal load in 10 days at 20 C', False),
    ('r_rep_max', R_REP_MAX, f'{PH}; full repair removes damage up to 10x as fast as it accrues', False),
    ('K_D', BIO['K_D'], f'{PH}; repair is half-maximal at a tenth of a lethal load', False),
    ('C_repair', C_REPAIR, f'{PH}; repairing a whole lethal load costs 10% of building a daughter', False),
    ('repair_level', 1.0, 'definitional; full repair by default', False),
    ('growth_level', 0.0, 'definitional; no growth until the player puts energy into it', False),
    ('J_ref', der(R_TOT, 'maintenance_coefficient x cell volume ^ maintenance_exponent', 'maintenance_coefficient',
                  'maintenance_exponent'), '', False),
    ('ATP_floor', par('build_ATP_floor'), '', False),
    ('M_build_in', 0.0, 'no pool connected (an open sum port is 0)', False),
    ('M_glyc_ATP', 0.0, 'definitional; no glycogen (used when unconnected)', False),
    ('D_init', D_REST, 'calibrated: steady damage at full repair', False),
    ('P_init', 0.0, 'definitional; no daughter under way', False),
]
bio_inv = [
    {'description': 'costs and the ATP port balance (exact): J_RNA, J_ribo, J_lipid from their counts; repair; growth',
     'expr': 'np.allclose(J_RNA, n_mRNA*L_mRNA*c_nt_recycle*np.log(2)/t_half_mRNA/N_A*f_ATP_bio, rtol=1e-12, atol=0) and '
             'np.allclose(J_lipid, 2*A_mem*surface_factor/area_per_lipid*c_lipid*k_lipid/N_A*f_ATP_bio, rtol=1e-12, '
             'atol=0) and np.allclose(J_repair, R_repair*C_repair, rtol=1e-12, atol=0) and '
             'np.allclose(j_ATP_bio, -J_ATP_bio, rtol=1e-12, atol=0) and np.allclose(j_ADP_bio, J_ATP_bio, rtol=1e-12, atol=0)'},
    {'description': 'growth is asked in multiples of the resting turnover and throttled: J_growth = growth_level J_ref '
                    'throttle; the daughter costs M_build_bio + M_build_protein; progress = P_growth / M_daughter, glycogen '
                    'progress = M_glyc_ATP / M_daughter (exact)',
     'expr': 'np.allclose(throttle, np.clip((ATP - ATP_floor)/(ATP_ref - ATP_floor), 0, 1), rtol=1e-12, atol=1e-15) and '
             'np.allclose(J_growth, growth_level*J_ref*throttle, rtol=1e-12, atol=0) and np.allclose(M_daughter, '
             'M_build_bio + M_build_protein, rtol=1e-12, atol=0) and np.allclose(growth_progress, P_growth/M_daughter, '
             'rtol=1e-12, atol=1e-15) and np.allclose(glycogen_progress, M_glyc_ATP/M_daughter, rtol=1e-12, atol=1e-15)'},
    {'description': 'the daughter receives exactly the growth ATP spent: P(t) - P(0) = int J_growth_used dt (relative 1e-6)',
     'expr': 'abs(P_growth[-1] - P_growth[0] - np.trapz(J_growth_used, t)) <= 1e-6*abs(np.trapz(J_growth_used, t)) + 1e-30'},
    {'description': 'damage balance: D(t) - D(0) = int (r_dmg f_T - R_repair) dt (relative 1e-6), damage never negative',
     'expr': 'abs(D[-1] - D[0] - np.trapz(r_dmg*f_T_dmg - R_repair, t)) <= 1e-6*np.trapz(np.abs(r_dmg*f_T_dmg) + '
             'np.abs(R_repair), t) + 1e-15 and np.all(D >= 0)'},
]
write_version('cell', 'biosynthesis', 'biosynthesis_fermenter', bio_vars, bio_eqs,
              {'required_citations': ['argus2026unpublished'], 'exit_ports': [], 'entrance_ports': [],
               'general_ports': [{'port_type': 'ATP_hydrolysis_port',
                                  'variables': ['ATP', 'ADP', 'Pi', 'j_ATP_bio', 'j_ADP_bio', 'j_Pi_bio']},
                                 {'port_type': 'membrane_area', 'variables': ['A_mem']},
                                 {'port_type': 'daughter_share', 'variables': ['M_build_in'], 'multi_port': 'sum'},
                                 {'port_type': 'glycogen_store', 'variables': ['M_glyc_ATP']}],
               'roles': bio_roles},
              bio_params,
              {'sim_time': 3600.0, 'pre_time': 0.0, 'dt': 1.0, 'solver': 'CVODE_myokit',
               'solver_info': {'rtol': 1e-10, 'atol': 1e-20},
               'outputs': ['D', 'P_growth', 'J_RNA', 'J_ribo', 'J_lipid', 'R_repair', 'J_repair', 'J_growth', 'J_growth_used',
                           'M_daughter', 'growth_progress', 'glycogen_progress', 'throttle', 'x_ATP', 'M_build_protein',
                           'J_ATP_bio', 'M_build_bio', 'j_ATP_bio', 'j_ADP_bio', 'f_ATP_bio', 'f_T_dmg', 'n_lipid'],
               'run_parameters': {'growth_level': 1.0, 'repair_level': 0.5}, 'invariants': bio_inv,
               'bc_sweep': {'ranges': {'T': [273.15, 313.15], 'repair_level': [0.0, 1.0], 'growth_level': [0.0, 20.0],
                                       'ATP': [0.5, 12.0]},
                            'rationale': 'temperatures 0-40 C; repair level and growth over the ranges the game offers'},
               # nearly linear: at small steps RK4 errors are round-off, so test convergence at steps where they are not
               'timestep': {'scheme': 'rk4', 'dts': [900.0, 450.0, 225.0], 't_end': 3600.0, 'min_order': 0.8, 'tol': 1e-3,
                            'cvode_tol': 1e-2},
               'stability': STAB},
              {'reviewed': False, 'notes': (
                  'The rest of a cell\'s energy budget besides pumps and protein turnover. Upkeep: mRNA turnover (n L '
                  'c_recycle ln2/t_half), ribosome rRNA upkeep (n L c_recycle k), membrane lipid turnover (n_lipid c_lipid '
                  'k, with n_lipid = 2 A_surface / area_per_lipid, so it scales with membrane area). Damage D (1 = lethal) '
                  'accrues at r_dmg (Q10 2) and is repaired at repair_level r_rep_max D/(K_D + D), costing C_repair per '
                  'unit: a cell can spend energy on repair or leave damage to accumulate and rely on replication. Growth: '
                  'the game sets growth_level, the ATP put into a daughter in multiples of the resting turnover J_ref, '
                  'throttled as ATP falls from rest to ATP_floor, like protein synthesis; P_growth is the ATP invested so far, and the cell can divide when growth_progress = '
                  'P_growth / M_daughter reaches 1, with M_daughter = M_build_bio (RNA, lipids, DNA) + M_build_protein '
                  '(every protein_pool\'s M_build, summed on daughter_share). glycogen_progress is the share of a daughter the glycogen store (M_glyc_ATP from '
                  'the cytosol) could finish. Everything slows as ATP runs out (f_ATP = 1 at the resting ATP). At rest the '
                  'cell builds a daughter for 1.8e11 ATP, close to Lynch & Marinov\'s growth cost 2.66e10 V^0.97 '
                  '(1.6e11 for 6.28 um3); proteins are 89% of it.')},
              ['argus2026unpublished', 'milo2010bionumbers', 'bernstein2002global', 'lynch2015bioenergetic',
               'bennett2009absolute'], V_NEW)

# ------------------------------------------------------------------ 7. the generic adaptive cell: a hierarchy of supermodules
# Module type cell/generic_adaptive_cell. Each version is a flat supermodule (libcuflynx has no nested supermodules): a
# version extends another by adding submodules and links, and the builder merges them. empty_membrane holds every
# general mechanism (membrane, ion pools, cytosol, biosynthesis, protein synthesis, the proteome and the glycolytic
# enzymes, leaks); a cell type adds its proteins (a protein_pool and the protein it scales). Later bases (an animal cell
# with the Na+/K+-ATPase, ...) extend empty_membrane the same way.
from supermodule_writer import module_config, supermodule_rows  # noqa: E402

GAC = 'generic_adaptive_cell'
GAC_DIR = os.path.join(REPO, 'modules', 'cell', GAC)
CELLS = {}


def S(name, module_type, version, instance='default'):
    return {'name': name, 'module_type': module_type, 'module_subtype': version, 'instance': instance}


def define_cell(name, extends=None, submodules=(), links=(), overrides=None, description=''):
    base = CELLS[extends] if extends else {'submodules': [], 'links': [], 'overrides': {}}
    subs = base['submodules'] + list(submodules)
    names = [s['name'] for s in subs]
    assert len(set(names)) == len(names), f'{name}: duplicate submodule names'
    lk = base['links'] + list(links)
    for a, b in lk:
        assert a in names and b in names, f'{name}: link {a} -> {b}'
    CELLS[name] = {'name': name, 'extends': extends, 'submodules': subs, 'links': lk,
                   'overrides': {**base['overrides'], **(overrides or {})}, 'description': description,
                   'added': [s['name'] for s in submodules]}
    return CELLS[name]


def with_links(cell):
    out = []
    for s in cell['submodules']:
        out.append({**s, 'inp_instances': [a for a, b in cell['links'] if b == s['name']],
                    'out_instances': [b for a, b in cell['links'] if a == s['name']]})
    return out


O = lambda value, ref, sourced=False: (value, ref, sourced)
LEAKS = [S('i_leak_Na', 'i_leak_Na', 'Argus2026_v01'), S('i_leak_Ca', 'i_leak_Ca', 'ArgusUNPUBLISHED_v01'),
         S('i_leak_H', 'i_leak_H', V_PROK)]
define_cell(
    'empty_membrane',
    submodules=[S('membrane', 'membrane_potential', V_NEW), S('Na_K', 'ion_concentrations', 'SN_soma_Argus2026_v01'),
                S('Ca', 'Ca_handling', V_NEW), S('H', 'H_handling', V_NEW), S('cytosol', 'energy_metabolism', V_NEW),
                S('bio', 'biosynthesis', V_NEW), S('protein', 'protein_expression', V_NEW),
                S('pool_proteome', 'protein_pool', POOL_V, 'proteome'), S('pool_Glyc', 'protein_pool', POOL_V, 'Glyc')]
               + LEAKS,
    links=[('membrane', x) for x in ('Na_K', 'Ca', 'H', 'cytosol', 'bio', 'protein', 'i_leak_Na', 'i_leak_Ca', 'i_leak_H')]
          + [('Na_K', 'i_leak_Na'), ('Ca', 'i_leak_Ca'), ('H', 'i_leak_H'),
             ('cytosol', 'H'), ('cytosol', 'protein'), ('cytosol', 'bio'),
             ('pool_proteome', 'protein'), ('pool_proteome', 'bio'),
             ('pool_Glyc', 'protein'), ('pool_Glyc', 'bio'), ('pool_Glyc', 'cytosol')],
    overrides={
        'T': O(T, 'operating point: 20 C, the temperature Lynch & Marinov normalise maintenance costs to'),
        'Nai_init_Na_K': O(NAI, par_ref('Nai')), 'Ki_init_Na_K': O(KI, par_ref('Ki'), True),
        'Nao_Na_K': O(NAO, par_ref('Nao'), True), 'Ko_Na_K': O(KO, par_ref('Ko'), True),
        'rho_leak_Na_i_leak_Na': O(RHO_NA, 'calibrated: Na+ entry balances NhaA less Na/Ca exchange at rest (the '
                                           'fermenter cell)'),
        'gamma_leak_Na_i_leak_Na': O(GAMMA_NA, par_ref('gamma_Na')),
    },
    description=(
        'The base of every AdaptOrDeity cell: the general mechanisms with no transport proteins. A ~2 um disc of '
        'membrane (membrane_potential adaptive: capacitor, disc geometry) around a cytosol (energy_metabolism adaptive: '
        'bond-graph glycolysis, maintenance ATPase and glycogen), its Na+/K+ (ion_concentrations SN_soma_Argus2026_v01), '
        'Ca2+ (Ca_handling adaptive, fixed buffer) and H+ (H_handling adaptive, buffered, pH and E_H) pools, the passive '
        'leaks of a lipid membrane (i_leak_Na, i_leak_Ca, i_leak_H), biosynthesis (RNA, ribosomes, lipids, damage and '
        'repair, growth towards a daughter), protein synthesis (protein_expression adaptive) and two protein_pools every '
        'cell has: the proteome (ribosomal, nucleoid and other unmodelled proteins) and the glycolytic enzymes (scaling '
        'V_gly). Its default values are the fermenter cell\'s, so a cell type that adds its proteins starts calibrated; '
        'alone it has no carriers or pumps, so it runs down (its glucose is spent and its gradients leak away). Cell types '
        'extend it: fermenter (generic_adaptive_cell fermenter); later bases add the proteins a class shares (an animal '
        'cell: the Na+/K+-ATPase, NCX, ...).'))

FERM_PROTEINS = [  # (pool key, the submodule it scales, (module_type, version))
    ('H_ATPase', 'i_H_ATPase', ('i_H_ATPase', 'FoF1_ArgusUNPUBLISHED_v01')),
    ('NhaA', 'i_NaH', ('i_NaH', 'NhaA_ArgusUNPUBLISHED_v01')),
    ('KefC', 'i_KH', ('i_KH', 'KefC_ArgusUNPUBLISHED_v01')),
    ('Trk', 'i_K_Trk', ('i_leak_K', V_NEW)),
    ('NCX', 'i_NaCa', ('i_NaCa', V_NEW)),
    ('Glf', 'membrane_face', ('cell_membrane_face', V_NEW)),
    ('LacT', 'membrane_face', None),   # the face's second protein_expression port (pools connect in submodule order)
    ('Flagellum', 'i_flagellum', ('i_H_flagellum', 'MotAB_ArgusUNPUBLISHED_v01')),
    # glucose carriers the game unlocks after level 1 (pools start at 0)
    ('PTS', 'glc_PTS', (GC, V_GC_ATP, 'PTS')),
    ('ABC', 'glc_ABC', (GC, V_GC_ATP, 'ABC')),
    ('GalP', 'glc_GalP', (GC, V_GC_H)),
]
GLC_CARRIERS = ['glc_PTS', 'glc_ABC', 'glc_GalP']
ATP_CARRIERS = {'glc_PTS': 'n_P_glc_PTS', 'glc_ABC': 'n_P_glc_ABC'}   # submodule: its n_P parameter
_targets = []
for key, target, mod in FERM_PROTEINS:
    if mod and target not in [t['name'] for t in _targets]:
        _targets.append(S(target, *mod))
define_cell(
    'fermenter', extends='empty_membrane',
    submodules=_targets + [S(f'pool_{key}', 'protein_pool', POOL_V, key) for key, _, _ in FERM_PROTEINS],
    links=[('membrane', t['name']) for t in _targets]
          + [('Na_K', x) for x in ('i_K_Trk', 'i_NaH', 'i_KH', 'i_NaCa')] + [('Ca', 'i_NaCa')]
          + [('H', x) for x in ('i_H_ATPase', 'i_NaH', 'i_KH', 'i_flagellum')]
          + [('cytosol', 'i_H_ATPase'), ('membrane_face', 'cytosol'), ('membrane_face', 'H')]
          + [(c, 'cytosol') for c in GLC_CARRIERS] + [('H', 'glc_GalP')]
          + [l for key, target, _ in FERM_PROTEINS
             for l in ((f'pool_{key}', 'protein'), (f'pool_{key}', 'bio'), (f'pool_{key}', target))],
    overrides={
        'rho_leak_K_i_K_Trk': O(RHO_K, 'calibrated: K+ entry through Trk balances KefC export at rest'),
        'gamma_leak_K_i_K_Trk': O(GAMMA_K, par_ref('gamma_K')),
        'Nao_i_NaCa': O(NAO, par_ref('Nao'), True), 'Cao_i_NaCa': O(CAO, par_ref('Cao'), True),
        'rho_NaCa_i_NaCa': O(RHO_NCX, 'calibrated: exchanger cycles = NCX_share of NhaA\'s at rest'),
        **{f'{k}_i_NaCa': O(NCX[k], par_ref(k)) for k in ('KmCa', 'KmNai', 'Ksat', 'alpha_NaCa', 'gamma_NaCa', 'r_NaCa')},
    },
    description=(
        'AdaptOrDeity\'s starting organism: generic_adaptive_cell empty_membrane plus the proteins of a ~2 um prokaryote '
        'in a Proterozoic sea (pH 7.3) that makes its ATP by homolactic fermentation and spends a quarter of it on an '
        'F1Fo H+-ATPase, which pumps protons out and so powers the membrane potential and the other ion gradients. Each '
        'protein is a protein_pool (its expression, costs, daughter share and membrane area) and the module it scales: '
        'i_flagellum (i_H_flagellum MotAB: the proton-driven flagellar motor; one flagellum at start), '
        'i_H_ATPase (FoF1, 10/3 H+ per ATP), i_NaH (NhaA, 2 H+ in per Na+ out), i_KH (KefC, K+/H+), i_K_Trk (i_leak_K '
        'adaptive as the Trk K+ uptake channel), i_NaCa (adaptive; archaeal NCX_Mj is the prokaryotic equivalent) and '
        'membrane_face (cell_membrane_face adaptive: the Glf glucose facilitator and a lactate-/H+ symporter). The '
        'default instance is an exact steady state at -105 mV and pH 7.6 with total ATP turnover equal to Lynch & '
        'Marinov\'s maintenance cost at 20 C. Three glucose carriers the game unlocks after level 1 sit beside Glf with '
        'their pools at 0 (m_init and m_target 0), so the level-1 cell is unchanged until the game sets a target: '
        'glc_PTS (glucose_carrier ATP, instance PTS: the phosphotransferase system, Km 5 uM, one PEP charged as one '
        'ATP per glucose), glc_ABC (instance ABC: binding-protein importer, Km 0.5 uM, 2 ATP per glucose) and glc_GalP '
        '(glucose_carrier H_symport: Km 0.1 mM, one H+ per glucose). Each takes its medium glucose from its own '
        'C_glc_e parameter, set to the same value as membrane_face\'s.'))

m = lambda sub_name, var: f'mod_{sub_name}__{var}'
D = lambda x: f'({x}[-1] - {x}[0])'


def cell_invariants(names, steady):
    CA_TOT = f'({m("Ca", "Cai")} + {m("Ca", "CaB_Cai")})'
    VOL0 = f'{m("membrane", "Vol")}[0]'
    # Protons that cross no membrane as charge: those released by ATP hydrolysis and those leaving with lactate-. By the
    # carbon and adenylate balances their total is exactly Vol (dC_lac - dC_ATP), so no flux integral is needed.
    CHARGE_IN = (f'96485*{VOL0}*({D(m("Na_K", "Nai"))} + {D(m("Na_K", "Ki"))} + 2*{D(CA_TOT)} + '
                 f'{D(m("H", "H_total_rel"))} - {D(m("cytosol", "C_lac"))} + {D(m("cytosol", "C_ATP"))})')
    inv = [
        {'description': 'charge conservation across modules: the charge on the membrane changes by exactly the cation '
                        'charge that crossed it, Cm dV = F Vol d(Na_i + K_i + 2 Ca_total + H_total - lactate + ATP): '
                        'protons from ATP hydrolysis and those leaving with lactate carry no membrane charge (relative '
                        '1e-6)',
         'expr': f'abs({m("membrane", "Cm")}[0]*1e-12*{D(m("membrane", "V"))}*1e-3 - {CHARGE_IN}) <= '
                 f'1e-6*96485*{VOL0}*(abs({D(m("Na_K", "Nai"))}) + abs({D(m("Na_K", "Ki"))}) + '
                 f'abs({D(m("H", "H_total_rel"))})) + 1e-21'},
        {'description': 'adenine and phosphoryl groups are conserved through every ATP user: ATP + ADP and ATP + Pi '
                        'constant (relative 1e-9)',
         'expr': f'np.ptp({m("cytosol", "C_ATP")} + {m("cytosol", "C_ADP")}) <= 1e-9*np.max({m("cytosol", "C_ATP")} + '
                 f'{m("cytosol", "C_ADP")}) and np.ptp({m("cytosol", "C_ATP")} + {m("cytosol", "C_Pi")}) <= '
                 f'1e-9*np.max({m("cytosol", "C_ATP")} + {m("cytosol", "C_Pi")})'},
        {'description': 'the ATP budget closes: the cytosol\'s external ATP use equals what the H+-ATPase (if present), '
                        'protein synthesis and biosynthesis take (relative 1e-9)',
         'expr': f'np.allclose({m("cytosol", "J_ATP_turnover")} - {m("cytosol", "J_ATP_maint")}, ' + ' + '.join(
             [f'{m("i_H_ATPase", "v_ATP_H")}'] * ('i_H_ATPase' in names) +
             [f'{m("protein", "J_ATP_total")}', f'{m("bio", "J_ATP_bio")}'] +
             [f'{m(c, "J_ATP_cost")}' for c in ATP_CARRIERS if c in names]) + ', rtol=1e-9, atol=1e-30)'},
        {'description': 'second law for the leaks: each runs down its electrochemical gradient, i (V - E) >= 0',
         'expr': ' and '.join(
             [f'np.all({m("i_leak_Na", "i_leak_Na")}*({m("membrane", "V")} - {m("Na_K", "E_Na")}) >= 0)',
              f'np.all({m("i_leak_Ca", "i_leak_Ca")}*({m("membrane", "V")} - {m("i_leak_Ca", "E_Ca")}) >= 0)',
              f'np.all({m("i_leak_H", "i_leak_H")}*({m("membrane", "V")} - {m("H", "E_H")}) >= 0)'] +
             [f'np.all({m("i_K_Trk", "i_leak_K")}*({m("membrane", "V")} - {m("Na_K", "E_K")}) >= 0)'] * ('i_K_Trk' in names))},
        {'description': 'concentrations stay non-negative: ATP, ADP, Pi, glucose, lactate and glycogen >= -1e-9 mM',
         'expr': ' and '.join(f'np.all({m("cytosol", "C_" + s_)} >= -1e-9)' for s_ in cyto_species + ['glyc'])},
        {'description': 'second law for the reactions: dissipation P_diss >= 0',
         'expr': f'np.all({m("cytosol", "P_diss")} >= -1e-30)'},
        {'description': 'protein synthesis never spends more than its capacity: J_build_total <= J_cap (relative 1e-9)',
         'expr': f'np.all({m("protein", "J_build_total")} <= {m("protein", "J_cap")}*(1 + 1e-9) + 1e-30)'},
    ]
    if 'i_H_ATPase' in names:
        inv += [
            {'description': 'second law for the H+-ATPase: it hydrolyses ATP only while the ATP affinity exceeds the '
                            'proton work, v (A_hyd - n_H F (E_H - V)) >= 0',
             'expr': f'np.all({m("i_H_ATPase", "v_ATP_H")}*({m("cytosol", "A_hyd")} - (10/3)*96485*({m("H", "E_H")} - '
                     f'{m("membrane", "V")})*1e-3) >= -1e-30)'},
            {'description': 'second law for the flagellar motor: H+ run down the proton motive force, j pmf >= 0',
             'expr': f'np.all({m("i_flagellum", "j_H_flag")}*{m("i_flagellum", "pmf")} >= 0)'},
            {'description': 'second law for the exchangers: NhaA and KefC run down their free energy, r drive >= 0',
             'expr': f'np.all({m("i_NaH", "r_cycle")}*-({m("membrane", "V")} + {m("Na_K", "E_Na")} - 2*{m("H", "E_H")}) '
                     f'>= 0) and np.all({m("i_KH", "r_cycle")}*({m("H", "E_H")} - {m("Na_K", "E_K")}) >= 0)'}]
    carriers = [c for c in GLC_CARRIERS if c in names]
    if carriers:
        glc_in = ' + '.join([m('membrane_face', 'J_glc')] + [m(c, 'J_glc') for c in carriers])
        inv += [
            {'description': 'carbon is conserved through every glucose carrier: Vol d(6 glc + 3 lac + 6 glycogen) = '
                            'int (6 sum J_glc + 3 J_lac) dt (relative 1e-5: trapezoid integration of the logged fluxes)',
             'expr': f'abs({VOL0}*(6*{D(m("cytosol", "C_glc"))} + 3*{D(m("cytosol", "C_lac"))} + '
                     f'6*{D(m("cytosol", "C_glyc"))}) - np.trapz(6*({glc_in}) + 3*{m("membrane_face", "J_lac")}, t)) '
                     f'<= 1e-5*np.trapz(np.abs(6*({glc_in})) + np.abs(3*{m("membrane_face", "J_lac")}), t) + 1e-30'},
            {'description': 'the carriers\' costs follow their stoichiometry (exact): n_P ATP per glucose for the '
                            'ATP-driven ones, one H+ per glucose for the symporter',
             'expr': ' and '.join([f'np.allclose({m(c, "J_ATP_cost")}, {n}*{m(c, "J_glc")}, rtol=1e-12, atol=0)'
                                   for c, n in ATP_CARRIERS.items() if c in names] +
                                  [f'np.allclose({m("glc_GalP", "j_H")}, n_H_glc_GalP*{m("glc_GalP", "J_glc")}, '
                                   'rtol=1e-12, atol=0)'] * ('glc_GalP' in names))},
            {'description': 'second law for the glucose carriers, each against the cytosol\'s own state: J (RT '
                            'ln(C_e/C_glc) + n_P A_hyd) >= 0 (ATP-driven), J (RT ln(C_e/C_glc) + n_H F (E_H - V)) >= 0 '
                            '(symporter)',
             'expr': ' and '.join(
                 [f'np.all({m(c, "J_glc")}*(R*T*np.log(C_glc_e_{c}/np.abs({m("cytosol", "C_glc")})) + '
                  f'{n}*{m("cytosol", "A_hyd")}) >= -1e-30)' for c, n in ATP_CARRIERS.items() if c in names] +
                 [f'np.all({m("glc_GalP", "J_glc")}*(R*T*np.log(C_glc_e_glc_GalP/np.abs({m("cytosol", "C_glc")})) + '
                  f'n_H_glc_GalP*F*({m("H", "E_H")} - {m("membrane", "V")})*1e-3) >= -1e-30)'] * ('glc_GalP' in names))},
        ]
    if steady:
        inv.append({'description': 'the default point is a steady state: V stays at -105 mV (1e-6 mV), pH at 7.6 (1e-9) '
                                   'and ATP at 9.6 mM (relative 1e-9)',
                    'expr': f'np.max(np.abs({m("membrane", "V")} + 105)) <= 1e-6 and np.max(np.abs({m("H", "pH_i")} - '
                            f'7.6)) <= 1e-9 and np.max(np.abs({m("cytosol", "C_ATP")} - 9.6)) <= 1e-8', 'applies': 'run'})
    return inv


def cell_outputs(names):
    base = [('membrane', 'V'), ('membrane', 'Cm'), ('membrane', 'Vol'), ('Na_K', 'Nai'), ('Na_K', 'Ki'),
            ('Na_K', 'E_Na'), ('Na_K', 'E_K'), ('Ca', 'Cai'), ('Ca', 'CaB_Cai'), ('H', 'H_i'), ('H', 'pH_i'),
            ('H', 'E_H'), ('H', 'H_total_rel'), ('cytosol', 'C_glc'), ('cytosol', 'C_lac'), ('cytosol', 'C_ATP'),
            ('cytosol', 'C_ADP'), ('cytosol', 'C_Pi'), ('cytosol', 'C_glyc'), ('cytosol', 'v_gly'), ('cytosol', 'v_hyd'),
            ('cytosol', 'v_store'), ('cytosol', 'v_mob'), ('cytosol', 'A_hyd'), ('cytosol', 'P_diss'),
            ('cytosol', 'Q_met'), ('cytosol', 'j_H_met'), ('cytosol', 'J_ATP_turnover'), ('cytosol', 'J_ATP_gly'),
            ('cytosol', 'J_ATP_maint'), ('i_leak_Na', 'i_leak_Na'), ('i_leak_Ca', 'i_leak_Ca'), ('i_leak_Ca', 'E_Ca'),
            ('i_leak_H', 'i_leak_H'), ('protein', 'J_ATP_total'), ('protein', 'J_keep_total'),
            ('protein', 'J_build_total'), ('protein', 'J_cap'), ('protein', 'M_build_protein'), ('protein', 'f_area'),
            ('bio', 'D'), ('bio', 'P_growth'), ('bio', 'J_ATP_bio'), ('bio', 'M_build_bio'), ('bio', 'M_daughter'),
            ('bio', 'growth_progress')]
    if 'i_H_ATPase' in names:
        base += [('i_K_Trk', 'i_leak_K'), ('i_H_ATPase', 'i_H_ATPase'), ('i_H_ATPase', 'v_ATP_H'),
                 ('i_H_ATPase', 'f_thermo'), ('i_NaH', 'r_cycle'), ('i_NaH', 'i_NaH'), ('i_KH', 'r_cycle'),
                 ('i_NaCa', 'i_NaCa'), ('membrane_face', 'J_glc'), ('membrane_face', 'J_lac'),
                 ('membrane_face', 'J_H_lac'), ('i_flagellum', 'j_H_flag'), ('i_flagellum', 'pmf'),
                 ('i_flagellum', 'v_swim')]
    for c in GLC_CARRIERS:
        if c in names:
            base += [(c, 'J_glc')] + ([(c, 'J_ATP_cost')] if c in ATP_CARRIERS else [(c, 'j_H')])
    base += [(n, 'm') for n in names if n.startswith('pool_')]
    return [m(s, v).replace('__', '/') for s, v in base]


def write_cell(cell, steady, sweep_extra, timestep=None):
    name = cell['name']
    d = os.path.join(GAC_DIR, 'versions', name)
    os.makedirs(os.path.join(d, 'instances', 'default'), exist_ok=True)
    subs = with_links(cell)
    names = [s['name'] for s in subs]
    stem = f'{GAC}_{name}'
    config = [{'module_type': GAC, 'module_subtype': name, 'module_format': 'supermodule',
               'default_instance': 'default', 'licence': 'CC0-1.0', 'creator': [],
               'required_citations': ['argus2026unpublished'], 'description': cell['description'],
               'submodules': subs}]
    json.dump(config, open(os.path.join(d, f'{stem}_modules_config.json'), 'w'), indent=2)
    rows = supermodule_rows(subs)
    have = {r[0] for r in rows}
    for k in cell['overrides']:
        assert k in have, f'{name}: override {k} is not a parameter of the supermodule'
    out = []
    for n, units, value, ref, sourced in rows:
        if n in cell['overrides']:
            value, ref, src_ok = cell['overrides'][n]
            sourced = 'yes' if src_ok else 'no'
        out.append([n, units, value, ref, sourced])
    with open(os.path.join(d, 'instances', 'default', 'default_parameters.csv'), 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['variable_name', 'units', 'value', 'data_reference', 'sourced'])
        for r in out:
            w.writerow([r[0], r[1], repr(float(r[2])), r[3], r[4]])
    vals = {r[0]: float(r[2]) for r in out}
    ranges = {'T': [273.15, 313.15], 'T_ref_cytosol': [283.15, 303.15],
              'growth_level_bio': [0.0, 20.0], 'repair_level_bio': [0.0, 1.0], 'store_alloc_cytosol': [0.0, 1.0],
              **{f'm_target_{n}': [0.25, 4.0] for n in names if n.startswith('pool_') and n != 'pool_proteome'},
              **sweep_extra}
    verif = {
        'module_type': GAC, 'version': name, 'sim_time': 60.0, 'pre_time': 0.0, 'dt': 0.01,
        'solver': 'CVODE_myokit', 'solver_info': {'rtol': 1e-10, 'atol': 1e-16},
        'outputs': cell_outputs(names), 'run_parameters': {}, 'invariants': cell_invariants(names, steady),
        'bc_sweep': {'ranges': ranges, 'extra_parameters': [k for k in ranges if k in vals and k != 'T'],
                     'exclude': ['R', 'F'],
                     'rationale': TEMPERATURE_SWEEP['rationale'] + '; R and F are physical constants, not swept; '
                                  'protein expression targets x0.25-4, growth 0-20x the resting turnover, repair and '
                                  'storage 0-1: the ranges AdaptOrDeity lets a player choose'},
        'timestep': timestep or {'scheme': 'rk4', 'dts': [2e-7, 1e-7, 5e-8], 't_end': 1e-4, 'min_order': 0.8,
                                 'tol': 1e-3, 'cvode_tol': 1e-2},
        # the Ca buffer (~1 us) makes the whole cell stiff: a short horizon and budget keep the solver matrix tractable
        'stability': {'supported': [{'solver': 'CVODE_myokit', 'rtol': 1e-8, 'atol': 1e-14}], 't_end': 1.0,
                      'time_budget': 10},
        'validation': {'default': {}},
        'supermodule': {'globals': ['R', 'F', 'T']},
    }
    json.dump(verif, open(os.path.join(d, f'{stem}_verification_config.json'), 'w'), indent=2)
    extends = f' It extends generic_adaptive_cell {cell["extends"]} with: {", ".join(cell["added"])}.' if cell['extends'] else ''
    yaml.safe_dump({'module_type': GAC, 'version': name, 'reviewed': False, 'description': cell['description'],
                    'extends': cell['extends'],
                    'notes': ('Assembled by tools/build_adaptive_cell.py from modules/cell/generic_adaptive_cell/'
                              'parameters.json: its parameter rows are its submodules\' default instances (every value '
                              'with its source there) with the overrides the cell sets.' + extends +
                              ' Protein expression: each protein_pool has a target m_target set by the game.')},
                   open(os.path.join(d, f'{stem}_tests.yaml'), 'w'), sort_keys=False, width=120)
    open(os.path.join(d, f'{stem}_references.bib'), 'w').write(bib_text(['argus2026unpublished']))
    return d


from module_writer import bib_text  # noqa: E402
# empty_membrane runs down from the fermenter's rest: at the fermenter's test steps its RK4 errors are already round-off,
# so convergence is tested at larger steps (still inside RK4's stability limit for the ~1 us Ca buffer, ~3 us)
CELL_DIRS = [write_cell(CELLS['empty_membrane'], False, {},
                        {'scheme': 'rk4', 'dts': [2e-6, 1e-6, 5e-7], 't_end': 1e-3, 'min_order': 0.8, 'tol': 1e-3,
                         'cvode_tol': 1e-2}),
             write_cell(CELLS['fermenter'], True, {'C_glc_e_membrane_face': [0.0, 50.0], 'a_motor_i_flagellum': [0.0, 1.0],
                                                   'm_target_pool_Flagellum': [1.0, 8.0]})]
_ferm_rows = {r[0]: r for r in csv.reader(open(os.path.join(CELL_DIRS[1], 'instances', 'default',
                                                            'default_parameters.csv')))}


# ------------------------------------------------------------------ every cited key in each version's references.bib
def existing_entry(key):
    """The entry text for key from any references.bib in the library (verbatim), or None."""
    import glob
    import re
    for path in glob.glob(os.path.join(REPO, 'modules', '**', '*_references.bib'), recursive=True):
        text = open(path).read()
        m_ = re.search(r'@\w+\{' + re.escape(key) + r',.*?\n\}', text, re.S)
        if m_:
            return m_.group(0)
    return None


def complete_bib(version_dir):
    import glob
    from cam_testing.bib import reference_key
    bib_path = glob.glob(os.path.join(version_dir, '*_references.bib'))[0]
    text = open(bib_path).read()
    keys = set()
    for csv_path in glob.glob(os.path.join(version_dir, 'instances', '*', '*_parameters.csv')):
        for row in csv.DictReader(open(csv_path)):
            k = reference_key(row['data_reference'])
            if k and k not in (PH, 'definitional'):
                keys.add(k)
    added = []
    for k in sorted(keys):
        if '{' + k + ',' in text:
            continue
        entry = bib_text([k]).strip() if k in BIB else existing_entry(k)
        assert entry, f'{version_dir}: no bib entry for {k}'
        text = text.rstrip('\n') + '\n\n' + entry + '\n'
        added.append(k)
    open(bib_path, 'w').write(text)
    return added


for d in [os.path.join(REPO, 'modules', *parts) for parts in [
        ('cell', 'energy_metabolism', 'versions', V_NEW), ('cell', 'membrane_potential', 'versions', V_NEW),
        ('cell', 'Ca_handling', 'versions', V_NEW), ('cell', 'ion_channels', 'i_leak_Ca', 'versions', 'ArgusUNPUBLISHED_v01'),
        ('cell', 'ion_channels', 'i_NaK', 'versions', V_NAK_VERSION), ('transport', 'cell_membrane_face', 'versions', V_NEW),
        ('cell', 'H_handling', 'versions', V_NEW), ('cell', 'ion_channels', 'i_leak_H', 'versions', V_PROK),
        ('cell', 'ion_channels', 'i_H_ATPase', 'versions', 'FoF1_ArgusUNPUBLISHED_v01'),
        ('cell', 'ion_channels', 'i_NaH', 'versions', 'NhaA_ArgusUNPUBLISHED_v01'),
        ('cell', 'ion_channels', 'i_KH', 'versions', 'KefC_ArgusUNPUBLISHED_v01'),
        ('cell', 'biosynthesis', 'versions', V_NEW), ('cell', 'protein_expression', 'versions', V_NEW),
        ('cell', 'protein_pool', 'versions', POOL_V), ('cell', 'ion_channels', 'i_leak_K', 'versions', V_NEW),
        ('cell', 'ion_channels', 'i_NaCa', 'versions', V_NEW),
        ('cell', 'ion_channels', 'i_H_flagellum', 'versions', 'MotAB_ArgusUNPUBLISHED_v01'),
        ('transport', GC, 'versions', V_GC_ATP), ('transport', GC, 'versions', V_GC_H)]] + CELL_DIRS:
    added = complete_bib(d)
    if added:
        print(f'{os.path.relpath(d, REPO)}: added {added}')
print('written')
if UNSOURCED:
    print(f'{len(UNSOURCED)} parameter values are not from parameters.json and not derived:')
    for u in UNSOURCED:
        print('  ' + u)
