"""
Import circulatory_autogen's system models (resources/<model>_vessel_array.csv and
_parameters.csv, plus obs_data / params_for_id when present) into

    modules/system/<category>/<model>/
        <model>_vessel_array.json, <model>_parameters.csv  the model, using the module library;
                                                           a heart is split into clock/chambers/valves
        reference/<model>_vessel_array.json, _parameters.csv  circulatory_autogen's original (vessel array as JSON),
                                                           for the equivalence test
        <model>_system.yaml                                test spec (never overwritten)

    python tools/import_systems.py --ca-dir ../circulatory_autogen [--only 3compartment ...]

Heart splitting (a vessel whose vessel_type starts with 'heart' and whose module has the
v_ivc/v_svc/v_pvn vessel ports):
  heart_clock (cardiac_clock nn, or nn_controlled for the controlled hearts)
  ra, rv, la, lv (chamber vv) and trv, puv, miv, aov (valve pp / pp_linear / pp_rmod)
  ra_inflow (flow_merge) when the right atrium has more than one inflow
  la_outflow (flow_split) + asd_shunt (resistor) for heart_ASD
It reproduces what circulatory_autogen does for a heart implicitly: inputs are ivc, svc,
pulmonary in that order (with two inputs: svc, pulmonary and a zero IVC), outputs are
aorta then pulmonary artery, and a heart with fewer than two outputs gets par/pvn vessels
(ModelParsers). Parameter names follow the new modules ({var}_{vessel}); valve names were
chosen so that valve parameters keep their original names (K_vo_aov, A_nn_trv, ...).
"""
import argparse
import csv
import glob
import json
import os
import shutil
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing import vessel_array  # noqa: E402
from cam_testing.library import read_config  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYSTEM_DIR = os.path.join(REPO, 'modules', 'system')

CATEGORIES = {
    'closed_loop_cvs': ['3compartment', '3compartment_extra_ops', '3compartment_nonstiff', 'simple_physiological',
                        'physiological', 'FinalModel', 'FTU_wCVS', 'cvs_model_0d', 'cvs_model_with_arm_0d',
                        'neonatal', 'new_valve_p_est', 'generic_junction_test_closed_loop',
                        'generic_junction_test2_closed_loop'],
    'controlled_cvs': ['control_phys', 'control_phys_asd', 'control_parasymp', 'elic', 'parasympathetic_model'],
    'open_loop_arterial': ['cerebral', 'cerebral_elic', 'behdad_test1', 'generic_junction_test_open_loop',
                           'generic_junction_test2_open_loop', 'test_fft'],
    'cellular': ['SN_simple', 'NKE_pump'],
    'benchmarks': ['Lotka_Volterra', 'VanDerPol', 'FitzHugh_Nagumo', 'TwoMinima', 'Simple_ODE_Benchmark',
                   'delay_test', 'pid_control', 'test_init_states', 'ports_test'],
}
# C++ / 1D models: not runnable by the CellML harness (to be brought in with a cpp test path)
EXCLUDED = ['aortic_bif_1d', 'aortic_bif_hybrid_V1', 'aortic_bif_hybrid_V2', 'cvs_model_with_arm_hybrid']

CHAMBERS = {'ra': 'a', 'la': 'a', 'rv': 'v', 'lv': 'v'}     # atrial / ventricular activation
VALVES = {'trv': ('ra', 'rv'), 'puv': ('rv', None), 'miv': ('la', 'lv'), 'aov': ('lv', None)}
VALVE_BC = {'heart_nonstiff': 'pp_linear', 'heart_new_valve': 'pp_rmod'}
CONTROLLED_BC = ('vp_wCont',)
HEART_OUTPUT_MAP = {   # heart variable -> (new vessel, variable), for the equivalence test
    'u_ra': ('ra', 'u'), 'u_rv': ('rv', 'u'), 'u_la': ('la', 'u'), 'u_lv': ('lv', 'u'),
    'q_ra': ('ra', 'q'), 'q_rv': ('rv', 'q'), 'q_la': ('la', 'q'), 'q_lv': ('lv', 'q'),
    'v_trv': ('trv', 'v'), 'v_puv': ('puv', 'v'), 'v_miv': ('miv', 'v'), 'v_aov': ('aov', 'v'),
    'zeta_trv': ('trv', 'zeta'), 'zeta_puv': ('puv', 'zeta'), 'zeta_miv': ('miv', 'zeta'), 'zeta_aov': ('aov', 'zeta'),
    'zeta_trv_pre': ('trv', 'zeta_pre'), 'zeta_puv_pre': ('puv', 'zeta_pre'),
    'zeta_miv_pre': ('miv', 'zeta_pre'), 'zeta_aov_pre': ('aov', 'zeta_pre'),
    'e_a': ('ra', 'e'), 'e_v': ('lv', 'e'), 's': ('heart_clock', 's'), 'mt': ('heart_clock', 'mt'),
    'chi_a': ('ra', 'chi'), 'chi_v': ('lv', 'chi'), 'v_asd': ('asd_shunt', 'v'),
    'T_wCont': ('heart_clock', 'T_eff'),
    'chi_afloor': ('ra', 'chi_floor'), 'chi_vfloor': ('lv', 'chi_floor'),
    'chi_afloor_final': ('ra', 'chi_floor_final'), 'chi_vfloor_final': ('lv', 'chi_floor_final'),
    'R_lin_trv': ('trv', 'R_lin'), 'R_lin_puv': ('puv', 'R_lin'), 'R_lin_miv': ('miv', 'R_lin'), 'R_lin_aov': ('aov', 'R_lin'),
    'vel_trv': ('trv', 'vel'), 'vel_puv': ('puv', 'vel'), 'vel_miv': ('miv', 'vel'), 'vel_aov': ('aov', 'vel'),
    'E_lv_A_wCont': ('lv', 'E_A_eff'), 'E_rv_A_wCont': ('rv', 'E_A_eff'),
    'T_ac_wCont': ('ra', 'T_c_eff'), 'T_ar_wCont': ('ra', 'T_r_eff'),
    'T_vc_wCont': ('lv', 'T_c_eff'), 'T_vr_wCont': ('lv', 'T_r_eff'),
}
for _v in ('trv', 'puv', 'miv', 'aov'):
    for _k in ('A_eff', 'B', 'L'):
        HEART_OUTPUT_MAP[f'{_k}_{_v}'] = (_v, _k)
# not single variables in the split heart, or not model outputs
EQUIVALENCE_IGNORE = {'q_heart': 'sum of the four chambers (no single variable after the split)',
                      'q': 'sum of the four chambers (no single variable after the split)'}
# the valves' opening state is stiff near closure: at rtol 1e-10 the split and original hearts differ by up to
# ~5e-5 in valve internals (solver error, converging to < 3e-7 at rtol 1e-12)
HEART_SOLVER_INFO = {'rtol': 1e-12, 'atol': 1e-14}
EQUIVALENCE_WRAPPED = {'mt': 1.0, 'chi_afloor': 1.0, 'chi_vfloor': 1.0,   # sawtooth outputs: compare modulo the period
                       'chi_afloor_final': 1.0, 'chi_vfloor_final': 1.0}
LVPROP = {'E_A_rv': ('E_lv_A', 0.8), 'E_A_ra': ('E_lv_A', 0.02), 'E_A_la': ('E_lv_A', 0.02),
          'E_B_la': ('E_lv_B', 0.9), 'E_B_ra': ('E_lv_B', 0.9), 'E_B_rv': ('E_lv_B', 0.6)}


def load_configs():
    configs = {}
    for path in glob.glob(os.path.join(REPO, 'modules', '**', '*_modules_config.json'), recursive=True):
        for e in read_config(path):   # either config format, as libcuflynx's names
            configs[(e['vessel_type'], e['BC_type'])] = e
    return configs


def read_rows(path):
    rows = []
    with open(path) as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        for r in reader:
            if not r or not r[0].strip():
                continue
            rows.append({h: (r[i].strip() if i < len(r) else '') for i, h in enumerate(header)})
    return header, rows


def split_list(s):
    return s.split() if s else []


def has_port(configs, row, side):
    '''Does this vessel have a vessel_port on the given side ('exit' or 'entrance')?'''
    cfg = configs.get((row['vessel_type'], row['BC_type']))
    return bool(cfg) and any(p['port_type'] == 'vessel_port' for p in cfg[f'{side}_ports'])


def port_types(configs, row, side):
    cfg = configs.get((row['vessel_type'], row['BC_type'])) or {}
    return [p['port_type'] for p in cfg.get(f'{side}_ports', []) + cfg.get('general_ports', [])]


def split_heart(rows, params, configs, notes):
    by = {r['name']: r for r in rows}
    hearts = [r for r in rows if r['vessel_type'].startswith('heart')
              and 'v_ivc' in [v[0] for v in (configs.get((r['vessel_type'], r['BC_type'])) or {}).get('variables_and_units', [])]]
    if not hearts:
        return rows, params, {}
    if len(hearts) > 1:
        raise ValueError('more than one heart')
    H = hearts[0]
    hname = H['name']
    inp, out = split_list(H['inp_vessels']), split_list(H['out_vessels'])

    # circulatory_autogen's ModelParsers adds par/pvn when the heart has fewer than two outputs
    if len(out) < 2:
        rows.append({'name': 'par', 'BC_type': 'vp', 'vessel_type': 'arterial_simple', 'inp_vessels': hname, 'out_vessels': 'pvn'})
        rows.append({'name': 'pvn', 'BC_type': 'vp', 'vessel_type': 'arterial_simple', 'inp_vessels': 'par', 'out_vessels': hname})
        out.append('par')
        inp.append('pvn')
        notes.append('the original heart had fewer than two outputs, so circulatory_autogen added par/pvn implicitly; they are explicit here')
        by = {r['name']: r for r in rows}

    vin = [v for v in inp if has_port(configs, by[v], 'exit')]
    eff = [v for v in inp if v not in vin]
    vout = [v for v in out if has_port(configs, by[v], 'entrance')]
    other_out = [v for v in out if v not in vout]
    if len(vin) == 3:
        systemic, pulmonary = vin[:2], vin[2]
    elif len(vin) == 2:
        systemic, pulmonary = vin[:1], vin[1]          # original: zero IVC
    else:
        raise ValueError(f'heart has {len(vin)} vessel inputs')
    aortic = vout[0]
    pulm_art = vout[1] if len(vout) > 1 else None

    asd = H['vessel_type'] == 'heart_ASD'
    controlled = H['BC_type'] in CONTROLLED_BC or H['vessel_type'] in ('heart_nonstiff', 'heart_ASD')
    valve_bc = VALVE_BC.get(H['vessel_type'], 'pp')

    # effectors: by the port type they export
    eff_for = {'rv': [], 'lv': [], 'heart_clock': []}
    for e in eff:
        types = port_types(configs, by[e], 'exit')
        if 'Elastance_delta_port_rv' in types:
            eff_for['rv'].append(e)
        elif 'Elastance_delta_port_lv' in types:
            eff_for['lv'].append(e)
        elif 'period_delta_port' in types:
            eff_for['heart_clock'].append(e)
        else:
            notes.append(f'heart input {e} (ports {types}) not mapped')

    ra_in = systemic + (['asd_shunt'] if asd else [])
    merge = len(ra_in) > 1
    new = []
    new.append({'name': 'heart_clock', 'BC_type': 'nn_controlled' if controlled else 'nn', 'vessel_type': 'cardiac_clock',
                'inp_vessels': ' '.join(eff_for['heart_clock']), 'out_vessels': 'ra rv la lv'})
    sums = other_out
    new.append({'name': 'ra', 'BC_type': 'vv', 'vessel_type': 'chamber',
                'inp_vessels': ' '.join((['ra_inflow'] if merge else ra_in) + ['heart_clock']),
                'out_vessels': ' '.join(['trv'] + sums)})
    new.append({'name': 'rv', 'BC_type': 'vv', 'vessel_type': 'chamber',
                'inp_vessels': ' '.join(['trv', 'heart_clock'] + eff_for['rv']), 'out_vessels': ' '.join(['puv'] + sums)})
    new.append({'name': 'la', 'BC_type': 'vv', 'vessel_type': 'chamber',
                'inp_vessels': ' '.join([pulmonary, 'heart_clock']),
                'out_vessels': ' '.join([('la_outflow' if asd else 'miv')] + sums)})
    new.append({'name': 'lv', 'BC_type': 'vv', 'vessel_type': 'chamber',
                'inp_vessels': ' '.join(['miv', 'heart_clock'] + eff_for['lv']), 'out_vessels': ' '.join(['aov'] + sums)})
    new.append({'name': 'trv', 'BC_type': valve_bc, 'vessel_type': 'valve', 'inp_vessels': 'ra', 'out_vessels': 'rv'})
    new.append({'name': 'puv', 'BC_type': valve_bc, 'vessel_type': 'valve', 'inp_vessels': 'rv', 'out_vessels': pulm_art or ''})
    new.append({'name': 'miv', 'BC_type': valve_bc, 'vessel_type': 'valve',
                'inp_vessels': 'la_outflow' if asd else 'la', 'out_vessels': 'lv'})
    new.append({'name': 'aov', 'BC_type': valve_bc, 'vessel_type': 'valve', 'inp_vessels': 'lv', 'out_vessels': aortic})
    if merge:
        new.append({'name': 'ra_inflow', 'BC_type': 'vp', 'vessel_type': 'flow_merge',
                    'inp_vessels': ' '.join(ra_in), 'out_vessels': 'ra'})
    if asd:
        new.append({'name': 'la_outflow', 'BC_type': 'pv', 'vessel_type': 'flow_split',
                    'inp_vessels': 'la', 'out_vessels': 'miv asd_shunt'})
        new.append({'name': 'asd_shunt', 'BC_type': 'pp', 'vessel_type': 'resistor',
                    'inp_vessels': 'la_outflow', 'out_vessels': 'ra_inflow' if merge else 'ra'})

    # rewire the neighbours
    def repl(lst, old, new_names):
        items = split_list(lst)
        outl = []
        for x in items:
            outl.extend(new_names if x == old else [x])
        return ' '.join(outl)
    for r in rows:
        if r is H:
            continue
        if r['name'] in systemic:
            r['out_vessels'] = repl(r['out_vessels'], hname, ['ra_inflow' if merge else 'ra'])
        elif r['name'] == pulmonary:
            r['out_vessels'] = repl(r['out_vessels'], hname, ['la'])
        if r['name'] == aortic:
            r['inp_vessels'] = repl(r['inp_vessels'], hname, ['aov'])
        if pulm_art and r['name'] == pulm_art:
            r['inp_vessels'] = repl(r['inp_vessels'], hname, ['puv'])
        if r['name'] in sums:
            r['inp_vessels'] = repl(r['inp_vessels'], hname, ['ra', 'rv', 'la', 'lv'])
        for target, effs in eff_for.items():
            if r['name'] in effs:
                r['out_vessels'] = repl(r['out_vessels'], hname, [target])
        if hname in split_list(r['inp_vessels']) + split_list(r['out_vessels']):
            notes.append(f'{r["name"]} still references {hname}')
    rows = [r for r in rows if r is not H] + new

    # ---- parameters -------------------------------------------------------------------------
    p = {row['variable_name']: row for row in params}
    val = lambda name: p[name]['value'] if name in p else None
    ref = lambda name: p[name].get('data_reference', '') if name in p else ''
    newp = []

    def add(name, units, value, reference):
        newp.append({'variable_name': name, 'units': units, 'value': value, 'data_reference': reference})

    for ch, kind in CHAMBERS.items():
        src = ch  # e.g. E_ra_A
        if H['vessel_type'] == 'heart_LVprop' and f'E_A_{ch}' in LVPROP:
            base, factor = LVPROP[f'E_A_{ch}']
            add(f'E_A_{ch}', 'J_per_m6', repr(float(val(base)) * factor), f'{ref(base)} x{factor} (heart_LVprop)')
        else:
            add(f'E_A_{ch}', 'J_per_m6', val(f'E_{src}_A'), ref(f'E_{src}_A'))
        if H['vessel_type'] == 'heart_LVprop' and f'E_B_{ch}' in LVPROP:
            base, factor = LVPROP[f'E_B_{ch}']
            add(f'E_B_{ch}', 'J_per_m6', repr(float(val(base)) * factor), f'{ref(base)} x{factor} (heart_LVprop)')
        else:
            add(f'E_B_{ch}', 'J_per_m6', val(f'E_{src}_B'), ref(f'E_{src}_B'))
        if val(f'q_{src}_us') is None and val(f'q_{src}_0') is not None:
            add(f'q_us_{ch}', 'm3', val(f'q_{src}_0'), ref(f'q_{src}_0') + ' (named q_' + src + '_0 in the original)')
        else:
            add(f'q_us_{ch}', 'm3', val(f'q_{src}_us'), ref(f'q_{src}_us'))
        add(f'q_init_{ch}', 'm3', val(f'q_{src}_init'), ref(f'q_{src}_init'))
        add(f'T_c_{ch}', 'second', val(f'T_{kind}c'), ref(f'T_{kind}c'))
        add(f'T_r_{ch}', 'second', val(f'T_{kind}r'), ref(f'T_{kind}r'))
        add(f't_start_{ch}', 'second', val(f't_{kind}start'), ref(f't_{kind}start'))
        add(f'u_ext_{ch}', 'J_per_m3', val(f'u_intraThor_{hname}') or '0', ref(f'u_intraThor_{hname}') or 'no intrathoracic pressure')
        if ch in ('ra', 'la') or not eff_for.get(ch):
            add(f'Delta_E_{ch}', 'J_per_m6', val(f'Delta_E_{ch}_{hname}') or '0', 'no elastance control')
    if not merge and not systemic:
        pass
    if controlled and not eff_for['heart_clock']:
        add('Delta_T_heart_clock', 'second', val(f'Delta_T_{hname}') or '0', 'no period control')
    if not pulm_art:
        add('u_out_puv', 'J_per_m3', val(f'u_par_{hname}'), ref(f'u_par_{hname}'))
    if asd:
        add('R_asd_shunt', 'Js_per_m6', val('R_ASD'), ref('R_ASD'))
    heart_globals = {'q_ra_us', 'q_rv_us', 'q_la_us', 'q_lv_us', 'q_ra_init', 'q_rv_init', 'q_la_init', 'q_lv_init',
                     'T_ac', 'T_ar', 't_astart', 'T_vc', 'T_vr', 't_vstart', 'E_ra_A', 'E_ra_B', 'E_rv_A', 'E_rv_B',
                     'E_la_A', 'E_la_B', 'E_lv_A', 'E_lv_B', 'R_ASD', 'q_ra_0', 'q_rv_0', 'q_la_0', 'q_lv_0'}
    kept = [row for row in params if row['variable_name'] not in heart_globals
            and not row['variable_name'].endswith('_' + hname)]
    missing = [r['variable_name'] for r in newp if r['value'] in (None, '')]
    if missing:
        notes.append(f'no original value for {missing}')
    mapping = {f'{hname}/{k}': f'{v[0]}/{v[1]}' for k, v in HEART_OUTPUT_MAP.items()}
    mapping['__ignore__'] = {f'{hname}/{k}': why for k, why in EQUIVALENCE_IGNORE.items()}
    mapping['__wrapped__'] = {f'{hname}/{k}': p for k, p in EQUIVALENCE_WRAPPED.items()}
    mapping['__reciprocal__'] = [f'{hname}/{k}_{v}' for k in ('L', 'B') for v in ('trv', 'puv', 'miv', 'aov')]
    return rows, kept + newp, mapping


# heart parameters (circulatory_autogen global constants) -> (vessels, parameter) after the split
HEART_PARAM_TARGETS = {}
for _c in ('ra', 'rv', 'la', 'lv'):
    HEART_PARAM_TARGETS.update({f'E_{_c}_A': (_c, 'E_A'), f'E_{_c}_B': (_c, 'E_B'), f'q_{_c}_us': (_c, 'q_us'),
                                f'q_{_c}_init': (_c, 'q_init'), f'q_{_c}_0': (_c, 'q_us')})
HEART_PARAM_TARGETS.update({'T_ac': ('ra la', 'T_c'), 'T_ar': ('ra la', 'T_r'), 't_astart': ('ra la', 't_start'),
                            'T_vc': ('rv lv', 'T_c'), 'T_vr': ('rv lv', 'T_r'), 't_vstart': ('rv lv', 't_start'),
                            'R_ASD': ('asd_shunt', 'R')})
for _v in ('trv', 'puv', 'miv', 'aov'):
    for _k in ('K_vo', 'K_vc', 'M_rg', 'M_st', 'A_nn', 'v_ref', 'Rmod_valve'):
        HEART_PARAM_TARGETS[f'{_k}_{_v}'] = (_v, _k)


def convert_params_for_id(src, dst, heart_name, notes):
    '''Rewrites a params_for_id CSV for the split heart (global/heart rows naming heart parameters).'''
    with open(src) as f:
        lines = [l for l in csv.reader(f)]
    header = [h.strip() for h in lines[0]]
    vi, pi = header.index('vessel_name'), header.index('param_name')
    out = [lines[0]]
    for row in lines[1:]:
        if len(row) > max(vi, pi) and row[vi].strip() in ('global', heart_name) and row[pi].strip() in HEART_PARAM_TARGETS:
            vessels, param = HEART_PARAM_TARGETS[row[pi].strip()]
            row = list(row)
            row[vi], row[pi] = vessels, param
        elif len(row) > vi and row[vi].strip() == heart_name:
            notes.append(f'params_for_id row {row} names a heart parameter with no split equivalent')
        out.append(row)
    with open(dst, 'w', newline='') as f:
        csv.writer(f).writerows(out)


def write_rows(path, header, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(h, '') for h in header])


def rename_operands(obj, mapping):
    if isinstance(obj, dict):
        return {k: rename_operands(v, mapping) for k, v in obj.items()}
    if isinstance(obj, list):
        return [rename_operands(v, mapping) for v in obj]
    if isinstance(obj, str) and obj in mapping:
        return mapping[obj]
    return obj


def import_model(ca_dir, model, category, configs):
    res = os.path.join(ca_dir, 'resources')
    dest = os.path.join(SYSTEM_DIR, category, model)
    ref_dir = os.path.join(dest, 'reference')
    os.makedirs(ref_dir, exist_ok=True)
    notes = []
    va = os.path.join(res, f'{model}_vessel_array.csv')
    pa = os.path.join(res, f'{model}_parameters.csv')
    vessel_array.write_records(os.path.join(ref_dir, f'{model}_vessel_array.json'), vessel_array.read_records(va))
    if os.path.isfile(pa):
        shutil.copyfile(pa, os.path.join(ref_dir, f'{model}_parameters.csv'))
    header, rows = read_rows(va)
    pheader, params = read_rows(pa) if os.path.isfile(pa) else (['variable_name', 'units', 'value', 'data_reference'], [])
    heart_name = next((r['name'] for r in rows if r['vessel_type'].startswith('heart')), 'heart')
    rows, params, mapping = split_heart(rows, params, configs, notes)
    ignore = mapping.pop('__ignore__', {})
    wrapped = mapping.pop('__wrapped__', {})
    reciprocal = mapping.pop('__reciprocal__', [])
    ignore.update({f'{r["name"]}/t': 'time' for r in rows})
    vessel_array.write_records(os.path.join(dest, f'{model}_vessel_array.json'), rows)
    if os.path.isfile(pa):
        write_rows(os.path.join(dest, f'{model}_parameters.csv'), ['variable_name', 'units', 'value', 'data_reference'], params)
    for extra in ('obs_data.json', 'params_for_id.csv'):
        src = os.path.join(res, f'{model}_{extra}')
        if not os.path.isfile(src):
            continue
        shutil.copyfile(src, os.path.join(ref_dir, f'{model}_{extra}'))
        if extra.endswith('.json'):
            data = rename_operands(json.load(open(src)), mapping)
            json.dump(data, open(os.path.join(dest, f'{model}_{extra}'), 'w'), indent=1)
        elif mapping:
            convert_params_for_id(src, os.path.join(dest, f'{model}_{extra}'), heart_name, notes)
        else:
            shutil.copyfile(src, os.path.join(dest, f'{model}_{extra}'))
    invariants = []
    if category in ('closed_loop_cvs', 'controlled_cvs') and any(r['vessel_type'] == 'volume_sum' and r['name'] == 'volume_sum' for r in rows):
        invariants.append({'description': 'total blood volume is conserved in the closed loop (relative drift < 1e-6)',
                           'expr': 'np.ptp(sum_blood_volume__q_volume_sum_sum) <= 1e-6 * np.max(np.abs(sum_blood_volume__q_volume_sum_sum))'})
    spec_path = os.path.join(dest, f'{model}_system.yaml')
    if not os.path.isfile(spec_path):
        spec = {'model': model, 'category': category, 'reviewed': False,
                'sim_time': 2.0, 'pre_time': 0.0, 'dt': 0.01,
                'equivalence': {'reference': 'reference', 'output_map': mapping, 'ignore': ignore,
                                'wrapped': wrapped, 'reciprocal': reciprocal, 'tol': 1e-6,
                                'solver_info': dict(HEART_SOLVER_INFO) if mapping else {'rtol': 1e-10, 'atol': 1e-12}},
                'invariants': invariants, 'notes': notes}
        yaml.safe_dump(spec, open(spec_path, 'w'), sort_keys=False, width=120, default_flow_style=None)
    return notes


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ca-dir', default=os.path.join(REPO, '..', 'circulatory_autogen'))
    parser.add_argument('--only', nargs='*')
    args = parser.parse_args(argv)
    configs = load_configs()
    for category, models in CATEGORIES.items():
        for model in models:
            if args.only and model not in args.only:
                continue
            if not os.path.isfile(os.path.join(args.ca_dir, 'resources', f'{model}_vessel_array.csv')):
                print(f'{model}: not in {args.ca_dir}/resources, skipped')
                continue
            notes = import_model(os.path.abspath(args.ca_dir), model, category, configs)
            print(f'{category}/{model}' + (': ' + '; '.join(notes) if notes else ''))


if __name__ == '__main__':
    main()
