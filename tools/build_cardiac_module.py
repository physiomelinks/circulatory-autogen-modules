"""
Build modules/cardiac/cardiac_modules.cellml from the monolithic heart module.

The heart (modules/heart) is split into reusable parts:

  cardiac_clock_type             s, mt = s - floor(s), T_eff = T            (one per heart)
  cardiac_clock_controlled_type  as above with T_eff = max(T + Delta_T, T_min) (period effector)
  chamber_type                   one time-varying-elastance chamber (LV, RV, LA or RA):
                                 its own activation timing (T_c, T_r, t_start, scaled by T_eff/T),
                                 u = (e (E_A + Delta_E) + E_B)(q - q_us) + u_ext,
                                 dq/dt = v_in - v_out
  valve_type                     Mynard valve: zeta dynamics, A_eff, B = rho/(2 A_eff^2), L, dv/dt
  valve_linear_type              as valve_type with the linearised loss R_lin = 2 B v_ref
  valve_rmod_type                as valve_type with B = Rmod rho / A_eff
  resistor_type                  v = (u_in - u_out)/R (e.g. an atrial septal defect)
  flow_merge_type / flow_split_type  algebraic nodes for several inflows / outflows

Equations are not retyped: they are copied from the heart's own MathML (heart_simple_wCont,
heart_nonstiff, heart_new_valve) with their variables renamed, so the split heart solves the
same equations as the original. Re-run after changing the heart module:

    python tools/build_cardiac_module.py

Pre-versions tool: it reads and writes the old per-module layout (modules/<name>/...), kept for
provenance. The library now uses modules/<category>/<module_type>/versions/<version>/ (see
modules/README.md); tools/restructure_modules.py moves an old-layout tree into it.

Pre-versions tool, kept for provenance: it reads and writes the old per-module layout
(modules/heart/heart_modules.cellml ...). Its output was moved into the versions layout by
tools/restructure_modules.py: cardiac/heart version Argus2026_v01 and the cardiac_clock, chamber and valve module_types. Edit those files directly now.
"""
import copy
import os
import sys

from lxml import etree

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEART = os.path.join(REPO, 'modules', 'heart', 'heart_modules.cellml')
OUT = os.path.join(REPO, 'modules', 'cardiac', 'cardiac_modules.cellml')

C = 'http://www.cellml.org/cellml/1.1#'
M = 'http://www.w3.org/1998/Math/MathML'
NS = {'c': C, 'm': M}


def _component(root, name):
    return root.xpath(f'//c:component[@name="{name}"]', namespaces=NS)[0]


def _lhs(apply):
    left = [child for child in apply if isinstance(child.tag, str)][1]
    if left.tag == f'{{{M}}}ci':
        return left.text.strip()
    if left.tag == f'{{{M}}}apply' and left[0].tag == f'{{{M}}}diff':
        return 'd/dt ' + left.xpath('m:ci', namespaces=NS)[0].text.strip()
    raise ValueError('unexpected equation form')


def equations(root, component):
    comp = _component(root, component)
    return {_lhs(a): a for m in comp.xpath('m:math', namespaces=NS) for a in m.xpath('m:apply', namespaces=NS)}


def renamed(apply, mapping):
    new = copy.deepcopy(apply)
    for ci in new.iter(f'{{{M}}}ci'):
        name = ci.text.strip()
        ci.text = mapping.get(name, name)
    return new


def parse_math(text):
    '''A MathML <apply> written by hand (for the few equations with no heart counterpart).'''
    return etree.fromstring(f'<apply xmlns="{M}" xmlns:cellml="{C}">{text}</apply>')


def build_component(name, variables, eqs):
    comp = etree.Element(f'{{{C}}}component', name=name)
    for var in variables:
        name_, units, iface = var[:3]
        attrs = {'name': name_, 'units': units}
        if iface:
            attrs['public_interface'] = iface
        if len(var) > 3 and var[3] is not None:
            attrs['initial_value'] = str(var[3])
        etree.SubElement(comp, f'{{{C}}}variable', **attrs)
    math = etree.SubElement(comp, f'{{{M}}}math', nsmap={None: M})
    for e in eqs:
        math.append(e)
    return comp


def main():
    root = etree.parse(HEART).getroot()
    wc = equations(root, 'heart_simple_wCont')
    ns = equations(root, 'heart_nonstiff')
    nv = equations(root, 'heart_new_valve')
    comps = []

    # ---- clocks -------------------------------------------------------------------------------
    clock_map = {'T_wCont': 'T_eff'}
    comps.append(build_component('cardiac_clock_type', [
        ('t', 'second', 'in'), ('T', 'second', 'in'), ('T_eff', 'second', 'out'),
        ('s', 'dimensionless', 'out', 0.0), ('mt', 'dimensionless', 'out'),
    ], [renamed(wc['d/dt s'], clock_map), renamed(wc['mt'], clock_map),
        parse_math('<eq/><ci>T_eff</ci><ci>T</ci>')]))
    comps.append(build_component('cardiac_clock_controlled_type', [
        ('t', 'second', 'in'), ('T', 'second', 'in'), ('T_min', 'second', 'in'),
        ('Delta_T', 'second', 'in'), ('T_eff', 'second', 'out'),
        ('s', 'dimensionless', 'out', 0.0), ('mt', 'dimensionless', 'out'),
    ], [renamed(wc['d/dt s'], clock_map), renamed(wc['mt'], clock_map), renamed(wc['T_wCont'], clock_map)]))

    # ---- chamber (from the right atrium of heart_simple_wCont; ventricles share the law) -----
    ch = {'chi_a': 'chi', 'chi_afloor': 'chi_floor', 'chi_afloor_final': 'chi_floor_final', 'e_a': 'e',
          'T_ac': 'T_c', 'T_ar': 'T_r', 'T_ac_wCont': 'T_c_eff', 'T_ar_wCont': 'T_r_eff', 'T_wCont': 'T_eff',
          't_astart': 't_start', 't_astart_norm': 't_start_norm',
          'E_ra_A': 'E_A_eff', 'E_ra_B': 'E_B', 'q_ra': 'q', 'q_ra_us': 'q_us', 'u_ra': 'u',
          'u_intraThor': 'u_ext',
          'E_lv_A_wCont': 'E_A_eff', 'E_lv_A': 'E_A', 'Delta_E_lv': 'Delta_E'}
    comps.append(build_component('chamber_type', [
        ('t', 'second', 'in'), ('T', 'second', 'in'), ('T_eff', 'second', 'in'), ('mt', 'dimensionless', 'in'),
        ('T_c', 'second', 'in'), ('T_r', 'second', 'in'), ('t_start', 'second', 'in'),
        ('T_c_eff', 'second', None), ('T_r_eff', 'second', None), ('t_start_norm', 'dimensionless', None),
        ('chi', 'dimensionless', 'out', 0.0), ('chi_floor', 'dimensionless', 'out'),
        ('chi_floor_final', 'dimensionless', None), ('e', 'dimensionless', 'out'),
        ('eps_1', 'dimensionless', None, 0.07), ('eps_2', 'dimensionless', None, 0.02),
        ('E_A', 'J_per_m6', 'in'), ('E_B', 'J_per_m6', 'in'), ('Delta_E', 'J_per_m6', 'in'),
        ('E_A_eff', 'J_per_m6', 'out'), ('q_us', 'm3', 'in'), ('q_init', 'm3', 'in'),
        ('u_ext', 'J_per_m3', 'in'), ('u', 'J_per_m3', 'out'),
        ('q', 'm3', 'out', 'q_init'), ('v_in', 'm3_per_s', 'in'), ('v_out', 'm3_per_s', 'in'),
    ], [renamed(wc['chi_afloor'], ch), renamed(wc['d/dt chi_a'], ch), renamed(wc['chi_afloor_final'], ch),
        renamed(wc['e_a'], ch), renamed(wc['T_ac_wCont'], ch), renamed(wc['T_ar_wCont'], ch),
        renamed(wc['t_astart_norm'], ch), renamed(wc['E_lv_A_wCont'], ch), renamed(wc['u_ra'], ch),
        parse_math('<eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>q</ci></apply>'
                   '<apply><minus/><ci>v_in</ci><ci>v_out</ci></apply>')]))

    # ---- valves (from the tricuspid valve) ------------------------------------------------------
    vm = {'zeta_trv_pre': 'zeta_pre', 'zeta_trv': 'zeta', 'A_eff_trv': 'A_eff', 'A_nn_trv': 'A_nn',
          'M_st_trv': 'M_st', 'M_rg_trv': 'M_rg', 'K_vo_trv': 'K_vo', 'K_vc_trv': 'K_vc', 'B_trv': 'B',
          'L_trv': 'L', 'v_trv': 'v', 'u_ra': 'u_in', 'u_rv': 'u_out', 'R_lin_trv': 'R_lin',
          'v_ref_trv': 'v_ref', 'Rmod_valve_trv': 'Rmod_valve', 'vel_trv': 'vel'}
    valve_vars = [
        ('t', 'second', 'in'), ('u_in', 'J_per_m3', 'in'), ('u_out', 'J_per_m3', 'in'),
        ('v', 'm3_per_s', 'out', 0), ('zeta_pre', 'dimensionless', 'out', 0), ('zeta', 'dimensionless', 'out'),
        ('A_eff', 'm2', 'out'), ('A_nn', 'm2', 'in'), ('M_st', 'dimensionless', 'in'),
        ('M_rg', 'dimensionless', 'in'), ('K_vo', 'm3_per_Js', 'in'), ('K_vc', 'm3_per_Js', 'in'),
        ('B', 'Js2_per_m9', 'out'), ('L', 'Js2_per_m6', 'out'), ('l_eff', 'metre', 'in'),
        ('rho', 'Js2_per_m5', 'in'), ('eps_m4', 'm4', None, '1e-14'), ('eps_m2', 'm2', None, '1e-14'),
        ('vel', 'm_per_s', 'out'),
    ]
    common = [renamed(wc[k], vm) for k in ('d/dt zeta_trv_pre', 'zeta_trv', 'A_eff_trv', 'L_trv')]
    vel = renamed(nv['vel_trv'], vm)
    comps.append(build_component('valve_type', valve_vars, common + [
        renamed(wc['B_trv'], vm), renamed(wc['d/dt v_trv'], vm), copy.deepcopy(vel)]))
    comps.append(build_component('valve_linear_type', valve_vars + [
        ('v_ref', 'm3_per_s', 'in'), ('R_lin', 'Js_per_m6', 'out')], [copy.deepcopy(e) for e in common] + [
        renamed(wc['B_trv'], vm), renamed(ns['R_lin_trv'], vm), renamed(ns['d/dt v_trv'], vm), copy.deepcopy(vel)]))
    comps.append(build_component('valve_rmod_type', valve_vars + [
        ('Rmod_valve', 'dimensionless', 'in')], [copy.deepcopy(e) for e in common] + [
        renamed(nv['B_trv'], vm), renamed(wc['d/dt v_trv'], vm), copy.deepcopy(vel)]))

    # ---- shunt resistor and flow nodes ----------------------------------------------------------
    comps.append(build_component('resistor_type', [
        ('u_in', 'J_per_m3', 'in'), ('u_out', 'J_per_m3', 'in'), ('R', 'Js_per_m6', 'in'), ('v', 'm3_per_s', 'out'),
    ], [parse_math('<eq/><ci>v</ci><apply><divide/><apply><minus/><ci>u_in</ci><ci>u_out</ci></apply><ci>R</ci></apply>')]))
    comps.append(build_component('flow_merge_type', [
        ('v_in', 'm3_per_s', 'in'), ('u_d', 'J_per_m3', 'in'), ('v_out', 'm3_per_s', 'out'), ('u', 'J_per_m3', 'out'),
    ], [parse_math('<eq/><ci>v_out</ci><ci>v_in</ci>'), parse_math('<eq/><ci>u</ci><ci>u_d</ci>')]))
    comps.append(build_component('flow_split_type', [
        ('u_in', 'J_per_m3', 'in'), ('v_out', 'm3_per_s', 'in'), ('v', 'm3_per_s', 'out'), ('u', 'J_per_m3', 'out'),
    ], [parse_math('<eq/><ci>v</ci><ci>v_out</ci>'), parse_math('<eq/><ci>u</ci><ci>u_in</ci>')]))

    model = etree.Element(f'{{{C}}}model', name='modules', nsmap={None: C, 'cellml': C})
    for comp in comps:
        model.append(comp)
    etree.indent(model, space='    ')
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'wb') as f:
        f.write(etree.tostring(model, xml_declaration=True, encoding='UTF-8', pretty_print=True))
    print(f'wrote {os.path.relpath(OUT, REPO)} ({len(comps)} components)')



OLD_LAYOUT_SOURCE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                 *'modules/heart/heart_modules.cellml'.split('/'))


if __name__ == '__main__':
    if not os.path.isfile(OLD_LAYOUT_SOURCE):
        sys.exit(f'{os.path.basename(__file__)} is a pre-versions tool: {OLD_LAYOUT_SOURCE} is gone '
                 '(see its docstring and modules/README.md)')
    main()
