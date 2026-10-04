'''
Vessel supermodules built from the lumped constitutive modules, one per monolithic vessel
version that splits into a compliance (P(A)), resistance and inertance, plus the empty
supermodules (templates) a vessel is assembled from in PhLynx.

For a monolithic version <module_type>/<version> this writes <module_type>/<version>_lumped: a
supermodule whose submodules are, in flow order,

    vp   C -> R(pv) -> I                 inlet compliance, outlet inertance
    pv   I -> R(vp) -> C
    vv   C_p -> R(pv) -> I -> C_d         the compliance in two halves
    pp   I_p -> R_p(vp) -> C -> R_d(pv) -> I_d
    noI  C -> R(pp), R(pp) -> C, C_p -> R(pp) -> C_d, R_p(pp) -> C -> R_d(pp), C alone, R alone

plus V (vessel_volume): the whole vessel's volume (and geometry), fed by every compliance, for
hosts that read the vessel's volume or radius and for resistance laws that need the whole volume.

Its default instance carries the monolithic default instance's values (TODOs filled with the
representative values in FILL), and its spec says it reproduces the monolithic version
(supermodule.equivalent_version, checked by cam_testing.supermodule).
'''
import copy
import csv
import os
from dataclasses import dataclass, field

from cam_testing.harness import VESSEL
from cam_testing.library import MODULES_DIR, load_version, write_spec

from .constitutive import ALL

CONSTITUTIVE = {(v.module_type, v.version): v for v in ALL}

# representative values for parameters the monolithic default instances leave TODO
FILL = {
    'R': 1e7, 'C': 1e-8, 'I': 1e5, 'q_0': 1e-5, 'u_0': 10600, 'u_ext': 0, 'q_C_init': 0,
    'R_T': 1e8, 'C_T': 1e-8, 'I_T': 1e5, 'I_T_1': 1e5, 'I_T_2': 1e5, 'I_T1': 1e5, 'I_T2': 1e5, 'R_T1': 2e7, 'R_T2': 8e7,
    'R_1_T': 2e7, 'R_2_T': 8e7, 'Zc': 1e7, 'frac_R_T_1_of_R_T': 0.2, 'q_us': 1e-4, 'q_us_0': 1e-4,
    'Rp': 0, 'fracR': 0.3, 'E': 1.6e6, 'l': 0.1, 'r_0': 0.01, 'r': 0.01, 'theta': 0, 'I_scale': 1, 'R_flag': 1,
    'm_tube': 0.5, 'n_tube': 0, 'R_local_multiplier': 1, 'blood_ctl_dir': 1,
    'rho': 1050, 'mu': 0.004, 'g': 9.81, 'beta_g': 1,
    'a_vessel': 0.2802, 'b_vessel': -505.3, 'c_vessel': 0.1324, 'd_vessel': -11.14,
}
# boundary values for the supermodule's open boundary conditions (and the monolithic version's, in
# the equivalence check)
BC_VALUES = {
    'v_in': 2e-6, 'v_out': 1e-6, 'u_in': 11000, 'u_out': 10000, 'u_ext': 0, 'K_tube': 1e5, 'K_tube_visco': 1e5,
    'Delta_R': 0.1, 'Delta_q_us': 0.1, 'Delta_C': 0.1, 'v_blood_intake': 1e-7,
}
# variables never shared over submodules (each submodule's own)
OWN = {'fraction', 'gravity_factor', 'q_C_init', 'v_init', 'u_in', 'u_out', 'v_in', 'v_out', 'R_other', 'q_ref'}
GLOBALS = {'rho', 'mu', 'g', 'beta_g', 'a_vessel', 'b_vessel', 'c_vessel', 'd_vessel'}


@dataclass
class Sub:
    name: str
    module_type: str
    version: str
    params: dict = field(default_factory=dict)     # constant -> value (others: the version's default)

    @property
    def cv(self):
        return CONSTITUTIVE[(self.module_type, self.version)]


@dataclass
class Layout:
    chain: list                                    # Subs in flow order (V excluded)
    volume: str = 'nn'                             # vessel_volume version
    outputs: dict = field(default_factory=dict)    # monolithic output -> '<sub>/<var>'
    bc_values: dict = field(default_factory=dict)  # open boundary condition -> value, beyond BC_VALUES
    notes: list = field(default_factory=list)


def C(name='C', version='vv_linear', **params):
    return Sub(name, 'compliance', version, params)


def R(name='R', version='pv_linear', **params):
    return Sub(name, 'resistance', version, params)


def I(name='I', version='pp_linear', **params):
    return Sub(name, 'inertance', version, params)


# ---- the monolithic families -------------------------------------------------------------------

def _fit(chain):
    '''Drops, from each submodule's parameters, the ones its version does not have (the family
    builders pass the monolithic version's parameters by name).'''
    for s in chain:
        have = {v.name for v in s.cv.variables}
        s.params = {k: v for k, v in s.params.items() if k in have}
    return chain


def _p(M, *names):
    return {n: M[n] for n in names if n in M}


def lumped_rci(bc, M, comp='vv_linear', r_law='linear', i_law='linear', geometric=False, extra_c=None,
               state='q_C', u0=None, compare_total=False):
    '''The R-C-I vessels, by inlet/outlet convention. Parameter names as the monolithic versions':
    C, R, I, q_0, u_0, u_ext (linear) or E, l, r_0, theta, I_scale, R_flag (geometric).'''
    cpar = _p(M, 'C', 'q_0', 'u_0', 'E', 'l', 'r_0', 'm_tube', 'n_tube', 'q_us_0', 'q_C_init')
    if u0 is not None:
        cpar['u_0'] = u0
    cpar.update(extra_c or {})
    rpar = _p(M, 'R', 'l', 'r_0', 'R_flag')
    ipar = _p(M, 'I', 'l', 'r_0', 'I_scale', 'theta')
    if i_law == 'geometric':
        ipar.setdefault('I_scale', 1)
        ipar.setdefault('theta', 0)
    if r_law in ('geometric', 'nonlinear'):
        rpar.setdefault('R_flag', 1)
    out = {'u': 'C/u', 'v': 'I/v', state: f'C/{state}'}
    if bc == 'vp':
        chain = [C(version=comp, **cpar), R(version=f'pv_{r_law}', **rpar), I(version=f'pp_{i_law}', **ipar)]
    elif bc == 'pv':
        chain = [I(version=f'pp_{i_law}', **ipar), R(version=f'vp_{r_law}', **rpar), C(version=comp, **cpar)]
    elif bc == 'vv':
        half = dict(cpar, fraction=0.5)
        chain = [C('C_p', comp, **half), R(version=f'pv_{r_law}', **rpar), I(version=f'pp_{i_law}', **ipar),
                 C('C_d', comp, **half)]
        out = {'u': 'C_p/u', 'u_d': 'C_d/u', 'v': 'I/v', state: f'C_p/{state}', f'{state}_d': f'C_d/{state}'}
    elif bc == 'pp':
        rh, ih = dict(rpar, fraction=0.5), dict(ipar, fraction=0.5)
        chain = [I('I_p', f'pp_{i_law}', **ih), R('R_p', f'vp_{r_law}', **rh), C(version=comp, **cpar),
                 R('R_d', f'pv_{r_law}', **rh), I('I_d', f'pp_{i_law}', **dict(ih, gravity_factor=0) if
                                                    i_law == 'geometric' else ih)]
        out = {'u': 'C/u', 'v': 'I_p/v', 'v_d': 'I_d/v', state: f'C/{state}'}
    else:
        raise ValueError(bc)
    if compare_total:
        # the stressed volume of a tiny compliance (the microvasculature: C ~ 1e-17 m^6/J) is a
        # difference of near-equal numbers; its total volume is compared instead
        out = {k: v for k, v in out.items() if not k.startswith(state)}
        out['q'] = 'V/q'
    return Layout(_fit(chain), volume='nn_geometric' if geometric else 'nn', outputs=out)


def lumped_noI(bc, M, comp='vv_linear_novisco', r_law='linear', geometric=False, q_state='q', cpar=None, rpar=None):
    '''Vessels without inertance: the resistance in flow form.'''
    cpar = cpar if cpar is not None else _p(M, 'C', 'q_0', 'u_0', 'E', 'l', 'r_0', 'q_C_init')
    rpar = rpar if rpar is not None else _p(M, 'R', 'l', 'r_0', 'R_flag')
    if r_law == 'geometric':
        rpar.setdefault('R_flag', 1)
    if bc == 'vp':
        chain, out = [C(version=comp, **cpar), R(version=f'pp_{r_law}', **rpar)], {'u': 'C/u', 'v': 'R/v'}
    elif bc == 'pv':
        chain, out = [R(version=f'pp_{r_law}', **rpar), C(version=comp, **cpar)], {'u': 'C/u', 'v': 'R/v'}
    elif bc == 'vv':
        if 'R' in rpar or r_law == 'geometric':
            half = dict(cpar, fraction=0.5)
            chain = [C('C_p', comp, **half), R(version=f'pp_{r_law}', **rpar), C('C_d', comp, **half)]
            out = {'u': 'C_p/u', 'u_d': 'C_d/u', 'v': 'R/v', 'q_C': 'C_p/q_C', 'q_C_d': 'C_d/q_C'}
        else:
            chain, out = [C(version=comp, **cpar)], {'u': 'C/u'}
    elif bc == 'pp':
        if cpar:
            rh = dict(rpar, fraction=0.5)
            chain = [R('R_p', f'pp_{r_law}', **rh), C(version=comp, **cpar), R('R_d', f'pp_{r_law}', **rh)]
            out = {'u': 'C/u', 'v': 'R_p/v', 'v_d': 'R_d/v'}
        else:
            chain, out = [R(version=f'pp_{r_law}', **rpar)], {'v': 'R/v'}
    else:
        raise ValueError(bc)
    if q_state == 'q' and any(s.module_type == 'compliance' for s in chain):
        out['q'] = 'V/q'
    elif q_state == 'q_C':
        out.update({'q_C': 'C/q_C'} if 'C' in [s.name for s in chain] else {})
    return Layout(_fit(chain), volume='nn_geometric' if geometric else 'nn', outputs=out)


def _q_init(M, q0_name):
    '''q_C_init for a monolithic version whose volume state q starts at q_init.'''
    return float(M['q_init']) - float(M[q0_name])


def simple_noI(bc):
    def build(M):
        cpar = _p(M, 'C', 'q_0', 'u_0')
        if bc == 'pp':
            cpar = {}
        else:
            cpar['u_ext'] = 0
            cpar['q_C_init'] = _q_init(M, 'q_0')
        layout = lumped_noI(bc, M, cpar=cpar)
        layout.bc_values['u_ext'] = 0
        return layout
    return build


def micro_noI(bc):
    def build(M):
        layout = lumped_noI(bc, M, comp='vv_linear_geometric_novisco', r_law='geometric', geometric=True,
                            q_state='q_C')
        if bc == 'pp' and 'u_0' in M:
            # the instance's u_in and u_out sit either side of u_0, which starts the vessel at its
            # equilibrium (nothing moves); drive it off it
            layout.bc_values['u_out'] = M['u_0'] - 2000
        return layout
    return build


def capillary_r(bc):
    '''capillary/constriction: geometric C and R from a radius parameter r, u = u_ext + q_C/C.'''
    def build(M):
        geo = {'E': M['E'], 'l': M['l'], 'r_0': M['r'], 'u_0': 0}
        if 'q_C_init' in M:
            geo['q_C_init'] = M['q_C_init']
        layout = lumped_noI(bc, M, comp='vv_linear_geometric_novisco', r_law='geometric', geometric=True,
                            q_state='q_C', cpar=geo, rpar={'l': M['l'], 'r_0': M['r'], 'R_flag': 1})
        layout.notes.append('the monolithic radius parameter r is the submodules\' r_0')
        return layout
    return build


def rci(bc, **kw):
    def build(M):
        return lumped_rci(bc, M, **kw)
    return build


def nonlinear(bc, r_law='nonlinear', visco=False):
    def build(M):
        comp = 'vv_nonlinear_geometric_visco' if visco else 'vv_nonlinear_geometric'
        return lumped_rci(bc, M, comp=comp, r_law=r_law, i_law='geometric', geometric=True)
    return build


def controlled(comp, u0=0):
    def build(M):
        par = {'q_us_0': M.get('q_us_0', M.get('q_0'))}
        layout = lumped_rci('vp', M, comp=comp, extra_c=par, state='q_C_change' if 'stressed' not in comp else 'q_C',
                            u0=u0)
        if 'stressed' not in comp:
            layout.outputs['q_C'] = 'C/q_C'
        return layout
    return build


def vol_flux(M):
    layout = lumped_rci('vv', M)
    layout.chain[0].version = 'vv_linear_volume_flux'
    layout.chain[0].params['blood_ctl_dir'] = M.get('blood_ctl_dir', 1)
    return layout


def arterial_simple_vp_nonlinear(M):
    cpar = _p(M, 'q_0', 'u_0', 'm_tube', 'n_tube')
    cpar['q_C_init'] = _q_init(M, 'q_0')
    chain = [C(version='vv_nonlinear', **cpar), R(version='pv_volume_scaled', R=M['R'], q_ref=M['q_init']),
             I(version='pp_linear', I=M['I'])]
    return Layout(chain, outputs={'u': 'C/u', 'v': 'I/v', 'q': 'V/q'})


def arterial_simple_novisco_pp(M):
    rp = float(M['Rp'])
    frac = rp / float(M['R']) if rp > 0 else (float(M['fracR']) if float(M['fracR']) >= 0 else 0.5)
    cpar = _p(M, 'C', 'q_0', 'u_0')
    chain = [I('I_p', 'pp_linear', I=M['I'], fraction=0.5), R('R_p', 'vp_linear', R=M['R'], fraction=frac),
             C(version='vv_linear_novisco', **cpar), R('R_d', 'pv_linear', R=M['R'], fraction=1 - frac),
             I('I_d', 'pp_linear', I=M['I'], fraction=0.5)]
    layout = Layout(chain, outputs={'u': 'C/u', 'v': 'I_p/v', 'v_d': 'I_d/v', 'q_C': 'C/q_C'})
    layout.notes.append('the monolithic Rp/fracR split of R is the two resistances\' fractions (R_p: Rp/R if Rp > 0, '
                        'else fracR if fracR >= 0, else 1/2); set fraction_R_p and fraction_R_d to change it')
    return layout


def venous_Rctl(M):
    layout = lumped_rci('vp', M, r_law='controlled_additive')
    return layout


def resistor(M):
    return Layout([R(version='pp_linear', R=M['R'])], outputs={'v': 'R/v'})


# ---- terminals ---------------------------------------------------------------------------------

def _terminal_c(M, comp='vv_linear', **extra):
    par = {'C': M.get('C_T'), 'q_0': M['q_us'], 'u_0': 0, 'q_C_init': _q_init(M, 'q_us')}
    par.update(extra)
    return _fit([C(version=comp, **{k: v for k, v in par.items() if v is not None})])[0]


def terminal_pp(r1, r2, i1, i2, comp='vv_linear', law='linear', r_extra=None, notes=()):
    '''I_p -> R_p(vp) -> C -> R_d(pv) -> I_d, values given as functions of the instance.'''
    def build(M):
        rx = r_extra(M) if r_extra else ({}, {})
        chain = [I('I_p', 'pp_linear', **i1(M)), R('R_p', f'vp_{law}', **r1(M), **rx[0]), _terminal_c(M, comp),
                 R('R_d', f'pv_{law}', **r2(M), **rx[1]), I('I_d', 'pp_linear', **i2(M))]
        return Layout(chain, outputs={'u': 'C/u', 'v': 'I_p/v', 'v_T': 'I_d/v', 'q': 'V/q'}, notes=list(notes))
    return build


def terminal_pp_noI(r1, r2, comp='vv_linear_novisco', law='linear', r_extra=None):
    def build(M):
        rx = r_extra(M) if r_extra else ({}, {})
        chain = [R('R_p', f'pp_{law}', **r1(M), **rx[0]), _terminal_c(M, comp), R('R_d', f'pp_{law}', **r2(M), **rx[1])]
        return Layout(chain, outputs={'u': 'C/u', 'v': 'R_p/v', 'v_T': 'R_d/v', 'q': 'V/q'})
    return build


def terminal_vp(r2, i=None, law='linear', comp='vv_linear_novisco', r_extra=None):
    '''C -> R (pv with an inertance after it, else pp).'''
    def build(M):
        rx = r_extra(M) if r_extra else {}
        if i:
            chain = [_terminal_c(M, comp), R(version=f'pv_{law}', **r2(M), **rx), I(version='pp_linear', **i(M))]
            out = {'u': 'C/u', 'v_T': 'I/v', 'q': 'V/q'}
        else:
            chain = [_terminal_c(M, comp), R(version=f'pp_{law}', **r2(M), **rx)]
            out = {'u': 'C/u', 'v_T': 'R/v', 'q': 'V/q'}
        return Layout(chain, outputs=out)
    return build


def _half(key):
    return lambda M: {'R': M[key], 'fraction': 0.5}


def _whole(key, var='R'):
    return lambda M: {var: M[key]}


def _frac(M, first):
    f = float(M['frac_R_T_1_of_R_T'])
    return {'R': M['R_T'], 'fraction': f if first else 1 - f}


def _local(M):
    return ({'R_local_multiplier': M.get('R_local_multiplier', 1)},) * 2


def _share(M):
    return ({'R_other': M['R_T2'], 'R_local_multiplier': M.get('R_local_multiplier', 1)},
            {'R_other': M['R_T1'], 'R_local_multiplier': M.get('R_local_multiplier', 1)})


def _terminal2_pp(M):
    M = dict(M, q_us=0, q_init=0)
    layout = terminal_pp(_half('R_T'), _half('R_T'), lambda M: {'I': 1}, lambda M: {'I': 1})(M)
    layout.outputs = {'u': 'C/u', 'v': 'I_p/v', 'v_T': 'I_d/v', 'q_T': 'C/q_C'}
    layout.bc_values['u_out'] = 0
    layout.notes.append('the monolithic version fixes u_out = 0 and I_T = 1: the supermodule\'s u_out and I '
                        'parameters, set so in its default instance')
    return layout


def _zcir(M):
    rt2 = float(M['R_T']) - float(M['Zc'])
    return terminal_vp(lambda M: {'R': rt2}, lambda M: {'I': M['I_T']}, law='controlled_share',
                       r_extra=lambda M: {'R_other': M['Zc'], 'R_local_multiplier': M.get('R_local_multiplier', 1)})(M)


FAMILIES = {
    **{('arterial_simple', bc): rci(bc) for bc in ('vp', 'pv', 'vv', 'pp')},
    ('arterial_simple', 'vp_nonlinear'): arterial_simple_vp_nonlinear,
    ('arterial_simple', 'vp_usV_C_ctl'): controlled('vv_linear_controlled_stressed', u0=None),
    ('arterial_simple_novisco', 'pp'): arterial_simple_novisco_pp,
    **{('capillary_simple', bc): rci(bc) for bc in ('vp', 'pv', 'vv', 'pp')},
    ('capillary_simple_vol_ctl_intake', 'vv'): vol_flux,
    ('capillary_simple_vol_ctl_uptake', 'vv'): vol_flux,
    **{(t, f'{bc}_micro'): rci(bc, compare_total=True) for t in ('arteriole', 'venule', 'vein')
       for bc in ('vp', 'pv', 'vv', 'pp')},
    **{('capillary', f'{bc}_micro'): rci(bc, compare_total=True) for bc in ('pv', 'vv')},
    **{(t, f'{bc}_micro_noI'): micro_noI(bc) for t in ('arteriole', 'venule', 'vein', 'capillary') for bc in ('vp', 'pv')},
    **{('artery', f'{bc}_micro'): micro_noI(bc) for bc in ('vp', 'pv', 'vv', 'pp')},
    **{('capillary', v): capillary_r(v[:2]) for v in ('vp', 'pp', 'vp_micro', 'pp_micro')},
    ('constriction', 'pv'): capillary_r('pv'),
    **{('arterial', bc): rci(bc, comp='vv_linear_geometric', r_law='geometric', i_law='geometric', geometric=True)
       for bc in ('vp', 'pv', 'vv', 'pp')},
    ('arterial', 'vv_novisco'): rci('vv', comp='vv_linear_geometric_novisco', r_law='geometric', i_law='geometric',
                                    geometric=True),
    **{(t, f'{bc}_nonlinear'): nonlinear(bc) for t in ('arterial',) for bc in ('vp', 'pv', 'vv', 'pp')},
    **{('arterial', f'{bc}_nonlinear_constR'): nonlinear(bc, r_law='geometric') for bc in ('vp', 'pv', 'vv', 'pp')},
    **{('arterial', f'{bc}_nonlinear_visco'): nonlinear(bc, visco=True) for bc in ('vp', 'pv', 'vv', 'pp')},
    ('venous', 'pp_nonlinear'): nonlinear('pp'),
    ('venous', 'vp_nonlinear'): nonlinear('vp'),
    **{('venous', bc): rci(bc) for bc in ('pp', 'pv', 'vv')},
    ('venous', 'vp_woCont'): rci('vp'),
    ('venous', 'vp'): controlled('vv_linear_controlled'),
    ('venous', 'vp_wCont'): controlled('vv_linear_controlled'),
    ('venous_novisco', 'vp'): controlled('vv_linear_controlled_novisco'),
    ('venous_novisco', 'vp_woCont'): rci('vp', comp='vv_linear_novisco'),
    ('venous_Rctl', 'vp_woCont'): venous_Rctl,
    **{(t, bc): simple_noI(bc[:2]) for t in ('simple', 'simple_2') for bc in ('vp_noI', 'pv_noI', 'vv_noI', 'pp_noI')},
    ('resistor', 'pp'): resistor,
    ('terminal', 'pp'): terminal_pp(_half('R_T'), _half('R_T'), lambda M: {'I': 1}, lambda M: {'I': 1},
                                    notes=['the monolithic version fixes I_T = 1: the inertances\' I, set so in the '
                                           'default instance']),
    ('terminal', 'pp_spec_I'): terminal_pp(_half('R_T'), _half('R_T'), _whole('I_T', 'I'), _whole('I_T', 'I')),
    ('terminal', 'pp_I'): terminal_pp(_whole('R_1_T'), _whole('R_2_T'), _whole('I_T', 'I'), _whole('I_T', 'I')),
    ('terminal', 'pp_R1ICR2I'): terminal_pp(_whole('R_T1'), _whole('R_T2'), _whole('I_T1', 'I'), _whole('I_T2', 'I'),
                                            comp='vv_linear_novisco'),
    ('terminal', 'pp_RICRI'): terminal_pp(lambda M: _frac(M, True), lambda M: _frac(M, False), _whole('I_T_1', 'I'),
                                          _whole('I_T_2', 'I'), comp='vv_linear_novisco'),
    ('terminal', 'pp_RICRI_nonlinear'): terminal_pp(lambda M: _frac(M, True), lambda M: _frac(M, False),
                                                    _whole('I_T_1', 'I'), _whole('I_T_2', 'I'), comp='vv_nonlinear'),
    **{('terminal', v): terminal_pp(_half('R_T'), _half('R_T'), _whole('I_T', 'I'), _whole('I_T', 'I'),
                                    law='controlled', r_extra=_local)
       for v in ('pp_wCont', 'pp_wCont_no_local', 'pp_wLocal')},
    ('terminal', 'pp_RCR_wCont'): terminal_pp_noI(_half('R_T'), _half('R_T'), law='controlled', r_extra=_local),
    ('terminal', 'pp_R1CR2_wCont'): terminal_pp_noI(_whole('R_T1'), _whole('R_T2'), law='controlled_share',
                                                    r_extra=_share),
    ('terminal2', 'pp'): _terminal2_pp,
    ('single_terminal', 'vp_wCont'): terminal_vp(_whole('R_T2'), law='controlled', r_extra=lambda M: _local(M)[0]),
    ('single_terminal', 'vp_R1CR2_wCont'): terminal_vp(_whole('R_T2'), law='controlled_share',
                                                       r_extra=lambda M: _share(M)[1]),
    ('single_terminal', 'vp_R1CR2_wI'): terminal_vp(_whole('R_T2'), _whole('I_T', 'I')),
    ('single_terminal', 'vp_R1CR2_wI_wCont'): terminal_vp(_whole('R_T2'), _whole('I_T', 'I'), law='controlled_share',
                                                          r_extra=lambda M: _share(M)[1]),
    ('single_terminal', 'vp_ZCIR_wCont'): _zcir,
}

PHLYNX_NO_SUPERMODULES = ('PhLynx has no supermodule support yet: it builds models from CellML module configs only '
                          '(see modules/README.md, "Lumped vessels", for what it needs)')

# vessels that do not split into R, C and I
NOT_SPLIT = {
    ('rigid_vessel', 'vv'): 'u = (v_in - v_out)/R: a pressure set by the net flow, which is neither a compliance nor '
                            'a resistance between two pressures',
    ('terminal', 'pp_controller'): 'a learning flow controller (aa*/cc* states) inside the terminal',
    ('terminal2', 'pp_controller'): 'a learning flow controller inside the terminal',
    ('terminalOL', 'pp_controller'): 'a learning flow controller inside the terminal',
    ('terminalOL', 'pp_wFB'): 'a learning flow controller inside the terminal',
}


# ---- building a supermodule version ------------------------------------------------------------

def mono_values(version):
    '''The monolithic default instance's values, TODOs filled from FILL (by name).'''
    values = {}
    for p in version.parameters():
        v = p.value
        if v in (None, '', 'TODO') or str(v).upper().startswith('TODO') or str(v).startswith('EMPTY'):
            if p.variable_name in FILL:
                v = FILL[p.variable_name]
            elif p.variable_name in BC_VALUES:
                v = BC_VALUES[p.variable_name]
            else:
                continue
        values[p.variable_name] = float(v)
    # constants and globals of the config that the instance has no row for
    for name, _units, _access, kind in (version.config.get('variables_and_units') or []):
        if name not in values and kind in ('constant', 'global_constant') and name in FILL:
            values[name] = FILL[name]
    return values


def edges(layout):
    '''(from, to) submodule connections: the chain, every compliance into V, V into the resistances
    whose law needs the whole volume.'''
    names = [s.name for s in layout.chain]
    out = list(zip(names, names[1:]))
    if layout.volume:
        out += [(s.name, 'V') for s in layout.chain if s.module_type == 'compliance']
        out += [('V', s.name) for s in layout.chain if s.module_type == 'resistance'
                and s.version.split('_', 1)[1] in ('nonlinear', 'volume_scaled')]
    return out


def all_subs(layout):
    subs = list(layout.chain)
    if layout.volume and any(s.module_type == 'compliance' for s in layout.chain):
        geometry = {}
        if layout.volume == 'nn_geometric':
            # the vessel's r_0 and l, as its compliances have them
            for s in layout.chain:
                geometry.update({k: s.params[k] for k in ('r_0', 'l') if k in s.params and k not in geometry})
        subs.append(Sub('V', 'vessel_volume', layout.volume, geometry))
    else:
        layout.volume = None
    return subs


def _ports(cv, side):
    return cv.entrance_ports if side == 'entrance' else cv.exit_ports


def open_boundary_conditions(layout):
    '''{sub: [boundary-condition variables no internal connection closes]}.'''
    subs = {s.name: s for s in all_subs(layout)}
    closed = {name: set() for name in subs}
    for a, b in edges(layout):
        out_ports = {p['port_type']: p['variables'] for p in subs[a].cv.exit_ports}
        for p in subs[b].cv.entrance_ports:
            if p['port_type'] in out_ports:
                closed[a].update(out_ports[p['port_type']])
                closed[b].update(p['variables'])
    return {name: [v.name for v in s.cv.variables if v.kind == 'boundary_condition' and v.name not in closed[name]]
            for name, s in subs.items()}


def routes(layout, mono_config):
    '''The supermodule's routes for the monolithic version's ports, and the ports no submodule offers.'''
    subs = all_subs(layout)
    chain = layout.chain
    result, dropped = {'inputs': {}, 'outputs': {}}, []
    for side, key, ports in (('entrance', 'inputs', mono_config.get('entrance_ports') or []),
                             ('exit', 'outputs', mono_config.get('exit_ports') or [])):
        for p in ports:
            t = p['port_type']
            if t in result[key]:
                continue
            if t == 'vessel_port':
                result[key][t] = chain[0].name if key == 'inputs' else chain[-1].name
                continue
            offering = [s.name for s in subs if t in [q['port_type'] for q in _ports(s.cv, side)]]
            if t == 'volume_port' and 'V' in offering:
                offering = ['V']
            if not offering:
                dropped.append(f'{key[:-1]} {t} {p.get("variables", [])}')
                continue
            result[key][t] = offering[0] if len(offering) == 1 else offering
    return {k: v for k, v in result.items() if v}, dropped


def instance_rows(layout, mono_rows, bc_values, M):
    '''The supermodule's default instance: shared rows for a constant every submodule that has it
    gives the same value, <var>_<sub> rows otherwise, globals by name, and the open boundary
    conditions (the monolithic instance's value where it has one, else a test value).'''
    subs = all_subs(layout)
    units = {}
    values = {}                                 # var -> {sub: value}
    for s in subs:
        for var, value in s.params.items():
            cv_var = next((v for v in s.cv.variables if v.name == var), None)
            if cv_var is None:
                raise KeyError(f'{s.module_type}/{s.version} has no variable {var}')
            units[var] = cv_var.units
            values.setdefault(var, {})[s.name] = value
    open_bcs = open_boundary_conditions(layout)
    boundary = set()
    for name, bcs in open_bcs.items():
        s = next(x for x in subs if x.name == name)
        for var in bcs:
            units[var] = s.cv.var(var).units
            values.setdefault(var, {})[name] = bc_values.get(var, M.get(var, BC_VALUES.get(var, s.cv.var(var).default)))
            boundary.add(var)
    shared, rows = [], []
    for var, by_sub in values.items():
        have = [s.name for s in subs if any(v.name == var for v in s.cv.variables)]
        same = len(set(map(float, by_sub.values()))) == 1 and set(by_sub) == set(have)
        if var not in OWN and same and len(have) > 1:
            shared.append(var)
            rows.append([var, units[var], next(iter(by_sub.values())), _ref(var, mono_rows, var in boundary), 'no'])
        else:
            for sub, value in by_sub.items():
                rows.append([f'{var}_{sub}', units[var], value, _ref(var, mono_rows, var in boundary), 'no'])
    globals_ = sorted({v.name for s in subs for v in s.cv.variables if v.kind == 'global_constant'})
    for g in globals_:
        if g in mono_rows:
            rows.append([g, mono_rows[g][0], mono_rows[g][1], mono_rows[g][2], 'no'])
    return rows, shared, globals_


def _ref(var, mono_rows, boundary=False):
    if boundary:
        return 'cam_testing: a test boundary value (connected in a model)'
    if var in mono_rows and mono_rows[var][2] and 'TODO' not in mono_rows[var][2].upper():
        return mono_rows[var][2]
    return 'lumped split of the monolithic version (tools/lumped/vessels.py)'


def supermodule_entry(module_type, version, layout, mono_config, description, template=False):
    subs = all_subs(layout)
    conns = edges(layout)
    records = []
    for s in subs:
        rec = {'name': s.name, 'module_type': s.module_type}
        if template:
            rec['choices'] = choices(s)
        else:
            rec['module_subtype'] = s.version
            rec['instance'] = 'default'
        rec['inp_instances'] = [a for a, b in conns if b == s.name]
        rec['out_instances'] = [b for a, b in conns if a == s.name]
        records.append(rec)
    route, dropped = routes(layout, mono_config)
    entry = {'module_type': module_type, 'module_subtype': version, 'module_format': 'supermodule',
             'default_instance': 'default', 'licence': 'CC0-1.0', 'creator': [], 'description': description}
    if template:
        entry['template'] = True
    entry.update({'routes': route, 'submodules': records})
    return entry, dropped


def choices(sub):
    '''The versions of sub's module_type that fit its place in the layout.'''
    if sub.module_type == 'resistance':
        form = sub.version.split('_', 1)[0]
        return sorted(v for (t, v) in CONSTITUTIVE if t == 'resistance' and v.startswith(form + '_'))
    return sorted(v for (t, v) in CONSTITUTIVE if t == sub.module_type)


def _mono_rows(version):
    return {p.variable_name: (p.units, p.value, p.data_reference) for p in version.parameters()}


def _description(mono, layout):
    names = ' -> '.join(f'{s.name} ({s.module_type} {s.version})' for s in layout.chain)
    extra = f'; V (vessel_volume {layout.volume}) gathers the compliances\' volumes' if layout.volume else ''
    return (f'{mono.vessel_type}/{mono.BC_type} built from lumped constitutive submodules '
            f'(modules/haemodynamics/lumped_constitutive): {names}{extra}. Reproduces the monolithic version '
            f'(supermodule.equivalent_version in its spec).')


def write_vessel(mono, layout, units_library, bib):
    from .write import fresh_version_dir, write_instance, write_json, GENERATED_NOTE
    module_dir = os.path.dirname(os.path.dirname(mono.dir))
    name = f'{mono.BC_type}_lumped'
    vdir = fresh_version_dir(module_dir, name)
    stem = f'{mono.vessel_type}_{name}'
    M = mono_values(mono)
    mono_rows = _mono_rows(mono)
    rows, shared, globals_ = instance_rows(layout, mono_rows, layout.bc_values, M)
    entry, dropped = supermodule_entry(mono.vessel_type, name, layout, mono.config, _description(mono, layout))
    if shared:
        entry['shared_parameters'] = shared
    # keep key order: shared_parameters before submodules
    entry = {k: entry[k] for k in [k for k in entry if k != 'submodules'] + ['submodules']}
    write_json(os.path.join(vdir, f'{stem}_modules_config.json'), [entry])
    write_instance(vdir, rows)
    mono_bcs = [v[0] for v in mono.config.get('variables_and_units') or [] if v[3] == 'boundary_condition']
    mono_params = dict(M)
    for b in mono_bcs:
        mono_params[b] = layout.bc_values.get(b, M.get(b, BC_VALUES.get(b, 0)))
    notes = [GENERATED_NOTE] + layout.notes
    if dropped:
        notes.append('ports of the monolithic version no submodule offers (connect them to the submodule that has '
                     'the variable instead): ' + '; '.join(dropped))
    spec = {
        'module_type': mono.vessel_type, 'version': name, 'reviewed': False,
        'description': entry['description'], 'notes': ' '.join(notes),
        'expected_failures': {'phlynx_export_test': PHLYNX_NO_SUPERMODULES},
        'sim_time': mono.spec.get('sim_time', 1.0), 'pre_time': 0.0, 'dt': mono.spec.get('dt', 0.001),
        'solver': 'CVODE_myokit', 'solver_info': {'rtol': 1e-8, 'atol': 1e-14},
        'reference_solver_info': {'rtol': 1e-10, 'atol': 1e-16},
        'outputs': sorted(f'{VESSEL}_{o}' for o in set(layout.outputs.values())), 'run_parameters': {},
        'supermodule': {'globals': globals_, 'equivalent_version': {
            'version': f'{mono.vessel_type}/{mono.BC_type}',
            'monolithic_parameters': {k: mono_params[k] for k in sorted(mono_params)},
            'output_map': layout.outputs, 'tol': 1e-5,
            'solver_info': {'rtol': 1e-10, 'atol': 1e-20}}},
    }
    write_spec(vdir, stem, spec)
    with open(os.path.join(vdir, f'{stem}_references.bib'), 'w') as f:
        f.write(bib)
    return vdir


def hide_from_phlynx(mono):
    '''Marks a monolithic version available_to_phlynx: false (PhLynx offers its _lumped supermodule
    instead); the library keeps it, and models that name it still build.'''
    import json
    with open(mono.config_path) as f:
        raw = json.load(f)
    entry = raw[0]
    if entry.get('available_to_phlynx') is False:
        return
    out = {}
    for k, v in entry.items():
        out[k] = v
        if k == 'creator':
            out['available_to_phlynx'] = False
    raw[0] = out
    with open(mono.config_path, 'w') as f:
        json.dump(raw, f, indent=2)
        f.write('\n')


def _add_system_equivalences():
    '''Builds the system models rebuilt with lumped vessels (tools/lumped/systems.py) and adds
    their supermodule.equivalent entries to the holding _lumped versions' specs.'''
    from cam_testing.library import merge_spec, read_spec_files
    from . import systems
    for (mt, v), entries in systems.equivalence_entries().items():
        version = load_version(mt, v)
        spec = merge_spec(*read_spec_files(version.dir, version.stem))
        spec['supermodule'] = {**spec.get('supermodule', {}), 'equivalent': entries}
        write_spec(version.dir, version.stem, spec)


# ---- templates ----------------------------------------------------------------------------------

TEMPLATE_DIR = os.path.join(MODULES_DIR, 'haemodynamics', 'lumped_vessel')
TEMPLATES = {
    'vp_empty': ([C(), R(version='pv_linear'), I()], 'inlet compliance, resistance, outlet inertance'),
    'pv_empty': ([I(), R(version='vp_linear'), C()], 'inlet inertance, resistance, outlet compliance'),
    'vv_empty': ([C('C_p'), R(version='pv_linear'), I(), C('C_d')], 'compliance in two halves around a resistance '
                 'and inertance (give the compliances fraction 0.5)'),
    'pp_empty': ([I('I_p'), R('R_p', 'vp_linear'), C(), R('R_d', 'pv_linear'), I('I_d')], 'resistance and inertance '
                 'in two halves around a compliance (give the halves fraction 0.5)'),
    'vp_noI_empty': ([C(), R(version='pp_linear')], 'compliance then resistance, no inertance'),
    'pv_noI_empty': ([R(version='pp_linear'), C()], 'resistance then compliance, no inertance'),
    'pp_noI_empty': ([R('R_p', 'pp_linear'), C(), R('R_d', 'pp_linear')], 'compliance between two resistances, '
                     'no inertance (give the resistances fraction 0.5)'),
}
_TEMPLATE_PORTS = {
    'entrance_ports': [{'port_type': 'vessel_port'}, {'port_type': 'external_pressure_port'},
                       {'port_type': 'resistance_delta_port'}, {'port_type': 'venous_us_volume_delta_port'},
                       {'port_type': 'venous_us_and_compliance_delta_port'}, {'port_type': 'blood_uptake_port'}],
    'exit_ports': [{'port_type': t} for t in ('vessel_port', 'volume_port', 'flow_port', 'flow_feedback_port',
                                              'stiffness_port', 'visco_coeff_port', 'radius_port', 'len_port',
                                              'radius_0_and_s_port', 'pressure_and_deriv_port')],
}


def write_templates(units_library, bib):
    from .write import fresh_version_dir, write_instance, write_json, GENERATED_NOTE
    for name, (chain, what) in TEMPLATES.items():
        layout = Layout(copy.deepcopy(chain))
        vdir = fresh_version_dir(TEMPLATE_DIR, name)
        stem = f'lumped_vessel_{name}'
        description = (f'Empty lumped vessel ({what}): drag it in, then choose a version for each submodule. '
                       'Each slot names the module_type that fits it and the versions that fit its place '
                       '(choices). V (vessel_volume) gathers the compliances\' volumes; choose nn_geometric for a '
                       'vessel with material-property modules or a nonlinear resistance.')
        entry, _ = supermodule_entry('lumped_vessel', name, layout, _TEMPLATE_PORTS, description, template=True)
        write_json(os.path.join(vdir, f'{stem}_modules_config.json'), [entry])
        write_instance(vdir, [])
        write_spec(vdir, stem, {'module_type': 'lumped_vessel', 'version': name, 'reviewed': False,
                                'description': description, 'notes': GENERATED_NOTE,
                                'skip': 'a template: it has no version for its submodules until one is chosen'})
        with open(os.path.join(vdir, f'{stem}_references.bib'), 'w') as f:
            f.write(bib)


def write_all(units_library, bib, only=None):
    written, failed = [], []
    for (mt, v), build in sorted(FAMILIES.items()):
        if only and f'{mt}/{v}' not in only:
            continue
        try:
            mono = load_version(mt, v)
        except Exception as e:  # noqa: BLE001
            failed.append(f'{mt}/{v}: {type(e).__name__}: {e}')
            continue
        try:
            layout = build(mono_values(mono))
            written.append(write_vessel(mono, layout, units_library, bib))
            hide_from_phlynx(mono)
        except Exception as e:  # noqa: BLE001
            failed.append(f'{mt}/{v}: {type(e).__name__}: {e}')
    if not only:
        write_templates(units_library, bib)
        _add_system_equivalences()
    print(f'{len(written)} vessel supermodules written; {len(TEMPLATES)} templates')
    for f in failed:
        print('  FAILED', f)
    return written, failed
