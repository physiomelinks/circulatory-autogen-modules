'''
The lumped constitutive modules a vessel supermodule is built from (modules/haemodynamics/
lumped_constitutive/): compliance (C, the vessel's pressure-volume law P(A)), resistance (R),
inertance (I), and vessel_volume (the vessel's total volume and geometry, for the hosts and laws
that need the whole vessel).

Every compliance, resistance and inertance version has a dimensionless ``fraction``: the share of
the vessel it stands for. A vessel whose compliance is split in two halves (vv) uses two
compliances with fraction 0.5; one whose resistance and inertance are split (pp) uses two of each.
fraction scales C, q_0, R and I and divides the viscous terms, which is what the monolithic
versions' 2*R_v, 2*K_tube_visco, R/2 and I/2 are.

Port conventions (BC letters): a compliance is vv (it takes the flow on both sides and sets the
pressure; both its vessel_ports sum over their node). A resistance in pressure-drop form is pv
(upstream of an inertance: u_out = u_in - f R v, its v passed back upstream as v_in) or vp
(downstream of one: u_in = u_out + f R v); in flow form it is pp (v = (u_in - u_out)/(f R), for
vessels without inertance). An inertance is pp (it takes both pressures and gives the flow).
Resistances and inertances never sum: the compliance next to them owns the node.
'''
from dataclasses import dataclass, field

from . import mathml as m

SUM_TRUE = ['sum', 'True']


@dataclass
class Variable:
    name: str
    units: str
    kind: str                  # constant | global_constant | boundary_condition | variable
    initial: str = None        # a state's initial value (a number or a constant's name)
    default: object = None     # the default instance's value (constants and boundary conditions)


@dataclass
class Version:
    module_type: str
    version: str
    description: str
    variables: list
    equations: list
    entrance_ports: list
    exit_ports: list
    outputs: list = field(default_factory=list)
    sim_time: float = 1.0
    dt: float = 0.001

    @property
    def component(self):
        return f'{self.module_type}_{self.version}'

    def var(self, name):
        return next(v for v in self.variables if v.name == name)


def port(port_type, variables, multi_port=None):
    p = {'port_type': port_type, 'variables': list(variables)}
    if multi_port is not None:
        p['multi_port'] = multi_port
    return p


def V(name, units, kind, default=None, initial=None):
    return Variable(name, units, kind, initial, default)


def _geometry_globals():
    return [V('a_vessel', 'dimensionless', 'global_constant', 0.2802), V('b_vessel', 'per_m', 'global_constant', -505.3),
            V('c_vessel', 'dimensionless', 'global_constant', 0.1324), V('d_vessel', 'per_m', 'global_constant', -11.14)]


def _fraction():
    return V('fraction', 'dimensionless', 'constant', 1)


PI = '<pi/>'


def _wall(r):
    '''h = r (a e^{b r} + c e^{d r}), the Avolio wall-thickness fit.'''
    return m.times(m.ci(r), m.plus(m.times(m.ci('a_vessel'), m.apply('exp', m.times(m.ci('b_vessel'), m.ci(r)))),
                                   m.times(m.ci('c_vessel'), m.apply('exp', m.times(m.ci('d_vessel'), m.ci(r))))))


# --------------------------------------------------------------------------------------------
# compliance: P(A)
# --------------------------------------------------------------------------------------------

def _compliance_ports(extra_entrance=(), extra_exit=()):
    entrance = [port('vessel_port', ['v_in', 'u'], SUM_TRUE), port('external_pressure_port', ['u_ext'])]
    exit_ = [port('vessel_port', ['v_out', 'u'], SUM_TRUE), port('volume_port', ['q']),
             port('pressure_and_deriv_port', ['u', 'du_C_dt']), port('pressure_feedback_port', ['u'])]
    return entrance + list(extra_entrance), exit_ + list(extra_exit)


def _compliance_common():
    return [_fraction(), V('u_0', 'J_per_m3', 'constant', 10600), V('q_C_init', 'm3', 'constant', 0),
            V('u_ext', 'J_per_m3', 'boundary_condition', 0), V('v_in', 'm3_per_s', 'boundary_condition', 1e-6),
            V('v_out', 'm3_per_s', 'boundary_condition', 1e-6)]


def _linear_compliance(version, description, geometric, visco):
    variables = _compliance_common() + [
        V('C_eff', 'm6_per_J', 'variable'), V('q_C', 'm3', 'variable', initial='q_C_init'),
        V('u_C', 'J_per_m3', 'variable'), V('u', 'J_per_m3', 'variable'), V('q', 'm3', 'variable'),
        V('du_C_dt', 'J_per_m3s', 'variable')]
    equations = []
    if geometric:
        variables = [V('r_0', 'metre', 'constant', 0.01), V('l', 'metre', 'constant', 0.1),
                     V('E', 'J_per_m3', 'constant', 1.6e6), V('h', 'metre', 'variable'), V('C', 'm6_per_J', 'variable'),
                     V('q_0', 'm3', 'variable')] + _geometry_globals() + variables
        equations += [m.eq(m.ci('h'), _wall('r_0')),
                      m.eq(m.ci('C'), m.divide(m.times(m.cn(2), PI, m.ci('l'), m.power(m.ci('r_0'), m.cn(3))),
                                               m.times(m.ci('E'), m.ci('h')))),
                      m.eq(m.ci('q_0'), m.times(PI, m.power(m.ci('r_0'), m.cn(2)), m.ci('l')))]
    else:
        variables = [V('C', 'm6_per_J', 'constant', 1e-8), V('q_0', 'm3', 'constant', 1e-5)] + variables
    net = m.minus(m.ci('v_in'), m.ci('v_out'))
    equations += [m.eq(m.ci('C_eff'), m.times(m.ci('fraction'), m.ci('C'))),
                  m.ode('q_C', net),
                  m.eq(m.ci('u_C'), m.divide(m.ci('q_C'), m.ci('C_eff'))),
                  m.eq(m.ci('q'), m.plus(m.times(m.ci('fraction'), m.ci('q_0')), m.ci('q_C'))),
                  m.eq(m.ci('du_C_dt'), m.divide(net, m.ci('C_eff')))]
    terms = [m.ci('u_0'), m.ci('u_C')]
    if visco:
        variables.append(V('R_v', 'Js_per_m6', 'variable'))
        equations.append(m.eq(m.ci('R_v'), m.divide(m.cn('0.01', 'second'), m.ci('C_eff'))))
        terms.append(m.times(m.ci('R_v'), net))
    equations.append(m.eq(m.ci('u'), m.plus(*terms, m.ci('u_ext'))))
    entrance, exit_ = _compliance_ports()
    return Version('compliance', version, description, variables, equations, entrance, exit_,
                   outputs=['u', 'q', 'q_C', 'u_C'])


def _nonlinear_compliance(version, description, geometric, visco):
    variables = _compliance_common() + [
        V('m_tube', 'dimensionless', 'constant', 0.5), V('n_tube', 'dimensionless', 'constant', 0),
        V('K_tube', 'J_per_m3', 'boundary_condition', 1e5),
        V('q_0f', 'm3', 'variable'), V('q_C', 'm3', 'variable', initial='q_C_init'), V('q', 'm3', 'variable'),
        V('u', 'J_per_m3', 'variable'), V('du_C_dt', 'J_per_m3s', 'variable')]
    equations = []
    extra_exit = [port('stiffness_port', ['K_tube'])]
    if geometric:
        variables = [V('r_0', 'metre', 'constant', 0.01), V('l', 'metre', 'constant', 0.1),
                     V('q_0', 'm3', 'variable')] + variables
        equations.append(m.eq(m.ci('q_0'), m.times(PI, m.power(m.ci('r_0'), m.cn(2)), m.ci('l'))))
    else:
        variables = [V('q_0', 'm3', 'constant', 1e-5)] + variables
    net = m.minus(m.ci('v_in'), m.ci('v_out'))
    ratio = m.divide(m.ci('q'), m.ci('q_0f'))
    equations += [m.eq(m.ci('q_0f'), m.times(m.ci('fraction'), m.ci('q_0'))),
                  m.ode('q_C', net),
                  m.eq(m.ci('q'), m.plus(m.ci('q_0f'), m.ci('q_C'))),
                  m.eq(m.ci('du_C_dt'), m.times(m.ci('K_tube'), net, m.minus(
                      m.divide(m.times(m.ci('m_tube'), m.power(m.ci('q'), m.minus(m.ci('m_tube'), m.cn(1)))),
                               m.power(m.ci('q_0f'), m.ci('m_tube'))),
                      m.divide(m.times(m.ci('n_tube'), m.power(m.ci('q'), m.minus(m.ci('n_tube'), m.cn(1)))),
                               m.power(m.ci('q_0f'), m.ci('n_tube'))))))]
    terms = [m.times(m.ci('K_tube'), m.minus(m.power(ratio, m.ci('m_tube')), m.power(ratio, m.ci('n_tube'))))]
    if visco:
        variables.append(V('K_tube_visco', 'Js_per_m6', 'boundary_condition', 1e6))
        extra_exit.append(port('visco_coeff_port', ['K_tube_visco']))
        terms.append(m.times(m.divide(m.ci('K_tube_visco'), m.ci('fraction')), net))
    equations.append(m.eq(m.ci('u'), m.plus(*terms, m.ci('u_0'), m.ci('u_ext'))))
    entrance, exit_ = _compliance_ports(extra_exit=extra_exit)
    return Version('compliance', version, description, variables, equations, entrance, exit_,
                   outputs=['u', 'q', 'q_C'])


def _controlled_compliance(version, description, visco, conserve_total):
    '''A linear compliance whose unstressed volume and compliance a controller scales:
    q_us = f q_us_0 (1 - Delta_q_us), C_wCont = f C (1 - Delta_C); the viscous term keeps f C.
    conserve_total: the state is the total volume less the base unstressed volume (venous: a change
    of Delta_q_us moves volume between stressed and unstressed); otherwise the state is the stressed
    volume and the unstressed volume adds to the total.'''
    variables = [V('C', 'm6_per_J', 'constant', 1e-8), V('q_us_0', 'm3', 'constant', 1e-5)] + _compliance_common() + [
        V('Delta_q_us', 'dimensionless', 'boundary_condition', 0), V('Delta_C', 'dimensionless', 'boundary_condition', 0),
        V('C_eff', 'm6_per_J', 'variable'), V('C_wCont', 'm6_per_J', 'variable'), V('q_us_wCont', 'm3', 'variable'),
        V('q_C', 'm3', 'variable'), V('u_C', 'J_per_m3', 'variable'), V('u', 'J_per_m3', 'variable'),
        V('q', 'm3', 'variable'), V('du_C_dt', 'J_per_m3s', 'variable')]
    net = m.minus(m.ci('v_in'), m.ci('v_out'))
    equations = [m.eq(m.ci('C_eff'), m.times(m.ci('fraction'), m.ci('C'))),
                 m.eq(m.ci('C_wCont'), m.times(m.ci('C_eff'), m.minus(m.cn(1), m.ci('Delta_C')))),
                 m.eq(m.ci('q_us_wCont'), m.times(m.ci('fraction'), m.ci('q_us_0'), m.minus(m.cn(1), m.ci('Delta_q_us'))))]
    if conserve_total:
        variables.append(V('q_C_change', 'm3', 'variable', initial='q_C_init'))
        equations += [m.ode('q_C_change', net),
                      m.eq(m.ci('q_C'), m.minus(m.plus(m.ci('q_C_change'), m.times(m.ci('fraction'), m.ci('q_us_0'))),
                                                m.ci('q_us_wCont')))]
    else:
        next(v for v in variables if v.name == 'q_C').initial = 'q_C_init'
        equations.append(m.ode('q_C', net))
    equations += [m.eq(m.ci('u_C'), m.divide(m.ci('q_C'), m.ci('C_wCont'))),
                  m.eq(m.ci('q'), m.plus(m.ci('q_C'), m.ci('q_us_wCont'))),
                  m.eq(m.ci('du_C_dt'), m.divide(net, m.ci('C_wCont')))]
    terms = [m.ci('u_0'), m.ci('u_C')]
    if visco:
        variables.append(V('R_v', 'Js_per_m6', 'variable'))
        equations.append(m.eq(m.ci('R_v'), m.divide(m.cn('0.01', 'second'), m.ci('C_eff'))))
        terms.append(m.times(m.ci('R_v'), net))
    equations.append(m.eq(m.ci('u'), m.plus(*terms, m.ci('u_ext'))))
    entrance, exit_ = _compliance_ports(extra_entrance=[
        port('venous_us_volume_delta_port', ['Delta_q_us']),
        port('venous_us_and_compliance_delta_port', ['Delta_q_us', 'Delta_C'])])
    return Version('compliance', version, description, variables, equations, entrance, exit_,
                   outputs=['u', 'q', 'q_C'])


def _volume_flux_compliance(version, description):
    '''A linear compliance with one more flux into it: blood_ctl_dir * v_blood_intake.'''
    v = _linear_compliance(version, description, geometric=False, visco=True)
    v.variables += [V('v_blood_intake', 'm3_per_s', 'boundary_condition', 0),
                    V('blood_ctl_dir', 'dimensionless', 'constant', 1)]
    net = m.plus(m.minus(m.ci('v_in'), m.ci('v_out')), m.times(m.ci('blood_ctl_dir'), m.ci('v_blood_intake')))
    v.equations = [m.ode('q_C', net) if e.startswith('<apply><eq/><apply><diff/>') else e for e in v.equations]
    # the organ volume controllers connect through either name
    v.entrance_ports += [port('blood_intake_port', ['v_blood_intake']), port('blood_uptake_port', ['v_blood_intake'])]
    return v


COMPLIANCE = [
    _linear_compliance('vv_linear', 'Linear compliance: u = u_0 + q_C/(f C) + R_v (v_in - v_out) + u_ext, with the '
                       'viscous resistance R_v = 0.01 s/(f C).', geometric=False, visco=True),
    _linear_compliance('vv_linear_novisco', 'Linear compliance without the viscous term: u = u_0 + q_C/(f C) + u_ext.',
                       geometric=False, visco=False),
    _linear_compliance('vv_linear_geometric', 'Linear compliance from the vessel geometry: C = 2 pi l r_0^3/(E h) '
                       '(h the Avolio wall-thickness fit), q_0 = pi l r_0^2, with the viscous term R_v = 0.01 s/(f C).',
                       geometric=True, visco=True),
    _linear_compliance('vv_linear_geometric_novisco', 'Linear compliance from the vessel geometry, without the viscous '
                       'term.', geometric=True, visco=False),
    _nonlinear_compliance('vv_nonlinear', 'Nonlinear tube law u = K_tube ((q/(f q_0))^m - (q/(f q_0))^n) + u_0 + u_ext, '
                          'q_0 a parameter.', geometric=False, visco=False),
    _nonlinear_compliance('vv_nonlinear_geometric', 'Nonlinear tube law with q_0 = pi l r_0^2 from the geometry; '
                          'K_tube from a material-property module.', geometric=True, visco=False),
    _nonlinear_compliance('vv_nonlinear_geometric_visco', 'Nonlinear tube law from the geometry, with the '
                          'viscoelastic term (K_tube_visco/f)(v_in - v_out).', geometric=True, visco=True),
    _controlled_compliance('vv_linear_controlled', 'Linear compliance whose unstressed volume and compliance a '
                           'controller scales (Delta_q_us, Delta_C); the total volume is the state, so a change in '
                           'Delta_q_us moves volume between stressed and unstressed.', visco=True, conserve_total=True),
    _controlled_compliance('vv_linear_controlled_novisco', 'vv_linear_controlled without the viscous term.',
                           visco=False, conserve_total=True),
    _controlled_compliance('vv_linear_controlled_stressed', 'Linear compliance whose unstressed volume and compliance '
                           'a controller scales (Delta_q_us, Delta_C); the stressed volume is the state, so a change '
                           'in Delta_q_us changes the total volume.', visco=True, conserve_total=False),
    _volume_flux_compliance('vv_linear_volume_flux', 'Linear compliance with a controlled extra flux '
                            'blood_ctl_dir * v_blood_intake.'),
]


# --------------------------------------------------------------------------------------------
# resistance
# --------------------------------------------------------------------------------------------

_LAWS = {
    'linear': 'R a parameter',
    'geometric': 'Poiseuille R = 8 mu l/(pi r_0^4) R_flag',
    'nonlinear': 'Poiseuille R = 8 mu l/(pi r^4) R_flag, r = sqrt(q_total/(pi l)) from the whole vessel\'s volume '
                 '(total_volume_port, from vessel_volume)',
    'volume_scaled': 'R (q_ref/q_total)^2, q_total the whole vessel\'s volume (total_volume_port, from vessel_volume)',
    'controlled': 'R R_local_multiplier (1 + Delta_R), Delta_R from a controller',
    'controlled_share': 'R_local_multiplier (R + Delta_R R/(R + R_other)): a resistance change Delta_R shared out '
                        'over this resistance and R_other in proportion',
    'controlled_additive': 'R + Delta_R, Delta_R from a controller',
}
_FORMS = {'pv': 'pressure-drop form, upstream of an inertance: u_out = u_in - f R v',
          'vp': 'pressure-drop form, downstream of an inertance: u_in = u_out + f R v',
          'pp': 'flow form, for a vessel without inertance: v = (u_in - u_out)/(f R)'}


def _resistance(form, law):
    variables = [_fraction(), V('R_eff', 'Js_per_m6', 'variable')]
    equations = []
    entrance_extra = []
    if law in ('linear', 'volume_scaled', 'controlled', 'controlled_share', 'controlled_additive'):
        variables.append(V('R', 'Js_per_m6', 'constant', 1e7))
    if law == 'linear':
        resistance = m.ci('R')
    elif law == 'volume_scaled':
        variables += [V('q_ref', 'm3', 'constant', 1e-5), V('q_total', 'm3', 'boundary_condition', 1e-5),
                      V('R_nonlinear', 'Js_per_m6', 'variable')]
        equations.append(m.eq(m.ci('R_nonlinear'), m.times(m.ci('R'), m.power(m.divide(m.ci('q_ref'), m.ci('q_total')),
                                                                                m.cn(2)))))
        entrance_extra.append(port('total_volume_port', ['q_total']))
        resistance = m.ci('R_nonlinear')
    elif law.startswith('controlled'):
        variables += [V('Delta_R', 'dimensionless' if law == 'controlled' else 'Js_per_m6', 'boundary_condition', 0),
                      V('R_wCont', 'Js_per_m6', 'variable')]
        if law == 'controlled':
            variables.append(V('R_local_multiplier', 'dimensionless', 'boundary_condition', 1))
            rhs = m.times(m.ci('R'), m.ci('R_local_multiplier'), m.plus(m.ci('Delta_R'), m.cn(1)))
        elif law == 'controlled_share':
            variables += [V('R_other', 'Js_per_m6', 'constant', 1e7),
                          V('R_local_multiplier', 'dimensionless', 'boundary_condition', 1)]
            rhs = m.times(m.ci('R_local_multiplier'), m.plus(m.ci('R'), m.divide(
                m.times(m.ci('Delta_R'), m.ci('R')), m.plus(m.ci('R'), m.ci('R_other')))))
        else:
            rhs = m.plus(m.ci('R'), m.ci('Delta_R'))
        equations.append(m.eq(m.ci('R_wCont'), rhs))
        entrance_extra.append(port('resistance_delta_port', ['Delta_R']))
        if law != 'controlled_additive':
            entrance_extra.append(port('resistance_multiplier_port', ['R_local_multiplier']))
        resistance = m.ci('R_wCont')
    else:
        variables += [V('l', 'metre', 'constant', 0.1), V('mu', 'Js_per_m3', 'global_constant', 0.004),
                      V('R_flag', 'dimensionless', 'constant', 1), V('R', 'Js_per_m6', 'variable')]
        if law == 'geometric':
            variables.append(V('r_0', 'metre', 'constant', 0.01))
            radius = 'r_0'
        else:
            variables += [V('q_total', 'm3', 'boundary_condition', 3.14e-5), V('r', 'metre', 'variable')]
            equations.append(m.eq(m.ci('r'), m.apply('root', m.divide(m.ci('q_total'), m.times(PI, m.ci('l'))))))
            entrance_extra.append(port('total_volume_port', ['q_total']))
            radius = 'r'
        equations.append(m.eq(m.ci('R'), m.divide(m.times(m.cn(8), m.ci('mu'), m.ci('l')),
                                                  m.times(PI, m.power(m.ci(radius), m.cn(4))))))
        resistance = m.times(m.ci('R'), m.ci('R_flag'))
    equations.append(m.eq(m.ci('R_eff'), m.times(m.ci('fraction'), resistance)))
    if form == 'pv':
        variables += [V('u_in', 'J_per_m3', 'boundary_condition', 11000), V('v', 'm3_per_s', 'boundary_condition', 1e-6),
                      V('v_in', 'm3_per_s', 'variable'), V('u_out', 'J_per_m3', 'variable')]
        equations += [m.eq(m.ci('v_in'), m.ci('v')),
                      m.eq(m.ci('u_out'), m.minus(m.ci('u_in'), m.times(m.ci('R_eff'), m.ci('v'))))]
        entrance = [port('vessel_port', ['v_in', 'u_in'])]
        exit_ = [port('vessel_port', ['v', 'u_out'])]
        outputs = ['u_out']
    elif form == 'vp':
        variables += [V('u_out', 'J_per_m3', 'boundary_condition', 10000), V('v', 'm3_per_s', 'boundary_condition', 1e-6),
                      V('v_out', 'm3_per_s', 'variable'), V('u_in', 'J_per_m3', 'variable')]
        equations += [m.eq(m.ci('v_out'), m.ci('v')),
                      m.eq(m.ci('u_in'), m.plus(m.ci('u_out'), m.times(m.ci('R_eff'), m.ci('v'))))]
        entrance = [port('vessel_port', ['v', 'u_in'])]
        exit_ = [port('vessel_port', ['v_out', 'u_out'])]
        outputs = ['u_in']
    else:
        variables += [V('u_in', 'J_per_m3', 'boundary_condition', 11000),
                      V('u_out', 'J_per_m3', 'boundary_condition', 10000), V('v', 'm3_per_s', 'variable')]
        equations.append(m.eq(m.ci('v'), m.divide(m.minus(m.ci('u_in'), m.ci('u_out')), m.ci('R_eff'))))
        entrance = [port('vessel_port', ['v', 'u_in'])]
        exit_ = [port('vessel_port', ['v', 'u_out']), port('flow_port', ['v']), port('flow_feedback_port', ['v']),
                 port('pressure_port', ['u_out'])]
        outputs = ['v']
    return Version('resistance', f'{form}_{law}', f'Resistance in {_FORMS[form]}; {_LAWS[law]}.', variables,
                   equations, entrance + entrance_extra, exit_, outputs=outputs)


RESISTANCE = [_resistance(form, law) for form in ('pv', 'vp', 'pp') for law in _LAWS]


# --------------------------------------------------------------------------------------------
# inertance
# --------------------------------------------------------------------------------------------

def _inertance(version, description, law):
    variables = [_fraction(), V('u_in', 'J_per_m3', 'boundary_condition', 11000),
                 V('u_out', 'J_per_m3', 'boundary_condition', 10000), V('v_init', 'm3_per_s', 'constant', 0),
                 V('v', 'm3_per_s', 'variable', initial='v_init'), V('I_eff', 'Js2_per_m6', 'variable')]
    equations = []
    drive = m.minus(m.ci('u_in'), m.ci('u_out'))
    if law == 'linear':
        variables.append(V('I', 'Js2_per_m6', 'constant', 1e5))
    else:
        variables += [V('r_0', 'metre', 'constant', 0.01), V('l', 'metre', 'constant', 0.1),
                      V('I_scale', 'dimensionless', 'constant', 1), V('theta', 'dimensionless', 'constant', 0),
                      V('gravity_factor', 'dimensionless', 'constant', 1),
                      V('rho', 'Js2_per_m5', 'global_constant', 1050), V('g', 'm_per_s2', 'global_constant', 9.81),
                      V('beta_g', 'dimensionless', 'global_constant', 1), V('I', 'Js2_per_m6', 'variable')]
        equations.append(m.eq(m.ci('I'), m.divide(m.times(m.ci('I_scale'), m.ci('l'), m.ci('rho')),
                                                  m.times(PI, m.power(m.ci('r_0'), m.cn(2))))))
        gravity = m.times(m.ci('gravity_factor'), m.ci('beta_g'), m.ci('g'), m.ci('l'), m.ci('rho'),
                          m.apply('cos', m.divide(m.times(PI, m.ci('theta')), m.cn(180))))
        drive = m.minus(drive, gravity)
    equations += [m.eq(m.ci('I_eff'), m.times(m.ci('fraction'), m.ci('I'))),
                  m.ode('v', m.divide(drive, m.ci('I_eff')))]
    entrance = [port('vessel_port', ['v', 'u_in'])]
    exit_ = [port('vessel_port', ['v', 'u_out']), port('flow_port', ['v']), port('flow_feedback_port', ['v']),
             port('pressure_port', ['u_out'])]
    return Version('inertance', version, description, variables, equations, entrance, exit_, outputs=['v'],
                   sim_time=0.1, dt=0.0001)


INERTANCE = [
    _inertance('pp_linear', 'Inertance: dv/dt = (u_in - u_out)/(f I), I a parameter.', 'linear'),
    _inertance('pp_geometric', 'Inertance from the geometry, I = I_scale rho l/(pi r_0^2), with the hydrostatic '
               'term gravity_factor beta_g g l rho cos(pi theta/180): dv/dt = (u_in - u_out - that)/(f I).', 'geometric'),
]


# --------------------------------------------------------------------------------------------
# vessel_volume: the whole vessel
# --------------------------------------------------------------------------------------------

def _vessel_volume(version, description, geometric):
    # q: the total as a variable of the vessel's own (q_total is the sum libcuflynx writes)
    variables = [V('q_total', 'm3', 'boundary_condition', 3.14e-5), V('q', 'm3', 'variable')]
    equations = [m.eq(m.ci('q'), m.ci('q_total'))]
    # multi_port True: one vessel's volume and geometry feed every resistance and host that reads them
    exit_ = [port('volume_port', ['q'], 'True'), port('total_volume_port', ['q'], 'True')]
    if geometric:
        variables += [V('r_0', 'metre', 'constant', 0.01), V('l', 'metre', 'constant', 0.1), V('r', 'metre', 'variable')]
        equations.append(m.eq(m.ci('r'), m.apply('root', m.divide(m.ci('q_total'), m.times(PI, m.ci('l'))))))
        exit_ += [port('radius_port', ['r_0'], 'True'), port('len_port', ['l'], 'True'),
                  port('radius_0_and_s_port', ['r_0', 'r'], 'True')]
    return Version('vessel_volume', version, description, variables, equations,
                   [port('volume_port', ['q_total'], ['sum'])], exit_, outputs=['q', 'r'] if geometric else ['q'])


VESSEL_VOLUME = [
    _vessel_volume('nn', 'The whole vessel\'s volume: the sum of its compliances\' volumes (volume_port in, summed), '
                   'passed to hosts (volume_port) and to resistance laws that need it (total_volume_port).', False),
    _vessel_volume('nn_geometric', 'The whole vessel\'s volume and geometry: q_total summed over its compliances, '
                   'r_0 and l, and the radius r = sqrt(q_total/(pi l)), for hosts (volume_port, radius_port, len_port, '
                   'radius_0_and_s_port: material-property modules) and resistance laws (total_volume_port).', True),
]

ALL = COMPLIANCE + RESISTANCE + INERTANCE + VESSEL_VOLUME
