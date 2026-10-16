"""Bond-graph descriptions of chemical modules: one description gives a module's equations, its BG-RDF (the CellDL
tools' bond-graph format, https://github.com/CellDL/BG-RDF) and a diagram (SVG carrying that RDF).

A description lists species (each a store, a C element on a 0-junction, with its element formula) and reactions (Re
elements between 1-junctions, with stoichiometry and a standard Gibbs energy of reaction at the module's pH). Species
standard potentials mu0 are found by least squares from the reactions' Gibbs energies (min-norm: they are defined only
up to the conserved moieties, which the rate constants absorb), so the reactions are thermodynamically consistent by
construction: around any cycle of reactions the Gibbs energies sum to zero (Wegscheider). The residual of that fit says
how far the textbook values were from consistent.

Kinetics: each reaction flows as v = kappa (Pf - Pr) / (1 + (Pf + Pr) / sigma), with Pf = prod_s (exp(mu0_s/RT) C_s/c_std_s)^n_s
over its substrates and Pr likewise over its products: exp(A/RT) of the forward and reverse affinities. With sigma
infinite (no saturation) this is the bond-graph chemical reaction, v = kappa (exp(Af/RT) - exp(Ar/RT)); a finite sigma
saturates it without moving its equilibrium, and allosteric activators multiply it by (1 + C_a/K_a), also without
moving it. Boundary species (chemostats) take their concentration from a port and
report their net production back on it. A reaction that moves charge q into a compartment with an electrical potential
psi (a carrier, a proton pump) does the electrical work q F psi, split evenly between its two directions: Pf x
exp(-q F psi / 2RT), Pr x exp(q F psi / 2RT). Its equilibrium is that of the species' electrochemical potentials
mu0 + RT ln(C/c_std) + z F psi, and a reaction inside one compartment does not feel the potential. A species with k_sat has
the activity a = C k_sat / (k_sat + C) in place of C (mu = mu0 + RT ln(a/c_std): a C element with a saturating
constitutive law), so a reaction that uses it stops as it runs out. With a capacitance (an electrical C element) the compartment's potential is
a state: C_m dpsi/dt = F sum_r q_r v_r, q_r the charge reaction r moves into it (sum of z n over its species there).

    from bondgraph import Species, Reaction, BondGraph
"""
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

R_GAS = 8.314462618        # J/(mol K)
FARADAY = 96485.33212      # C/mol

ELEMENTS = ('C', 'H', 'O', 'N', 'P', 'S')


@dataclass
class Species:
    name: str                       # variable stem: C_<name>
    label: str                      # e.g. 'glucose 6-phosphate'
    formula: Dict[str, int]         # element counts of the species as it is balanced (e.g. {'C': 6, 'H': 11, 'O': 9, 'P': 1})
    boundary: bool = False          # concentration from a port (a chemostat), not a state of this module
    charge: int = 0
    c_std: float = 1000.0           # standard concentration [mM]: 1 M; H+ at 1e-4 mM (pH 7) and water at its own, so
                                    # transformed standard Gibbs energies (pH 7, 1 M) apply with H+ and water explicit
    compartment: str = ''           # a compartment with an electrical potential (BondGraph.potentials) adds z F psi
    k_sat: Optional[float] = None   # [mM] a store whose activity saturates, a = C k_sat / (k_sat + C) (a polymer's chain
                                    # ends, say): ideal while C << k_sat; None: ideal, a = C


@dataclass
class Reaction:
    name: str                       # variable stem: v_<name>
    enzyme: str
    stoich: Dict[str, int]          # species -> coefficient: negative consumed, positive made
    dG0: float                      # standard Gibbs energy of reaction [J/mol] at the module's pH
    ref: str = ''                   # where dG0 comes from
    sigma: Optional[float] = None   # saturation (None: mass action)
    activators: Dict[str, float] = field(default_factory=dict)   # species -> K_a [mM]: rate x (1 + C/K_a), equilibrium kept
    inhibitors: Dict[str, float] = field(default_factory=dict)   # species -> K_i [mM]: rate x K_i/(K_i + C), equilibrium kept


@dataclass
class BondGraph:
    name: str
    label: str
    species: List[Species]
    reactions: List[Reaction]
    temperature: float = 293.15     # K
    gas_constant: float = R_GAS     # J/(mol K): a lumped module's own value, to describe it exactly
    potentials: Dict[str, float] = field(default_factory=dict)   # compartment -> electrical potential psi [V]
    capacitance: Dict[str, float] = field(default_factory=dict)  # compartment -> membrane capacitance [F/m3]: its psi a
                                                                 # state, charged by the reactions that move charge in
    mu0: Dict[str, float] = field(default_factory=dict)
    fit_residual: float = 0.0       # J/mol: largest change the consistency fit made to a reaction's dG0

    def __post_init__(self):
        names = [s.name for s in self.species]
        assert len(set(names)) == len(names), f'{self.name}: duplicate species'
        for r in self.reactions:
            for s in r.stoich:
                assert s in names, f'{self.name}: {r.name} uses unknown species {s}'
        self.fit_potentials()

    # ---------------------------------------------------------------- structure
    def S(self):
        """Stoichiometric matrix, species x reactions."""
        idx = {s.name: i for i, s in enumerate(self.species)}
        m = np.zeros((len(self.species), len(self.reactions)))
        for j, r in enumerate(self.reactions):
            for s, n in r.stoich.items():
                m[idx[s], j] = n
        return m

    def element_imbalance(self):
        """Per reaction, the elements not balanced (products minus substrates); empty when balanced."""
        out = {}
        f = {s.name: s.formula for s in self.species}
        for r in self.reactions:
            bad = {}
            for e in ELEMENTS:
                d = sum(n * f[s].get(e, 0) for s, n in r.stoich.items())
                if d:
                    bad[e] = d
            if bad:
                out[r.name] = bad
        return out

    def charge_imbalance(self):
        c = {s.name: s.charge for s in self.species}
        return {r.name: sum(n * c[s] for s, n in r.stoich.items()) for r in self.reactions
                if sum(n * c[s] for s, n in r.stoich.items())}

    def fit_potentials(self):
        """mu0 from the reactions' dG0 by least squares (min-norm): dG0_r = sum_s n_sr mu0_s."""
        S = self.S()
        g = np.array([r.dG0 for r in self.reactions])
        mu, *_ = np.linalg.lstsq(S.T, g, rcond=None)
        self.mu0 = {s.name: float(m) for s, m in zip(self.species, mu)}
        self.fit_residual = float(np.max(np.abs(S.T @ mu - g))) if len(g) else 0.0

    def centre_potentials(self, conc):
        """Moves mu0 along the conserved moieties (the null space of S^T, which changes no reaction's Gibbs energy) so the
        species' chemical potentials at conc are as near zero as they can be: Pf and Pr near 1 there, and the rate
        constants near the fluxes (well scaled numbers, in the module and its BG-RDF)."""
        rt = self.gas_constant * self.temperature
        S = self.S()
        mu = np.array([self.mu0[s.name] for s in self.species])
        y = np.array([rt * math.log(self.activity(s.name, conc[s.name]) / s.c_std) for s in self.species])
        _, sv, vh = np.linalg.svd(S.T)
        rank = int(np.sum(sv > 1e-9 * sv.max()))
        N = vh[rank:].T
        z, *_ = np.linalg.lstsq(N, -(mu + y), rcond=None)
        mu = mu + N @ z
        g = np.array([r.dG0 for r in self.reactions])
        assert np.max(np.abs(S.T @ mu - g)) < 1e-6 + self.fit_residual, 'centring moved a reaction\'s Gibbs energy'
        self.mu0 = {s.name: float(m) for s, m in zip(self.species, mu)}

    def dG0_fitted(self, r):
        return sum(n * self.mu0[s] for s, n in r.stoich.items())

    # ---------------------------------------------------------------- kinetics
    def affinity_terms(self, r, conc):
        """(Pf, Pr): exp(Af/RT) and exp(Ar/RT) at concentrations conc [mM]."""
        rt = self.gas_constant * self.temperature
        std = {sp.name: sp.c_std for sp in self.species}
        pf = pr = 1.0
        for s, n in r.stoich.items():
            x = (math.exp(self.mu0[s] / rt) * self.activity(s, conc[s]) / std[s]) ** abs(n)
            if n < 0:
                pf *= x
            else:
                pr *= x
        work = self.electrical_work(r)
        return pf * math.exp(-0.5 * work / rt), pr * math.exp(0.5 * work / rt)

    def electrical_work(self, r):
        """q F psi summed over the compartments with a potential [J/mol]: the electrical part of r's Gibbs energy."""
        return sum(self.charge_into(r, c) * FARADAY * psi for c, psi in self.potentials.items())

    def electrical(self, s):
        """z F psi of species s [J/mol] (0 outside a compartment with a potential)."""
        sp = self.by_name()[s]
        return sp.charge * FARADAY * self.potentials[sp.compartment] if sp.compartment in self.potentials else 0.0

    def activity(self, s, c):
        """The species' activity at concentration c [mM]: c, or c k_sat / (k_sat + c) for a saturating store."""
        k = self.by_name()[s].k_sat
        return c if k is None else c * k / (k + c)

    def by_name(self):
        return {sp.name: sp for sp in self.species}

    @staticmethod
    def activation(r, conc):
        """The allosteric factor on r's rate: activators and inhibitors (they change its rate, not its equilibrium)."""
        a = 1.0
        for s, k in r.activators.items():
            a *= 1.0 + conc[s] / k
        for s, k in r.inhibitors.items():
            a *= k / (k + conc[s])
        return a

    def rate(self, r, kappa, conc):
        pf, pr = self.affinity_terms(r, conc)
        sat = 1.0 + ((pf + pr) / r.sigma if r.sigma else 0.0)
        return kappa * self.activation(r, conc) * (pf - pr) / sat

    def calibrate(self, conc, flux, forward=None):
        """kappa per reaction so it carries flux[r] [mM/s] at conc; refuses a flux against the reaction's affinity. A
        reaction with no net flux (at equilibrium at conc) takes its kappa from its forward rate there, forward[r]
        [mM/s] (0 without one)."""
        kappa = {}
        for r in self.reactions:
            pf, pr = self.affinity_terms(r, conc)
            want = flux[r.name]
            if want == 0.0:
                fwd = (forward or {}).get(r.name, 0.0)
                sat = 1.0 + ((pf + pr) / r.sigma if r.sigma else 0.0)
                kappa[r.name] = fwd * sat / pf / self.activation(r, conc) if fwd else 0.0
                continue
            if (pf - pr) * want <= 0.0:
                gamma = -self.gas_constant * self.temperature * math.log(pf / pr)
                raise ValueError(f'{self.name}: {r.name} cannot carry {want:+g} mM/s at these concentrations '
                                 f'(Gibbs energy {gamma / 1e3:+.2f} kJ/mol)')
            sat = 1.0 + ((pf + pr) / r.sigma if r.sigma else 0.0)
            kappa[r.name] = want * sat / (pf - pr) / self.activation(r, conc)
        return kappa

    def near_equilibrium(self, r, conc, unknown, ratio):
        """The concentration of `unknown` (a species of r with coefficient +-1) at which r's forward rate is `ratio` times
        its net rate (Pr / Pf = 1 - 1/ratio: Gibbs energy RT ln(1 - 1/ratio)), the others as in conc."""
        n = r.stoich[unknown]
        assert abs(n) == 1, f'{r.name}: {unknown} has coefficient {n}'
        trial = dict(conc, **{unknown: 1.0})
        pf, pr = self.affinity_terms(r, trial)
        target = 1.0 - 1.0 / ratio   # Pr / Pf
        # Pr/Pf scales as x for a product (n = +1) and as 1/x for a substrate (n = -1)
        return target * pf / pr if n > 0 else pr / (target * pf)

    def solve_near_equilibrium(self, names, conc, unknowns, ratio):
        """Concentrations of `unknowns` at which each reaction in `names` runs forwards at `ratio` times its net rate (Pr/Pf
        = 1 - 1/ratio; a reaction with no net flux takes ratio infinite, at equilibrium: pass it in names with ratio
        None in `equilibrium`). Linear in the logarithms of the concentrations: one least-squares solve, which must be
        exact (as many independent conditions as unknowns)."""
        rt = self.gas_constant * self.temperature
        std = {sp.name: sp.c_std for sp in self.species}
        rx = {r.name: r for r in self.reactions}
        idx = {u: i for i, u in enumerate(unknowns)}
        A = np.zeros((len(names), len(unknowns)))
        b = np.zeros(len(names))
        for k, (name, rho) in enumerate(names.items() if isinstance(names, dict) else ((n, ratio) for n in names)):
            target = math.log(1.0 - 1.0 / rho) if rho else 0.0   # ln(Pr/Pf)
            rhs = target
            for sp, n in rx[name].stoich.items():
                term = (self.mu0[sp] + self.electrical(sp)) / rt - math.log(std[sp])
                if sp in idx:
                    A[k, idx[sp]] += n
                    rhs -= n * term
                else:
                    rhs -= n * (term + math.log(self.activity(sp, conc[sp])))
            b[k] = rhs
        assert all(self.by_name()[u].k_sat is None for u in unknowns), 'solve for ideal species only'
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        resid = float(np.max(np.abs(A @ x - b))) if len(b) else 0.0
        assert resid < 1e-9 and np.linalg.matrix_rank(A) == len(unknowns), \
            f'{self.name}: near-equilibrium conditions do not fix {unknowns} (rank {np.linalg.matrix_rank(A)}, residual {resid:g})'
        return dict(conc, **{u: float(math.exp(x[idx[u]])) for u in unknowns})

    def charge_into(self, r, comp):
        """The charge (in units of e per turnover) reaction r moves into compartment comp."""
        by = self.by_name()
        return sum(n * by[s].charge for s, n in r.stoich.items() if by[s].compartment == comp)

    def current(self, comp, flux):
        """Net charge into comp at these fluxes [mM/s of elementary charge]."""
        return sum(self.charge_into(r, comp) * flux[r.name] for r in self.reactions)

    def steady(self, flux):
        """Net production of each species at these fluxes [mM/s]."""
        S = self.S()
        v = np.array([flux[r.name] for r in self.reactions])
        return {s.name: float(x) for s, x in zip(self.species, S @ v)}

    def reaction_gibbs(self, r, conc):
        """The reaction's Gibbs energy at conc [J/mol] (negative: it runs forwards)."""
        pf, pr = self.affinity_terms(r, conc)
        return -self.gas_constant * self.temperature * math.log(pf / pr)

    # ---------------------------------------------------------------- CellML (module_writer's infix equations)
    def factor(self, s):
        if self.by_name()[s].k_sat is not None:
            return f'exp(mu0_{s} / (R * T)) * C_{s} * k_sat_{s} / (k_sat_{s} + C_{s}) / c_std_{s}'
        return f'exp(mu0_{s} / (R * T)) * C_{s} / c_std_{s}'

    def _work(self, r):
        """r's electrical work as infix (empty when it moves no charge across a membrane)."""
        terms = [f'{float(self.charge_into(r, c))} * F * psi_{c}' for c in self.potentials if self.charge_into(r, c)]
        return ' + '.join(terms)

    def equations(self):
        eqs = []
        for r in self.reactions:
            sub = ' * '.join(self._power(s, -n) for s, n in r.stoich.items() if n < 0) or '1.0'
            prod = ' * '.join(self._power(s, n) for s, n in r.stoich.items() if n > 0) or '1.0'
            work = self._work(r)
            if work:
                sub = f'{sub} * exp(-0.5 * ({work}) / (R * T))'
                prod = f'{prod} * exp(0.5 * ({work}) / (R * T))'
            eqs.append((f'Pf_{r.name}', sub))
            eqs.append((f'Pr_{r.name}', prod))
            sat = f' / (1.0 + (Pf_{r.name} + Pr_{r.name}) / sigma_{r.name})' if r.sigma else ''
            act = ''.join(f' * (1.0 + C_{s} / Ka_{s}_{r.name})' for s in r.activators) + \
                ''.join(f' * Ki_{s}_{r.name} / (Ki_{s}_{r.name} + C_{s})' for s in r.inhibitors)
            eqs.append((f'v_{r.name}', f'kappa_{r.name}{act} * (Pf_{r.name} - Pr_{r.name}){sat}'))
            eqs.append((f'A_{r.name}', f'R * T * (ln(Pf_{r.name}) - ln(Pr_{r.name}))'))
        for s in self.species:
            net = self.net(s.name)
            if s.boundary:
                eqs.append((f'j_{s.name}', net))   # its net production, reported on the port
            else:
                eqs.append((f'd/dt C_{s.name}', net))
        for comp in self.capacitance:
            terms = [(self.charge_into(r, comp), r.name) for r in self.reactions]
            cur = ' + '.join(f'{float(q)} * v_{name}' for q, name in terms if q)
            eqs.append((f'I_{comp}', cur or '0.0'))
            eqs.append((f'd/dt psi_{comp}', f'F * I_{comp} / C_m_{comp}'))
        eqs.append(('P_diss', ' + '.join(f'v_{r.name} * A_{r.name}' for r in self.reactions)))
        return eqs

    def net(self, s):
        """The species' net production as infix: a signed sum of its reactions' rates."""
        out = ''
        for r in self.reactions:
            n = r.stoich.get(s, 0)
            if not n:
                continue
            term = f'v_{r.name}' if abs(n) == 1 else f'{float(abs(n))} * v_{r.name}'
            out += (' - ' if n < 0 else ' + ') + term if out else ('-' if n < 0 else '') + term
        return out or '0.0'

    def _power(self, s, n):
        f = f'({self.factor(s)})'
        return f if n == 1 else f'{f} ** {float(n)}'

    # ---------------------------------------------------------------- BG-RDF (Turtle)
    def bg_rdf(self, conc, kappa, volume):
        """The bond graph in BG-RDF (bgf: framework, chemical template), in its reduced form (as the CellDL editor writes
        it): each species a ChemicalStorageNode (a store that is its own 0-junction), each reaction a
        ChemicalReactionNode (with its 1-junctions), bonds from substrates to the reaction and from it to products, a
        stoichiometric coefficient n > 1 as the bond's bgf:bondCount (the reaction sees n mu, the store n v); for a compartment of volume [m3] at concentrations conc, stores
        in mol with K [1/mol] and reactions with kappa [mol/s]. The chemical template is mass action: a reaction's
        saturation (sigma) and activators are folded into its kappa at conc (the same flow there; mass action away from
        it), and noted. Chemostats are exported as stores (the template has no held potential for them yet), noted."""
        rt = self.gas_constant * self.temperature
        std = {s.name: s.c_std for s in self.species}
        out = ['@prefix bgf: <https://bg-rdf.org/ontologies/bondgraph-framework#> .',
               '@prefix cdt: <https://w3id.org/cdt/> .', '@prefix tpl: <https://bg-rdf.org/templates/> .', '',
               '@prefix : <#> .', '',
               f'# {self.label}: written by tools/bondgraph.py from the module\'s description (circulatory-autogen-modules)',
               f'# T = {self.temperature} K; volume {volume:g} m3; stores q = C V, K = exp(mu0/RT) / (c_std V)',
               *([f'# electrical potentials folded into K: ' + ', '.join(f'{c} {v * 1e3:g} mV' for c, v in self.potentials.items())]
                 if self.potentials else []), '']
        elems = [f':q_{s.name}' for s in self.species] + [f':Re_{r.name}' for r in self.reactions]
        out += [f':{self.name}', '    a bgf:BondgraphModel ;', '    bgf:usesTemplate tpl:chemical.ttl ;',
                '    bgf:hasBondElement ' + ', '.join(elems) + ' .', '']
        for s in self.species:
            # mM = mol/m3; a saturating store is exported linear with its activity coefficient a/C at the export point
            k = math.exp((self.mu0[s.name] + self.electrical(s.name)) / rt) / (std[s.name] * volume)
            if s.k_sat is not None:
                k *= self.activity(s.name, conc[s.name]) / conc[s.name]
            note = ('  # a chemostat of the module (held by its port), here a store' if s.boundary else
                    f'  # activity C k_sat / (k_sat + C), k_sat {s.k_sat:g} mM: linear here, at the export point'
                    if s.k_sat is not None else '')
            out += [f':q_{s.name}{note}', '    a bgf:ChemicalStorageNode ;',
                    f'    bgf:hasValue "{conc[s.name] * volume:.15e} mol"^^cdt:ucum ;',
                    f'    bgf:parameterValue [ bgf:varName "K" ; bgf:hasValue "{k:.15e} mol-1"^^cdt:ucum ] .', '']
        for r in self.reactions:
            # the stores' K carry z F psi; kappa takes back the substrates' share and the even split of q F psi
            sub_el = sum(-n * self.electrical(sp) for sp, n in r.stoich.items() if n < 0)
            factor, note = math.exp(-(0.5 * self.electrical_work(r) + sub_el) / rt), ''
            if r.sigma or r.activators or r.inhibitors:
                pf, pr = self.affinity_terms(r, conc)
                factor *= self.activation(r, conc) / (1.0 + ((pf + pr) / r.sigma if r.sigma else 0.0))
                parts = ([f'saturation sigma {r.sigma:.4g}'] if r.sigma else []) + \
                    [f'activation by {a} (K_a {k:g} mM)' for a, k in r.activators.items()] + \
                    [f'inhibition by {a} (K_i {k:g} mM)' for a, k in r.inhibitors.items()]
                note = '  # ' + ' and '.join(parts) + ': folded into kappa at the export\'s concentrations'
            out += [f':Re_{r.name}{note}', '    a bgf:ChemicalReactionNode ;',
                    f'    bgf:parameterValue [ bgf:varName "kappa" ; bgf:hasSymbol "kappa_{r.name}" ; '
                    f'bgf:hasValue "{kappa[r.name] * factor * volume:.15e} mol/s"^^cdt:ucum ] .', '']
        for r in self.reactions:
            for sp, n in r.stoich.items():
                a, b = (f'q_{sp}', f'Re_{r.name}') if n < 0 else (f'Re_{r.name}', f'q_{sp}')
                count = f' ; bgf:bondCount {abs(n)}' if abs(n) > 1 else ''
                out.append(f':{a}.{b} bgf:hasSource :{a} ; bgf:hasTarget :{b}{count} .')
        return '\n'.join(out) + '\n'

    # ---------------------------------------------------------------- diagram
    def svg(self, rdf):
        """A diagram: reactions down the middle in their order; a species that links two consecutive reactions sits
        between them, one used by a single reaction beside it, and the rest (cofactors, species joining distant steps)
        in columns at the sides, each species once (a 0-junction with its store); chemostats dashed. Bonds point from
        substrates to reactions to products, with the coefficient when it is not 1. The BG-RDF is in the SVG's metadata
        where CellDL keeps it (id celldl-rdf-metadata, Turtle), so the CellDL tools read the diagram as a bond-graph
        model (bg2cellml makes CellML from it), and each element's id is its bond-graph name."""
        rr = self.reactions
        dy, xr, bh, rad = 110.0, 0.0, 36.0, 21.0
        bw = max(170.0, max(5.6 * len(r.enzyme) + 20.0 for r in rr))   # wide enough for the enzyme names (10 px)
        order = {r.name: i for i, r in enumerate(rr)}
        ypos = {r.name: 60.0 + i * dy for i, r in enumerate(rr)}
        uses = {sp.name: sorted((r.name for r in rr if sp.name in r.stoich), key=order.get) for sp in self.species}
        place, sides, hubs, gaps = {}, {}, [], {}
        for sp in self.species:
            u = uses[sp.name]
            if len(u) == 2 and order[u[1]] - order[u[0]] == 1:
                k = gaps.get(u[0], 0)   # a second species between the same two steps goes beside the first
                gaps[u[0]] = k + 1
                place[sp.name] = (xr + (0.0 if k == 0 else (-1 if k % 2 else 1) * 2.4 * rad * ((k + 1) // 2)),
                                  0.5 * (ypos[u[0]] + ypos[u[1]]))
            elif len(u) == 1:
                k = sides.get(u[0], 0)
                sides[u[0]] = k + 1
                place[sp.name] = (xr + (-1 if k % 2 == 0 else 1) * (bw / 2 + 70.0 + 55.0 * (k // 2)), ypos[u[0]])
            else:
                hubs.append(sp.name)
        # side columns: by mean height, alternating sides; a hub too close to another on its side moves out a column
        hubs.sort(key=lambda n: sum(ypos[r] for r in uses[n]) / len(uses[n]))
        taken = {-1: [], 1: []}
        for i, n in enumerate(hubs):
            side = -1 if i % 2 == 0 else 1
            y = sum(ypos[r] for r in uses[n]) / len(uses[n])
            col = 0
            while any(c == col and abs(y - yy) < 2.5 * rad for c, yy in taken[side]):
                col += 1
            taken[side].append((col, y))
            place[n] = (xr + side * (bw / 2 + 220.0 + 80.0 * col), y)
        xs = [x for x, _ in place.values()] + [xr - bw / 2, xr + bw / 2]
        x0, x1 = min(xs) - 60.0, max(xs) + 60.0
        height = max(ypos.values()) + 80.0
        bnd = {sp.name for sp in self.species if sp.boundary}
        lab = {sp.name: sp.label for sp in self.species}

        def edge_of_box(cx, cy, tx, ty):
            dx, dy_ = tx - cx, ty - cy
            if dx == 0 and dy_ == 0:
                return cx, cy
            t = min((bw / 2) / abs(dx) if dx else float('inf'), (bh / 2) / abs(dy_) if dy_ else float('inf'))
            return cx + dx * t, cy + dy_ * t

        def edge_of_circle(cx, cy, tx, ty):
            d = math.hypot(tx - cx, ty - cy) or 1.0
            return cx + (tx - cx) * rad / d, cy + (ty - cy) * rad / d

        out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0:.0f} 0 {x1 - x0:.0f} {height:.0f}" '
               f'width="{x1 - x0:.0f}" height="{height:.0f}" font-family="sans-serif" font-size="11">',
               f'<title>{_x(self.label)}: bond graph</title>',
               '<metadata id="celldl-rdf-metadata" data-content-type="text/turtle"><![CDATA[\n' + rdf + ']]></metadata>',
               '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
               'orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#666"/></marker></defs>']
        for r in rr:
            y = ypos[r.name]
            for sp, n in r.stoich.items():
                sx, sy = place[sp]
                a = edge_of_circle(sx, sy, xr, y)
                b = edge_of_box(xr, y, sx, sy)
                (x1_, y1_), (x2_, y2_) = (a, b) if n < 0 else (b, a)
                out.append(f'<line class="bond" id="bond_{sp}_{r.name}" x1="{x1_:.1f}" y1="{y1_:.1f}" x2="{x2_:.1f}" '
                           f'y2="{y2_:.1f}" stroke="#999" stroke-width="1.2" marker-end="url(#arrow)"/>')
                if abs(n) > 1:
                    out.append(f'<text x="{(x1_ + x2_) / 2 + 4:.1f}" y="{(y1_ + y2_) / 2 - 4:.1f}" fill="#444" '
                               f'font-weight="bold">{abs(n)}</text>')
        for r in rr:
            y = ypos[r.name]
            out.append(f'<g id="Re_{r.name}" class="reaction"><rect x="{xr - bw / 2:.1f}" y="{y - bh / 2:.1f}" width="{bw}" '
                       f'height="{bh}" rx="6" fill="#fff4d6" stroke="#b08a2e"/>'
                       f'<text x="{xr:.1f}" y="{y - 3:.1f}" text-anchor="middle" font-weight="bold">Re: {_x(r.name)}</text>'
                       f'<text x="{xr:.1f}" y="{y + 11:.1f}" text-anchor="middle" font-size="10">{_x(r.enzyme)}</text></g>')
        for sp in self.species:
            x, y = place[sp.name]
            dash = ' stroke-dasharray="4 3"' if sp.name in bnd else ''
            out.append(f'<g id="q_{sp.name}" class="species"><circle cx="{x:.1f}" cy="{y:.1f}" r="{rad}" fill="#e3f0fb" '
                       f'stroke="#3b6e9e"{dash}/><text x="{x:.1f}" y="{y + 4:.1f}" text-anchor="middle">{_x(sp.name)}</text>'
                       f'<title>{_x(lab[sp.name])}{" (chemostat)" if sp.name in bnd else ""}</title></g>')
        out.append(f'<text x="{x0 + 10:.0f}" y="18" fill="#555">{_x(self.label)}: species (0-junction + C), '
                   'reactions (Re between 1-junctions); dashed: chemostats</text>')
        out.append('</svg>')
        return '\n'.join(out) + '\n'

def _x(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
