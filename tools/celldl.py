"""CellDL documents of a module's bond graph (tools/bondgraph.py): its BG-RDF, with each element the BG-RDF template of
its constitutive relation, laid out for the CellDL editor.

    from celldl import celldl_spec

celldl_spec(bg, conc, kappa, volume) gives the document's Turtle and its elements and connections with their positions.
The CellDL editor's render tool (CellDLEditorComponent `bun run render-celldl spec.json out.celldl.svg`) draws each
element as the editor does and writes the CellDL file, which CellDL's bg2cellml (BG-Tools) turns into CellML.

The elements:
  a species in the module     a store: bgf:ChemicalStorageNode, or bgf:SaturatingChemicalStorageNode for a species
                              with a saturating activity (k_sat); q = C V [mol], K [1/mol] so that K q is the
                              module's exp(mu0/RT) a/c_std (the template's RT then only scales u, not the rates)
  a reaction                  bgf:ChemicalReactionNode, or bgf:SaturatingChemicalReactionNode when it saturates
                              (sigma); kappa [mol/s]
  a chemostat                 bgf:ChemicalPotential, u = RT ln(K C V), one beside each reaction that uses it (each
                              copy holds the same potential, so the model is unchanged and its bonds stay short)
Not representable in the templates, so folded into kappa at `conc` (exact while they stay at it, noted in the result):
allosteric activators and inhibitors, and an electrical potential (BondGraph.potentials, applied to the stores' K, the
electrical work's even split moved into kappa and sigma).
"""
import math

from bondgraph import R_GAS

TEMPLATE_RT = 2578.73058   # J/mol: the BG-RDF chemical template's RT (37 C)
BGF = 'https://bg-rdf.org/ontologies/bondgraph-framework#'

DY = 95.0          # between consecutive reactions down the column
BESIDE = 95.0      # from a reaction to the species beside it
SIDE = 125.0       # spacing of the side columns


def _ucum(value, units):
    return f'"{value:.15g} {units}"^^cdt:ucum'


def _cycle_layout(bg):
    """Reaction positions: the reactions before the cycle (BondGraph.cycle) in a column, the cycle a ring, the rest in a
    column below. Returns (positions, ring, ring centre y, ring radius)."""
    rr = [r.name for r in bg.reactions]
    order = {n: i for i, n in enumerate(rr)}
    ring = bg.cycle()
    pos, y = {}, 70.0
    first = min((order[c] for c in ring), default=len(rr))
    for n in [n for n in rr if n not in ring and order[n] < first]:
        pos[n], y = (0.0, y), y + DY
    yc = R = 0.0
    if ring:
        R = max(150.0, 34.0 * len(ring))
        yc = y + R - DY / 2
        for k, n in enumerate(ring):
            a = -0.5 * math.pi + 2.0 * math.pi * k / len(ring)
            pos[n] = (R * math.cos(a), yc + R * math.sin(a))
        y = yc + R + DY
    for n in [n for n in rr if n not in ring and order[n] > first]:
        pos[n], y = (0.0, y), y + DY
    return pos, ring, yc, R


def _outward(pos, n, ring, yc, distance, k):
    """A point beside reaction n: outside the ring along its radius for a ring reaction, else to alternate sides."""
    x, y = pos[n]
    if n in ring:
        # copies fan out around the radius, alternately either side of it
        a = math.atan2(y - yc, x) + (1 if k % 2 else -1) * 0.38 * ((k + 1) // 2)
        return x + math.cos(a) * distance, y + math.sin(a) * distance
    side = -1 if k % 2 == 0 else 1
    return x + side * distance * (1 + 0.55 * (k // 2)), y + 28.0 * ((k // 2) % 2)


def celldl_spec(bg, conc, kappa, volume, location='c'):
    """The CellDL document of bg at concentrations conc [mM] with rate constants kappa [mM/s], for a compartment of
    volume [m3]: {'label', 'turtle', 'elements': [...], 'connections': [...], 'folded': [...]}."""
    rt = bg.gas_constant * bg.temperature
    by = bg.by_name()
    pos, ring, yc, R = _cycle_layout(bg)

    def factor(s):
        sp = by[s]
        return math.exp((bg.mu0[s] + bg.electrical(s)) / rt) * bg.activity(s, conc[s]) / sp.c_std

    elements, connections, folded, turtle = [], [], [], []
    bond_elements = []

    def element(eid, rdf_type, symbol, species, x, y, values, where=None):
        where = where or location
        elements.append({'id': eid, 'type': BGF + rdf_type, 'symbol': symbol, 'species': species,
                         'location': where, 'x': round(x, 2), 'y': round(y, 2)})
        bond_elements.append(f':{eid}')
        lines = [f':{eid} a celldl:Component, bgf:{rdf_type} ;', f'  bgf:hasSymbol "{symbol}" ;',
                 f'  bgf:hasSpecies "{species}" ;', f'  bgf:hasLocation "{where}"']
        for v in values:
            lines[-1] += ' ;'
            lines.append(f'  {v}')
        turtle.append('\n'.join(lines) + ' .')

    def connect(a, b, n):
        cid = f'{a}.{b}'
        connections.append({'id': cid, 'source': a, 'target': b, 'count': n})
        count = f' ;\n  bgf:bondCount {n}' if n > 1 else ''
        turtle.append(f':{cid} a celldl:Connection ;\n  celldl:hasSource :{a} ;\n  celldl:hasTarget :{b} ;\n'
                      f'  bgf:hasSource :{a} ;\n  bgf:hasTarget :{b}{count} .')

    # stores: the module's species, placed between the reactions they link, beside the one that uses them, or in the
    # side columns
    internal = [s for s in bg.species if not s.boundary]
    uses = {s.name: [r.name for r in bg.reactions if s.name in r.stoich] for s in internal}
    rx = {r.name: r for r in bg.reactions}
    place, besides, hubs, between = {}, {}, [], {}
    for s in internal:
        u = uses[s.name]
        on_ring = [k for k in range(len(ring)) if rx[ring[k]].stoich.get(s.name, 0) > 0
                   and rx[ring[(k + 1) % len(ring)]].stoich.get(s.name, 0) < 0 and s.name not in bg.cofactors]
        if on_ring:
            k = on_ring[0]
            j = between.get(('ring', k), 0)
            between[('ring', k)] = j + 1
            a = -0.5 * math.pi + 2.0 * math.pi * (k + 0.5) / len(ring)
            r_ = R - 50.0 * j
            place[s.name] = (r_ * math.cos(a), yc + r_ * math.sin(a))
        elif len(u) == 2 and math.dist(pos[u[0]], pos[u[1]]) < 1.6 * DY:
            j = between.get(tuple(u), 0)
            between[tuple(u)] = j + 1
            (xa, ya), (xb, yb) = pos[u[0]], pos[u[1]]
            place[s.name] = (0.5 * (xa + xb) + (0.0 if j == 0 else (-1) ** j * 55.0 * ((j + 1) // 2)), 0.5 * (ya + yb))
        elif len(u) == 1:
            k = besides.get(u[0], 0)
            besides[u[0]] = k + 1
            place[s.name] = _outward(pos, u[0], ring, yc, BESIDE, k)
        else:
            hubs.append(s.name)
    edge = R + 2.4 * BESIDE
    taken = {-1: [], 1: []}
    for i, n in enumerate(sorted(hubs, key=lambda n: sum(pos[r][1] for r in uses[n]) / len(uses[n]))):
        side = -1 if i % 2 == 0 else 1
        y = sum(pos[r][1] for r in uses[n]) / len(uses[n])
        col = 0
        while any(c == col and abs(y - y2) < 50.0 for c, y2 in taken[side]):
            col += 1
        taken[side].append((col, y))
        place[n] = (side * (edge + SIDE * col), y)
    for s in internal:
        x, y = place[s.name]
        k = math.exp((bg.mu0[s.name] + bg.electrical(s.name)) / rt) / (s.c_std * volume)
        values = [f'bgf:hasValue {_ucum(conc[s.name] * volume, "mol")}']
        if s.k_sat is not None:
            # K q k_sat/(k_sat + q) = K V a with a = C k_sat/(k_sat + C): the module's exp(mu0/RT) a/c_std
            values.append(f'bgf:parameterValue [ bgf:varName "K" ; bgf:hasValue {_ucum(k, "mol-1")} ] , '
                          f'[ bgf:varName "k_sat" ; bgf:hasValue {_ucum(s.k_sat * volume, "mol")} ]')
            element(f'q_{s.name}', 'SaturatingChemicalStorageNode', 'q', s.name, x, y, values, s.compartment)
        else:
            values.append(f'bgf:parameterValue [ bgf:varName "K" ; bgf:hasValue {_ucum(k, "mol-1")} ]')
            element(f'q_{s.name}', 'ChemicalStorageNode', 'q', s.name, x, y, values, s.compartment)

    # reactions, with a copy of each chemostat they use beside them
    for r in bg.reactions:
        x, y = pos[r.name]
        sub_el = sum(-n * bg.electrical(sp) for sp, n in r.stoich.items() if n < 0)
        c = math.exp((0.5 * bg.electrical_work(r) + sub_el) / rt)   # the templates' Pf and Pr over the module's
        regulation = bg.activation(r, conc)
        if r.activators or r.inhibitors:
            folded.append(f'{r.name}: ' + ', '.join([f'activation by {a}' for a in r.activators] +
                                                    [f'inhibition by {i}' for i in r.inhibitors]) +
                          f' folded into kappa (x{regulation:.4g}) at the export concentrations')
        if bg.electrical_work(r):
            folded.append(f'{r.name}: the electrical work at psi = ' +
                          ', '.join(f'{v * 1e3:g} mV' for v in bg.potentials.values()) + ' folded into K, kappa, sigma')
        values = [f'bgf:parameterValue [ bgf:varName "kappa" ; bgf:hasValue '
                  f'{_ucum(kappa[r.name] * regulation * volume / c, "mol/s")} ]']
        if r.sigma:
            values[-1] += f' , [ bgf:varName "sigma" ; bgf:hasValue {_ucum(r.sigma * c, "1")} ]'
            element(f'Re_{r.name}', 'SaturatingChemicalReactionNode', 'Re', r.name, x, y, values)
        else:
            element(f'Re_{r.name}', 'ChemicalReactionNode', 'Re', r.name, x, y, values)
        k = besides.get(r.name, 0)
        for sp, n in r.stoich.items():
            if not by[sp].boundary:
                continue
            sx, sy = _outward(pos, r.name, ring, yc, BESIDE, k)
            k += 1
            se = f'Se_{sp}_{r.name}'
            element(se, 'ChemicalPotential', 'Se', sp, sx, sy,
                    [f'bgf:hasValue {_ucum(TEMPLATE_RT * math.log(factor(sp)), "J/mol")}'], by[sp].compartment)
            connect(se, f'Re_{r.name}', abs(n)) if n < 0 else connect(f'Re_{r.name}', se, n)
        for sp, n in r.stoich.items():
            if not by[sp].boundary:
                connect(f'q_{sp}', f'Re_{r.name}', abs(n)) if n < 0 else connect(f'Re_{r.name}', f'q_{sp}', n)
    head = ['@prefix : <#> .', '@prefix bgf: <https://bg-rdf.org/ontologies/bondgraph-framework#> .',
            '@prefix cdt: <https://w3id.org/cdt/> .', '@prefix celldl: <http://celldl.org/ontologies/celldl#> .',
            '@prefix owl: <http://www.w3.org/2002/07/owl#> .', '@prefix tpl: <https://bg-rdf.org/templates/> .', '',
            f'# {bg.label}: written by tools/celldl.py (circulatory-autogen-modules) from the module\'s bond graph',
            *[f'# {f}' for f in folded], '',
            '<> a celldl:Document, bgf:BondgraphModel ;', '  owl:versionInfo "1.0" ;',
            '  bgf:usesTemplate tpl:chemical.ttl ;', '  bgf:hasBondElement ' + ', '.join(bond_elements) + ' .', '']
    return {'label': bg.label, 'turtle': '\n'.join(head) + '\n\n'.join(turtle) + '\n',
            'elements': elements, 'connections': connections, 'folded': folded}
