"""
The checks of a version whose bond graph is a CellDL document, <stem>_bondgraph.svg in its directory. A
version without one is not a bond graph, and both checks are N/A.

bondgraph_celldl_test     The document opens in the CellDL editor: an SVG whose celldl-rdf-metadata is Turtle
                          with one bgf:BondgraphModel; every bond element drawn as a celldl-Component and every
                          connection as a celldl-Connection between two of them; nothing drawn that the RDF
                          doesn't have. <stem>_bondgraph.ttl, when there is one, is the same graph. With BG_TOOLS
                          set (a BG-Tools checkout with its .venv), bg2cellml also makes CellML from it.
bondgraph_reactions_test  Every reaction annotated with a Rhea reaction (bqbiol:isVersionOf rhea:N, or
                          bqbiol:hasPart rhea:N for each reaction a lumped one combines) has Rhea's participants
                          and stoichiometry. Its participants are the species elements bonded to it (bqbiol:is
                          CHEBI:N on each), on the side the bond gives, the bond count their coefficient. The
                          direction may be reversed and one factor may scale every coefficient; a difference
                          only in H+ or H2O is a warning; any other fails. Rhea's reactions are read from the
                          package's cache (cam_testing/rhea.py).

The annotation format is BG-RDF's: see the CellDL editor's Annotations panel (FinbarArgus/CellDLEditorComponent,
src/plugins/bondgraph/annotations.ts) and bg2cellml's OMEX export of it.
"""
import itertools
import os
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from fractions import Fraction

from cam_testing import checks, rhea

BONDGRAPH_TESTS = ['bondgraph_celldl_test', 'bondgraph_reactions_test']
NOT_A_BONDGRAPH = 'not a bond graph: the version has no CellDL document (<version>_bondgraph.svg)'
NO_REACTIONS = 'the bond graph has no reactions'

BGF = 'https://bg-rdf.org/ontologies/bondgraph-framework#'
CELLDL = 'http://celldl.org/ontologies/celldl#'
BQBIOL = 'http://biomodels.net/biology-qualifiers/'
RDF_TYPE = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type'
IDENTIFIERS = 'https://identifiers.org/'
JUNCTIONS = {BGF + 'ZeroNode', BGF + 'OneNode'}
SOFT = {rhea.H_PLUS, rhea.WATER}


def document_path(version):
    return version.path('bondgraph.svg')


def is_bondgraph(version):
    return os.path.isfile(document_path(version))


# ---- reading a CellDL document -----------------------------------------------------------------

def _local(iri):
    return re.split(r'[#/]', str(iri))[-1]


def _tag(element):
    return element.tag.rsplit('}', 1)[-1]


def read_document(path):
    '''(svg root, rdflib graph, the model's URI, problems): the document's SVG and its BG-RDF.'''
    import rdflib
    problems = []
    root = ET.parse(path).getroot()
    if _tag(root) != 'svg':
        return root, None, None, [f'the document is not SVG (its root is <{_tag(root)}>)']
    metadata = next((e for e in root.iter() if e.get('id') == 'celldl-rdf-metadata'), None)
    if metadata is None:
        return root, None, None, ['no celldl-rdf-metadata: the SVG is not a CellDL document']
    if metadata.get('data-content-type', 'text/turtle') != 'text/turtle':
        problems.append(f'celldl-rdf-metadata is {metadata.get("data-content-type")}, not text/turtle')
    graph = rdflib.Graph()
    try:
        graph.parse(data=''.join(metadata.itertext()), format='turtle', publicID=_base(path))
    except Exception as e:
        return root, None, None, problems + [f'its Turtle does not parse: {e}']
    models = list(graph.subjects(rdflib.RDF.type, rdflib.URIRef(BGF + 'BondgraphModel')))
    if len(models) != 1:
        return root, graph, None, problems + [f'{len(models)} bgf:BondgraphModel in the document, not one']
    return root, graph, models[0], problems


def _base(path):
    return 'file://' + os.path.abspath(path)


def _ids(root, cls):
    return {e.get('id') for e in root.iter() if cls in (e.get('class') or '').split() and e.get('id')}


def celldl_problems(path, ttl_path=None):
    '''What would keep the CellDL editor from opening the document as a bond graph ([] when nothing).'''
    import rdflib
    from rdflib.compare import isomorphic
    root, graph, model, problems = read_document(path)
    if graph is None or model is None:
        return problems
    uri = rdflib.URIRef
    elements = set(graph.objects(model, uri(BGF + 'hasBondElement')))
    junctions = set(graph.objects(model, uri(BGF + 'hasJunctionStructure')))
    components, connections = _ids(root, 'celldl-Component'), _ids(root, 'celldl-Connection')
    for e in sorted(elements):
        if not any(str(t).startswith(BGF) for t in graph.objects(e, rdflib.RDF.type)):
            problems.append(f'{_local(e)} has no BG-RDF type')
        if _local(e) not in components:
            problems.append(f'{_local(e)} is not drawn (no celldl-Component with its id)')
    drawn = elements | junctions
    rdf_connections = set(graph.subjects(rdflib.RDF.type, uri(CELLDL + 'Connection')))
    for c in sorted(rdf_connections):
        if _local(c) not in connections:
            problems.append(f'connection {_local(c)} is not drawn')
        for end in ('hasSource', 'hasTarget'):
            node = graph.value(c, uri(BGF + end)) or graph.value(c, uri(CELLDL + end))
            if node is None:
                problems.append(f'connection {_local(c)} has no {end}')
            elif node not in drawn:
                problems.append(f'connection {_local(c)}: its {end[3:].lower()} {_local(node)} is not a bond element')
    known = {_local(e) for e in drawn} | {_local(c) for c in rdf_connections}
    for drawn_id in sorted((components | connections) - known):
        problems.append(f'{drawn_id} is drawn but not in the RDF')
    if ttl_path and os.path.isfile(ttl_path):
        ttl = rdflib.Graph()
        ttl.parse(ttl_path, format='turtle', publicID=_base(path))
        if not isomorphic(graph, ttl):
            problems.append(f'{os.path.basename(ttl_path)} is not the document\'s RDF')
    return problems


def bg2cellml_problem(path, bg_tools):
    '''bg2cellml's error making CellML from the document, or None when it makes it.'''
    python = os.path.join(bg_tools, '.venv', 'bin', 'python')
    with tempfile.TemporaryDirectory() as out:
        result = subprocess.run([python, '-m', 'bgtool', '--output', out, path], cwd=bg_tools,
                                env=dict(os.environ, PYTHONPATH=os.path.join(bg_tools, 'cli')),
                                capture_output=True, text=True, timeout=300)
        made = [f for f in os.listdir(out) if f.endswith('.cellml')]
    if result.returncode != 0 or not made:
        return (result.stderr or result.stdout).strip()[-1500:] or 'bg2cellml made no CellML'
    return None


# ---- reactions and their participants ------------------------------------------------------------

@dataclass
class Participant:
    element: str
    species: str
    location: str
    chebi: str | None
    coefficient: Fraction        # negative on the left (a substrate), positive on the right (a product)


@dataclass
class Reaction:
    element: str
    rhea: str | None = None
    parts: list = field(default_factory=list)
    participants: list = field(default_factory=list)

    @property
    def annotated(self):
        return bool(self.rhea or self.parts)


def _identifier(iri, prefix):
    s = str(iri)
    for start in (f'{IDENTIFIERS}{prefix}:', f'http://identifiers.org/{prefix}:', f'{IDENTIFIERS}{prefix.lower()}:'):
        if s.startswith(start):
            return s[len(start):]
    return None


def reactions(graph, model):
    '''The model's reactions (elements whose class is a ...ReactionNode, or annotated with a Rhea reaction),
    each with its Rhea terms and participants.'''
    import rdflib
    uri = rdflib.URIRef
    elements = set(graph.objects(model, uri(BGF + 'hasBondElement')))
    types = {e: {str(t) for t in graph.objects(e, rdflib.RDF.type)} for e in elements}
    rhea_of = {e: [i for i in (_identifier(o, 'rhea') for o in graph.objects(e, uri(BQBIOL + 'isVersionOf'))) if i]
               for e in elements}
    parts_of = {e: [i for i in (_identifier(o, 'rhea') for o in graph.objects(e, uri(BQBIOL + 'hasPart'))) if i]
                for e in elements}
    is_reaction = {e for e in elements
                   if any(t.endswith('ReactionNode') for t in types[e]) or rhea_of[e] or parts_of[e]}
    species = {e for e in elements - is_reaction if graph.value(e, uri(BGF + 'hasSpecies')) is not None
               and not types[e] & JUNCTIONS}
    junctions = set(graph.objects(model, uri(BGF + 'hasJunctionStructure'))) | {e for e in elements if types[e] & JUNCTIONS}
    bonds = []
    for c in graph.subjects(uri(BGF + 'hasSource'), None):
        source, target = graph.value(c, uri(BGF + 'hasSource')), graph.value(c, uri(BGF + 'hasTarget'))
        count = graph.value(c, uri(BGF + 'bondCount'))
        if source is not None and target is not None:
            bonds.append((source, target, Fraction(str(count)) if count is not None else Fraction(1)))
    out = []
    for r in sorted(is_reaction, key=str):
        reaction = Reaction(_local(r), (rhea_of[r] or [None])[0], parts_of[r])
        # walk out from the reaction through junctions; a bond into the reaction makes the far side its left
        stack = [(r, None, Fraction(1), frozenset([r]))]
        while stack:
            node, sign, scale, seen = stack.pop()
            for source, target, count in bonds:
                if node not in (source, target):
                    continue
                other = target if source == node else source
                if other in seen:
                    continue
                into = target == node
                side = sign if sign is not None else (-1 if into else 1)
                if other in species:
                    chebis = [_identifier(o, 'CHEBI') for o in graph.objects(other, uri(BQBIOL + 'is'))]
                    reaction.participants.append(Participant(
                        _local(other), str(graph.value(other, uri(BGF + 'hasSpecies'))),
                        str(graph.value(other, uri(BGF + 'hasLocation')) or ''),
                        next((f'CHEBI:{c}' for c in chebis if c), None), side * scale * count))
                elif other in junctions:
                    stack.append((other, side, scale * count, seen | {other}))
        out.append(reaction)
    return out


def annotated_rhea_ids(version):
    '''Every Rhea reaction the version's bond graph annotates.'''
    if not is_bondgraph(version):
        return []
    _, graph, model, _ = read_document(document_path(version))
    if model is None:
        return []
    return sorted({i for r in reactions(graph, model) for i in ([r.rhea] if r.rhea else []) + r.parts})


# ---- comparing a reaction with Rhea --------------------------------------------------------------

def _coefficient(text):
    '''A Rhea coefficient as (multiple of n, constant): "2" (0, 2), "n" (1, 0), "n+1" (1, 1).'''
    t = str(text).replace(' ', '')
    m = re.fullmatch(r'(\d*)n([+-]\d+)?', t)
    if m:
        return Fraction(int(m.group(1) or 1)), Fraction(int(m.group(2) or 0))
    return Fraction(0), Fraction(t)


def rhea_vector(entry):
    '''{(ChEBI, in|out|None): signed net coefficient} of a Rhea reaction: left negative. A polymer, its
    repeating unit's ChEBI term on both sides with indices n and n+1, nets to the units it gains; unlocated
    participants of a transport reaction are on the in side.'''
    nets = {}
    for p in entry['participants']:
        n, c = _coefficient(p.get('polymer_index') or '1')
        times = Fraction(str(p['coefficient']))
        n, c = n * times, c * times
        location = p['location'] or ('in' if entry.get('is_transport') else None)
        sign = -1 if p['side'] == 'L' else 1
        key = (p['chebi'], location)
        a, b = nets.get(key, (Fraction(0), Fraction(0)))
        nets[key] = (a + sign * n, b + sign * c)
    return {k: b for k, (a, b) in nets.items() if b != 0 or a != 0}


@dataclass
class Comparison:
    status: str                  # match, warn, fail
    text: str
    warnings: list = field(default_factory=list)


def _assignments(model_locations, rhea_locations):
    '''Ways to read Rhea's locations (in, out or None) as the model's.'''
    rhea_locations = sorted({l for l in rhea_locations}, key=str)
    if rhea_locations == [None]:
        yield {None: None}
        return
    for chosen in itertools.permutations(sorted(model_locations), len(rhea_locations)):
        yield dict(zip(rhea_locations, chosen))


def compare(reaction, entries):
    '''Compare a reaction with its Rhea reaction (one entry) or the reactions it lumps (several).'''
    names = {p['chebi']: p['name'] for e in entries for p in e['participants'] if p['chebi']}
    names.update({p.chebi: p.species for p in reaction.participants if p.chebi and p.chebi not in names})
    missing = [p.element for p in reaction.participants if not p.chebi]
    if missing:
        return Comparison('fail', f'no ChEBI term on {", ".join(sorted(set(missing)))}')
    model = {}
    for p in reaction.participants:
        model[(p.chebi, p.location)] = model.get((p.chebi, p.location), 0) + p.coefficient
    model = {k: v for k, v in model.items() if v != 0}
    vectors = [rhea_vector(e) for e in entries]
    model_locations = {l for _, l in model}
    best = None
    # a lumped reaction: each part forwards or backwards
    for signs in itertools.product((1, -1), repeat=len(vectors)):
        net = {}
        for s, v in zip(signs, vectors):
            for k, c in v.items():
                net[k] = net.get(k, 0) + s * c
        net = {k: c for k, c in net.items() if c != 0}
        for assignment in _assignments(model_locations, {l for _, l in net}):
            if None in assignment and len(model_locations) > 1:
                continue                              # Rhea's reaction isn't a transport but this one spans locations
            ref = {}
            for (chebi, l), c in net.items():
                key = (chebi, assignment[l] if l is not None else next(iter(model_locations), None))
                ref[key] = ref.get(key, 0) + c
            result = _match(model, ref, names)
            if best is None or result[0] < best[0]:
                best = result
    if best is None:
        spans = 'spans ' + ' and '.join(sorted(model_locations)) if len(model_locations) > 1 else 'is in one location'
        kind = 'a transport' if any(l for v in vectors for _, l in v) else 'not a transport'
        return Comparison('fail', f'Rhea\'s reaction is {kind}, but this one {spans}')
    _, comparison = best
    return comparison


def _match(model, ref, names):
    '''(cost, Comparison): the model's coefficients against k times Rhea's, for the best k.'''
    hard = lambda d: {k: v for k, v in d.items() if k[0] not in SOFT}
    m, r = hard(model), hard(ref)
    shared = [k for k in m if k in r]
    if not shared:
        return 10 ** 6, Comparison('fail', 'no participant in common with Rhea\'s reaction')
    # the factor most participants agree on (1 or -1 first on a tie), so the report names the odd ones out
    ratios = [m[key] / r[key] for key in shared]
    k = max(ratios, key=lambda q: (ratios.count(q), abs(q) == 1, q > 0))
    label = lambda key: f'{names.get(key[0], key[0])} ({key[0]}{", " + key[1] if key[1] else ""})'
    problems = []
    for key in sorted(set(m) | set(r), key=str):
        want, have = k * r.get(key, 0), m.get(key, 0)
        if want == have:
            continue
        if key not in m:
            problems.append(f'missing {label(key)}: Rhea has {abs(want)}')
        elif key not in r:
            problems.append(f'extra {label(key)}')
        else:
            problems.append(f'{label(key)}: {abs(have)} where Rhea has {abs(want)}')
    warnings = []
    for key in sorted(set(model) | set(ref), key=str):
        if key[0] in SOFT and k * ref.get(key, 0) != model.get(key, 0):
            name = 'H+' if key[0] == rhea.H_PLUS else 'H2O'
            warnings.append(f'{name}: {abs(model.get(key, 0))} where Rhea has {abs(k * ref.get(key, 0))}')
    how = ('reversed' if k < 0 else 'as written') + (f', scaled by {abs(k)}' if abs(k) != 1 else '')
    if problems:
        return len(problems), Comparison('fail', '; '.join(problems), warnings)
    if warnings:
        return 0.5, Comparison('warn', f'matches ({how}) but for {"; ".join(warnings)}', warnings)
    return 0, Comparison('match', f'matches ({how})')


# ---- the checks ----------------------------------------------------------------------------------

def celldl_check(version):
    test = 'bondgraph_celldl_test'
    if not is_bondgraph(version):
        return checks.Result(test, checks.NOT_APPLICABLE, NOT_A_BONDGRAPH)
    path = document_path(version)
    problems = celldl_problems(path, version.path('bondgraph.ttl'))
    bg_tools = os.environ.get('BG_TOOLS')
    made = None
    if not problems and bg_tools:
        problem = bg2cellml_problem(path, bg_tools)
        made = problem is None
        if problem:
            problems.append(f'bg2cellml: {problem}')
    metrics = {'document': os.path.basename(path), 'bg2cellml': made}
    if problems:
        return checks.Result(test, checks.FAILED, f'{len(problems)} problem(s): {problems[0]}', metrics, details=problems)
    return checks.Result(test, checks.PASSED, 'a CellDL bond graph the editor opens'
                         + ('; bg2cellml makes its CellML' if made else ' (set BG_TOOLS to run bg2cellml too)'), metrics)


def reactions_check(version, load=rhea.load):
    test = 'bondgraph_reactions_test'
    if not is_bondgraph(version):
        return checks.Result(test, checks.NOT_APPLICABLE, NOT_A_BONDGRAPH)
    _, graph, model, problems = read_document(document_path(version))
    if model is None:
        return checks.Result(test, checks.FAILED, f'the document does not read: {problems[0]}', details=problems)
    found = reactions(graph, model)
    if not found:
        return checks.Result(test, checks.NOT_APPLICABLE, NO_REACTIONS)
    lines, counts = [], {'match': 0, 'warn': 0, 'fail': 0, 'not annotated': 0}
    for r in found:
        if not r.annotated:
            counts['not annotated'] += 1
            lines.append(f'{r.element}: not annotated')
            continue
        ids = [r.rhea] if r.rhea else r.parts
        entries = [load(i) for i in ids]
        uncached = [i for i, e in zip(ids, entries) if e is None]
        if uncached:
            counts['fail'] += 1
            lines.append(f'{r.element}: Rhea {", ".join(uncached)} not in the cache '
                         f'(python -m cam_testing.rhea fetch {" ".join(uncached)})')
            continue
        c = compare(r, entries)
        counts[c.status] += 1
        equation = ' + '.join(e['equation'] for e in entries)
        lines.append(f'{r.element} ~ Rhea {"+".join(ids)} ({equation}): {c.text}')
    metrics = dict(counts, reactions=len(found))
    summary = (f'{counts["match"] + counts["warn"]} of {len(found)} reactions match Rhea'
               + (f' ({counts["warn"]} with an H+ or H2O difference)' if counts['warn'] else '')
               + (f', {counts["fail"]} do not' if counts['fail'] else '')
               + (f', {counts["not annotated"]} not annotated' if counts['not annotated'] else ''))
    if counts['fail']:
        status = checks.FAILED
    elif counts['not annotated']:
        status = checks.PENDING
    else:
        status = checks.PASSED
    return checks.Result(test, status, summary, metrics, details=lines)
