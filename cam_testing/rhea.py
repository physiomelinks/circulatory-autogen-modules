"""
Rhea's reactions, as the bond-graph reaction check compares a module's reactions with them.

The check reads a cache in the package (cam_testing/data/rhea/<id>.json), so CI never goes online. Entries come
from Rhea's SPARQL endpoint, or when that fails from its REST API (whose entries have no charges or formulas).
Refresh or add entries with

    python -m cam_testing.rhea fetch 17825 18353 ...     # those reactions
    python -m cam_testing.rhea fetch --annotated         # every Rhea term the library's bond graphs annotate
                                                         # (and the ChEBI ancestors of their species and Rhea's)
    python -m cam_testing.rhea find CHEBI:4167 CHEBI:30616 ...   # Rhea reactions with all these participants

An entry: {"id", "equation", "is_transport", "participants": [{"chebi", "name", "side" (L or R),
"coefficient", "location" (in, out or null), "charge", "formula", "polymer_index"}]}. A polymer participant has
the ChEBI term of its repeating unit (Rhea's underlyingChebi) and its polymerisation index ("n", "n+1"), and no
charge. Rhea's IDs here are master
(undirected) reactions; its participants are ChEBI's major microspecies at pH 7.3.
"""
import argparse
import json
import os
import re
import sys
import urllib.parse
import urllib.request

CACHE_DIR = os.path.join(os.path.dirname(__file__), 'data', 'rhea')
# {ChEBI term: its is_a ancestors}, for a model's species written more or less specifically than Rhea's
ANCESTORS_PATH = os.path.join(CACHE_DIR, 'chebi_ancestors.json')
OLS = 'https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms'
SPARQL = 'https://sparql.rhea-db.org/sparql'
REST = 'https://www.rhea-db.org/rhea'
H_PLUS, WATER = 'CHEBI:15378', 'CHEBI:15377'
# Rhea's REST API refuses Python's default user agent
HEADERS = {'User-Agent': 'cam_testing (circulatory-autogen-modules bond-graph check)'}

PARTICIPANTS = '''PREFIX rh: <http://rdf.rhea-db.org/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
SELECT ?side ?coef ?chebi ?polymer ?index ?name ?charge ?formula ?loc ?transport WHERE {{
  rh:{id} rh:side ?side ; rh:isTransport ?transport . ?side ?contains ?p .
  ?contains rdfs:subPropertyOf rh:contains ; rh:coefficient ?coef .
  ?p rh:compound ?c . ?c rh:name ?name .
  OPTIONAL {{ ?c rh:chebi ?chebi }} OPTIONAL {{ ?c rh:charge ?charge }} OPTIONAL {{ ?c rh:formula ?formula }}
  OPTIONAL {{ ?c rh:underlyingChebi ?polymer ; rh:polymerizationIndex ?index }}
  OPTIONAL {{ ?p rh:location ?loc }} }}'''

class RheaError(Exception):
    pass


def cache_path(rhea_id):
    return os.path.join(CACHE_DIR, f'{rhea_id}.json')


def load(rhea_id):
    '''The cached reaction, or None when it isn't cached.'''
    path = cache_path(rhea_id)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


_ancestors = None


def ancestors(chebi):
    '''A ChEBI term's is_a ancestors from the cache (empty when it isn't cached).'''
    global _ancestors
    if _ancestors is None:
        _ancestors = json.load(open(ANCESTORS_PATH)) if os.path.isfile(ANCESTORS_PATH) else {}
    return set(_ancestors.get(chebi, []))


def related(a, b):
    '''Whether one ChEBI term is the other or is_a the other.'''
    return a == b or b in ancestors(a) or a in ancestors(b)


def fetch_ancestors(chebi):
    '''A ChEBI term's is_a ancestors from the EBI's Ontology Lookup Service.'''
    number = chebi.split(':')[1]
    iri = urllib.parse.quote(urllib.parse.quote(f'http://purl.obolibrary.org/obo/CHEBI_{number}', safe=''), safe='')
    found, url = [], f'{OLS}/{iri}/ancestors?size=500'
    while url:
        request = urllib.request.Request(url, headers={**HEADERS, 'Accept': 'application/json'})
        with urllib.request.urlopen(request, timeout=120) as response:
            data = json.loads(response.read().decode())
        found += [t['obo_id'] for t in data.get('_embedded', {}).get('terms', []) if (t.get('obo_id') or '').startswith('CHEBI:')]
        url = data.get('_links', {}).get('next', {}).get('href')
    return sorted(set(found))


def save_ancestors(chebis):
    '''Add the ancestors of ``chebis`` (those not cached yet) to the cache.'''
    global _ancestors
    cached = json.load(open(ANCESTORS_PATH)) if os.path.isfile(ANCESTORS_PATH) else {}
    for chebi in sorted(set(chebis) - set(cached)):
        cached[chebi] = fetch_ancestors(chebi)
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(ANCESTORS_PATH, 'w') as f:
        json.dump(cached, f, indent=1, sort_keys=True)
        f.write('\n')
    _ancestors = cached
    return ANCESTORS_PATH


def _sparql(query):
    url = SPARQL + '?' + urllib.parse.urlencode({'query': query, 'format': 'json'})
    request = urllib.request.Request(url, headers={**HEADERS, 'Accept': 'application/sparql-results+json'})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode())['results']['bindings']


def _rest(query, columns='rhea-id,equation,chebi-id', limit=None):
    params = {'query': query, 'columns': columns, 'format': 'tsv'}
    if limit:
        params['limit'] = limit
    url = REST + '?' + urllib.parse.urlencode(params)
    with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=120) as response:
        lines = response.read().decode().strip().splitlines()
    return [line.split('\t') for line in lines[1:]]


TERM = re.compile(r'^(?:(\d+|n) )?(.*?)(?:\((n(?:[+-]\d+)?)\))?(?:\((in|out)\))?$')


def parse_equation(equation, chebis):
    '''A cache entry's participants from Rhea's equation and its ChEBI list (distinct compounds, in the order the
    equation first names them): coefficients, (in)/(out) and a polymer's (n) or (n+1) read from the text.'''
    parts = []
    for side, text in zip('LR', equation.split(' = ')):
        for term in text.split(' + '):
            coefficient, name, index, location = TERM.match(term.strip()).groups()
            parts.append({'side': side, 'coefficient': coefficient if coefficient and coefficient != 'n' else '1',
                          'name': name, 'polymer_index': index, 'location': location})
    names = list(dict.fromkeys(p['name'] for p in parts))
    if len(names) != len(chebis):
        raise RheaError(f'{len(names)} compounds in "{equation}" but {len(chebis)} ChEBI terms')
    chebi_of = dict(zip(names, chebis))
    for p in parts:
        p.update(chebi=chebi_of[p['name']], charge=None, formula=None)
    return sorted(parts, key=lambda p: (p['side'], p['chebi'] or '', p['location'] or ''))


def fetch_rest(rhea_id):
    '''Rhea's reaction from its REST API, as a cache entry (no charges or formulas).'''
    rows = _rest(f'RHEA:{rhea_id}')
    row = next((r for r in rows if r[0] == f'RHEA:{rhea_id}'), None)
    if row is None:
        raise RheaError(f'Rhea has no reaction {rhea_id}')
    participants = parse_equation(row[1], row[2].split(';'))
    return {'id': str(rhea_id), 'equation': row[1], 'is_transport': any(p['location'] for p in participants),
            'participants': participants, 'source': 'rest'}


def _equation(rhea_id):
    url = REST + '?' + urllib.parse.urlencode({'query': f'RHEA:{rhea_id}', 'columns': 'rhea-id,equation',
                                               'format': 'tsv'})
    with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=120) as response:
        lines = response.read().decode().strip().splitlines()
    return lines[1].split('\t', 1)[1] if len(lines) > 1 else ''


def _chebi(iri):
    return 'CHEBI:' + iri.rsplit('CHEBI_', 1)[1] if 'CHEBI_' in iri else iri


def _charge(text):
    '''A participant's charge, or None for one that isn't a number (a polymer's is written per unit).'''
    try:
        return int(text)
    except (TypeError, ValueError):
        return None


def fetch(rhea_id):
    '''Rhea's reaction from its SPARQL endpoint and REST API, or from the REST API alone when SPARQL fails.'''
    try:
        return fetch_sparql(rhea_id)
    except (OSError, TimeoutError) as e:
        print(f'Rhea SPARQL failed for {rhea_id} ({e}); using the REST API', file=sys.stderr)
        return fetch_rest(rhea_id)


def fetch_sparql(rhea_id):
    '''Rhea's reaction from its SPARQL endpoint (participants, charges, formulas) and REST API (its equation).'''
    rows = _sparql(PARTICIPANTS.format(id=rhea_id))
    if not rows:
        raise RheaError(f'Rhea has no reaction {rhea_id} (or it is not a master reaction)')
    value = lambda row, key: row[key]['value'] if key in row else None
    participants = []
    for row in rows:
        loc = value(row, 'loc')
        chebi = value(row, 'chebi') or value(row, 'polymer')
        participants.append({
            'chebi': _chebi(chebi) if chebi else None, 'polymer_index': value(row, 'index'),
            'name': value(row, 'name'), 'side': value(row, 'side').rsplit('_', 1)[1],
            'coefficient': value(row, 'coef'),
            'location': loc.rsplit('/', 1)[1].lower() if loc else None,
            'charge': _charge(value(row, 'charge')),
            'formula': value(row, 'formula')})
    participants.sort(key=lambda p: (p['side'], p['chebi'] or '', p['location'] or ''))
    return {'id': str(rhea_id), 'equation': _equation(rhea_id),
            'is_transport': value(rows[0], 'transport') == 'true', 'participants': participants}


def save(entry):
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(cache_path(entry['id']), 'w') as f:
        json.dump(entry, f, indent=1, sort_keys=True)
        f.write('\n')
    return cache_path(entry['id'])


def find(chebis):
    '''Rhea reactions with every one of ``chebis`` as a participant: [(id, equation)], from its REST search.'''
    query = ' AND '.join(f'chebi:{c.split(":")[1]}' for c in chebis)
    rows = _rest(query, columns='rhea-id,equation', limit=50)
    return sorted({(r[0].split(':', 1)[1], r[1]) for r in rows}, key=lambda r: int(r[0]))


def main(argv=None):
    parser = argparse.ArgumentParser(description='The Rhea reactions the bond-graph check compares with.')
    sub = parser.add_subparsers(dest='command', required=True)
    get = sub.add_parser('fetch', help='fetch reactions into the cache')
    get.add_argument('ids', nargs='*')
    get.add_argument('--annotated', action='store_true', help="every Rhea term the library's bond graphs annotate")
    look = sub.add_parser('find', help='Rhea reactions with all these ChEBI participants')
    look.add_argument('chebis', nargs='+')
    args = parser.parse_args(argv)
    if args.command == 'find':
        for rid, equation in find(args.chebis):
            print(f'{rid}\t{equation}')
        return 0
    ids = list(args.ids)
    if args.annotated:
        from cam_testing import bondgraph
        from cam_testing.library import all_versions
        for version in all_versions():
            ids += bondgraph.annotated_rhea_ids(version)
    chebis = set()
    if args.annotated:
        from cam_testing import bondgraph
        from cam_testing.library import all_versions
        for version in all_versions():
            chebis |= set(bondgraph.annotated_chebis(version))
    for rid in sorted(set(ids), key=int):
        entry = fetch(rid)
        print(f'{rid}: {save(entry)}')
        chebis |= {p['chebi'] for p in entry['participants'] if p['chebi']}
    if chebis:
        print(f'ChEBI ancestors: {save_ancestors(chebis)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
