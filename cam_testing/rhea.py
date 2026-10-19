"""
Rhea's reactions, as the bond-graph reaction check compares a module's reactions with them.

The check reads a cache in the package (cam_testing/data/rhea/<id>.json), so CI never goes online. Refresh
or add entries with

    python -m cam_testing.rhea fetch 17825 18353 ...     # those reactions
    python -m cam_testing.rhea fetch --annotated         # every Rhea term the library's bond graphs annotate
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
import sys
import urllib.parse
import urllib.request

CACHE_DIR = os.path.join(os.path.dirname(__file__), 'data', 'rhea')
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

FIND = '''PREFIX rh: <http://rdf.rhea-db.org/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX CHEBI: <http://purl.obolibrary.org/obo/CHEBI_>
SELECT DISTINCT ?acc ?equation WHERE {{
  ?r rdfs:subClassOf rh:Reaction ; rh:accession ?acc ; rh:equation ?equation .
  {participants}
}} LIMIT 50'''


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


def _sparql(query):
    url = SPARQL + '?' + urllib.parse.urlencode({'query': query, 'format': 'json'})
    request = urllib.request.Request(url, headers={**HEADERS, 'Accept': 'application/sparql-results+json'})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode())['results']['bindings']


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
    '''Rhea's reaction from its SPARQL endpoint and REST API, as a cache entry.'''
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
    '''Approved Rhea reactions with every one of ``chebis`` as a participant: [(id, equation)].'''
    clauses = '\n  '.join(f'?r rh:side ?s{k} . ?s{k} rh:contains ?p{k} . ?p{k} rh:compound/rh:chebi CHEBI:{c.split(":")[1]} .'
                          for k, c in enumerate(chebis))
    rows = _sparql(FIND.format(participants=clauses))
    found = {(r['acc']['value'].split(':', 1)[1], r['equation']['value']) for r in rows}
    return sorted(found, key=lambda r: int(r[0]))


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
    for rid in sorted(set(ids), key=int):
        print(f'{rid}: {save(fetch(rid))}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
