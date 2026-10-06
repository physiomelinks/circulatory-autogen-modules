"""
Rekey every BibTeX key in the library to the 2026-10-06 convention, <first author's surname><year><first
significant title word>, all lowercase ASCII (belluzzi1986quantitative; cam_testing.bib.convention_key),
and merge entries that cite the same paper under two keys.

    venv/bin/python tools/rekey_bib.py --dry-run    # print the map, the merges and every change, write nothing
    venv/bin/python tools/rekey_bib.py              # rewrite the files and tools/bib_rekey_map.json

What it rewrites (tracked files only, so other sessions' uncommitted work is never touched):
- the entry keys of every *_references.bib and *_references_proposed.bib. Where two old keys become one
  key in the same file (a duplicate pair), the entry with more fields is kept and the other dropped;
- a key cited as a key in the library's data and records: under modules/ and system_models/ the
  parameter CSVs, tests.yaml (data_reference / reference_proposals reference strings and notes,
  also_cites), verification_config.json (validation source strings, notes), obs_data.json
  (reference_key, comments), source_figures.json and SOURCES.md; the review records reviews/*.yaml;
  and the tools that write such files (tools/*.py, this script and the historical proposal generator
  excepted).

How a key is recognised outside the .bib files, so that nothing that merely looks like a key changes:
- a key is matched case-sensitively, as a whole token: not inside a longer word or identifier
  (Paci2013_v01, the version, and davis2020_wistar, the instance, stay), not a file name or part of one
  (no '.', ':', '@', '#', '\\' or '-' before it, no '.<extension>' after it), and not in a whitespace-delimited
  token that is a URL or a DOI. "Key-derived" and "Key1/Key2" are citations and are rewritten;
- the same rule rewrites keys cited inside .bib entries (notes such as "as tabulated by Albanese2014");
- the software citation's old key circulatory_autogen is also the software's name in prose, so it is
  replaced only where it is a reference string's leading key ("circulatory_autogen; ..." in a CSV
  data_reference, a YAML reference/source value or a JSON string) or a YAML list item;
- README files, CellML comments, generated HTML and the historical review pages are left alone.

Keys:
- each old key's entry (its first occurrence; every occurrence of a key must give the same new key) gives
  the new key with cam_testing.bib.convention_key;
- two old keys whose entries are the same paper (same DOI, or same title and year) share the key (merged);
- two different papers with the same key both take further title words until they differ
  (deterministic; listed as collisions in the map);
- OVERRIDES: the software citation circulatory_autogen (@misc, no year) becomes argus2026circulatory, the
  key the required_citations proposal uses, and its entries get year = {2026} and a placeholder note.

The map (old -> new, merges, collisions, per-file change counts) is written to tools/bib_rekey_map.json.
The script can be re-run: old keys no longer in any .bib are taken from that map, so a later run rewrites
citations an earlier run missed (the counts accumulate).
"""
import argparse
import collections
import json
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from cam_testing import bib  # noqa: E402

MAP_PATH = os.path.join(REPO, 'tools', 'bib_rekey_map.json')
EXCLUDED_DIRS = ('modules/system/', 'modules/poiseuille_transport/', 'handovers/')
SKIP_FILES = {'tools/rekey_bib.py', 'tools/library_restructure_proposal.py', 'tools/bib_rekey_map.json'}

# keys that are also ordinary words: replaced only where they are a reference's leading key
PROSE_KEYS = {'circulatory_autogen'}
OVERRIDES = {'circulatory_autogen': 'argus2026circulatory'}
OVERRIDE_FIELDS = {
    'circulatory_autogen': [
        ('year', '2026'),
        ('note', 'Placeholder software citation for circulatory_autogen (libcuflynx) and this module library, '
                 'until a release or paper is cited; year = the version cited'),
    ],
}

BEFORE = r'(?<![A-Za-z0-9_.:@#\\-])'
AFTER = r'(?![A-Za-z0-9_]|\.[A-Za-z0-9])'


def git_files(*patterns):
    out = subprocess.run(['git', 'ls-files', '-z', '--', *patterns], cwd=REPO, capture_output=True, text=True,
                         check=True).stdout
    return sorted(p for p in out.split('\0') if p and not p.startswith(EXCLUDED_DIRS))


def bib_files():
    return [p for p in git_files('*.bib') if p.endswith(('_references.bib', '_references_proposed.bib'))]


def text_files():
    '''The non-.bib files a key may be cited in (see the module docstring).'''
    out = []
    for p in git_files('modules', 'system_models'):
        base = os.path.basename(p)
        if base.endswith(('.csv', '_tests.yaml', '_verification_config.json', '_obs_data.json')) \
                or base in ('source_figures.json', 'SOURCES.md'):
            out.append(p)
    out += [p for p in git_files('reviews/*.yaml', 'tools/*.py') if p not in SKIP_FILES]
    return out


def read(p):
    with open(os.path.join(REPO, p), encoding='utf-8', newline='') as f:    # keep CRLF files CRLF
        return f.read()


def same_paper(a, b):
    doi = lambda e: (e.get('doi') or '').lower().strip()          # noqa: E731
    if doi(a) and doi(a) == doi(b):
        return True
    norm = lambda e: re.sub(r'[^a-z0-9]', '', bib._ascii(e.get('title') or '').lower())   # noqa: E731
    return bool(norm(a)) and norm(a) == norm(b) and (a.get('year') or '') == (b.get('year') or '')


def build_map(files):
    '''old key -> new key, merges {new: [old, ...]}, collisions {new: [old, ...]}.'''
    first, new_of = {}, {}
    for p in files:
        for kind, key, _, _, raw in bib.entries(read(p)):
            fields = dict(raw)
            for name, value in OVERRIDE_FIELDS.get(key, []):
                fields.setdefault(name, value)
            nk = OVERRIDES.get(key) or bib.convention_key(fields)
            if key in new_of and new_of[key] != nk:
                sys.exit(f'{key}: its entries give different keys ({new_of[key]} in an earlier file, {nk} in {p})')
            first.setdefault(key, fields)
            new_of[key] = nk
    by_new = collections.defaultdict(list)
    for k, nk in new_of.items():
        by_new[nk].append(k)
    merges, collisions = {}, {}
    for nk, olds in sorted(by_new.items()):
        if len(olds) == 1:
            continue
        # group the old keys into papers
        papers = []
        for k in sorted(olds):
            for g in papers:
                if same_paper(first[g[0]], first[k]):
                    g.append(k)
                    break
            else:
                papers.append([k])
        if len(papers) == 1:
            merges[nk] = sorted(olds)
            continue
        collisions[nk] = sorted(olds)
        n = 2
        while True:
            keys = [OVERRIDES.get(g[0]) or bib.convention_key(first[g[0]], n) for g in papers]
            if len(set(keys)) == len(keys) or n > 12:
                break
            n += 1
        for g, k2 in zip(papers, keys):
            for k in g:
                new_of[k] = k2
            if len(g) > 1:
                merges[k2] = sorted(g)
    taken = collections.Counter(new_of.values())
    clash = [k for k, n in taken.items() if n > 1 and k not in merges]
    if clash:
        sys.exit(f'unresolved key clashes: {clash}')
    return new_of, merges, collisions


def rekey_bib_text(text, new_of):
    '''Rename every entry key; drop the poorer of two entries that now share a key.'''
    ents = bib.entries(text)
    keep = {}
    for i, (kind, key, s, e, raw) in enumerate(ents):
        nk = new_of[key]
        if nk in keep:
            j = keep[nk]
            if len(raw) > len(ents[j][4]):
                keep[nk] = i
        else:
            keep[nk] = i
    drop = {i for i in range(len(ents))} - set(keep.values())
    out, pos, n_renamed, n_dropped = [], 0, 0, 0
    for i, (kind, key, s, e, raw) in enumerate(ents):
        out.append(text[pos:s])
        entry = text[s:e]
        pos = e
        if i in drop:
            n_dropped += 1
            # also drop the blank line that separated it
            if text[pos:pos + 2] == '\n\n':
                pos += 1
            continue
        nk = new_of[key]
        if nk != key:
            n_renamed += 1
            entry = re.sub(r'^(@\w+\s*\{\s*)' + re.escape(key) + r'(\s*,)', lambda m: m[1] + nk + m[2], entry, count=1)
        for name, value in OVERRIDE_FIELDS.get(key, []):
            if name not in raw:
                indent = re.search(r'\n(\s+)\w+\s*=', entry)
                indent = indent[1] if indent else '  '
                body = entry.rstrip()
                assert body.endswith('}')
                body = body[:-1].rstrip()
                if not body.endswith(','):
                    body += ','
                entry = f'{body}\n{indent}{name} = {{{value}}}\n}}'
        out.append(entry)
    out.append(text[pos:])
    return ''.join(out), n_renamed, n_dropped


def token_pattern(keys):
    alts = '|'.join(sorted((re.escape(k) for k in keys), key=len, reverse=True))
    return re.compile(BEFORE + '(' + alts + ')' + AFTER)


def in_url(text, start, end):
    s = max(text.rfind(c, 0, start) for c in ' \t\n"\'(<[') + 1
    tok_end = min((i for i in (text.find(c, end) for c in ' \t\n"\')>]') if i >= 0), default=len(text))
    tok = text[s:tok_end].lower()
    return '://' in tok or 'doi.org' in tok or tok.startswith(('10.', 'doi:', 'www.'))


def prose_key_patterns(p, key):
    k = re.escape(key)
    if p.endswith('.csv'):
        return [re.compile(r'(?:(?<=^)|(?<=,)|(?<=,"))' + k + r'(?=;|"|,|$)', re.M)]
    if p.endswith('.yaml'):
        return [re.compile(r"(?:(?<=reference: )|(?<=source: )|(?<=source_short: ))'?" + k + r"(?=;|'|$)", re.M),
                re.compile(r'(?:(?<=^-\s)|(?<=^\s\s-\s)|(?<=^\s\s\s\s-\s))' + k + r'(?=\s*$)', re.M)]
    if p.endswith('.json'):
        return [re.compile(r'(?<=")' + k + r'(?=;|")')]
    return []


def rekey_text(p, text, new_of, pattern):
    count = collections.Counter()

    def sub(m):
        if in_url(text, m.start(), m.end()):
            return m[0]
        count[m[1]] += 1
        return new_of[m[1]]
    text2 = pattern.sub(sub, text)
    for key in PROSE_KEYS:
        if key not in new_of or new_of[key] == key:
            continue
        for rx in prose_key_patterns(p, key):
            def sub2(m, key=key):
                count[key] += 1
                return m[0].replace(key, new_of[key])
            text2 = rx.sub(sub2, text2)
    return text2, count


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dry-run', action='store_true', help='print the map and the changes; write nothing')
    ap.add_argument('--verbose', action='store_true', help='with --dry-run: list every file changed')
    args = ap.parse_args()

    bibs = bib_files()
    new_of, merges, collisions = build_map(bibs)
    prev = None
    if os.path.isfile(MAP_PATH):
        with open(MAP_PATH) as f:
            prev = json.load(f)
        for k, v in prev['map'].items():
            if new_of.get(k, v) != v:
                sys.exit(f'{k}: maps to {new_of[k]} now but to {v} in {os.path.relpath(MAP_PATH, REPO)}')
            new_of[k] = v
        merges = {**prev['merged'], **merges}
        collisions = {**{k: v['old'] for k, v in prev['collisions'].items()}, **collisions}
    changed = {k: v for k, v in new_of.items() if k != v}
    pattern = token_pattern([k for k in changed if k not in PROSE_KEYS])

    per_file, per_key = {}, collections.Counter()
    n_ren = n_drop = 0
    writes = {}
    for p in bibs:
        t = read(p)
        t2, r, d = rekey_bib_text(t, new_of)
        t2, c = rekey_text(p, t2, new_of, pattern)
        per_key.update(c)
        n_ren += r
        n_drop += d
        if t2 != t:
            writes[p] = t2
            per_file[p] = {'entries_renamed': r, 'entries_dropped': d, 'citations': sum(c.values())}
    for p in text_files():
        t = read(p)
        t2, c = rekey_text(p, t, new_of, pattern)
        if t2 != t:
            writes[p] = t2
            per_file[p] = {'citations': sum(c.values())}
            per_key.update(c)

    if prev is not None:     # accumulate over runs
        per_key.update(prev['citations_rewritten_per_key'])
        n_ren += prev['summary']['n_bib_entries_renamed']
        n_drop += prev['summary']['n_bib_entries_dropped']
        files = dict(prev['files'])
        for f, c in per_file.items():
            old = files.get(f, {})
            files[f] = {k: old.get(k, 0) + c.get(k, 0) for k in set(old) | set(c)}
        per_file = files
    unused = sorted(k for k in changed if per_key[k] == 0)
    summary = {
        'convention': '<first author surname><year><first significant title word>, lowercase ASCII (cam_testing.bib.convention_key)',
        'n_keys_rekeyed': len(changed), 'n_keys_after': len(set(new_of.values())),
        'n_merged_pairs': sum(len(v) - 1 for v in merges.values()), 'n_collisions': len(collisions),
        'n_bib_files': len(bibs), 'n_bib_entries_renamed': n_ren, 'n_bib_entries_dropped': n_drop,
        'n_files_changed': len(per_file), 'n_citations_rewritten': sum(per_key.values()),
    }
    out = {'summary': summary,
           'map': dict(sorted(new_of.items(), key=lambda kv: kv[0].lower())),
           'merged': merges, 'collisions': {k: {'old': v, 'new': sorted({new_of[o] for o in v})} for k, v in collisions.items()},
           'overrides': OVERRIDES,
           'citations_rewritten_per_key': dict(sorted(per_key.items(), key=lambda kv: kv[0].lower())),
           'keys_cited_only_in_bib_files': unused,
           'files': dict(sorted(per_file.items()))}
    print(json.dumps(summary, indent=1))
    print('merged:', json.dumps(merges))
    print('collisions:', json.dumps(out['collisions']))
    if args.dry_run:
        for k, v in out['map'].items():
            print(f'  {k:32s} -> {v}' + ('   (merged)' if any(k in o for o in merges.values()) else ''))
        if args.verbose:
            for p, c in out['files'].items():
                print(f'  {p}: {c}')
        print('dry run: nothing written')
        return
    for p, t in writes.items():
        with open(os.path.join(REPO, p), 'w', encoding='utf-8', newline='') as f:
            f.write(t)
    with open(MAP_PATH, 'w') as f:
        json.dump(out, f, indent=1)
        f.write('\n')
    print(f'wrote {len(writes)} files and {os.path.relpath(MAP_PATH, REPO)}')


if __name__ == '__main__':
    main()
