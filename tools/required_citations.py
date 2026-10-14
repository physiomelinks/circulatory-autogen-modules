"""
Write "required_citations" (and "required_citations_uncertain") into every version's config entry, from
the proposal approved on 2026-10-06 (reviews/library_restructure_proposal.html, section 4), and make
sure each cited key has its entry in the version's references.bib.

    venv/bin/python tools/required_citations.py --dry-run    # list every version's citations, write nothing
    venv/bin/python tools/required_citations.py              # write the configs and the .bib entries

- The proposal's rules (tools/library_restructure_proposal.py CITATION_RULES, first match wins) are
  regexes on the proposal's post-rename "<type path>/<version>". Today's versions are mapped back to
  that name through tools/restructure_modules_renames.json (today's path -> its pre-move path) and the
  proposal's CHANGES table (pre-move path -> the proposal's name), so the rules apply as approved even
  where the move named a version differently (nn_ prefixes dropped, source-named versions).
- "required_citations": the rule's keys, in its order. "required_citations_uncertain": the keys whose
  attribution the rule flags as a best guess (UNCERTAIN_KEYS below names them where the rule's note
  makes that clear; otherwise every key of an uncertain rule). It is left out when nothing is uncertain.
- A key not yet in the version's references.bib is added: moved from its references_proposed.bib when
  it is proposed there (a model citation is a decision, not a proposal; an emptied proposed .bib is
  removed), else copied from another version's .bib (without its note, which is about that version's
  values), else written from PLACEHOLDERS (the owner-original placeholders and the papers the proposal
  added, NEW_BIB there).
- Existing required_citations are kept: the script only fills entries that have none (--force rewrites).
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tools'))
from cam_testing import bib  # noqa: E402
from cam_testing.library import with_record_keys  # noqa: E402
import library_restructure_proposal as proposal  # noqa: E402

MODULES = os.path.join(REPO, 'modules')
EXCLUDED = ('modules/system/', 'modules/poiseuille_transport/')
S, U, G, C = proposal.ARGUS_SN, proposal.ARGUS_UNP, proposal.GEE_ARGUS, proposal.CA_SOFT

# For an uncertain rule: the keys that are the guess, where the rule's note names them (regex as in
# CITATION_RULES). Rules not listed mark all their keys uncertain.
UNCERTAIN_KEYS = {
    r'cell/ion_channels/i_A/Tao2011_v01$': ['belluzzi1985fast'],
    r'cell/ion_channels/i_CaN/Argus2026_v01$': ['belluzzi1989calcium'],
    r'cell/ion_channels/i_Kv1_5/Argus2026_v01$': ['philipson1991sequence', 'ranjanchannelpedia'],
    r'cell/ion_channels/i_SK/Argus2026_v01$': ['hirschberg1998gating'],
    r'cell/ion_channels/i_M/Argus2026_v01$': ['martin2023modelling'],
    r'cell/ion_channels/i_NaK/Argus2026_v01$': ['nygren1998mathematical'],
    r'cell/ion_channels/i_NaK/ArgusUNPUBLISHED_v01$': ['smith2004development', 'terkildsen2007balance', 'pan2020cardiac'],
    r'cell/ion_channels/IP3R/Argus2026_v01$': ['deyoung1992single'],
    r'cell/Ca_handling/SN_soma_Argus2026_v01$': ['shannon2004mathematical'],
    r'cell/Ca_handling/SN_varicosity_Argus2026_v01$': ['shannon2004mathematical'],
    r'cell/neurotransmitter_release/SN_varicosity_Tao2011_v01$': ['schneggenburger2000intracellular'],
    r'control/NTS/ArgusUNPUBLISHED_v01$': ['park2020investigating'],
    r'control/afferent_to_(sympathetic|vagal)_efferent/': ['ursino2000acute'],
    r'vessels/properties/material_prop_visco': ['alastruey2011pulse'],
    r'respiratory/lung_mechanics/ventilated_ArgusUNPUBLISHED_v01$': ['albanese2016integrated'],
}

PLACEHOLDER_NOTE = 'Placeholder citation (required_citations, 2026-10-06): replace with the real reference when it exists'
# BibTeX for the keys no .bib had: the owner-original placeholders and the papers the proposal added.
PLACEHOLDERS = {
    S: ('unpublished', {'author': 'Argus, Finbar and others', 'title': 'Sympathetic neuron model (title to be decided)',
                        'year': '2026', 'note': PLACEHOLDER_NOTE + ': the planned sympathetic-neuron paper, in preparation'}),
    U: ('unpublished', {'author': 'Argus, Finbar', 'title': 'Unpublished modules of the circulatory-autogen module library',
                        'year': '2026', 'note': PLACEHOLDER_NOTE + ": the owner's unpublished or modified models, no publication"}),
    G: ('unpublished', {'author': 'Gee, Michelle M. and Argus, Finbar', 'title': 'Unpublished respiratory control work (title to be decided)',
                        'year': '2026', 'note': PLACEHOLDER_NOTE + ': the Gee and Argus respiratory-control work, in preparation'}),
    C: ('misc', {'author': 'Argus, Finbar and Maso Talou, Gonzalo D. and others', 'title': 'circulatory\\_autogen (libcuflynx)',
                 'howpublished': 'GitHub', 'url': 'https://github.com/physiomelinks/circulatory_autogen', 'year': '2026',
                 'note': 'Placeholder software citation for circulatory_autogen (libcuflynx) and this module library, '
                         'until a release or paper is cited; year = the version cited'}),
    'park2020investigating': ('article', {
        'author': 'Park, James H. and Gorky, Jonathan and Ogunnaike, Babatunde and Vadigepalli, Rajanikanth and Schwaber, James S.',
        'title': 'Investigating the effects of brainstem neuronal adaptation on cardiovascular homeostasis',
        'journal': 'Frontiers in Neuroscience', 'year': '2020', 'volume': '14', 'pages': '470',
        'note': 'Added for required_citations (2026-10-06); details to verify'}),
    'westerhof2009arterial': ('article', {
        'author': 'Westerhof, Nico and Lankhaar, Jan-Willem and Westerhof, Berend E.', 'title': 'The arterial Windkessel',
        'journal': 'Medical \\& Biological Engineering \\& Computing', 'year': '2009', 'volume': '47', 'number': '2',
        'pages': '131--141', 'doi': '10.1007/s11517-008-0359-2', 'note': 'Added for required_citations (2026-10-06); details to verify'}),
    'fitzhugh1961impulses': ('article', {
        'author': 'FitzHugh, Richard', 'title': 'Impulses and physiological states in theoretical models of nerve membrane',
        'journal': 'Biophysical Journal', 'year': '1961', 'volume': '1', 'number': '6', 'pages': '445--466',
        'doi': '10.1016/S0006-3495(61)86902-6'}),
    'nagumo1962active': ('article', {
        'author': 'Nagumo, J. and Arimoto, S. and Yoshizawa, S.', 'title': 'An active pulse transmission line simulating nerve axon',
        'journal': 'Proceedings of the IRE', 'year': '1962', 'volume': '50', 'number': '10', 'pages': '2061--2070',
        'doi': '10.1109/JRPROC.1962.288235'}),
    'lotka1925elements': ('book', {
        'author': 'Lotka, Alfred J.', 'title': 'Elements of Physical Biology', 'publisher': 'Williams \\& Wilkins',
        'address': 'Baltimore', 'year': '1925'}),
    'volterra1926fluctuations': ('article', {
        'author': 'Volterra, Vito', 'title': 'Fluctuations in the abundance of a species considered mathematically',
        'journal': 'Nature', 'year': '1926', 'volume': '118', 'pages': '558--560', 'doi': '10.1038/118558a0'}),
    'vanderpol1926relaxation': ('article', {
        'author': 'van der Pol, Balthasar', 'title': 'On ``relaxation-oscillations\'\'',
        'journal': 'The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science', 'year': '1926',
        'volume': '2', 'number': '11', 'pages': '978--992', 'doi': '10.1080/14786442608564127'}),
}


def rel(p):
    return os.path.relpath(p, REPO)


def versions():
    '''[(type path under modules/, module_type, version, version dir, config path)].'''
    out = []
    for cfg in sorted(glob.glob(os.path.join(MODULES, '**', 'versions', '*', '*_modules_config.json'), recursive=True)):
        if rel(cfg).startswith(EXCLUDED):
            continue
        vdir = os.path.dirname(cfg)
        tdir = os.path.dirname(os.path.dirname(vdir))
        out.append((os.path.relpath(tdir, MODULES), os.path.basename(tdir), os.path.basename(vdir), vdir, cfg))
    return out


def proposal_names():
    '''today's "<type path>/<version>" -> the name the proposal's rules were written for.'''
    to_old = {r['to_path']: r['from_path'] for r in json.load(open(os.path.join(REPO, 'tools', 'restructure_modules_renames.json')))}
    prop = {f'{o}/{ov}': f'{n}/{nv}' for o, ov, n, nv, *_ in proposal.CHANGES}
    return lambda cur: prop.get(to_old.get(cur, cur), to_old.get(cur, cur))


def citations(name):
    '''(keys, uncertain keys, basis) for the proposal name.'''
    for rx, keys, unc, basis in proposal.CITATION_RULES:
        if re.search(rx, name):
            if not unc:
                return list(keys), [], basis
            return list(keys), [k for k in UNCERTAIN_KEYS.get(rx, keys) if k in keys], basis
    raise SystemExit(f'{name}: no citation rule')


def library_entries():
    '''key -> the entry's BibTeX, from the first .bib (references before proposed) that has it, without its
    note (a note says what that version took from the paper, which need not hold for another version).'''
    files = subprocess.run(['git', 'ls-files', '-z', '--', '*.bib'], cwd=REPO, capture_output=True, text=True).stdout.split('\0')
    files = sorted((f for f in files if f and not f.startswith(EXCLUDED)), key=lambda f: ('proposed' in f, f))
    out = {}
    for f in files:
        text = open(os.path.join(REPO, f), encoding='utf-8').read()
        for kind, key, s, e, raw in bib.entries(text):
            if key not in out:
                fields = {k: v for k, v in raw.items() if k not in ('note', 'annote', 'comment')}
                out[key] = format_bibtex(key, kind.lower(), fields) if len(fields) < len(raw) else text[s:e]
    return out


def format_bibtex(key, kind, fields):
    width = max(len(k) for k in fields)
    lines = [f'  {k.ljust(width)} = {{{v}}}' for k, v in fields.items()]
    return f'@{kind}{{{key},\n' + ',\n'.join(lines) + '\n}'


def append_entry(path, entry):
    text = open(path, encoding='utf-8').read() if os.path.isfile(path) else ''
    text = text.rstrip('\n')
    with open(path, 'w', encoding='utf-8') as f:
        f.write((text + '\n\n' if text else '') + entry + '\n')


def remove_entry(path, key):
    text = open(path, encoding='utf-8').read()
    for kind, k, s, e, raw in bib.entries(text):
        if k == key:
            entry = text[s:e]
            text = (text[:s].rstrip('\n') + ('\n\n' if text[:s].strip() else '') + text[e:].lstrip('\n')).strip('\n')
            if text.strip():
                with open(path, 'w', encoding='utf-8') as f:
                    f.write(text + '\n')
            else:
                os.remove(path)
            return entry
    return None


def with_citations(entry, keys, uncertain):
    '''The config entry with the two keys right after "creator" (with_record_keys places those).'''
    entry = with_record_keys(entry)
    out = {}
    for k, v in entry.items():
        if k in ('required_citations', 'required_citations_uncertain'):
            continue
        out[k] = v
        if k == 'creator':
            out['required_citations'] = keys
            if uncertain:
                out['required_citations_uncertain'] = uncertain
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--force', action='store_true', help='rewrite required_citations that are already set')
    args = ap.parse_args()

    name_of = proposal_names()
    lib = library_entries()
    stats = dict(versions=0, uncertain_versions=0, uncertain_keys=0, moved=0, copied=0, written=0, kept=0)
    for tpath, mtype, version, vdir, cfg in versions():
        name = name_of(f'{tpath}/{version}')
        keys, uncertain, basis = citations(name)
        with open(cfg) as f:
            entries = json.load(f)
        if any(e.get('required_citations') for e in entries) and not args.force:
            stats['kept'] += 1
            continue
        stats['versions'] += 1
        stats['uncertain_versions'] += bool(uncertain)
        stats['uncertain_keys'] += len(uncertain)
        refs = os.path.join(vdir, f'{mtype}_{version}_references.bib')
        proposed = os.path.join(vdir, f'{mtype}_{version}_references_proposed.bib')
        have = set(bib.read(refs))
        actions = []
        for k in keys:
            if k in have:
                continue
            if k in bib.read(proposed):
                actions.append(('moved', k))
                if not args.dry_run:
                    append_entry(refs, remove_entry(proposed, k))
            elif k in lib:
                actions.append(('copied', k))
                if not args.dry_run:
                    append_entry(refs, lib[k])
            elif k in PLACEHOLDERS:
                kind, fields = PLACEHOLDERS[k]
                actions.append(('written', k))
                text = format_bibtex(k, kind, fields)
                lib[k] = text
                if not args.dry_run:
                    append_entry(refs, text)
            else:
                raise SystemExit(f'{tpath}/{version}: no BibTeX for {k}')
            stats[actions[-1][0]] += 1
        print(f'{tpath}/{version}: {keys}' + (f' uncertain {uncertain}' if uncertain else '')
              + (f'  [{", ".join(a + " " + k for a, k in actions)}]' if actions else ''))
        if not args.dry_run:
            entries = [with_citations(e, keys, uncertain) for e in entries]
            with open(cfg, 'w') as f:
                json.dump(entries, f, indent=2)
                f.write('\n')
    print(json.dumps(stats))
    if args.dry_run:
        print('dry run: nothing written')


if __name__ == '__main__':
    main()
