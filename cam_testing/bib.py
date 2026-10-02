"""
Per-version BibTeX: versions/<version>/<module_type>_<version>_references.bib holds the sources
cited by the version's parameters (all its instances). A parameter's data_reference is "<bibkey>; <note>" (the note is optional), or
starts with "definitional" / names a test fixture when there is nothing to cite.

A small reader, enough for the entries this repo writes (@type{key, field = {...}, ...}).
"""
import os
import re

_ENTRY_RE = re.compile(r'@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)\n\}', re.S)
_FIELD_RE = re.compile(r'(\w+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|"[^"]*"|\w+)\s*,?', re.S)


def bib_path(version):
    return version.path('references.bib')


def proposed_bib_path(version):
    '''Citations proposed for review, not yet confirmed (<module_type>_<version>_references_proposed.bib).'''
    return version.path('references_proposed.bib')


def read(path):
    '''{key: {'type': ..., field: value}} (braces and surrounding quotes removed).'''
    if not os.path.isfile(path):
        return {}
    with open(path) as f:
        text = f.read()
    entries = {}
    for kind, key, body in _ENTRY_RE.findall(text):
        fields = {'type': kind.lower()}
        for name, value in _FIELD_RE.findall(body):
            value = value.strip()
            if value[:1] in '{"':
                value = value[1:-1]
            fields[name.lower()] = re.sub(r'\s+', ' ', value.replace('{', '').replace('}', '')).strip()
        entries[key] = fields
    return entries


def reference_key(reference):
    '''The bib key a data_reference cites ("Key; note" -> "Key"), or None.'''
    head = (reference or '').split(';', 1)[0].strip()
    return head if re.fullmatch(r'[A-Za-z][\w:-]*', head) else None


def format_entry(fields):
    '''A one-line citation: Authors (year). Title. Venue volume(number):pages.'''
    authors = fields.get('author', '')
    names = [a.strip() for a in authors.split(' and ') if a.strip()]
    if len(names) > 3:
        authors = f'{names[0]} et al.'
    else:
        authors = ', '.join(names)
    venue = fields.get('journal') or fields.get('booktitle') or fields.get('publisher') or fields.get('howpublished', '')
    vol = fields.get('volume', '')
    num = f"({fields['number']})" if fields.get('number') else ''
    pages = f":{fields['pages'].replace('--', '–')}" if fields.get('pages') else ''
    parts = [f"{authors} ({fields.get('year', 'n.d.')}).", f"{fields.get('title', '')}."]
    if venue:
        parts.append(f'{venue} {vol}{num}{pages}'.strip() + '.')
    return ' '.join(p for p in parts if p.strip('. '))


def link(fields):
    if fields.get('doi'):
        return f"https://doi.org/{fields['doi']}"
    return fields.get('url')
