"""
Per-version BibTeX: versions/<version>/<module_type>_<version>_references.bib holds the sources
cited by the version's parameters (all its instances). A parameter's data_reference is "<bibkey>; <note>" (the note is optional), or
starts with "definitional" / names a test fixture when there is nothing to cite.

A small reader, enough for the entries this repo writes (@type{key, field = {...}, ...}).
"""
import os
import re
import unicodedata

_ENTRY_HEAD_RE = re.compile(r'@(\w+)\s*\{\s*([^,\s]+)\s*,')
_FIELD_NAME_RE = re.compile(r'(\w+)\s*=\s*')


def bib_path(version):
    return version.path('references.bib')


def proposed_bib_path(version):
    '''Citations proposed for review, not yet confirmed (<module_type>_<version>_references_proposed.bib).'''
    return version.path('references_proposed.bib')


def _balanced(text, i):
    '''Index just past the brace that closes the one opened before text[i].'''
    depth = 1
    while i < len(text) and depth:
        depth += {'{': 1, '}': -1}.get(text[i], 0)
        i += 1
    return i


def entries(text):
    '''[(kind, key, start, end, raw fields)] for every entry in a .bib text (start/end: the entry's
    span, '@' to its closing brace). Brace-balanced, so values like {Bing{\\"u}l} parse.'''
    out = []
    for m in _ENTRY_HEAD_RE.finditer(text):
        end = _balanced(text, m.end())
        body = text[m.end():end - 1]
        fields, j = {}, 0
        while True:
            f = _FIELD_NAME_RE.search(body, j)
            if not f:
                break
            k = f.end()
            if k < len(body) and body[k] == '{':
                stop = _balanced(body, k + 1)
                value = body[k + 1:stop - 1]
            elif k < len(body) and body[k] == '"':
                stop = body.find('"', k + 1)
                stop = len(body) if stop < 0 else stop + 1
                value = body[k + 1:stop - 1]
            else:
                stop = k + len(re.match(r'[^,\n]*', body[k:])[0])
                value = body[k:stop].strip()
            fields.setdefault(f[1].lower(), value)
            j = stop
        out.append((m[1], m[2], m.start(), end, fields))
    return out


def read(path):
    '''{key: {'type': ..., field: value}} (braces and surrounding quotes removed).'''
    if not os.path.isfile(path):
        return {}
    with open(path, encoding='utf-8') as f:
        text = f.read()
    out = {}
    for kind, key, _, _, raw in entries(text):
        fields = {'type': kind.lower()}
        for name, value in raw.items():
            fields[name] = re.sub(r'\s+', ' ', value.replace('{', '').replace('}', '')).strip()
        out[key] = fields
    return out


# ---------------------------------------------------------------------------------------------- keys
# The key convention (2026-10-06): <first author's surname><year><first significant title word>, all
# lowercase ASCII (belluzzi1986quantitative). Particles join the surname ("van der Pol" -> vanderpol,
# "Hernandez-Cruz" -> hernandezcruz); a missing author or year is left out (web pages without a date).
# Two different papers with the same key take further title words (tools/rekey_bib.py).
STOPWORDS = {'a', 'an', 'the', 'on', 'of', 'in', 'for', 'to', 'and', 'at', 'by', 'with', 'from', 'is', 'are'}
KEY_RE = re.compile(r'^[a-z]+[0-9]{4}[a-z]+$')


def _ascii(s):
    s = re.sub(r'\\[\'"`^~=.uvHcdbtk]\s*\{?([A-Za-z])\}?', r'\1', s)      # \'{a} -> a, \v{z} -> z
    s = s.replace('{', '').replace('}', '').replace('\\', '')
    return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode()


def first_surname(author):
    first = re.split(r'\s+and\s+', (author or '').strip())[0]
    if first.startswith('{') and first.endswith('}'):
        sur = first                                  # a corporate author, {U.S. National Library of Medicine}
    elif ',' in first:
        sur = first.split(',')[0]
    else:
        parts = first.split()
        sur = parts[-1] if parts else first
    return re.sub(r'[^a-z]', '', _ascii(sur).lower())


def title_words(title):
    '''The significant title words (stopwords skipped), lowercase ASCII.'''
    words = []
    for w in _ascii(title or '').split():
        m = re.match(r'[`\'"]*([A-Za-z]+)', w)       # opening quotes skipped; a word in brackets ([Na]) is not a word
        if m and m[1].lower() not in STOPWORDS:
            words.append(m[1].lower())
    return words


def convention_key(fields, n_words=1):
    '''The entry's key under the convention, from its raw or read fields; n_words > 1 adds further
    title words (to tell two papers with the same key apart).'''
    year = re.sub(r'[^0-9]', '', fields.get('year') or '')[:4]
    return first_surname(fields.get('author')) + year + ''.join(title_words(fields.get('title'))[:n_words])


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
