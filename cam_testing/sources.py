"""
Where each calibration and validation data item comes from: the rows of the review pages' "Data sources"
table and the data_sources_test (modules/README.md, "Data sources").

An item's source is, in this order:

  1. its ``data_reference`` key (an obs_data data_item or prediction_item, or a baseline target): a string
     (the citation) or an object with any of the keys in REFERENCE_KEYS;
  2. its section in the SOURCES.md beside the data (the parameterisation's or system model's directory):
     a Markdown heading that is the item's ``data_item_name`` (a trailing " — <value>" or "(note)" on the
     heading is ignored), then optional labelled lines ``- <Label>: <text>`` (LABELS), then free prose.

Value and sigma always come from the data file, never from SOURCES.md. Kind (measured / assumed /
derived) is read only from the labels; role (target / prediction / held out / validation) is the item's
place in the data: a data_item with weight > 0 is a target, a prediction_item with a value is held out,
one without is a prediction, a baseline target or data column is validation. A data_item with weight 0
whose sigma is >= TERM_SIGMA, or that another item's operation_kwargs names, is a term of another item
(x1/x2 of a subtraction or division): its source must say "derived: term of <item>".
"""
import json
import os
import re

SOURCES_MD = 'SOURCES.md'
TERM_SIGMA = 1e10
MISSING = 'source missing'

# data_reference object keys (and SOURCES.md fields)
REFERENCE_KEYS = ('source', 'figure', 'obtained', 'species', 'preparation', 'conditions', 'kind', 'role',
                  'sigma_reason', 'primary_read', 'note')
# SOURCES.md label (lower case) -> field
LABELS = {
    'source': 'source', 'citation': 'source', 'reference': 'source',
    'figure': 'figure', 'table': 'figure', 'figure/table': 'figure',
    'obtained': 'obtained', 'how obtained': 'obtained',
    'species': 'species', 'preparation': 'preparation', 'conditions': 'conditions',
    'species/preparation': 'species', 'species/preparation/conditions': 'species',
    'kind': 'kind', 'role': 'role',
    'sigma': 'sigma_reason', 'sigma reason': 'sigma_reason',
    'primary read': 'primary_read',
    'note': 'note', 'notes': 'note',
}
_HEADING = re.compile(r'^(#{1,6})\s+(.*?)\s*#*\s*$')
_LABELLED = re.compile(r'^\s*[-*]\s+(?:\*\*)?([A-Za-z][A-Za-z /]*?)(?::\*\*|\*\*:|:)\s*(.*)$')
_TERM_OF = re.compile(r'term of\s+`?([^\s`,;:()]+)', re.I)
_DASH_SPLIT = re.compile(r'\s+[—–-]\s+')   # " — ", " – ", " - "


def heading_name(text):
    '''The data_item_name a SOURCES.md heading names: without backticks, a trailing " — <value>" or a
    trailing "(note)".'''
    name = text.replace('`', '').strip()
    name = _DASH_SPLIT.split(name, 1)[0]
    name = re.sub(r'\s*\([^)]*\)\s*$', '', name)
    return name.strip()


def parse_sources_md(text):
    '''{heading name: {field: text, 'prose': text}} of a SOURCES.md, plus the order of its headings.
    Every heading is a section; its labelled lines (LABELS) become fields, the rest prose. A section ends
    at the next heading.'''
    sections, order = {}, []
    current, last_field = None, None
    for line in text.splitlines():
        m = _HEADING.match(line)
        if m:
            name = heading_name(m.group(2))
            current, last_field = sections.setdefault(name, {'prose': []}), None
            current.setdefault('_level', len(m.group(1)))
            order.append(name)
            continue
        if current is None:
            continue
        lm = _LABELLED.match(line)
        field = LABELS.get(lm.group(1).strip().lower()) if lm else None
        if field:
            current[field] = ' '.join(filter(None, [current.get(field), lm.group(2).strip()]))
            last_field = field
        elif last_field and line.startswith((' ', '\t')) and line.strip():
            current[last_field] += ' ' + line.strip()       # a labelled line's continuation
        else:
            last_field = None
            current['prose'].append(line.strip())
    for s in sections.values():
        s['prose'] = re.sub(r'\n{3,}', '\n\n', '\n'.join(s['prose'])).strip()
    return sections, order


def read_sources_md(directory):
    '''(sections, headings, path) of directory/SOURCES.md; ({}, [], None) when there is none.'''
    path = os.path.join(directory, SOURCES_MD)
    if not os.path.isfile(path):
        return {}, [], None
    with open(path, encoding='utf-8') as f:
        sections, order = parse_sources_md(f.read())
    return sections, order, path


def read_reference(ref):
    '''A data_reference as {field: text}: a string is the source; an object's REFERENCE_KEYS are read
    (other keys are kept under note).'''
    if ref is None or ref == '' or ref == {}:
        return {}
    if not isinstance(ref, dict):
        return {'source': str(ref)}
    out, extra = {}, []
    for k, v in ref.items():
        key = LABELS.get(str(k).lower().replace('_', ' '), k) if k not in REFERENCE_KEYS else k
        if key in REFERENCE_KEYS and v not in (None, ''):
            out[key] = v if isinstance(v, str) else json.dumps(v)
        elif v not in (None, ''):
            extra.append(f'{k}: {v}')
    if extra:
        out['note'] = '; '.join(filter(None, [out.get('note'), *extra]))
    return out


def classify_kind(text):
    '''(kind, detail) from a Kind label: measured, assumed, derived, term (detail: the item it is a term
    of) or other; (None, None) without one.'''
    if not text:
        return None, None
    t = text.strip()
    m = _TERM_OF.search(t)
    if m:
        return 'term', m.group(1).rstrip('.')
    low = t.lower()
    for kind in ('measured', 'assumed', 'derived'):
        if low.startswith(kind):
            return kind, t[len(kind):].strip(' :;,-–—()') or None
    return 'other', t


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _referenced_by(items):
    '''{item name: the item whose operation_kwargs names it}.'''
    names = {i.get('data_item_name') for i in items}
    out = {}
    for i in items:
        for v in (i.get('operation_kwargs') or {}).values():
            if isinstance(v, str) and v in names and v != i.get('data_item_name'):
                out.setdefault(v, i.get('data_item_name'))
    return out


def _row(name, *, variable='', operation='', value=None, sigma=None, unit='', role='', comment='',
         reference=None, section=None, term_of=None, structural_term=False):
    ref = read_reference(reference)
    sec = {k: v for k, v in (section or {}).items() if not k.startswith('_')}
    fields = {**{k: v for k, v in sec.items() if k != 'prose'}, **ref}
    origin = 'data_reference' if ref else ('SOURCES.md' if section is not None else None)
    kind, detail = classify_kind(fields.get('kind'))
    problem = None
    if structural_term and kind != 'term':
        problem = ('weight-0 term of ' + (term_of or 'another item') + ': its source must say "derived: term of '
                   + (term_of or '<item>') + '"')
    elif origin is None:
        problem = MISSING
    if kind == 'term':
        kind_label = f'term of {detail}, not a measurement'
    elif structural_term:
        kind_label = f'term of {term_of or "another item"}, not a measurement'
    elif kind == 'other':
        kind_label = detail
    elif kind:
        kind_label = kind + (f' ({detail})' if detail else '')
    else:
        kind_label = 'not stated'
    stated_role = (fields.get('role') or '').strip()
    species = '; '.join(fields[k] for k in ('species', 'preparation', 'conditions') if fields.get(k))
    return {
        'name': name, 'variable': variable, 'operation': operation, 'value': value, 'sigma': sigma, 'unit': unit,
        'source': fields.get('source'), 'figure': fields.get('figure'), 'obtained': fields.get('obtained'),
        'primary_read': fields.get('primary_read'), 'sigma_reason': fields.get('sigma_reason'),
        'note': fields.get('note'), 'prose': sec.get('prose') or '', 'comment': comment or '',
        'species': species, 'kind': 'term' if structural_term else kind, 'kind_label': kind_label,
        'role': role, 'stated_role': stated_role if stated_role and stated_role.lower() != role.lower() else '',
        'origin': origin, 'missing': origin is None, 'problem': problem,
    }


def obs_data_rows(obs, sections):
    '''One row per data_item and prediction_item of an obs_data dict, with sections its SOURCES.md's.'''
    data, preds = obs.get('data_items') or [], obs.get('prediction_items') or []
    referenced = _referenced_by(data + preds)
    rows = []
    for where, items in (('data', data), ('prediction', preds)):
        for it in items:
            name = it.get('data_item_name') or it.get('variable') or '?'
            weight, std = _float(it.get('weight', 1.0)), _float(it.get('std'))
            term = where == 'data' and weight == 0 and ((std is not None and std >= TERM_SIGMA) or name in referenced)
            if where == 'prediction':
                role = 'held out' if it.get('value') is not None else 'prediction'
            elif term:
                role = 'none (term)'
            else:
                role = 'target' if weight is None or weight > 0 else 'held out'
            ops = it.get('operands') or ([it['variable']] if it.get('variable') else [])
            rows.append(_row(name, variable=', '.join(map(str, ops)), operation=it.get('operation') or '',
                             value=it.get('value'), sigma=it.get('std'), unit=it.get('unit') or '', role=role,
                             comment=it.get('comment') or '', reference=it.get('data_reference'),
                             section=sections.get(name), term_of=referenced.get(name), structural_term=term))
    return rows


def baseline_rows(block, sections):
    '''Rows for a validation baseline (validation.<p>.baseline): one per scalar target, or one per data
    column of its CSV (whose source may also be the block's own data_reference).'''
    rows = []
    for tg in block.get('targets') or []:
        rng = tg.get('range')
        rows.append(_row(tg.get('name', '?'), variable=tg.get('expr', ''), operation='expression',
                         value=f'[{rng[0]}, {rng[1]}]' if rng else tg.get('value'), sigma=tg.get('std'),
                         unit=tg.get('units', ''), role='validation', comment=tg.get('note') or '',
                         reference=tg.get('data_reference'), section=sections.get(tg.get('name'))))
    if block.get('data') and block.get('variables'):
        for model_var, column in block['variables'].items():
            rows.append(_row(column, variable=model_var, operation='series', value=f'series ({block["data"]})',
                             role='validation', comment=block.get('source') or '',
                             reference=block.get('data_reference'), section=sections.get(column)))
    return rows


def _active(block):
    return bool(block) and block.get('status', 'pending') not in ('pending', 'skipped', 'not_applicable')


def data_table(directory, obs=None, baseline=None):
    '''The "Data sources" table of one data set: {'rows', 'n_missing', 'n_problems', 'extra_headings',
    'sources_md'} for the obs_data dict obs and/or the validation baseline block, with the SOURCES.md in
    directory. None when there are no items.'''
    sections, order, path = read_sources_md(directory)
    rows = obs_data_rows(obs, sections) if obs else []
    if baseline and _active(baseline):
        rows += baseline_rows(baseline, sections)
    if not rows:
        return None
    names = {r['name'] for r in rows}
    # headings at the level the items use that name no item (e.g. proposed data not in the obs_data yet)
    levels = {sections[n]['_level'] for n in order if n in names}
    extra = [n for n in dict.fromkeys(order) if n not in names and sections[n]['_level'] in levels]
    return {'rows': rows, 'n_missing': sum(r['missing'] for r in rows),
            'n_problems': sum(bool(r['problem']) for r in rows), 'extra_headings': extra,
            'sources_md': path}


def read_obs(path):
    if not path or not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def parameterisation_table(inst):
    '''The data_table of a parameterisation: its obs_data and its validation baseline.'''
    return data_table(inst.dir, obs=read_obs(inst.obs_data_path), baseline=inst.validation.get('baseline'))


def version_tables(version):
    '''[(parameterisation, table)] for each parameterisation with data items.'''
    out = []
    for inst in version.parameterisations():
        t = parameterisation_table(inst)
        if t:
            out.append((inst, t))
    return out


def version_problems(version):
    '''Every unsourced data item of the version: ["<parameterisation>: <item>: <problem> (add ...)"].'''
    problems = []
    for inst, t in version_tables(version):
        rel = os.path.relpath(inst.dir, version.dir)
        for r in t['rows']:
            if r['problem']:
                problems.append(f'{inst.name}: {r["name"]}: {r["problem"]} (add a data_reference to the item, or '
                                f'a "## {r["name"]}" section to {rel}/{SOURCES_MD})')
    return problems


def check(version):
    '''The data_sources_test result: passed when every calibration and validation data item of the version
    has a source (a data_reference or its SOURCES.md section), weight-0 terms saying "derived: term of
    <item>"; not applicable without data items.'''
    from cam_testing import checks
    tables = version_tables(version)
    n = sum(len(t['rows']) for _, t in tables)
    extra = [f'{i.name}: SOURCES.md heading "{h}" names no data item' for i, t in tables for h in t['extra_headings']]
    if not n:
        return checks.Result('data_sources_test', checks.NOT_APPLICABLE, 'no calibration or validation data items')
    problems = version_problems(version)
    if problems:
        return checks.Result('data_sources_test', checks.FAILED,
                             f'{len(problems)} of {n} data items have no source', {'n_items': n, 'n_unsourced': len(problems)},
                             details=problems + extra)
    return checks.Result('data_sources_test', checks.PASSED, f'all {n} data items have a source',
                         {'n_items': n, 'n_unsourced': 0}, details=extra)

