"""
Where each calibration and validation data item comes from: the rows of the review pages' "Data sources"
table and the data_sources_test (modules/README.md, "Data sources").

An item's source combines:

  1. its ``data_reference`` key (an obs_data data_item or prediction_item, or a baseline target): a string
     (``<BibTeX key>; <short note>``) or an object with any of the keys in REFERENCE_KEYS;
  2. its section in the SOURCES.md beside the data (the parameterisation's or system model's directory):
     a Markdown heading that is the item's ``data_item_name`` (a trailing " — <value>" or "(note)" on the
     heading is ignored), then optional labelled lines ``- <Label>: <text>`` (LABELS), then free prose.

The section extends the data_reference: a field in both is taken from the data_reference. When an item has
both, the section must have a ``- Data reference: <text>`` line repeating the item's data_reference string
exactly (for an object, its ``source``); a missing or different line is a problem. A section with nothing
under its heading is not a source.

Value and sigma always come from the data file, never from SOURCES.md. Kind (measured / assumed /
derived / synthetic) is read only from the labels; role (target / prediction / held out / term / validation) is the
item's place in the data: a data_item with weight > 0 is a target, a prediction_item with a value is held
out, one without is a prediction, a data_item with weight 0 is a term (never a target), a baseline target or
data column is validation. Weight 0 is the only way an item is used without being scored (a large sigma
never marks one): a weight-0 data_item is a term of another item, usually one whose operation_kwargs names
it (x1/x2 of a subtraction or division), and its source must say "derived: term of <item>", naming that item
when one does.

Each problem has a kind (PROBLEMS): no source, a Data reference line missing or different, or a term not labelled.
"""
import json
import os
import re

SOURCES_MD = 'SOURCES.md'
MISSING = 'source missing'
DATA_REFERENCE = '_data_reference'   # a section's "Data reference" line as written
# problem kind -> how a summary counts it
PROBLEMS = {'missing': 'no source', 'reference': 'Data reference mismatch/missing', 'term': 'term not labelled'}
KINDS = ('measured', 'assumed', 'derived', 'synthetic')

# data_reference object keys (and SOURCES.md fields)
REFERENCE_KEYS = ('source', 'figure', 'obtained', 'species', 'preparation', 'conditions', 'kind', 'role',
                  'sigma_reason', 'primary_read', 'note')
# SOURCES.md label (lower case) -> field
LABELS = {
    'source': 'source', 'citation': 'source', 'reference': 'source', 'data reference': 'source',
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
_DERIVED_TERM = re.compile(r'^\s*derived\s*:\s*term of\b', re.I)
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
    current, last_fields = None, ()
    for line in text.splitlines():
        m = _HEADING.match(line)
        if m:
            name = heading_name(m.group(2))
            current, last_fields = sections.setdefault(name, {'prose': []}), ()
            current.setdefault('_level', len(m.group(1)))
            order.append(name)
            continue
        if current is None:
            continue
        lm = _LABELLED.match(line)
        label = lm.group(1).strip().lower() if lm else None
        field = LABELS.get(label)
        if field:
            # a Data reference line is also kept as written, to compare with the data file's data_reference
            last_fields = (field, DATA_REFERENCE) if label == 'data reference' else (field,)
            for f in last_fields:
                current[f] = ' '.join(filter(None, [current.get(f), lm.group(2).strip()]))
        elif last_fields and line.startswith((' ', '\t')) and line.strip():
            for f in last_fields:
                current[f] += ' ' + line.strip()            # a labelled line's continuation
        else:
            last_fields = ()
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


def _detail(rest):
    '''A Kind's detail: the text after the kind, without separators or the parentheses around all of it.'''
    rest = rest.strip(' :;,-–—')
    if rest.startswith('(') and rest.endswith(')'):
        depth = 0
        for i, c in enumerate(rest):
            depth += {'(': 1, ')': -1}.get(c, 0)
            if depth == 0:
                break
        if i == len(rest) - 1:          # the first parenthesis closes at the end
            rest = rest[1:-1].strip()
    return rest or None


def classify_kind(text):
    '''(kind, detail) from a Kind label: one of KINDS, term (detail: the item it is a term of) or other;
    (None, None) without one.'''
    if not text:
        return None, None
    t = text.strip()
    m = _TERM_OF.search(t)
    if m:
        return 'term', m.group(1).rstrip('.')
    low = t.lower()
    for kind in KINDS:
        if low.startswith(kind):
            return kind, _detail(t[len(kind):])
    return 'other', t


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _referenced_by(items):
    '''{item name: [the items whose operation_kwargs name it]}.'''
    names = {i.get('data_item_name') for i in items}
    out = {}
    for i in items:
        for v in (i.get('operation_kwargs') or {}).values():
            if isinstance(v, str) and v in names and v != i.get('data_item_name'):
                users = out.setdefault(v, [])
                if i.get('data_item_name') not in users:
                    users.append(i.get('data_item_name'))
    return out


def _reference_string(reference):
    '''The data_reference string a SOURCES.md "Data reference" line repeats: the string, or an object's source.'''
    if isinstance(reference, dict):
        reference = reference.get('source')
    return str(reference).strip() if reference not in (None, '') else None


def _term_problem(kind_text, kind, detail, users):
    '''The problem with a weight-0 term's Kind, or None: it must say "derived: term of <item>", the item one
    whose operation_kwargs name it (users) when there is one.'''
    term_of = users[0] if users else None
    want = 'derived: term of ' + (term_of or '<item>')
    where = 'weight-0 term of ' + (' / '.join(users) if users else 'another item')
    if kind != 'term':
        return f'{where}: its source must say "{want}"'
    if not _DERIVED_TERM.match(kind_text):
        return f'{where}: its Kind "{kind_text.strip()}" must say "{want}"'
    if users and detail not in users:
        return (f'{where}: its Kind names {detail}, but {" / ".join(users)} uses it (operation_kwargs): '
                f'it must say "{want}"')
    return None


def _row(name, *, variable='', operation='', value=None, sigma=None, unit='', role='', comment='',
         reference=None, section=None, term_users=(), structural_term=False, data_file=None, sources_md=None):
    ref = read_reference(reference)
    sec = {k: v for k, v in (section or {}).items() if not k.startswith('_')}
    fields = {**{k: v for k, v in sec.items() if k != 'prose'}, **ref}
    from_md = section is not None and any(sec.values())     # a heading with nothing under it is no source
    origin = ' + '.join(o for o, used in (('data_reference', bool(ref)), ('SOURCES.md', from_md)) if used) or None
    term_of = term_users[0] if term_users else None
    kind_text = fields.get('kind') or ''
    kind, detail = classify_kind(kind_text)
    in_data, in_md = _reference_string(reference), (section or {}).get(DATA_REFERENCE)
    mismatch = bool(in_data and in_md and in_data != in_md.strip())
    files = f'{data_file or "the data file"}', f'{sources_md or SOURCES_MD}'
    problem, problem_kind = None, None
    if mismatch:
        problem, problem_kind = (f'data_reference "{in_data}" in {files[0]} differs from the Data reference '
                                 f'"{in_md.strip()}" in {files[1]}: they must be identical'), 'reference'
    elif in_data and section is not None and in_md is None:
        problem, problem_kind = (f'data_reference "{in_data}" in {files[0]} has no Data reference line in its '
                                 f'section of {files[1]}: add "- Data reference: {in_data}"'), 'reference'
    elif structural_term and _term_problem(kind_text, kind, detail, list(term_users)):
        problem, problem_kind = _term_problem(kind_text, kind, detail, list(term_users)), 'term'
    elif origin is None:
        problem, problem_kind = MISSING, 'missing'
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
        'origin': origin, 'missing': origin is None, 'problem': problem, 'problem_kind': problem_kind,
        'mismatch': mismatch,
    }


def obs_data_rows(obs, sections, data_file=None, sources_md=None):
    '''One row per data_item and prediction_item of an obs_data dict, with sections its SOURCES.md's
    (data_file and sources_md: the two files as problems name them).'''
    data, preds = obs.get('data_items') or [], obs.get('prediction_items') or []
    referenced = _referenced_by(data + preds)
    rows = []
    for where, items in (('data', data), ('prediction', preds)):
        for it in items:
            name = it.get('data_item_name') or it.get('variable') or '?'
            weight = _float(it.get('weight', 1.0))
            term = where == 'data' and weight == 0
            if where == 'prediction':
                role = 'held out' if it.get('value') is not None else 'prediction'
            else:
                role = 'term' if term else 'target'
            ops = it.get('operands') or ([it['variable']] if it.get('variable') else [])
            rows.append(_row(name, variable=', '.join(map(str, ops)), operation=it.get('operation') or '',
                             value=it.get('value'), sigma=it.get('std'), unit=it.get('unit') or '', role=role,
                             comment=it.get('comment') or '', reference=it.get('data_reference'),
                             section=sections.get(name), term_users=referenced.get(name, ()), structural_term=term,
                             data_file=data_file, sources_md=sources_md))
    return rows


def baseline_rows(block, sections, data_file=None, sources_md=None):
    '''Rows for a validation baseline (validation.<p>.baseline): one per scalar target, or one per data
    column of its CSV (whose source may also be the block's own data_reference).'''
    files = {'data_file': data_file, 'sources_md': sources_md}
    rows = []
    for tg in block.get('targets') or []:
        rng = tg.get('range')
        rows.append(_row(tg.get('name', '?'), variable=tg.get('expr', ''), operation='expression',
                         value=f'[{rng[0]}, {rng[1]}]' if rng else tg.get('value'), sigma=tg.get('std'),
                         unit=tg.get('units', ''), role='validation', comment=tg.get('note') or '',
                         reference=tg.get('data_reference'), section=sections.get(tg.get('name')), **files))
    if block.get('data') and block.get('variables'):
        for model_var, column in block['variables'].items():
            rows.append(_row(column, variable=model_var, operation='series', value=f'series ({block["data"]})',
                             role='validation', comment=block.get('source') or '',
                             reference=block.get('data_reference'), section=sections.get(column), **files))
    return rows


def problem_counts(rows):
    '''{problem kind: number of rows with it}, over PROBLEMS.'''
    return {k: sum(r['problem_kind'] == k for r in rows) for k in PROBLEMS}


def add_counts(counts):
    '''The sum of several problem_counts (any iterable).'''
    counts = list(counts)
    return {k: sum((c or {}).get(k, 0) for c in counts) for k in PROBLEMS}


def summary(counts):
    '''"1 no source, 2 term not labelled": the problems of problem_counts, each kind counted.'''
    return ', '.join(f'{counts[k]} {label}' for k, label in PROBLEMS.items() if counts.get(k))


def _active(block):
    return bool(block) and block.get('status', 'pending') not in ('pending', 'skipped', 'not_applicable')


def data_table(directory, obs=None, baseline=None, obs_path=None, baseline_path=None, base_dir=None):
    '''The "Data sources" table of one data set: {'rows', 'n_missing', 'n_problems', 'problems' (problem_counts),
    'extra_headings', 'sources_md'} for the obs_data dict obs and/or the validation baseline block, with the SOURCES.md in
    directory. obs_path and baseline_path are the files they were read from, named (relative to base_dir,
    default directory) in a data reference mismatch. None when there are no items.'''
    sections, order, path = read_sources_md(directory)
    base = base_dir or directory

    def shown(p):
        return os.path.relpath(p, base).replace(os.sep, '/') if p else None
    md = shown(path or os.path.join(directory, SOURCES_MD))
    rows = obs_data_rows(obs, sections, shown(obs_path), md) if obs else []
    if baseline and _active(baseline):
        rows += baseline_rows(baseline, sections, shown(baseline_path), md)
    if not rows:
        return None
    names = {r['name'] for r in rows}
    # headings at the level the items use that name no item (e.g. proposed data not in the obs_data yet)
    levels = {sections[n]['_level'] for n in order if n in names}
    extra = [n for n in dict.fromkeys(order) if n not in names and sections[n]['_level'] in levels]
    return {'rows': rows, 'n_missing': sum(r['missing'] for r in rows),
            'n_problems': sum(bool(r['problem']) for r in rows), 'problems': problem_counts(rows),
            'extra_headings': extra,
            'sources_md': path}


def read_obs(path):
    if not path or not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def parameterisation_table(inst):
    '''The data_table of a parameterisation: its obs_data and its validation baseline.'''
    return data_table(inst.dir, obs=read_obs(inst.obs_data_path), baseline=inst.validation.get('baseline'),
                      obs_path=inst.obs_data_path, baseline_path=inst.version.verification_config_path,
                      base_dir=inst.version.dir)


def version_tables(version):
    '''[(parameterisation, table)] for each parameterisation with data items.'''
    out = []
    for inst in version.parameterisations():
        t = parameterisation_table(inst)
        if t:
            out.append((inst, t))
    return out


def table_problems(table, prefix, md_dir):
    '''["<prefix><item>: <problem>"] of a data_table; an unsourced item with where to add its source
    (md_dir: the SOURCES.md's directory as shown).'''
    out = []
    for r in table['rows']:
        if r['problem_kind'] == 'missing':
            out.append(f'{prefix}{r["name"]}: {r["problem"]} (add a data_reference to the item, or a '
                       f'"## {r["name"]}" section to {md_dir}{SOURCES_MD})')
        elif r['problem']:
            out.append(f'{prefix}{r["name"]}: {r["problem"]}')
    return out


def version_problems(version):
    '''Every problem of the version's data items: ["<parameterisation>: <item>: <problem>"].'''
    problems = []
    for inst, t in version_tables(version):
        rel = os.path.relpath(inst.dir, version.dir).replace(os.sep, '/')
        problems += table_problems(t, f'{inst.name}: ', f'{rel}/')
    return problems


def _result(tables, problems, extra):
    '''The data_sources_test result of data tables, their problems and extra-heading notes.'''
    from cam_testing import checks
    n = sum(len(t['rows']) for t in tables)
    if not n:
        return checks.Result('data_sources_test', checks.NOT_APPLICABLE, 'no calibration or validation data items')
    counts = add_counts(t['problems'] for t in tables)
    metrics = {'n_items': n, 'n_problems': len(problems), **{f'n_{k}': v for k, v in counts.items()}}
    if problems:
        return checks.Result('data_sources_test', checks.FAILED,
                             f'{len(problems)} of {n} data items have a problem: {summary(counts)}', metrics,
                             details=problems + extra)
    return checks.Result('data_sources_test', checks.PASSED, f'all {n} data items have a source', metrics,
                         details=extra)


def check(version):
    '''The data_sources_test result: passed when every calibration and validation data item of the version
    has a source (a data_reference and/or its SOURCES.md section, their data reference strings identical),
    weight-0 terms saying "derived: term of <item>"; not applicable without data items.'''
    tables = version_tables(version)
    extra = [f'{i.name}: SOURCES.md heading "{h}" names no data item' for i, t in tables for h in t['extra_headings']]
    return _result([t for _, t in tables], version_problems(version), extra)


def system_table(system):
    '''The data_table of a system model's obs_data (<model>_obs_data.json beside it).'''
    obs_path = os.path.join(system.dir, f'{system.name}_obs_data.json')
    return data_table(system.dir, obs=read_obs(obs_path), obs_path=obs_path)


def check_system(system):
    '''data_sources_test's result for a system model's obs_data (system_data_sources_test).'''
    t = system_table(system)
    if not t:
        return _result([], [], [])
    extra = [f'SOURCES.md heading "{h}" names no data item' for h in t['extra_headings']]
    return _result([t], table_problems(t, '', ''), extra)

