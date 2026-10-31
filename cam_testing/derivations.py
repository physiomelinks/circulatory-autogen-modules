"""
Parameter derivations of a parameterisation or system model: the checks the data_sources_test (a version)
and system_derivations_test (a system model) make, and the "Derived parameters" table of their pages
(modules/README.md, "Derived parameters").

A ``<name>_derivations.json`` beside ``<name>_parameters.json`` says how parameters are calculated from
experimental data (libcuflynx.derivations; its inputs are obs_data ``input_items``, ``data_items`` or other
parameters). Its checks: every input resolves, every CellML derivation model validates, and every parameter
a derivation sets holds the derived value as its cache (libcuflynx's ``check_derivations``, rtol
``CACHE_RTOL``). The input_items' sources are checked with the data items' (cam_testing.sources). A
reviewed version's or system model's parameter whose data_reference starts with "derived" but names no
``derivation`` is a warning, not a failure.

Without a libcuflynx that has ``libcuflynx.derivations``, a derivations file is not checked: the check is
skipped and says why.
"""
import inspect
import json
import math
import os

DERIVATIONS_SUFFIX = '_derivations.json'
PARAMETERS_SUFFIX = '_parameters.json'
DERIVED_PREFIX = 'derived'


def unavailable():
    '''Why derivations cannot be checked with the installed libcuflynx, or None when they can.'''
    try:
        import libcuflynx.derivations  # noqa: F401
    except ImportError as e:
        return (f'libcuflynx has no libcuflynx.derivations ({e}); derivations need a libcuflynx with '
                f'physiomelinks/circulatory_autogen#621')
    return None


def derivations_path(parameters_path):
    '''``<name>_derivations.json`` beside ``<name>_parameters.json``.'''
    base = str(parameters_path)
    stem = base[:-len(PARAMETERS_SUFFIX)] if base.endswith(PARAMETERS_SUFFIX) else os.path.splitext(base)[0]
    return stem + DERIVATIONS_SUFFIX


def _obs(obs_data_path):
    return [obs_data_path] if obs_data_path and os.path.isfile(obs_data_path) else []


def _error(e):
    return f'{type(e).__name__}: {e}'


def problems(parameters_path, obs_data_path):
    '''(results, problems) of the derivations file beside ``parameters_path``, its inputs read from the obs_data
    at ``obs_data_path``: the libcuflynx DerivationResults (None when they could not be evaluated) and every
    input that does not resolve, CellML model that does not validate, and cached value that is not the derived
    one.'''
    from libcuflynx import derivations as dv
    from libcuflynx.derivations.cellml_method import load_derivation_model
    path = derivations_path(parameters_path)
    try:
        found = dv.load_derivations(path)
        doc = dv.read_parameters(parameters_path)
        obs_items = dv.load_obs_items(_obs(obs_data_path))
    except Exception as e:
        return None, [_error(e)]
    out = []
    for d in found:
        try:
            dv.resolve_inputs(d, doc, obs_items)
        except Exception as e:
            out.append(_error(e))
        if d.method == 'cellml':
            try:
                load_derivation_model(d.files()[0])
            except Exception as e:
                out.append(f"derivation '{d.name}': {_error(e)}")
    if out:
        return None, out
    try:
        results = dv.evaluate_derivations(found, doc, _obs(obs_data_path))
    except Exception as e:
        return None, [_error(e)]
    return results, dv.check_derivations(results, doc)


def _leaves(node, path=()):
    for key, child in node.items():
        if isinstance(child, dict):
            if 'value' in child:
                yield path + (key,), child
            else:
                yield from _leaves(child, path + (key,))


def derived_without_derivation(parameters_path):
    '''The parameters of a parameters file whose data_reference starts with "derived" but that name no
    ``derivation``, as "instance/variable" paths.'''
    if not os.path.isfile(parameters_path):
        return []
    with open(parameters_path, encoding='utf-8-sig') as f:
        doc = json.load(f)
    return ['/'.join(p) for p, leaf in _leaves(doc)
            if isinstance(leaf.get('data_reference'), str)
            and leaf['data_reference'].strip().lower().startswith(DERIVED_PREFIX) and not leaf.get('derivation')]


def check(sets, reviewed):
    '''{'problems', 'warnings', 'n_files', 'skipped'} over ``sets``: [(label, parameters_path, obs_data_path)],
    each label prefixing its messages. ``skipped`` says why derivations files were not checked (None when they
    were, or there are none); warnings only when ``reviewed``.'''
    out = {'problems': [], 'warnings': [], 'n_files': 0, 'skipped': None}
    for label, parameters_path, obs_data_path in sets:
        if reviewed:
            out['warnings'] += [f'{label}: {p}: data_reference says "derived" but it names no derivation (add it to '
                                f'{os.path.basename(derivations_path(parameters_path))})'
                                for p in derived_without_derivation(parameters_path)]
        if not os.path.isfile(derivations_path(parameters_path)):
            continue
        out['n_files'] += 1
        out['skipped'] = out['skipped'] or unavailable()
        if out['skipped']:
            continue
        out['problems'] += [f'{label}: {p}' for p in problems(parameters_path, obs_data_path)[1]]
    return out


def version_sets(version):
    '''[(parameterisation name, parameters path, obs_data path)] of a version.'''
    return [(i.name, i.parameters_path, i.obs_data_path) for i in version.parameterisations()]


def system_sets(system):
    '''[(model name, parameters path, obs_data path)] of a system model.'''
    return [(system.name, os.path.join(system.dir, f'{system.name}{PARAMETERS_SUFFIX}'),
             os.path.join(system.dir, f'{system.name}_obs_data.json'))]


def system_check(system):
    '''The system_derivations_test result of a system model, as (status, message, details): its derivations
    check and its input_items have sources (cam_testing.sources); not_applicable without a derivations file,
    skipped without a libcuflynx that has derivations. Warnings (a reviewed model's "derived" parameter
    without a derivation) are in details and never fail it.'''
    from cam_testing import sources
    der = check(system_sets(system), bool(system.spec.get('reviewed')))
    warns = [f'warning: {w}' for w in der['warnings']]
    if not der['n_files']:
        return 'not_applicable', f'no {system.name}{DERIVATIONS_SUFFIX}', warns
    if der['skipped']:
        return 'skipped', f'derivations not checked: {der["skipped"]}', warns
    table = sources.data_table(system.dir, obs=sources.read_obs(system_sets(system)[0][2])) or {'rows': []}
    unsourced = [f'{r["name"]}: input_item {r["problem"]} (add a data_reference to the item, or a '
                 f'"## {r["name"]}" section to {sources.SOURCES_MD})'
                 for r in table['rows'] if r['role'] == 'input' and r['problem']]
    fails = unsourced + der['problems']
    if fails:
        return 'failed', f'{len(fails)} derivation problems', fails + warns
    return 'passed', 'every derivation input resolves and is sourced; every derived value matches its cache', warns


def _python_doc(spec, base_dir):
    from libcuflynx.derivations.python_method import derivation_function
    try:
        doc = inspect.getdoc(derivation_function(spec, base_dir)) or ''
    except Exception as e:
        return _error(e)
    return doc.strip().splitlines()[0] if doc.strip() else ''


def _input_source(record, source_rows):
    '''(source text, missing) of a resolved input: an obs_data item's from its Data sources row, a parameter's
    from its data_reference.'''
    if record['kind'] == 'parameter':
        ref = record.get('data_reference')
        return (str(ref) if ref else 'no data_reference'), not ref
    row = source_rows.get(record['reference'])
    if row is None or row['missing']:
        return 'source missing', True
    return row['source'] or row['kind_label'], False


def table(parameters_path, obs_data_path, source_table, base_dir):
    '''The "Derived parameters" table of the derivations file beside ``parameters_path``: {'file', 'rows',
    'problems', 'skipped'}, a row per derivation with its equation (LaTeX of its CellML model) or Python
    function (name and first docstring line: not declarative), its inputs (value, std, units, source) and its
    outputs (derived and cached value). ``source_table``: the Data sources table (cam_testing.sources) of the
    same obs_data, for the inputs' sources. None without a derivations file.'''
    path = derivations_path(parameters_path)
    if not os.path.isfile(path):
        return None
    out = {'file': os.path.relpath(path, base_dir).replace(os.sep, '/'), 'rows': [], 'problems': [],
           'skipped': unavailable()}
    if out['skipped']:
        return out
    from libcuflynx import derivations as dv
    from cam_testing.mathml import model_equations
    results, out['problems'] = problems(parameters_path, obs_data_path)
    by_name = {r.derivation.name: r for r in results or []}
    rows_by_name = {r['name']: r for r in (source_table or {}).get('rows') or []}
    try:
        found = dv.load_derivations(path)
        doc = dv.read_parameters(parameters_path)
        obs_items = dv.load_obs_items(_obs(obs_data_path))
    except Exception:
        return out
    for d in found:
        row = {'name': d.name, 'method': d.method, 'spec': d.spec, 'comment': d.comment, 'equations': [],
               'function': None, 'doc': None, 'inputs': [], 'outputs': [], 'error': None}
        if d.method == 'cellml':
            try:
                row['equations'] = model_equations(d.files()[0])[0]
            except Exception as e:
                row['error'] = _error(e)
        else:
            row['function'] = d.spec.rsplit(':', 1)[-1]
            row['doc'] = _python_doc(d.spec, d.base_dir)
        result = by_name.get(d.name)
        if result is not None:
            inputs = result.inputs
        else:
            try:
                inputs = dv.resolve_inputs(d, doc, obs_items)
            except Exception as e:
                inputs, row['error'] = {}, row['error'] or _error(e)
        for symbol, rec in inputs.items():
            source, missing = _input_source(rec, rows_by_name)
            row['inputs'].append({'symbol': symbol, 'kind': rec['kind'], 'reference': rec['reference'],
                                  'value': rec['value'], 'std': rec['std'], 'units': rec['units'],
                                  'source': source, 'missing': missing})
        for symbol, name in d.outputs.items():
            o = (result.outputs.get(symbol) if result else None) or {}
            try:
                cached = float(o.get('cached'))
            except (TypeError, ValueError):
                cached = None
            value = o.get('value')
            row['outputs'].append({'symbol': symbol, 'parameter': name, 'value': value, 'units': o.get('units'),
                                   'cached': o.get('cached'),
                                   'stale': value is not None and (cached is None or not math.isclose(
                                       cached, value, rel_tol=dv.CACHE_RTOL, abs_tol=0.0))})
        out['rows'].append(row)
    return out
