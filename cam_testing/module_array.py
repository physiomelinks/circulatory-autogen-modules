"""
Module arrays: <prefix>_module_array.json, a JSON list of instance records in PhLynx's key names

    {"name": "venous_svc", "module_type": "venous", "module_subtype": "vp",
     "inp_instances": ["systemic_T"], "out_instances": ["heart", "volume_sum"]}

(module_subtype is the version). Each record is an instance: one use of a version in the
network. It may name the version's parameterisation it uses, "parameterisation": "default" (its
parameters are the defaults; the model's own parameters file wins; several instances can use one
parameterisation). A supermodule instance has no inp/out lists of its own; it gives
per_submodule_inputs / per_submodule_outputs instead (see modules/README.md). libcuflynx reads
these files, and the CSV layout with the same columns (name, module_subtype, module_type,
inp_instances, out_instances, lists space-separated), which read_records here reads as the same
records.
"""
import csv
import json
import os

LIST_KEYS = ('inp_instances', 'out_instances')
KEY_ORDER = ('name', 'module_type', 'module_subtype', 'parameterisation', 'inp_instances', 'out_instances')
PARAMETERISATION_KEY = 'parameterisation'


def parameterisation_of(record, default=None):
    '''The parameterisation a module-array record or supermodule submodule names, else ``default``.'''
    return record.get(PARAMETERISATION_KEY) or default


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, str):
        return v.split()
    return [str(x) for x in v]


def normalise_record(rec):
    '''A record (or a CSV row) -> lists as lists, key order kept tidy.'''
    out = {}
    for k, v in rec.items():
        k = k.strip() if isinstance(k, str) else k
        out[k] = v.strip() if isinstance(v, str) and k not in LIST_KEYS else v
    for k in LIST_KEYS:
        if k in out:
            out[k] = _as_list(out[k])
    if 'per_submodule_inputs' not in out and 'per_submodule_outputs' not in out:
        for k in LIST_KEYS:
            out.setdefault(k, [])
    head = {k: out[k] for k in KEY_ORDER if k in out}
    return {**head, **{k: v for k, v in out.items() if k not in head}}


def read_records(path):
    '''A module array (.json or .csv) as a list of records.'''
    if path.endswith('.json'):
        with open(path) as f:
            data = json.load(f)
        if not isinstance(data, list):
            raise ValueError(f'{path}: a module array is a JSON list of instance records')
        return [normalise_record(r) for r in data]
    with open(path, newline='') as f:
        rows = [{(k or '').strip(): (v or '').strip() for k, v in r.items()}
                for r in csv.DictReader(f, skipinitialspace=True)]
    return [normalise_record(r) for r in rows if r.get('name') and not r['name'].startswith('#')]


def write_records(path, records):
    '''Writes records as <...>_module_array.json, one record per line.'''
    recs = [normalise_record(r) for r in records]
    with open(path, 'w') as f:
        f.write('[\n' + ',\n'.join(' ' + json.dumps(r) for r in recs) + '\n]\n')
    return path


def find(directory, prefix):
    '''<prefix>_module_array.json, else .csv; None if neither exists.'''
    for ext in ('.json', '.csv'):
        p = os.path.join(directory, f'{prefix}_module_array{ext}')
        if os.path.isfile(p):
            return p
    return None


def harness_rows(harness):
    '''A verification_config harness's network rows (its "module_array").'''
    return (harness or {}).get('module_array')


def from_rows(rows):
    '''Spec harness rows [name, module_subtype, module_type, inp, out(, parameterisation)] (inp/out
    space-separated) -> records. A row may instead be a module-array record (a dict), e.g. a
    supermodule instance with per_submodule_inputs / per_submodule_outputs.'''
    out = []
    for r in rows:
        if isinstance(r, dict):
            out.append(normalise_record(r))
            continue
        n, bc, vt, i, o = r[:5]
        rec = {'name': str(n), 'module_type': str(vt), 'module_subtype': str(bc)}
        if len(r) > 5 and r[5]:
            rec[PARAMETERISATION_KEY] = str(r[5])
        rec.update({'inp_instances': str(i or ''), 'out_instances': str(o or '')})
        out.append(normalise_record(rec))
    return out


def to_library_versions(records, parameterisation='default'):
    '''Records naming circulatory_autogen's (or this library's pre-versions) module pairs -> this
    library's (module_type, version), renamed where the move to versions changed them
    (library.legacy_renames), each naming ``parameterisation`` unless it names one already.'''
    from cam_testing.library import legacy_renames
    table = legacy_renames()
    out = []
    for rec in records:
        rec = normalise_record(rec)
        new = table.get((rec.get('module_type'), rec.get('module_subtype')))
        if new:
            rec['module_type'], rec['module_subtype'] = new
        if not parameterisation_of(rec):
            rec[PARAMETERISATION_KEY] = parameterisation
        out.append(normalise_record(rec))
    return out

