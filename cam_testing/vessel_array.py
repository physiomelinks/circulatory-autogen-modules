"""
Vessel arrays: <prefix>_vessel_array.json, a JSON list of instance records in PhLynx's key names

    {"name": "venous_svc", "module_type": "venous", "module_subtype": "vp",
     "inp_instances": ["systemic_T"], "out_instances": ["heart", "volume_sum"]}

(module_type is libcuflynx's vessel_type, module_subtype its BC_type: the version). A record may
name an instance of that version, "instance": "default" (its parameters are the defaults; the
model's own parameters file wins). A supermodule instance has no inp/out lists of its own; it
gives per_submodule_inputs / per_submodule_outputs instead (see modules/README.md). libcuflynx reads these files, and still reads the older CSV layout
(name, BC_type, vessel_type, inp_vessels, out_vessels with space-separated lists), which
read_records here converts to the same records.
"""
import csv
import json
import os

# libcuflynx name -> PhLynx name
TO_PHLYNX = {'vessel_type': 'module_type', 'BC_type': 'module_subtype',
             'inp_vessels': 'inp_instances', 'out_vessels': 'out_instances'}
LIST_KEYS = ('inp_instances', 'out_instances')
KEY_ORDER = ('name', 'module_type', 'module_subtype', 'instance', 'inp_instances', 'out_instances')


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, str):
        return v.split()
    return [str(x) for x in v]


def normalise_record(rec):
    '''A record in either key style (or a CSV row) -> PhLynx keys, lists as lists, key order kept tidy.'''
    if {'vessel_type', 'BC_type'} & set(rec) and {'module_subtype', 'inp_instances', 'out_instances'} & set(rec):
        raise ValueError(f'vessel-array record mixes libcuflynx and PhLynx keys: {rec}')
    out = {}
    for k, v in rec.items():
        k = TO_PHLYNX.get(k.strip() if isinstance(k, str) else k, k)
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
    '''A vessel array (.json, or legacy .csv) as a list of PhLynx-key records.'''
    if path.endswith('.json'):
        with open(path) as f:
            data = json.load(f)
        if not isinstance(data, list):
            raise ValueError(f'{path}: a vessel array is a JSON list of instance records')
        return [normalise_record(r) for r in data]
    with open(path, newline='') as f:
        rows = [{(k or '').strip(): (v or '').strip() for k, v in r.items()}
                for r in csv.DictReader(f, skipinitialspace=True)]
    return [normalise_record(r) for r in rows if r.get('name') and not r['name'].startswith('#')]


def write_records(path, records):
    '''Writes records as <...>_vessel_array.json, one record per line.'''
    recs = [normalise_record(r) for r in records]
    with open(path, 'w') as f:
        f.write('[\n' + ',\n'.join(' ' + json.dumps(r) for r in recs) + '\n]\n')
    return path


def find(directory, prefix):
    '''<prefix>_vessel_array.json, else the legacy .csv; None if neither exists.'''
    for ext in ('.json', '.csv'):
        p = os.path.join(directory, f'{prefix}_vessel_array{ext}')
        if os.path.isfile(p):
            return p
    return None


def from_rows(rows):
    '''Spec harness rows [name, module_subtype, module_type, inp, out(, instance)] (inp/out
    space-separated) -> records.'''
    out = []
    for r in rows:
        n, bc, vt, i, o = r[:5]
        rec = {'name': str(n), 'vessel_type': str(vt), 'BC_type': str(bc)}
        if len(r) > 5 and r[5]:
            rec['instance'] = str(r[5])
        rec.update({'inp_vessels': str(i or ''), 'out_vessels': str(o or '')})
        out.append(normalise_record(rec))
    return out


def libcuflynx_reads_json():
    '''Whether the installed libcuflynx reads JSON vessel arrays (its config_schemas.read_module_array_records,
    called read_vessel_array_records before circulatory_autogen #549).'''
    try:
        from libcuflynx.utilities import config_schemas
    except ImportError:
        return False
    return hasattr(config_schemas, 'read_module_array_records') or hasattr(config_schemas, 'read_vessel_array_records')


def to_library_versions(records, instance='default'):
    '''Records naming circulatory_autogen's (or this library's pre-versions) module pairs -> this
    library's (module_type, version), renamed where the move to versions changed them
    (library.legacy_renames), each naming ``instance`` unless it names one already.'''
    from cam_testing.library import legacy_renames
    table = legacy_renames()
    out = []
    for rec in records:
        rec = normalise_record(rec)
        new = table.get((rec.get('module_type'), rec.get('module_subtype')))
        if new:
            rec['module_type'], rec['module_subtype'] = new
        rec.setdefault('instance', instance)
        out.append(normalise_record(rec))
    return out

