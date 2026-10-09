"""
Module arrays: <prefix>_module_array.json, a JSON list of instance records in PhLynx's key names

    {"name": "venous_svc", "module_type": "venous", "module_subtype": "vp",
     "inp_instances": ["systemic_T"], "out_instances": ["heart", "volume_sum"]}

(module_type is libcuflynx's vessel_type, module_subtype its BC_type: the version). Each record is
an instance: one use of a version in the network. It may name the version's parameterisation it
uses, "parameterisation": "default" (its parameters are the defaults; the model's own parameters
file wins; several instances can use one parameterisation). The key was "instance" before
2026-10-09 and is still read (parameterisation_of). A supermodule instance has no inp/out lists of its own; it
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
KEY_ORDER = ('name', 'module_type', 'module_subtype', 'parameterisation', 'inp_instances', 'out_instances')
PARAMETERISATION_KEY, LEGACY_PARAMETERISATION_KEY = 'parameterisation', 'instance'


def parameterisation_of(record, default=None):
    '''The parameterisation a module-array record or supermodule submodule names (its
    "parameterisation", or the former key "instance"), else ``default``.'''
    return record.get(PARAMETERISATION_KEY) or record.get(LEGACY_PARAMETERISATION_KEY) or default


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, str):
        return v.split()
    return [str(x) for x in v]


def normalise_record(rec):
    '''A record in either key style (or a CSV row) -> PhLynx keys, lists as lists, key order kept tidy.'''
    if {'vessel_type', 'BC_type'} & set(rec) and {'module_subtype', 'inp_instances', 'out_instances'} & set(rec):
        raise ValueError(f'module-array record mixes libcuflynx and PhLynx keys: {rec}')
    out = {}
    for k, v in rec.items():
        k = TO_PHLYNX.get(k.strip() if isinstance(k, str) else k, k)
        out[k] = v.strip() if isinstance(v, str) and k not in LIST_KEYS else v
    for k in LIST_KEYS:
        if k in out:
            out[k] = _as_list(out[k])
    if LEGACY_PARAMETERISATION_KEY in out:   # the former key
        legacy = out.pop(LEGACY_PARAMETERISATION_KEY)
        if out.get(PARAMETERISATION_KEY, legacy) != legacy:
            raise ValueError(f'module-array record names two parameterisations ("{PARAMETERISATION_KEY}" and '
                             f'the former "{LEGACY_PARAMETERISATION_KEY}"): {rec}')
        out.setdefault(PARAMETERISATION_KEY, legacy)
    if 'per_submodule_inputs' not in out and 'per_submodule_outputs' not in out:
        for k in LIST_KEYS:
            out.setdefault(k, [])
    head = {k: out[k] for k in KEY_ORDER if k in out}
    return {**head, **{k: v for k, v in out.items() if k not in head}}


def read_records(path):
    '''A module array (.json, or legacy .csv) as a list of PhLynx-key records.'''
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


# Module arrays were called vessel arrays (circulatory_autogen #549); the old file names are read too,
# after the new ones: system_models/**/reference/ keeps circulatory_autogen's originals under them.
NAMES = ('module_array', 'vessel_array')


def find(directory, prefix):
    '''<prefix>_module_array.json, else .csv, else the same under the old name vessel_array; None if
    none exists.'''
    for name in NAMES:
        for ext in ('.json', '.csv'):
            p = os.path.join(directory, f'{prefix}_{name}{ext}')
            if os.path.isfile(p):
                return p
    return None


def harness_rows(harness):
    '''A verification_config harness's network rows: its "module_array", or the old key "vessel_array".'''
    harness = harness or {}
    return harness.get('module_array') or harness.get('vessel_array')


def from_rows(rows):
    '''Spec harness rows [name, module_subtype, module_type, inp, out(, parameterisation)] (inp/out
    space-separated) -> records.'''
    out = []
    for r in rows:
        n, bc, vt, i, o = r[:5]
        rec = {'name': str(n), 'vessel_type': str(vt), 'BC_type': str(bc)}
        if len(r) > 5 and r[5]:
            rec[PARAMETERISATION_KEY] = str(r[5])
        rec.update({'inp_vessels': str(i or ''), 'out_vessels': str(o or '')})
        out.append(normalise_record(rec))
    return out


def libcuflynx_reads_json():
    '''Whether the installed libcuflynx reads JSON module arrays (its config_schemas.read_module_array_records,
    called read_vessel_array_records before circulatory_autogen #549).'''
    try:
        from libcuflynx.utilities import config_schemas
    except ImportError:
        return False
    return hasattr(config_schemas, 'read_module_array_records') or hasattr(config_schemas, 'read_vessel_array_records')


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

