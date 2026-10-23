"""Replace a system model's written-out heart with one heart supermodule record.

    python tools/heart_to_supermodule.py system_models/<category>/<model> --version Argus2026_v01

The written-out heart is the model's cardiac_clock, chamber and valve records. They must be wired
as the heart version's own module array wires its submodules (the clock is matched by module_type,
whatever its name). The tool:

- writes one "heart" record of the version where the first heart record was, with the heart's host
  neighbours in per_submodule_inputs / per_submodule_outputs, and names "heart" in place of the
  submodules in the host records;
- folds in a flow_merge record that feeds only a chamber, and a flow_split record fed only by a
  chamber: the chamber's vessel ports sum their flows, so the junction's neighbours connect to the
  chamber directly;
- writes the heart's values into the parameters file with ``cuflynx-upgrade supermodule-instances``
  (the version's template), then the model's own values of the submodules over them, under
  "heart/<submodule>";
- renames the submodules in params_for_id targets, obs_data operands and the system spec's
  output_map and ignore (``ra/E_A`` -> ``heart/ra/E_A``, ``heart_clock/s`` -> ``heart/clock/s``).

The tool does not change reference/ (circulatory_autogen's flat-heart original). Once the conversion
is checked (the model's outputs are unchanged), remove reference/ and set the spec's equivalence to
{status: not_applicable, reason: ...}; closed_loop_cvs/obsolete_3compartment is the flat-heart check.
"""
import argparse
import glob
import json
import os
import re
import sys

HEART_TYPES = ('cardiac_clock', 'chamber', 'valve')
# a written-out submodule version -> the supermodule's version of the same maths
SAME_MATHS = {('chamber', 'Liang2009NoVisco_v01'): 'Liang2009NoVisco_v02'}
INSTANCE = 'heart'


def _load(path):
    with open(path) as f:
        return json.load(f)


def _write_array(records, path):
    with open(path, 'w') as f:
        f.write('[\n' + ',\n'.join(' ' + json.dumps(r) for r in records) + '\n]\n')


def _dedupe(names):
    out = []
    for n in names:
        if n not in out:
            out.append(n)
    return out


def version_array(library, version):
    path = os.path.join(library, 'heart', 'versions', version, f'heart_{version}_module_array.json')
    return _load(path)


def convert_array(records, sub_array, version):
    """(new records, {old name: submodule name}, folded junction names)."""
    heart = [r for r in records if r['module_type'] in HEART_TYPES]
    if not heart:
        raise ValueError('no cardiac_clock, chamber or valve records')
    by_type = {}
    for s in sub_array:
        by_type.setdefault(s['module_type'], []).append(s['name'])
    local = {}
    for r in heart:
        if r['module_type'] == 'cardiac_clock':
            local[r['name']] = by_type['cardiac_clock'][0]
        elif r['name'] in {s['name'] for s in sub_array}:
            local[r['name']] = r['name']
        else:
            raise ValueError(f'"{r["name"]}" is not a submodule of heart {version}')
    names = set(local)
    by_name = {r['name']: r for r in records}
    junction_into, junction_from = {}, {}   # junction -> chamber
    for r in records:
        if r['module_type'] == 'flow_merge' and len(r['out_instances']) == 1 and r['out_instances'][0] in names:
            junction_into[r['name']] = r['out_instances'][0]
        if r['module_type'] == 'flow_split' and len(r['inp_instances']) == 1 and r['inp_instances'][0] in names:
            junction_from[r['name']] = r['inp_instances'][0]
    folded = set(junction_into) | set(junction_from)
    for j in folded:
        if any(n in folded for n in by_name[j]['inp_instances'] + by_name[j]['out_instances']):
            raise ValueError(f'junction "{j}" is connected to another junction')

    subs = {s['name']: s for s in sub_array}
    for r in heart:
        sub = subs[local[r['name']]]
        if SAME_MATHS.get((r['module_type'], r['module_subtype']), r['module_subtype']) != sub['module_subtype']:
            raise ValueError(f'"{r["name"]}" is {r["module_type"]} {r["module_subtype"]}; heart {version} '
                             f'has {sub["module_subtype"]}')
        inner_out = [local[n] for m in r['out_instances']
                     for n in (by_name[m]['out_instances'] if m in junction_from else [m]) if n in names]
        if sorted(inner_out) != sorted(sub['out_instances']):
            raise ValueError(f'"{r["name"]}" outputs {inner_out} in the heart; heart {version} wires '
                             f'{sub["out_instances"]}')

    inputs, outputs = {}, {}
    for r in heart:
        ins = []
        for n in r['inp_instances']:
            if n in junction_into:
                ins += [m for m in by_name[n]['inp_instances'] if m not in names]
            elif n not in names and n not in junction_from:
                ins.append(n)
        outs = []
        for n in r['out_instances']:
            if n in junction_from:
                outs += [m for m in by_name[n]['out_instances'] if m not in names]
            elif n not in names and n not in junction_into:
                outs.append(n)
        if ins:
            inputs[local[r['name']]] = ins
        if outs:
            outputs[local[r['name']]] = outs

    replaced = names | folded
    record = {'name': INSTANCE, 'module_type': 'heart', 'module_subtype': version, 'parameterisation': 'default',
              'per_submodule_inputs': inputs, 'per_submodule_outputs': outputs}
    out = []
    for r in records:
        if r['name'] in replaced:
            if r is heart[0]:
                out.append(record)
            continue
        r = dict(r)
        for key in ('inp_instances', 'out_instances'):
            r[key] = _dedupe([INSTANCE if n in replaced else n for n in r[key]])
        out.append(r)
    return out, local, folded


def rename_path(text_name, local):
    head, sep, rest = text_name.partition('/')
    if sep and head in local:
        return f'{INSTANCE}/{local[head]}/{rest}'
    return text_name


def convert_parameters(model_dir, prefix, local, folded, library):
    from libcuflynx import upgrade
    from libcuflynx.utilities.parameters_file import read_parameters_file, write_parameters_file
    path = os.path.join(model_dir, f'{prefix}_parameters.json')
    doc = read_parameters_file(path) if os.path.exists(path) else {'schema_version': 1}
    own = {local[k]: doc.pop(k) for k in list(doc) if k in local}
    for j in folded:
        leftover = doc.pop(j, None)
        if leftover:
            raise ValueError(f'junction "{j}" has parameters {sorted(leftover)}')
    doc.pop(INSTANCE, None)
    write_parameters_file(path, doc)
    upgrade.main(['supermodule-instances', os.path.join(model_dir, f'{prefix}_module_array.json'),
                  '--module-library-dir', library, '--no-builtin-modules', '--replace'])
    doc = read_parameters_file(path)
    for sub, values in own.items():
        doc[INSTANCE].setdefault(sub, {}).update(values)
    write_parameters_file(path, doc)


def convert_params_for_id(path, local):
    text = open(path).read()
    targets = {t for p in _load(path).get('params', []) for key in ('targets', 'modifies') for t in p.get(key, [])}
    for old in sorted(targets, key=len, reverse=True):
        new = rename_path(old, local)
        if new != old:
            text = text.replace(f'"{old}"', f'"{new}"')
    with open(path, 'w') as f:
        f.write(text)


def convert_obs_data(path, local):
    text = open(path).read()
    pattern = re.compile(r'"(%s)/([^"/]+)"' % '|'.join(re.escape(n) for n in local))
    text = pattern.sub(lambda m: f'"{INSTANCE}/{local[m.group(1)]}/{m.group(2)}"', text)
    with open(path, 'w') as f:
        f.write(text)


def convert_spec(path, local, folded, version):
    text = open(path).read()
    pattern = re.compile(r'(?<![\w/])(%s)/(?=\w)' % '|'.join(re.escape(n) for n in local))
    text = pattern.sub(lambda m: f'{INSTANCE}/{local[m.group(1)]}/', text)
    for j in folded:
        text = re.sub(r'\s*%s/t: time,' % re.escape(j), '', text)
    with open(path, 'w') as f:
        f.write(text)


def convert(model_dir, version, library):
    prefix = os.path.basename(os.path.normpath(model_dir))
    array_path = os.path.join(model_dir, f'{prefix}_module_array.json')
    records, local, folded = convert_array(_load(array_path), version_array(library, version), version)
    _write_array(records, array_path)
    convert_parameters(model_dir, prefix, local, folded, library)
    for p in glob.glob(os.path.join(model_dir, f'{prefix}_params_for_id.json')):
        convert_params_for_id(p, local)
    for p in glob.glob(os.path.join(model_dir, f'{prefix}_obs_data.json')):
        convert_obs_data(p, local)
    convert_spec(os.path.join(model_dir, f'{prefix}_system.yaml'), local, folded, version)
    return local, folded


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('model_dirs', nargs='+')
    parser.add_argument('--version', required=True, help='the heart supermodule version')
    parser.add_argument('--module-library-dir', default=os.path.join(os.path.dirname(__file__), '..', 'modules'))
    args = parser.parse_args(argv)
    for d in args.model_dirs:
        local, folded = convert(d, args.version, os.path.abspath(args.module_library_dir))
        print(f'{d}: heart {args.version}; {len(local)} submodules'
              + (f'; folded {", ".join(sorted(folded))}' if folded else ''))


if __name__ == '__main__':
    main()
