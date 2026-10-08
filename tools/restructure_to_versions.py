"""
Moves the module library to the versions/instances layout, driven by tools/restructure_map.yaml:

    modules/<category path>/<module_type>/versions/<version>/instances/<instance>/

and modules/system/ to system_models/. See modules/README.md for the layout and its rules.

    python tools/restructure_to_versions.py                     # from the old layout in modules/
    python tools/restructure_to_versions.py --source <old modules/ copy>   # rebuild from a backup

Each version directory gets its one config entry (module_subtype == version, "default_instance":
"default"), its CellML component (shared components copied and renamed
<component>__<module_type>_<version>), the units it uses, its spec (<module_type>_<version>_tests.yaml:
the old component spec with the module defaults merged in, the module's review block and reviewed
flag, and a validation block per instance), the bibliography entries it cites, its risk files, and
its instances: default (the old parameter rows, vessel suffix dropped) plus one per validation data
set (named as in the map; obs_data files get "obs_data_name"). Supermodules become versions of the
module_type they replace, with their submodule records rewritten to the new (module_type,
module_subtype) pairs. System-model vessel records get the new pairs plus "instance": "default".
Earlier test results and plots are carried over, so the reports keep what the last run found.

Idempotent: the new tree is built from the source in a staging directory and moved into place, so a
run from the same source gives the same files. With the old layout gone and no --source, there is
nothing to do. Another session's uncommitted work is left alone: modules/poiseuille_transport/ and
modules/system/poiseuille/ stay where they are (see LEFT_IN_PLACE).

Historical: this made the first versions layout (categories cardiac and cell/neurons). The later move
to nested module_types (cell/neuron and heart as module_types holding their parts, the alternative
hearts as versions of heart) was done by hand and is not in the map; modules/README.md describes the
current layout, and cam_testing.library.HEART_VERSION_RENAMES the heart renames.
"""
import argparse
import copy
import csv
import glob
import json
import os
import re
import shutil
import sys
import tempfile

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tools'))
from cam_testing import bib as bibmod  # noqa: E402
from cam_testing.library import DEFAULT_LICENCE, write_spec  # noqa: E402
from import_from_libcuflynx import read_units_blocks, units_closure, write_units_file  # noqa: E402

MAP = os.path.join(REPO, 'tools', 'restructure_map.yaml')
MODULES = os.path.join(REPO, 'modules')
SYSTEM_MODELS = os.path.join(REPO, 'system_models')
# another session's uncommitted work: not moved, not edited
LEFT_IN_PLACE = ('poiseuille_transport', os.path.join('system', 'poiseuille'))
TEMPLATE_DEST = os.path.join(REPO, 'manifests', 'template_modules.cellml')
REVIEWS = os.path.join(REPO, 'reviews')

CELLML_HEAD = ("<?xml version='1.0' encoding='UTF-8'?>\n"
               '<model name="{name}" xmlns="http://www.cellml.org/cellml/1.1#" '
               'xmlns:cellml="http://www.cellml.org/cellml/1.1#">\n')
INSTANCE_COLUMNS = ['variable_name', 'units', 'value', 'data_reference', 'sourced']
# validation files that are not part of a data instance: which versions' default instance gets them
EXTRA_VALIDATION_FILES = {
    ('gas_exchange', 'severinghaus1979_odc.csv'): [('gas_transport_simple', 'nn')],
    ('gas_exchange', 'SOURCES.md'): [('gas_transport_simple', 'nn')],
    ('lung', 'SOURCES.md'): [('simple_lung_bg', 'nn')],
}
# system-model specs' text that named the old places
SYSTEM_TEXT = [
    ('the heart supermodule (modules/supermodules/heart)', 'the heart supermodule version (cardiac/heart Argus2026_v01)'),
    ('(modules/supermodules/heart/heart_parameters.csv)',
     '(modules/cardiac/heart/versions/Argus2026_v01/instances/default/default_parameters.csv)'),
    ('the supermodule defaults', 'the default instance of the supermodule version'),
    ('tests/test_supermodules.py', 'tests/test_modules.py supermodule_equivalence_test'),
    ('one sympathetic_neuron supermodule instance "SN" (modules/supermodules/sympathetic_neuron), which\n  nests the SN_soma and SN_varicosity supermodules.',
     'one record "SN" of the supermodule version cell/neurons/neuron sympathetic, which\n  nests the soma and varicosity sympathetic supermodule versions.'),
    ("the supermodules'' defaults", "the supermodule versions'' default instances"),
]
# supermodule descriptions' text that named the old files and places
DESCRIPTION_TEXT = [
    ('modules/ion_channel, *_SN', 'cell/ion_channels, versions Argus2026_v01'),
    ('(modules/cardiac)', '(cardiac: cardiac_clock, chamber, valve)'),
    ('heart_modules_config.json', 'heart_Argus2026_v01_modules_config.json'),
    ('in heart_parameters.csv', 'in its default instance'),
    ('(SN_soma_parameters.csv)', '(its default instance)'),
    ('(SN_varicosity_parameters.csv)', '(its default instance)'),
    ('sympathetic_neuron_parameters.csv holds', 'its default instance holds'),
    ("the SN_soma / SN_varicosity supermodules'", "the soma / varicosity sympathetic versions'"),
    ('Tested through the sympathetic_neuron supermodule.', 'Tested through the neuron sympathetic version.'),
]
YAML_KW = dict(sort_keys=False, width=120, default_flow_style=None, allow_unicode=True)

notes = []


def note(msg):
    notes.append(msg)
    print('NOTE:', msg)


def deep_merge(base, override):
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def looks_sourced(reference):
    from cam_testing.library import looks_sourced as ls
    return ls(reference)


# ----------------------------------------------------------------------------------------------
# the old layout
# ----------------------------------------------------------------------------------------------

class OldModule:
    def __init__(self, d):
        self.dir = d
        self.name = os.path.basename(d)
        n = self.name
        self.config = json.load(open(os.path.join(d, f'{n}_modules_config.json')))
        sp = os.path.join(d, f'{n}_tests.yaml')
        self.spec = yaml.safe_load(open(sp)) if os.path.isfile(sp) else {}
        self.cellml = open(os.path.join(d, f'{n}_modules.cellml')).read()
        up = os.path.join(d, f'{n}_units.cellml')
        self.units = read_units_blocks([up]) if os.path.isfile(up) else {}
        pp = os.path.join(d, f'{n}_parameters.csv')
        self.params = list(csv.DictReader(open(pp, newline=''))) if os.path.isfile(pp) else []
        self.bib = self._bib_blocks(os.path.join(d, f'{n}_references.bib'))
        self.bib_proposed = self._bib_blocks(os.path.join(d, f'{n}_references_proposed.bib'))
        self.comp_specs = {(c['vessel_type'], c['BC_type']): c for c in self.spec.get('components', [])}

    @staticmethod
    def _bib_blocks(path):
        if not os.path.isfile(path):
            return {}
        text = open(path).read()
        out = {}
        for e in re.split(r'(?m)^(?=@)', text):
            e = e.strip()
            m = re.match(r'@\w+\s*\{\s*([^,\s]+)\s*,', e)
            if m:
                out[m.group(1)] = e
        return out


def find_old_modules(src):
    found = {}
    for root, dirs, files in os.walk(src):
        rel = os.path.relpath(root, src)
        if rel.split(os.sep)[0] in ('system', 'supermodules', 'template', 'poiseuille_transport'):
            dirs[:] = []
            continue
        dirs[:] = sorted(d for d in dirs if d not in ('plots', 'results', 'risk', 'validation', 'versions'))
        n = os.path.basename(root)
        if f'{n}_modules_config.json' in files and f'{n}_modules.cellml' in files:
            found[n] = OldModule(root)
    return found


def component_block(text, name):
    m = re.search(r'<component name="%s"[\s>].*?</component>' % re.escape(name), text, re.S)
    return m.group(0) if m else None


# ----------------------------------------------------------------------------------------------
# the map
# ----------------------------------------------------------------------------------------------

FROM_RE = re.compile(r'^(?P<module>[\w/]+): \((?P<mt>[^,]+), (?P<st>[^)]+)\)$')


def load_map():
    m = yaml.safe_load(open(MAP))
    versions = []          # (category, module_type, version, info)
    for cat, mts in m['categories'].items():
        for mt, node in mts.items():
            for v, info in node['versions'].items():
                versions.append((cat, mt, v, info))
    return versions


def rename_table(versions, supermodule_dirs):
    '''(old module_type, old subtype) -> (new module_type, version), for every entry and supermodule.'''
    table = {}
    for cat, mt, v, info in versions:
        m = FROM_RE.match(info['from'])
        if m:
            table[(m['mt'], m['st'])] = (mt, v)
        else:
            sdir = re.match(r'supermodules/(\w+)', info['from']).group(1)
            old_type = supermodule_dirs[sdir]['module_type']
            old_sub = supermodule_dirs[sdir]['module_subtype']
            table[(old_type, old_sub)] = (mt, v)
    return table


# ----------------------------------------------------------------------------------------------
# writing
# ----------------------------------------------------------------------------------------------

def write_json(path, data):
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')


def write_yaml(path, data, header=None):
    with open(path, 'w') as f:
        if header:
            f.write(header)
        yaml.safe_dump(data, f, **YAML_KW)


def write_instance_parameters(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=INSTANCE_COLUMNS, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def config_head(entry, mt, v, component_file, component_type, extra):
    '''The version's config entry in PhLynx key order.'''
    head = {'module_type': mt, 'module_subtype': v}
    if 'module_format' in entry:
        head['module_format'] = entry['module_format']
    if component_file is not None:
        head['component_file'] = component_file
        head['component_type'] = component_type
    head['default_instance'] = 'default'
    head['licence'], head['creator'] = entry.get('licence', DEFAULT_LICENCE), list(entry.get('creator', []))
    head.update(extra)
    skip = {'module_type', 'module_subtype', 'module_format', 'component_file', 'component_type',
            'default_parameters', 'default_instance', 'vessel_type', 'BC_type', 'module_file', 'licence', 'creator'}
    return {**head, **{k: val for k, val in entry.items() if k not in skip}}


def rename_rows(rows, table, where):
    '''Harness vessel_array rows [name, BC_type, vessel_type, inp, out(, instance)] -> new names + instance.'''
    out = []
    for r in rows:
        r = list(r)
        new = table.get((str(r[2]), str(r[1])))
        if new is None:
            note(f'{where}: harness row {r[0]} ({r[2]}, {r[1]}) is not in the library map; left as is')
        else:
            r[2], r[1] = new
        r = r[:5] + ['default']
        out.append(r)
    return out


def cited_keys(rows, spec):
    keys = {bibmod.reference_key(r.get('data_reference', '')) for r in rows}
    for prop in (spec.get('reference_proposals') or {}).values():
        if isinstance(prop, dict):
            keys.add(bibmod.reference_key(str(prop.get('reference', ''))))
            keys.update(prop.get('also_cites') or [])
    return {k for k in keys if k}


def write_bibs(dest, stem, rows, spec, own, library_bibs):
    '''<stem>_references.bib / _references_proposed.bib: only the entries this version cites.'''
    text = json.dumps(spec, ensure_ascii=False) + json.dumps(rows)
    confirmed, proposed = {}, {}
    all_keys = set(own[0]) | set(own[1]) | set(library_bibs[0]) | set(library_bibs[1])
    cites = cited_keys(rows, spec) | {k for k in all_keys if re.search(r'(?<![\w-])%s(?![\w-])' % re.escape(k), text)}
    for k in sorted(cites):
        if k in own[0]:
            confirmed[k] = own[0][k]
        elif k in own[1]:
            proposed[k] = own[1][k]
        elif k in library_bibs[0]:
            confirmed[k] = library_bibs[0][k]
        elif k in library_bibs[1]:
            proposed[k] = library_bibs[1][k]
    for suffix, entries in (('_references.bib', confirmed), ('_references_proposed.bib', proposed)):
        path = os.path.join(dest, stem + suffix)
        if entries or suffix == '_references.bib':
            with open(path, 'w') as f:
                f.write('\n\n'.join(entries.values()) + ('\n' if entries else ''))


def _relocate_plots(result, old_prefix, new_prefix, src_dir, dest_dir):
    out = []
    for p in result.get('plots') or []:
        base = os.path.basename(p)
        new_base = new_prefix + (base[len(old_prefix):] if base.startswith(old_prefix) else base)
        src = os.path.join(src_dir, p)
        if os.path.isfile(src):
            os.makedirs(os.path.join(dest_dir, 'plots'), exist_ok=True)
            shutil.copy2(src, os.path.join(dest_dir, 'plots', new_base))
        out.append(f'plots/{new_base}')
    result['plots'] = out
    return result


# ----------------------------------------------------------------------------------------------
# one version
# ----------------------------------------------------------------------------------------------

def validation_split(old_val, instance_files):
    '''The old validation block -> {instance: {kind: block}}: blocks whose data files moved to the
    data instance go there (file keys dropped / made instance-relative), the rest stay on default.'''
    out = {'default': {}}
    data_inst = instance_files[0] if instance_files else None
    files = set(instance_files[1]) if instance_files else set()
    for kind in ('baseline', 'calibrate'):
        blk = copy.deepcopy((old_val or {}).get(kind))
        if blk is None:
            continue
        moved = any(isinstance(blk.get(k), str) and blk[k] in files for k in ('data', 'obs_data', 'params_for_id',
                                                                              'validation_obs_data'))
        if moved:
            # file references point at the instance's files, relative to the version directory
            idir = f'instances/{data_inst}/'
            if 'data' in blk:
                blk['data'] = idir + os.path.basename(blk['data'])
            if kind == 'calibrate':
                same = blk.get('obs_data') == blk.get('validation_obs_data')
                blk['obs_data'] = f'{idir}{data_inst}_obs_data.json'
                blk['params_for_id'] = f'{idir}{data_inst}_params_for_id.csv'
                if same:
                    blk.pop('validation_obs_data', None)
                else:
                    blk['validation_obs_data'] = f'{idir}{data_inst}_validation_obs_data.json'
            out.setdefault(data_inst, {})[kind] = blk
        elif kind == 'baseline':
            out['default'][kind] = blk
        else:
            # calibration needs obs_data in the instance; keep the old reason as a note
            if blk.get('reason'):
                out['default'][kind] = {'note': f"{blk.get('status', 'pending')}: {blk['reason']}"}
    return out


def known_issue_applies(issue, entry, comp, block, undefined_units, mt, comp_spec):
    """Whether one of the old module's known_issues (substrings of structure problems, or notes)
    is about this version."""
    m = re.match(r'port variable (\S+) not in variables_and_units', issue)
    if m:
        names = {x[0] for x in entry.get('variables_and_units', [])}
        ports = [x for side in ('entrance_ports', 'exit_ports', 'general_ports') for p in entry.get(side, [])
                 for x in p.get('variables', [])]
        return m[1] in ports and m[1] not in names
    m = re.match(r'undefined units: (\S+)', issue)
    if m:
        return m[1] in undefined_units
    m = re.match(r"variable (\S+) has kind", issue)
    if m:
        return any(x[0] == m[1] and x[3] != x[3].strip() for x in entry.get('variables_and_units', []))
    m = re.match(r'module_type (\S+) is not a component', issue)
    if m:
        return comp == m[1] and block is None
    if re.search(r'(?<!\w)%s(?!\w)' % re.escape(mt), issue):
        return True
    words = [w for w in re.findall(r'[A-Za-z_]\w{3,}', issue) if '_' in w or any(c.isdigit() for c in w)]
    text = json.dumps(comp_spec)
    return any(w in text for w in words if w not in ('variables_and_units',))


def build_version(stage, cat, mt, v, info, old, table, library_bibs, all_blocks, all_modules):
    m = FROM_RE.match(info['from'])
    omt, ost = m['mt'], m['st']
    entry = next(e for e in old.config if e['module_type'] == omt and e['module_subtype'] == ost)
    comp_spec = old.comp_specs.get((omt, ost))
    if comp_spec is None:
        note(f'{mt}/{v}: no component spec in {old.name}_tests.yaml; an empty one is used')
        comp_spec = {}
    dest = os.path.join(stage, cat, mt, 'versions', v)
    os.makedirs(dest, exist_ok=True)
    stem = f'{mt}_{v}'
    fmt = entry.get('module_format', 'cellml')

    # CellML
    comp = entry['component_type']
    new_comp = info.get('component_renamed', comp)
    block = component_block(old.cellml, comp)
    if block is None and fmt == 'cellml':
        # a component defined in another module's file (libcuflynx merges the whole library)
        owner = next((o for o in all_modules.values() if component_block(o.cellml, comp)), None)
        if owner is None:
            note(f'{mt}/{v}: component {comp} is not defined in any module (a known defect); '
                 'its CellML file is left without it')
        else:
            block = component_block(owner.cellml, comp)
            note(f'{mt}/{v}: component {comp} copied from {owner.name}_modules.cellml (it was defined there)')
    body = ''
    if block is not None:
        body = '    ' + block.replace(f'<component name="{comp}"', f'<component name="{new_comp}"', 1) + '\n'
        tail = '</model>\n'
    else:
        # no CellML component; the comment shares the closing line so libcuflynx, which copies
        # every module file's lines from the first mentioning "component" up to the last line,
        # copies nothing from it
        why = ('C++ 1D-solver module (module_format cpp): its component is the 1D solver, not CellML'
               if fmt == 'cpp' else f'{comp} is not defined anywhere (a known defect)')
        tail = f'    <!-- no CellML component: {why} --></model>\n'
    with open(os.path.join(dest, f'{stem}_modules.cellml'), 'w') as f:
        f.write(CELLML_HEAD.format(name=stem) + body + tail)

    # config
    extra = {'notes': info['notes']} if info.get('notes') else {}
    new_entry = config_head(entry, mt, v, f'{stem}_modules.cellml', new_comp, extra)
    write_json(os.path.join(dest, f'{stem}_modules_config.json'), [new_entry])

    # units
    used = set(re.findall(r'\bunits="([^"]+)"', body)) | {x[1] for x in entry.get('variables_and_units', [])}
    needed, missing = units_closure(used, old.units)
    undefined_units = set()
    if missing:
        needed2, missing2 = units_closure(missing, all_blocks)
        needed |= needed2
        undefined_units = set(missing2)
        if missing2:
            note(f'{mt}/{v}: units not defined anywhere: {sorted(missing2)}')
        else:
            note(f'{mt}/{v}: units {sorted(missing)} taken from other modules\' units files')
    blocks = dict(old.units)
    for k, b in all_blocks.items():
        blocks.setdefault(k, b)
    write_units_file(os.path.join(dest, f'{stem}_units.cellml'), needed, blocks)

    # parameters -> default instance
    kinds = {x[0]: x[3].strip() for x in entry.get('variables_and_units', [])}
    globals_ = {n for n, k in kinds.items() if k == 'global_constant'}
    own = [r for r in old.params if r['vessel_type'] == omt and r['BC_type'] == ost]
    glob_rows = [r for r in old.params if r['vessel_type'] == 'global' and r['variable_name'] in globals_]
    for r in own:
        k = (r.get('kind') or '').strip()
        if r['variable_name'] not in kinds:
            note(f'{mt}/{v}: parameter row {r["variable_name"]} is not a variable of the config entry')
        elif k and k != kinds[r['variable_name']]:
            note(f'{mt}/{v}: {r["variable_name"]} has kind {k} in the parameters file but '
                 f'{kinds[r["variable_name"]]} in the config (the config\'s kind is used)')
    rows = [{c: (r.get(c) or '').strip() for c in INSTANCE_COLUMNS} for r in own + glob_rows]
    # the old loader filled a TODO global from a reference proposal of any component of the module;
    # each version now has its own copy of the global, so it carries that proposal itself
    inherited = {}
    for r in glob_rows:
        if (r.get('value') or '').strip().upper() != 'TODO' or r['variable_name'] in (comp_spec.get('reference_proposals') or {}):
            continue
        for other in old.spec.get('components', []):
            prop = (other.get('reference_proposals') or {}).get(r['variable_name'])
            if isinstance(prop, dict) and 'value' in prop:
                prop = dict(prop)
                prop['note'] = ((prop.get('note') + ' ') if prop.get('note') else '') + (
                    f"(proposed for the {old.name} global {r['variable_name']} under {other['vessel_type']}/{other['BC_type']})")
                inherited[r['variable_name']] = prop
                break

    # spec
    spec = {'module_type': mt, 'version': v, 'reviewed': bool(old.spec.get('reviewed', False))}
    known = [k for k in (old.spec.get('known_issues') or [])
             if known_issue_applies(k, entry, comp, block, undefined_units, mt, comp_spec)]
    if known:
        spec['known_issues'] = known
    merged = deep_merge(old.spec.get('defaults') or {}, comp_spec)
    for k in ('vessel_type', 'BC_type', 'module_type', 'validation'):
        merged.pop(k, None)
    if isinstance((merged.get('harness') or {}).get('vessel_array'), list):
        merged['harness']['vessel_array'] = rename_rows(merged['harness']['vessel_array'], table, f'{mt}/{v}')
    if inherited:
        merged['reference_proposals'] = {**(merged.get('reference_proposals') or {}), **inherited}
    spec.update(merged)

    # instances
    inst_infos = [i for i in info['instances'] if i['name'] != 'default']
    data = None
    if inst_infos:
        if len(inst_infos) > 1:
            raise SystemExit(f'{mt}/{v}: more than one data instance is not handled')
        data = (inst_infos[0]['name'], [f for f in inst_infos[0]['from'] if f.endswith(('.json', '.csv'))])
    val = validation_split(comp_spec.get('validation'), data)
    inst_root = os.path.join(dest, 'instances')
    os.makedirs(os.path.join(inst_root, 'default'), exist_ok=True)
    write_instance_parameters(os.path.join(inst_root, 'default', 'default_parameters.csv'), rows)
    if data:
        name, files = data
        idir = os.path.join(inst_root, name)
        os.makedirs(idir, exist_ok=True)
        write_instance_parameters(os.path.join(idir, f'{name}_parameters.csv'), rows)
        cal = (comp_spec.get('validation') or {}).get('calibrate') or {}
        for f in inst_infos[0]['from']:
            src = os.path.join(old.dir, f)
            base = os.path.basename(f)
            if f == cal.get('obs_data'):
                target = f'{name}_obs_data.json'
            elif f == cal.get('validation_obs_data'):
                target = f'{name}_validation_obs_data.json'
            elif f == cal.get('params_for_id'):
                target = f'{name}_params_for_id.csv'
            else:
                target = base
            if target.endswith('_obs_data.json'):
                obs = json.load(open(src))
                obs = {'obs_data_name': name, **{k: x for k, x in obs.items() if k != 'obs_data_name'}}
                write_json(os.path.join(idir, target), obs)
            else:
                shutil.copy2(src, os.path.join(idir, target))
        if cal.get('obs_data') and cal.get('obs_data') == cal.get('validation_obs_data'):
            note(f'{mt}/{v}: calibration and validation obs_data are one file; the instance has only {name}_obs_data.json '
                 '(calibration validates on it)')
    for (omod, fname), targets in EXTRA_VALIDATION_FILES.items():
        if omod == old.name and (mt, v) in targets:
            shutil.copy2(os.path.join(old.dir, 'validation', fname), os.path.join(inst_root, 'default', fname))
    # a SOURCES.md for targets that live in the spec (no data files) goes with the default instance
    srcmd = os.path.join(old.dir, 'validation', 'SOURCES.md')
    if os.path.isfile(srcmd) and not data and any((b or {}).get('status') in ('active', 'proposed')
                                                  for b in val['default'].values()):
        shutil.copy2(srcmd, os.path.join(inst_root, 'default', 'SOURCES.md'))
    # text mentions of validation/<file> -> the file's new home
    spec['validation'] = {k: val[k] for k in val if val[k] or k == 'default'}
    txt = yaml.safe_dump(spec, **YAML_KW)
    for fname in re.findall(r'validation/([\w.]+)', txt):
        home = (f'instances/{data[0]}/{fname}' if data and any(os.path.basename(f) == fname for f in data[1])
                else f'instances/default/{fname}' if os.path.isfile(os.path.join(inst_root, 'default', fname))
                else f'{old.name} validation/{fname}')
        txt = txt.replace(f'validation/{fname}', home)
    spec = yaml.safe_load(txt)
    if old.spec.get('review'):
        spec['review'] = copy.deepcopy(old.spec['review'])
        if len(old.config) > 1:
            # one copy of a shared review, in reviews/; the versions point to it
            os.makedirs(REVIEWS, exist_ok=True)
            write_yaml(os.path.join(REVIEWS, f'{old.name}_review.yaml'), spec['review'])
            spec['review'] = f'reviews/{old.name}_review.yaml'
            spec['review_scope'] = (f'The review of the whole {old.name} module ({len(old.config)} entries) before the move to '
                                    'versions; only the parts about this version apply.')
    write_spec(dest, stem, spec)

    # bibliographies
    write_bibs(dest, stem, rows, spec, (old.bib, old.bib_proposed), library_bibs)

    # risk files and earlier results / plots
    cid = f'{omt}__{ost}'
    for p in glob.glob(os.path.join(old.dir, 'risk', f'{cid}_risk*')):
        os.makedirs(os.path.join(dest, 'risk'), exist_ok=True)
        suffix = os.path.basename(p)[len(cid):]
        target = os.path.join(dest, 'risk', f'{stem}{suffix}')
        if p.endswith('.json'):
            r = json.load(open(p))
            _relocate_plots(r, f'{cid}__', '', old.dir, dest)
            write_json(target, r)
        else:
            shutil.copy2(p, target)
    rdir = os.path.join(old.dir, 'results', cid)
    for p in sorted(glob.glob(os.path.join(rdir, '*.json'))):
        r = json.load(open(p))
        test = r.get('test') or os.path.splitext(os.path.basename(p))[0]
        if test.startswith('validation_test_'):
            kind = test.rsplit('_', 1)[1]
            inst = next((i for i, blocks in val.items() if kind in blocks and i != 'default'), 'default')
            _relocate_plots(r, f'{cid}__', f'{inst}__', old.dir, dest)
            target = os.path.join(dest, 'results', 'instances', inst, f'{test}.json')
        else:
            _relocate_plots(r, f'{cid}__', '', old.dir, dest)
            target = os.path.join(dest, 'results', f'{test}.json')
        os.makedirs(os.path.dirname(target), exist_ok=True)
        write_json(target, r)
    return dest


def build_zero_flow(stage, olds):
    '''libcuflynx's zero_flow helper (it feeds v_zero = 0 to a monolithic heart named "heart" that has one
    vena cava input; the generator writes the component by name) had no config entry in the old layout,
    so it isn't in the map. It becomes boundary_conditions/zero_flow version nn, keeping its name.'''
    owner = next((o for o in olds.values() if component_block(o.cellml, 'zero_flow')), None)
    if owner is None:
        return None
    mt, v = 'zero_flow', 'nn'
    dest = os.path.join(stage, 'boundary_conditions', mt, 'versions', v)
    os.makedirs(os.path.join(dest, 'instances', 'default'), exist_ok=True)
    stem = f'{mt}_{v}'
    block = component_block(owner.cellml, 'zero_flow')
    with open(os.path.join(dest, f'{stem}_modules.cellml'), 'w') as f:
        f.write(CELLML_HEAD.format(name=stem) + '    ' + block + '\n</model>\n')
    entry = {'module_type': mt, 'module_subtype': v, 'module_format': 'cellml', 'component_file': f'{stem}_modules.cellml',
             'component_type': 'zero_flow', 'default_instance': 'default', 'licence': DEFAULT_LICENCE, 'creator': [],
             'notes': ('libcuflynx helper: the generator connects v_zero to the v_ivc input of a monolithic heart '
                       'module named "heart" that has only one vena cava input, and finds the component by its name '
                       '(zero_flow), so the name must not change.'),
             'entrance_ports': [], 'exit_ports': [], 'general_ports': [],
             'variables_and_units': [['v_zero', 'm3_per_s', 'access', 'variable']]}
    write_json(os.path.join(dest, f'{stem}_modules_config.json'), [entry])
    needed, _ = units_closure({'m3_per_s'} | set(re.findall(r'\bunits="([^"]+)"', block)), owner.units)
    write_units_file(os.path.join(dest, f'{stem}_units.cellml'), needed, owner.units)
    write_instance_parameters(os.path.join(dest, 'instances', 'default', 'default_parameters.csv'), [])
    spec = {'module_type': mt, 'version': v, 'reviewed': False, 'sim_time': 1.0, 'pre_time': 0.0, 'dt': 0.01,
            'solver': 'CVODE_myokit', 'outputs': ['v_zero'],
            'invariants': [{'description': 'the flow is zero', 'expr': 'v_zero == 0'}],
            'notes': (f'Moved from {owner.name}_modules.cellml, where it had no config entry (libcuflynx used it by name). '
                      'Not in tools/restructure_map.yaml, which maps config entries.'),
            'validation': {'default': {}}}
    write_spec(dest, stem, spec)
    write_bibs(dest, stem, [], spec, ({}, {}), ({}, {}))
    note(f'zero_flow (a libcuflynx helper component with no config entry in {owner.name}) -> boundary_conditions/zero_flow/nn')
    return dest


def build_supermodule(stage, cat, mt, v, sdir, table, library_bibs):
    name = os.path.basename(sdir)
    raw = json.load(open(os.path.join(sdir, f'{name}_modules_config.json')))
    entry = next(e for e in raw if e.get('module_format') == 'supermodule')
    dest = os.path.join(stage, cat, mt, 'versions', v)
    os.makedirs(dest, exist_ok=True)
    stem = f'{mt}_{v}'
    subs = []
    for s in entry['submodules']:
        s = dict(s)
        new = table.get((s['module_type'], s['module_subtype']))
        if new is None:
            raise SystemExit(f'{mt}/{v}: submodule {s["name"]} ({s["module_type"]}, {s["module_subtype"]}) not in the map')
        s['module_type'], s['module_subtype'] = new
        head = {k: s[k] for k in ('name', 'module_type', 'module_subtype')}
        head['instance'] = 'default'
        subs.append({**head, **{k: x for k, x in s.items() if k not in head}})
    new_entry = config_head(entry, mt, v, None, None, {})
    new_entry['submodules'] = subs
    write_json(os.path.join(dest, f'{stem}_modules_config.json'), [new_entry])

    rows = []
    defaults = entry.get('default_parameters')
    if defaults:
        for r in csv.DictReader(open(os.path.join(sdir, defaults), newline='')):
            r = {k.strip(): (x or '').strip() for k, x in r.items()}
            if not r.get('sourced'):
                r['sourced'] = 'yes' if looks_sourced(r.get('data_reference', '')) else 'no'
            rows.append(r)
    os.makedirs(os.path.join(dest, 'instances', 'default'), exist_ok=True)
    write_instance_parameters(os.path.join(dest, 'instances', 'default', 'default_parameters.csv'), rows)

    old_spec = yaml.safe_load(open(os.path.join(sdir, f'{name}_supermodule.yaml'))) or {}
    description = old_spec.get('description', entry.get('description', ''))
    for old_text, new_text in DESCRIPTION_TEXT:
        description = description.replace(old_text, new_text)
    spec = {'module_type': mt, 'version': v, 'reviewed': bool(old_spec.get('reviewed', False)),
            'description': description}
    sm = {}
    if old_spec.get('globals'):
        sm['globals'] = old_spec['globals']
    if (old_spec.get('tests') or {}).get('equivalent'):
        sm['equivalent'] = old_spec['tests']['equivalent']
    spec['supermodule'] = sm
    spec['validation'] = {'default': {}}
    write_spec(dest, stem, spec)
    write_bibs(dest, stem, rows, spec, ({}, {}), library_bibs)
    return dest


# ----------------------------------------------------------------------------------------------
# system models
# ----------------------------------------------------------------------------------------------

def move_system_models(src_system, dest_root, table):
    count = 0
    for cat in sorted(os.listdir(src_system)):
        cdir = os.path.join(src_system, cat)
        if not os.path.isdir(cdir) or os.path.join('system', cat) in LEFT_IN_PLACE:
            continue
        for model in sorted(os.listdir(cdir)):
            mdir = os.path.join(cdir, model)
            if not os.path.isdir(mdir):
                continue
            target = os.path.join(dest_root, cat, model)
            shutil.copytree(mdir, target)
            sy = os.path.join(target, f'{model}_system.yaml')
            if os.path.isfile(sy):
                text = open(sy).read()
                for old_text, new_text in SYSTEM_TEXT:
                    text = text.replace(old_text, new_text)
                with open(sy, 'w') as f:
                    f.write(text)
            va = os.path.join(target, f'{model}_vessel_array.json')
            if os.path.isfile(va):
                recs = json.load(open(va))
                out = []
                for r in recs:
                    new = table.get((r.get('module_type'), r.get('module_subtype')))
                    if new is None:
                        note(f'system_models/{cat}/{model}: {r["name"]} ({r.get("module_type")}, '
                             f'{r.get("module_subtype")}) not in the library map; left as is')
                        out.append(r)
                        continue
                    r = dict(r)
                    r['module_type'], r['module_subtype'] = new
                    head = {k: r[k] for k in ('name', 'module_type', 'module_subtype') if k in r}
                    head['instance'] = 'default'
                    out.append({**head, **{k: x for k, x in r.items() if k not in head}})
                with open(va, 'w') as f:
                    f.write('[\n' + ',\n'.join(' ' + json.dumps(r) for r in out) + '\n]\n')
            count += 1
    return count


# ----------------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--source', help='the old modules/ directory (default: modules/ when it still has the old layout)')
    ap.add_argument('--out', help='write modules/ and system_models/ under this directory instead of the repo (a dry run)')
    args = ap.parse_args(argv)
    global MODULES, SYSTEM_MODELS, TEMPLATE_DEST, REVIEWS
    if args.out:
        MODULES = os.path.join(os.path.abspath(args.out), 'modules')
        SYSTEM_MODELS = os.path.join(os.path.abspath(args.out), 'system_models')
        TEMPLATE_DEST = os.path.join(os.path.abspath(args.out), 'template_modules.cellml')
        REVIEWS = os.path.join(os.path.abspath(args.out), 'reviews')
        os.makedirs(MODULES, exist_ok=True)
    src = os.path.abspath(args.source or os.path.join(REPO, 'modules'))
    olds = find_old_modules(src)
    if not olds:
        print(f'{src} has no modules in the old layout: nothing to do')
        return 0
    super_dirs = {}
    for cfg in glob.glob(os.path.join(src, 'supermodules', '*', '*_modules_config.json')):
        e = next(x for x in json.load(open(cfg)) if x.get('module_format') == 'supermodule')
        super_dirs[os.path.basename(os.path.dirname(cfg))] = {'dir': os.path.dirname(cfg), 'module_type': e['module_type'],
                                                             'module_subtype': e['module_subtype']}
    versions = load_map()
    table = rename_table(versions, super_dirs)
    all_blocks, lib_bib, lib_prop = {}, {}, {}
    component_users = {}
    for o in olds.values():
        for k, b in o.units.items():
            all_blocks.setdefault(k, b)
        for k, b in o.bib.items():
            lib_bib.setdefault(k, b)
        for k, b in o.bib_proposed.items():
            lib_prop.setdefault(k, b)

    stage = tempfile.mkdtemp(prefix='restructure_', dir=REPO)
    try:
        stage_modules = os.path.join(stage, 'modules')
        built = 0
        for cat, mt, v, info in versions:
            m = FROM_RE.match(info['from'])
            if m:
                old = olds[m['module'].split('/')[-1]]
                build_version(stage_modules, cat, mt, v, info, old, table, (lib_bib, lib_prop), all_blocks, olds)
            else:
                sdir = re.match(r'supermodules/(\w+)', info['from']).group(1)
                build_supermodule(stage_modules, cat, mt, v, super_dirs[sdir]['dir'], table, (lib_bib, lib_prop))
            built += 1
        if build_zero_flow(stage_modules, olds):
            built += 1
        # components with no config entry are not moved (libcuflynx can't use them without one)
        used = {e.get('component_type') for o in olds.values() for e in o.config} | {'zero_flow'}
        for o in olds.values():
            for c in re.findall(r'<component name="([^"]+)"', o.cellml):
                if c not in used:
                    note(f'{o.name}: CellML component {c} has no config entry: not moved (it stays in git history)')
        # every old config entry and parameter row accounted for
        mapped = {FROM_RE.match(i['from']).group('mt', 'st') for _, _, _, i in versions if FROM_RE.match(i['from'])}
        for o in olds.values():
            for e in o.config:
                if e.get('module_format') != 'supermodule' and (e['module_type'], e['module_subtype']) not in mapped:
                    note(f'{o.name}: config entry ({e["module_type"]}, {e["module_subtype"]}) is not in the map: dropped')
            for r in o.params:
                if r['vessel_type'] != 'global' and (r['vessel_type'], r['BC_type']) not in mapped:
                    note(f'{o.name}: parameter row {r["vessel_type"]}/{r["BC_type"]}/{r["variable_name"]} belongs to no version: dropped')
        n_sys = move_system_models(os.path.join(src, 'system'), os.path.join(stage, 'system_models'), table)

        # remove the old layout (only when rebuilding in place) and move the new tree in
        if src == MODULES:
            for n, o in olds.items():
                if os.path.isdir(o.dir):
                    shutil.rmtree(o.dir)
            for d in ('supermodules',):
                shutil.rmtree(os.path.join(MODULES, d), ignore_errors=True)
            tmpl = os.path.join(MODULES, 'template', 'template_modules.cellml')
            if os.path.isfile(tmpl):
                shutil.move(tmpl, TEMPLATE_DEST)
                shutil.rmtree(os.path.join(MODULES, 'template'), ignore_errors=True)
            for cat in os.listdir(os.path.join(MODULES, 'system')) if os.path.isdir(os.path.join(MODULES, 'system')) else []:
                if os.path.join('system', cat) not in LEFT_IN_PLACE:
                    shutil.rmtree(os.path.join(MODULES, 'system', cat))
            # emptied category dirs (e.g. modules/cell/ after its old module dirs went)
            for root, dirs, files in os.walk(MODULES, topdown=False):
                if root != MODULES and not os.listdir(root):
                    os.rmdir(root)
        else:
            if os.path.isfile(os.path.join(src, 'template', 'template_modules.cellml')):
                shutil.copy2(os.path.join(src, 'template', 'template_modules.cellml'), TEMPLATE_DEST)
        for top in os.listdir(stage_modules):
            target = os.path.join(MODULES, top)
            if os.path.isdir(target):
                _merge_tree(os.path.join(stage_modules, top), target)
            else:
                shutil.move(os.path.join(stage_modules, top), target)
        if os.path.isdir(SYSTEM_MODELS):
            shutil.rmtree(SYSTEM_MODELS)
        shutil.move(os.path.join(stage, 'system_models'), SYSTEM_MODELS)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
    print(f'{built} versions built; {n_sys} system models moved to system_models/')
    print(f'left in place (another session\'s work): {", ".join("modules/" + p for p in LEFT_IN_PLACE)}')
    if notes:
        print(f'{len(notes)} notes (above)')
    return 0


def _merge_tree(src, dst):
    '''Moves src into dst, replacing whatever version directories dst already has.'''
    for name in os.listdir(src):
        s, d = os.path.join(src, name), os.path.join(dst, name)
        if os.path.isdir(s) and os.path.isdir(d):
            if name == 'versions' or os.path.basename(dst) == 'versions':
                shutil.rmtree(d)
                shutil.move(s, d)
            else:
                _merge_tree(s, d)
        else:
            if os.path.isdir(d):
                shutil.rmtree(d)
            elif os.path.exists(d):
                os.remove(d)
            shutil.move(s, d)


if __name__ == '__main__':
    sys.exit(main())
