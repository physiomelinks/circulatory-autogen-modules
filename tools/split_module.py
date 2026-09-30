"""
Splits one module into several, e.g. modules/cell into modules/cell/{neurons,cardiomyocytes,myocytes}/.

    python tools/split_module.py cell --into neurons=SN_soma,SN_axon,... --into myocytes=myocyte_membrane_voltage \\
        [--parent cell]      # the new modules go in modules/<parent>/<new>/ (default: the old module's directory)

Each group (config module_types) becomes a module modules/<parent>/<new>/ with its own CellML
(the components of those entries, verbatim), config, units (only those it uses), parameters
(the group's rows; globals to every group that uses them), test spec (defaults, the group's
components, and the old module's review block), bibliographies (the entries its references
cite) and risk files. The old module's files are removed; its plots/results are regenerated
by the tests. Component ids (<module_type>__<subtype>) are unchanged, so risk files move as is.

Pre-versions tool: it reads and writes the old per-module layout (modules/<name>/...), kept for
provenance. The library now uses modules/<category>/<module_type>/versions/<version>/ (see
modules/README.md); tools/restructure_modules.py moves an old-layout tree into it.
"""
import argparse
import csv
import glob
import json
import os
import re
import shutil
import sys

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'tools'))
from cam_testing import bib  # noqa: E402
from cam_testing.library import load_module, module_dir, normalise_config_entry  # noqa: E402
from import_from_libcuflynx import read_units_blocks, units_closure, write_units_file  # noqa: E402

CELLML_HEAD = ("<?xml version='1.0' encoding='UTF-8'?>\n"
               '<model name="modules" xmlns="http://www.cellml.org/cellml/1.1#" '
               'xmlns:cellml="http://www.cellml.org/cellml/1.1#">\n')


def component_block(text, name):
    m = re.search(r'(\s*<!--(?:(?!-->).)*-->)*\s*<component name="%s"[\s>].*?</component>' % re.escape(name), text, re.S)
    if not m:
        raise KeyError(f'component {name} not found')
    block = m.group(0)
    # keep a comment directly above the component only if it isn't a file-level banner
    return block.lstrip('\n')


def split(old, groups, parent_dir):
    m = load_module(old)
    src = m.dir
    with open(os.path.join(src, f'{old}_modules.cellml')) as f:
        cellml = f.read()
    with open(os.path.join(src, f'{old}_modules_config.json')) as f:
        raw_config = json.load(f)
    blocks = read_units_blocks([os.path.join(src, f'{old}_units.cellml')])
    with open(os.path.join(src, f'{old}_tests.yaml')) as f:
        spec = yaml.safe_load(f)
    with open(os.path.join(src, f'{old}_parameters.csv'), newline='') as f:
        reader = csv.DictReader(f)
        pheader, prows = reader.fieldnames, list(reader)
    bib_all = {**bib.read(bib.bib_path(m)), **bib.read(bib.proposed_bib_path(m))}   # key -> fields

    assigned = [t for types in groups.values() for t in types]
    unassigned = sorted({e['module_type'] for e in raw_config} - set(assigned))
    if unassigned:
        raise SystemExit(f'module_types not assigned to a group: {unassigned}')

    for new, types in groups.items():
        dest = os.path.join(parent_dir, new)
        os.makedirs(dest, exist_ok=True)
        entries = [e for e in raw_config if e['module_type'] in types]
        comps = []
        for e in entries:
            if e['component_type'] not in comps:
                comps.append(e['component_type'])
        with open(os.path.join(dest, f'{new}_modules.cellml'), 'w') as f:
            f.write(CELLML_HEAD + '\n'.join('    ' + component_block(cellml, c).strip() for c in comps) + '\n</model>\n')
        for e in entries:
            e['component_file'] = f'{new}_modules.cellml'
        with open(os.path.join(dest, f'{new}_modules_config.json'), 'w') as f:
            json.dump(entries, f, indent=2)
            f.write('\n')
        used = set()
        for c in comps:
            used |= set(re.findall(r'\bunits="([^"]+)"', component_block(cellml, c)))
        for e in entries:
            used |= {v[1] for v in e['variables_and_units']}
        needed, missing = units_closure(used, blocks)
        if missing:
            print(f'{new}: units not defined in {old}_units.cellml: {sorted(missing)}')
        write_units_file(os.path.join(dest, f'{new}_units.cellml'), needed, blocks)

        normal = [normalise_config_entry(e) for e in entries]
        globals_used = {v[0] for e in normal for v in e['variables_and_units'] if v[3].strip() == 'global_constant'}
        rows = [r for r in prows if r['vessel_type'] in types
                or (r['vessel_type'] == 'global' and r['variable_name'] in globals_used)]
        with open(os.path.join(dest, f'{new}_parameters.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=pheader)
            w.writeheader()
            w.writerows(rows)

        new_spec = {k: v for k, v in spec.items() if k != 'components'}
        new_spec['module'] = new
        new_spec['components'] = [c for c in spec['components'] if c['vessel_type'] in types]
        if spec.get('review'):
            new_spec['review'] = dict(spec['review'])
            new_spec['review']['comment'] = (f'[Split from the {old} module on review; this block is the whole of {old}\'s review, '
                                             f'of which only the parts about this module\'s components apply.] '
                                             + str(spec['review'].get('comment', '')))
        with open(os.path.join(dest, f'{new}_tests.yaml'), 'w') as f:
            yaml.safe_dump(new_spec, f, sort_keys=False, width=120, default_flow_style=None, allow_unicode=True)

        cited = {bib.reference_key(r.get('data_reference', '')) for r in rows}
        for c in new_spec['components']:
            for p in (c.get('reference_proposals') or {}).values():
                cited.add(bib.reference_key(str(p.get('reference', ''))))
        text = json.dumps(new_spec)
        cited |= {k for k in bib_all if k in text}
        for path_fn in (bib.bib_path, bib.proposed_bib_path):
            path = path_fn(m)
            text_bib = open(path).read() if os.path.isfile(path) else ''
            entries_raw = [e.strip() for e in re.split(r'(?m)^(?=@)', text_bib) if e.strip().startswith('@')]
            keep = [e for e in entries_raw if re.match(r'@\w+\{([^,]+),', e) and re.match(r'@\w+\{([^,]+),', e).group(1).strip() in cited]
            target = os.path.join(dest, os.path.basename(path).replace(old, new, 1))
            with open(target, 'w') as f:
                f.write('\n\n'.join(keep) + ('\n' if keep else ''))

        os.makedirs(os.path.join(dest, 'risk'), exist_ok=True)
        for c in new_spec['components']:
            cid = f"{c['vessel_type']}__{c['BC_type']}"
            for p in glob.glob(os.path.join(src, 'risk', f'{cid}_risk*')):
                shutil.copy2(p, os.path.join(dest, 'risk', os.path.basename(p)))
        print(f'{os.path.relpath(dest, REPO)}: {len(entries)} config entries, {len(comps)} components, '
              f'{len(rows)} parameters, {len(needed)} units')
    return src


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('module')
    ap.add_argument('--into', action='append', required=True, help='new=module_type,module_type,... (repeatable)')
    ap.add_argument('--parent', help='directory under modules/ for the new modules (default: the old module\'s)')
    ap.add_argument('--keep-old', action='store_true', help='keep the old module\'s files')
    args = ap.parse_args(argv)
    groups = {}
    for g in args.into:
        new, types = g.split('=', 1)
        groups[new] = [t for t in types.split(',') if t]
    parent = os.path.join(REPO, 'modules', args.parent) if args.parent else module_dir(args.module)
    src = split(args.module, groups, parent)
    if not args.keep_old:
        for p in glob.glob(os.path.join(src, f'{args.module}_*')) + glob.glob(os.path.join(src, f'{args.module}.html')):
            os.remove(p)
        for d in ('plots', 'results', 'risk'):
            shutil.rmtree(os.path.join(src, d), ignore_errors=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
