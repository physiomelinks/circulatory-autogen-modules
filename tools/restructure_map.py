"""
Proposes the versions/instances layout: writes tools/restructure_map.yaml, mapping every config
entry (module_type, module_subtype) of the current modules to

    modules/<category path>/<module_type>/versions/<version>/instances/<instance>/

for review before tools/restructure_modules.py moves anything. Placement rule: a module goes in
the most specific category that covers every context it belongs to (soma only in neurons ->
cell/neurons; something several cell types share -> cell). No directory name may appear twice in
the tree. The version is the entry's module_subtype (new versions: <source>_vXX).

    python tools/restructure_map.py [--source <old modules/ copy>]

The map is made from the old per-module layout; once modules/ has moved to versions, pass --source
with a copy of the old tree (the script refuses to read the new layout).

Historical: this made the first versions layout (categories cardiac and cell/neurons). The later move
to nested module_types (cell/neuron and heart as module_types holding their parts, the alternative
hearts as versions of heart) was done by hand and is not in the map; modules/README.md describes the
current layout, and cam_testing.library.HEART_VERSION_RENAMES the heart renames.
"""
import collections
import glob
import json
import os
import re
import sys

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(REPO, 'modules')
OUT = os.path.join(REPO, 'tools', 'restructure_map.yaml')
SKIP = ('system', 'supermodules', 'poiseuille_transport')   # poiseuille: another session's uncommitted work


def category(module, mt):
    '''The category path for a module_type from current module `module`.'''
    micro_only = re.compile(r'^(artery|arteriole|capillary|venule|vein)_(Minlet|Nout|Min|MinNout|Noutlet|inlet|outlet)$')
    if module == 'microvasculature_network':
        if mt in ('artery', 'arteriole', 'venule', 'vein', 'P_outlet') or micro_only.match(mt):
            return 'vessels/microvasculature'
        if mt == 'capillary':
            return 'vessels/compartments'              # also BG's capillary: shared by both contexts
    if module in ('BG', 'cvs_mvp'):
        if mt in ('terminal', 'terminal2', 'single_terminal'):
            return 'vessels/terminals'
        if 'junction' in mt or mt in ('capillary_split', 'capillary_merge', 'flow_sum_2'):
            return 'vessels/junctions'
        if mt.startswith('heart_simple'):
            return 'cardiac'
        if mt.startswith('sum_blood_vol'):
            return 'coupling'
        return 'vessels/compartments'
    if module == 'elic':
        return 'vessels/terminals' if mt.startswith('terminal') else 'control'
    if module == 'open_loop':
        return 'boundary_conditions' if mt == 'inlet_flow' else 'vessels/terminals'
    if module == 'vessel_properties':
        return 'vessels/properties'
    if module == 'cardiac':
        if mt in ('flow_merge', 'flow_split'):
            return 'vessels/junctions'
        if mt == 'resistor':
            return 'vessels/compartments'
        return 'cardiac'
    if module == 'heart':
        return 'cardiac'
    if module in ('neurons',):
        return 'cell/neurons'
    if module == 'cardiomyocytes':
        return 'cell/cardiomyocytes'
    if module == 'myocytes':
        return 'cell/cardiomyocytes'                   # myocyte_membrane_voltage: the Paci cell's integrator (to confirm)
    if module == 'ion_channel':
        return 'cell/ion_channels'                     # every channel; the cell type is in its version's notes
    if module == 'NKE_pump':
        return 'cell'
    if module in ('lung', 'gas_exchange'):
        return 'respiratory'
    if module in ('control', 'parasympathetic', 'control_system'):
        return 'control'
    if module in ('boundary_condition', 'input_stimulation'):
        return 'boundary_conditions'
    if module in ('coupling',):
        return 'coupling'
    if module == 'diffusion_volume':
        return 'transport'
    if module.startswith('BVC_') or module == 'BG_volume_control':
        return 'organs'
    if module in ('Lotka_Volterra', 'VanDerPol', 'FitzHugh_Nagumo', 'TwoMinima', 'delay'):
        return 'benchmarks'
    return 'UNMAPPED'


# module_type renames (the monolithic sympathetic neuron pieces become versions of soma/varicosity/axon)
RENAME = {('SN_soma', 'nn'): ('soma', 'sympathetic_monolithic_v01'),
          ('SN_varicosity', 'nn'): ('varicosity', 'sympathetic_monolithic_v01'),
          ('SN_axon', 'nn'): ('axon', 'sympathetic_monolithic_v01')}
# ion channels: one module_type per channel kind, one version per source; the cell type each
# version was built for goes in the config entry's "notes"
CHANNEL_SOURCES = [('_Paci_2013', 'Paci2013_v01', 'hiPSC-derived ventricular-like cardiomyocyte (Paci et al. 2013)'),
                   ('_Tao_2011', 'Tao2011_v01', 'sympathetic neuron (Tao et al. 2011)'),
                   ('_SN', 'Argus2026_v01', 'sympathetic neuron, soma or varicosity (Argus sympathetic-neuron model, '
                                           'decomposed from SN_soma / SN_varicosity)'),
                   ('_Argus_unpublished', 'Argus_unpublished_v01', 'stretch-activated Na current (Piezo), not cell-specific')]


def channel_rename(mt, st):
    for suffix, version, cell in CHANNEL_SOURCES:
        if mt.endswith(suffix):
            return mt[:-len(suffix)], version, cell
    return mt, st, None


# supermodules (modules/supermodules/<dir>) -> (category, module_type, version)
SUPERMODULES = {'heart': ('cardiac', 'heart', 'Argus2026_v01'),
                'SN_soma': ('cell/neurons', 'soma', 'sympathetic'),
                'SN_varicosity': ('cell/neurons', 'varicosity', 'sympathetic'),
                'sympathetic_neuron': ('cell/neurons', 'neuron', 'sympathetic')}


def instances_for(module_dir, cid, comp_spec):
    '''Proposed instances of a component: default (its parameter rows) + one per validation data set
    named in its spec (baseline data, calibration obs_data / params_for_id / validation_obs_data).'''
    out = [{'name': 'default', 'from': 'parameter rows in the module parameters file'}]
    v = (comp_spec or {}).get('validation') or {}
    files = []
    for kind in ('baseline', 'calibrate'):
        blk = v.get(kind) or {}
        for key in ('data', 'obs_data', 'params_for_id', 'validation_obs_data'):
            f = blk.get(key)
            if isinstance(f, str) and f not in files:
                files.append(f)
    files = [f for f in files if os.path.isfile(os.path.join(module_dir, f))]
    if files:
        stems = [os.path.splitext(os.path.basename(f))[0] for f in files]
        named = next((s for s in stems if not s.startswith(cid)), None)
        name = named or cid.split('__')[0] + '_data'
        src = os.path.join(module_dir, 'validation', 'SOURCES.md')
        out.append({'name': name, 'from': files + (['validation/SOURCES.md'] if os.path.isfile(src) else []),
                    'status': {k: (v.get(k) or {}).get('status') for k in ('baseline', 'calibrate')}})
    return out


def main(argv=None):
    import argparse
    global MODULES
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--source', help='the old modules/ directory (default: modules/)')
    MODULES = os.path.abspath(ap.parse_args(argv).source or MODULES)
    if glob.glob(os.path.join(MODULES, '**', 'versions', '*', '*_modules_config.json'), recursive=True):
        print(f'{MODULES} is in the versions layout already; the map is made from the old layout (--source)')
        return 2
    entries, component_users = [], collections.defaultdict(list)
    for cfg in sorted(glob.glob(os.path.join(MODULES, '**', '*_modules_config.json'), recursive=True)):
        rel = os.path.relpath(os.path.dirname(cfg), MODULES)
        if rel.split(os.sep)[0] in SKIP:
            continue
        module = os.path.basename(os.path.dirname(cfg))
        spec_path = os.path.join(os.path.dirname(cfg), f'{module}_tests.yaml')
        specs = {}
        if os.path.isfile(spec_path):
            for c in (yaml.safe_load(open(spec_path)) or {}).get('components', []):
                specs[(c.get('vessel_type'), c.get('BC_type'))] = c
        for e in json.load(open(cfg)):
            if e.get('module_format') == 'supermodule':
                continue
            entries.append((module, os.path.dirname(cfg), e, specs.get((e['module_type'], e['module_subtype']))))
            component_users[e['component_type']].append((e['module_type'], e['module_subtype']))

    tree = collections.defaultdict(lambda: {'category': None, 'from_modules': set(), 'versions': {}})
    problems = []
    for module, mdir, e, comp_spec in entries:
        mt, st = e['module_type'], e['module_subtype']
        new_mt, version = RENAME.get((mt, st), (mt, st))
        cell_note = None
        if module == 'ion_channel':
            new_mt, version, cell_note = channel_rename(mt, st)
        cat = category(module, mt)
        node = tree[new_mt]
        if node['category'] and node['category'] != cat:
            problems.append(f'{new_mt}: categories {node["category"]} and {cat} (from {sorted(node["from_modules"])} + {module})')
        node['category'] = node['category'] or cat
        node['from_modules'].add(module)
        if version in node['versions']:
            problems.append(f'{new_mt}: version {version} defined twice ({node["versions"][version]["from"]} and {module})')
        comp = e['component_type']
        shared = len(component_users[comp]) > 1
        node['versions'][version] = {
            'from': f'{module}: ({mt}, {st})', 'component': comp,
            **({'notes': f'Built for: {cell_note}'} if cell_note else {}),
            **({'component_renamed': f'{comp}__{new_mt}_{version}'} if shared else {}),
            'instances': instances_for(mdir, f'{mt}__{st}', comp_spec)}
    for sdir, (cat, mt, version) in SUPERMODULES.items():
        node = tree[mt]
        node['category'] = node['category'] or cat
        if node['category'] != cat:
            problems.append(f'{mt}: supermodule category {cat} vs {node["category"]}')
        node['from_modules'].add(f'supermodules/{sdir}')
        node['versions'][version] = {'from': f'supermodules/{sdir} (module_format supermodule)',
                                     'instances': [{'name': 'default', 'from': f'supermodules/{sdir} default_parameters'}]}

    # no directory name twice: categories vs module_types
    cat_parts = {p for n in tree.values() for p in n['category'].split('/')}
    problems += [f'module_type {mt} has the same name as a category' for mt in tree if mt in cat_parts]
    problems += [f'{mt}: unmapped (from {sorted(n["from_modules"])})' for mt, n in tree.items() if n['category'] == 'UNMAPPED']

    by_cat = collections.defaultdict(dict)
    for mt, n in sorted(tree.items(), key=lambda x: (x[1]['category'], x[0].lower())):
        by_cat[n['category']][mt] = {'from_modules': sorted(n['from_modules']),
                                     'versions': dict(sorted(n['versions'].items()))}
    doc = {'summary': {'module_types': len(tree), 'versions': sum(len(n['versions']) for n in tree.values()),
                       'categories': {c: len(v) for c, v in sorted(by_cat.items())},
                       'shared_components_renamed': sum(1 for c, u in component_users.items() if len(u) > 1),
                       'problems': problems},
           'questions': ['myocyte_membrane_voltage placed in cell/cardiomyocytes (the Paci cell\'s charge integrator); '
                         'cell/myocytes is then empty and not created',
                         'ion channel versions: Paci2013_v01, Tao2011_v01, Argus2026_v01 (the _SN channels), Argus_unpublished_v01 (piezo); '
                         'module_subtype changes from nn to the version name',
                         'instance names are proposed from each data set\'s file name: confirm or rename (they become obs_data_name)'],
           'categories': dict(sorted(by_cat.items()))}
    with open(OUT, 'w') as f:
        yaml.safe_dump(doc, f, sort_keys=False, width=140, allow_unicode=True)
    print(json.dumps(doc['summary'], indent=1))
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
