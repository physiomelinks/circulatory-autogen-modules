"""
The 2026-10-06 library restructure: module_types named by mechanism, versions named by source, no nn
prefixes. Approved on reviews/library_restructure_proposal.html (tools/library_restructure_proposal.py),
with the owner's decisions of 2026-10-06:

  1. the page's rename/merge table as is (phase-2 merges deferred);
  2. ports may differ between versions of one type (documented in the versions' notes);
  3. no nn_ prefixes: nn_<x> -> <x>; a plain nn version is named by its source, <Source><Year>_vXX
     (ArgusUNPUBLISHED_v01 for owner-modified or software-only modules), from the page's
     required_citations proposal; vessel modules keep their BC prefix (vv/vp/pv/pp);
  4. a new type cell/reversal_potentials for Paci's reversal potentials.

    python tools/restructure_modules.py --dry-run      # print the plan, change nothing
    python tools/restructure_modules.py                # do it (git mv, then rewrite the references)

What it does, per renamed version (RENAMES below):
  - git mv modules/<old type>/versions/<old version> -> modules/<new type>/versions/<new version>
    (untracked results/, plots/ move with it); the files <old type>_<old version>_* inside are
    renamed to <new type>_<new version>_* (tracked ones with git mv; generated .html renamed, stale
    .omex archives deleted);
  - its config entry: module_type, module_subtype, component_file (the CellML component keeps its
    name, so generated models are unchanged); its CellML <model name>; tests.yaml and
    verification_config.json module_type / version;
  - emptied old type directories (and the category cell/cardiomyocytes) are removed.
Then, across the library:
  - supermodule submodule records and verification_config harness vessel_array rows;
  - system_models/*/*/*_vessel_array.json records (reference/ copies are circulatory_autogen's
    originals, left as they are); vessel (record) names do not change, so no parameter row does;
  - text that names a renamed version: "<type>/<version>" keys, "<type> <version>" in prose,
    "(<type>, <version>)" pairs and
    "<type path>/versions/<version>" paths with their "<type>_<version>_" file stems, in tests.yaml,
    config descriptions, verification configs, SOURCES.md, source_figures.json, reviews/*.yaml and
    system_models/*/*/*_system.yaml;
  - tools/restructure_modules_renames.json: old pair -> new pair, read by
    cam_testing.library.legacy_renames() so older vessel arrays keep resolving.

Idempotent: a rename whose source is gone and whose target exists is skipped. Another session's
uncommitted work (modules/system, modules/poiseuille_transport) is never touched.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from cam_testing.library import dumps_json  # noqa: E402

MODULES = os.path.join(REPO, 'modules')
SYSTEM_MODELS = os.path.join(REPO, 'system_models')
REVIEWS = os.path.join(REPO, 'reviews')
RENAMES_JSON = os.path.join(REPO, 'tools', 'restructure_modules_renames.json')
EXCLUDED = ('modules/system', 'modules/poiseuille_transport')   # another session's uncommitted work
GENERATED = ('results', 'plots')

ARGUS_UNP = 'ArgusUNPUBLISHED_v01'

# ================================================================================== the rename table
# (old type path, old version, new type path, new version, why). Type paths are under modules/.
PAGE_TABLE = [
    # ---- cell: one type per mechanism, versions by maths (page table, "merge")
    ('cell/cardiomyocytes/Ca_dynamics_Paci_2013', 'nn', 'cell/Ca_handling', 'cardiomyocyte_Paci2013_v01', 'page: merge'),
    ('cell/neuron/soma/SN_Ca_handling_soma', 'nn', 'cell/Ca_handling', 'SN_soma_Argus2026_v01', 'page: merge'),
    ('cell/neuron/varicosity/SN_Ca_handling_varicosity', 'nn', 'cell/Ca_handling', 'SN_varicosity_Argus2026_v01', 'page: merge'),
    ('cell/cardiomyocytes/myocyte_membrane_voltage', 'nn', 'cell/membrane_potential', 'cardiomyocyte_Paci2013_v01', 'page: merge'),
    ('cell/neuron/sympathetic_neuron_membrane_voltage', 'nn', 'cell/membrane_potential', 'SN_Tao2011_v01', 'page: merge'),
    ('cell/neuron/soma/SN_membrane_soma', 'nn', 'cell/membrane_potential', 'SN_soma_Argus2026_v01', 'page: merge'),
    ('cell/neuron/varicosity/SN_varicosity_membrane', 'nn', 'cell/membrane_potential', 'SN_varicosity_Argus2026_v01', 'page: merge'),
    ('cell/cardiomyocytes/Na_dynamics_Paci_2013', 'nn', 'cell/ion_concentrations', 'cardiomyocyte_Paci2013_v01', 'page: merge'),
    ('cell/neuron/soma/SN_Na_K_concentrations_soma', 'nn', 'cell/ion_concentrations', 'SN_soma_Argus2026_v01', 'page: merge'),
    ('cell/cardiomyocytes/electric_potentials_Paci_2013', 'nn', 'cell/reversal_potentials', 'cardiomyocyte_Paci2013_v01',
     'page: merge (new type, decision 4)'),
    ('cell/neuron/NE_release_Tao_2011', 'nn', 'cell/neurotransmitter_release', 'SN_varicosity_Tao2011_v01', 'page: merge'),
    ('cell/neuron/varicosity/SN_NE_release', 'nn', 'cell/neurotransmitter_release', 'SN_varicosity_Argus2026_v01', 'page: merge'),
    ('cell/NKE_pump', 'nn', 'cell/ion_channels/i_NaK', ARGUS_UNP, 'page: merge'),
    ('cell/ion_channels/i_piezo_Na', 'Argus_unpublished_v01', 'cell/ion_channels/i_piezo_Na', ARGUS_UNP, 'page: version'),
    # ---- boundary conditions
    ('boundary_conditions/electrical_stim_input_Paci_2013', 'nn', 'boundary_conditions/i_stim_periodic',
     'cardiomyocyte_Paci2013_v01', 'page: merge'),
    # ---- control: surname-free types; all Gee 2023 modules are ArgusUNPUBLISHED for now
    *[(f'control/{n}_Gee2023', 'nn', f'control/{n}', ARGUS_UNP, 'page: rename')
      for n in ('DMV', 'LCN', 'NA', 'NActr', 'NTS', 'PNDMV', 'PNNA', 'cardiopulmonary_receptor', 'lung_stretch_receptor',
                'peripheral_resistance_effector', 'sympathetic_afferent', 'sympathetic_efferent',
                'unstressed_volume_effector')],
    ('control/lung_stretch_receptor_with_threshold_OLD', 'nn', 'control/lung_stretch_receptor', 'ArgusUNPUBLISHED_OLD',
     'page: _OLD type folded in'),
    # the page's nn_lv_/nn_rv_ spelling without nn_ (decision 3; libcuflynx reads [:2] only for vessels)
    ('control/heart_effector_Gee2023', 'nn_lv', 'control/heart_effector', 'lv_ArgusUNPUBLISHED_v01', 'page: rename (nn_ dropped)'),
    ('control/heart_effector_Gee2023', 'nn_rv', 'control/heart_effector', 'rv_ArgusUNPUBLISHED_v01', 'page: rename (nn_ dropped)'),
    ('control/heart_effector_Gee2023_OLD', 'nn_rv', 'control/heart_effector', 'rv_ArgusUNPUBLISHED_OLD',
     'page: _OLD type folded in (nn_ dropped)'),
    ('control/efferent_resistance_effector_Ursino', 'nn', 'control/efferent_resistance_effector', 'Ursino1998_v01', 'page: merge'),
    ('control/efferent_respiratory_effector_Gee_Argus_2026', 'nn', 'control/efferent_respiratory_effector', 'GeeArgus2026_v01',
     'page: merge'),
    # ---- respiratory
    ('respiratory/simple_lung_bg', 'nn', 'respiratory/lung_mechanics', ARGUS_UNP, 'page: merge'),
    ('respiratory/lung_bg_Argus_Chase', 'nn', 'respiratory/lung_mechanics', 'ventilated_ArgusUNPUBLISHED_v01', 'page: merge'),
    # ---- _2 (SI duplicates)
    ('coupling/sum_blood_vol_2', 'nn', 'coupling/sum_blood_vol', ARGUS_UNP + '_SI',
     'page: _2 type folded in as _SI (page nn_SI; the base nn is now ArgusUNPUBLISHED_v01)'),
    *[('vessels/compartments/simple_2', v, 'vessels/compartments/simple', f'{v}_SI', 'page: _2 type folded in as _SI')
      for v in ('pp_noI', 'pv_noI', 'vp_noI', 'vv_noI')],
    *[('vessels/junctions/Min_junction_2', v, 'vessels/junctions/Min_junction', f'{v}_SI', 'page: _2 type folded in as _SI')
      for v in ('vp_noI', 'vv_noI')],
    ('vessels/junctions/MinNout_junction_2', 'vv_noI', 'vessels/junctions/MinNout_junction', 'vv_noI_SI',
     'page: _2 type folded in as _SI'),
    *[('vessels/junctions/Nout_junction_2', v, 'vessels/junctions/Nout_junction', f'{v}_SI', 'page: _2 type folded in as _SI')
      for v in ('pv_noI', 'vv_noI')],
    ('heart', 'vp_simple_2', 'heart', 'vp_simple_SI', 'page: version (_2 -> _SI)'),
]

# Decision 3, nn_<x> -> <x> (type unchanged).
NN_PREFIX = [
    ('boundary_conditions/inlet_flow', v) for v in
    ('nn_adan', 'nn_adan_2', 'nn_aorticbif', 'nn_constant', 'nn_constant_2', 'nn_constant_flow_port', 'nn_fourier_8')
] + [
    ('boundary_conditions/inlet_pressure', 'nn_constant'),
    ('boundary_conditions/muscle_pressure', 'nn_constant'),
    ('boundary_conditions/outlet_flow', 'nn_constant'),
    ('boundary_conditions/outlet_flow', 'nn_constant_flow_port'),
    ('boundary_conditions/outlet_pressure', 'nn_constant'),
    ('control/efferent_heart_elastance_effector', 'nn_lv'),
    ('control/efferent_heart_elastance_effector', 'nn_rv'),
    ('heart/cardiac_clock', 'nn_controlled'),
    ('respiratory/pulmonary_GE', 'nn_5_lobe'),
    ('vessels/microvasculature/P_outlet', 'nn_micro'),
]

# Decision 3, a plain nn version named by its source: the first key of the page's required_citations
# proposal (CITATION_RULES in tools/library_restructure_proposal.py); ArgusUNPUBLISHED where the
# proposal is the owner's placeholder (argus2026unpublished: owner-modified) or only the
# circulatory_autogen software citation (generic modules with no model paper).
NN_SOURCE = {
    'benchmarks/delay_test': (ARGUS_UNP, 'software citation only (test fixture)'),
    'benchmarks/FitzHugh_Nagumo': ('FitzHugh1961_v01', 'fitzhugh1961impulses'),
    'benchmarks/Lotka_Volterra': ('Lotka1925_v01', 'lotka1925elements'),
    'benchmarks/TwoMinima': (ARGUS_UNP, 'software citation only (test fixture)'),
    'benchmarks/VanDerPol': ('VanDerPol1926_v01', 'vanderpol1926relaxation'),
    'boundary_conditions/CO2_O2_generator': (ARGUS_UNP, 'software citation only'),
    'boundary_conditions/i_stim_periodic': (ARGUS_UNP, 'software citation only'),
    'boundary_conditions/lung_vol_generator': (ARGUS_UNP, 'software citation only'),
    'boundary_conditions/voltage_change': (ARGUS_UNP, 'software citation only'),
    'boundary_conditions/zero_flow': (ARGUS_UNP, 'software citation only'),
    'control/afferent_to_respiratory_efferent': ('Albanese2016_v01', 'albanese2016integrated'),
    'control/afferent_to_sympathetic_efferent': ('Ursino1998_v01', 'ursino1998interaction'),
    'control/afferent_to_vagal_efferent': ('Ursino1998_v01', 'ursino1998interaction'),
    'control/baroreceptor': ('Ursino1998_v01', 'ursino1998interaction'),
    'control/chemoreceptor': ('Ursino2000_v01', 'ursino2000acute'),
    'control/efferent_heart_period_effector': ('Ursino1998_v01', 'ursino1998interaction'),
    'control/efferent_resistance_effector': (ARGUS_UNP, "owner's saturating form (ursino1998interaction, argus2026unpublished)"),
    'control/efferent_respiratory_effector': ('Albanese2016_v01', 'albanese2016integrated'),
    'control/efferent_venous_us_volume_and_compliance_effector': (ARGUS_UNP, "owner's saturating form"),
    'control/flow_observer': (ARGUS_UNP, 'software citation only'),
    'control/local_flow_control_brain': ('Ursino2000_v01', 'ursino2000acute'),
    'control/local_flow_control_muscle': ('Ursino2000_v01', 'ursino2000acute'),
    'control/pid_control': (ARGUS_UNP, 'software citation only'),
    'control/pressure_observer': (ARGUS_UNP, 'software citation only'),
    'coupling/imposter_1D': (ARGUS_UNP, 'software citation only'),
    'coupling/sum_blood_vol': (ARGUS_UNP, 'software citation only'),
    'coupling/volume_sum': (ARGUS_UNP, 'software citation only'),
    'heart/cardiac_clock': ('Liang2009_v01', 'liang2009multi'),
    'organs/BG_volume_control': (ARGUS_UNP, 'argus2026unpublished'),
    'organs/BVC_GI': (ARGUS_UNP, 'argus2026unpublished'),
    'organs/BVC_Kidney_left': (ARGUS_UNP, 'argus2026unpublished'),
    'organs/BVC_Kidney_right': (ARGUS_UNP, 'argus2026unpublished'),
    'organs/BVC_Na_circ': (ARGUS_UNP, 'argus2026unpublished'),
    'respiratory/gas_transport_simple': ('Albanese2016_v01', 'albanese2016integrated'),
    'respiratory/GE_capillary': ('Albanese2016_v01', 'albanese2016integrated'),
    'respiratory/GE_capillary_merge': ('Albanese2016_v01', 'albanese2016integrated'),
    'respiratory/pulmonary_GE': ('Albanese2016_v01', 'albanese2016integrated'),
    'respiratory/respiratory_gas_transport': ('Albanese2016_v01', 'albanese2016integrated'),
    'respiratory/thoracic_cavity': ('Albanese2016_v01', 'albanese2016integrated'),
    'respiratory/tissue_GE_simple': ('Ursino2001_v01', 'ursino2001integrated'),
    'transport/tissue_diffusion_face': ('Fang2008_v01', 'fang2008oxygen'),
    'transport/tissue_diffusion_volume': ('Fang2008_v01', 'fang2008oxygen'),
    'vessels/junctions/flow_sum_2': ('Safaei2016_v01', 'safaei2016roadmap'),
    **{f'vessels/properties/{n}': ('Olufsen2000_v01', 'olufsen2000numerical') for n in (
        'material_prop_const', 'material_prop_const_ven', 'material_prop_fromPWV', 'material_prop_fromPWV_ven',
        'material_prop_visco_const', 'material_prop_visco_const_ven', 'material_prop_visco_fromPWV',
        'material_prop_visco_fromPWV_ven')},
}

# nn versions kept: libcuflynx writes these pairs itself (ModelParsers adds a ('FV1D_volume_sum', 'nn')
# volume_sum_1D record to every 1D-coupled model; convert_0d_to_1d writes ('FV1D_vessel', 'nn')).
KEEP_NN = {
    'coupling/FV1D_vessel/nn': "libcuflynx's 1D coupling writes FV1D_vessel records with BC_type nn",
    'coupling/FV1D_volume_sum/nn': "libcuflynx's ModelParsers adds a ('FV1D_volume_sum', 'nn') record itself",
}

# Port tables for the merged types go into each version's notes (decision 2).
MERGED_TYPES = ('Ca_handling', 'membrane_potential', 'ion_concentrations', 'neurotransmitter_release', 'i_NaK',
                'i_stim_periodic', 'heart_effector', 'lung_stretch_receptor', 'lung_mechanics',
                'efferent_resistance_effector', 'efferent_respiratory_effector')
PORTS_NOTE_MARK = 'Ports of the versions of '


def rename_table():
    '''[(old type path, old version, new type path, new version, why)], the whole move.'''
    rows = list(PAGE_TABLE)
    rows += [(t, v, t, v[3:], 'decision 3: nn_ prefix dropped') for t, v in NN_PREFIX]
    rows += [(t, 'nn', t, new, f'decision 3: plain nn named by source ({basis})') for t, (new, basis) in NN_SOURCE.items()]
    return rows


# ================================================================================== helpers
def rel(p):
    return os.path.relpath(p, REPO)


def excluded(p):
    r = rel(p)
    return any(r == e or r.startswith(e + '/') for e in EXCLUDED)


def vdir(type_path, version):
    return os.path.join(MODULES, type_path, 'versions', version)


def tname(type_path):
    return type_path.split('/')[-1]


def all_version_dirs():
    out = []
    for root, dirs, _ in os.walk(MODULES):
        if excluded(root):
            dirs[:] = []
            continue
        dirs[:] = sorted(d for d in dirs if d not in GENERATED and not d.startswith('.'))
        if os.path.basename(root) == 'versions':
            # nested module_types sit beside versions/, never inside it
            out += [(os.path.relpath(os.path.dirname(root), MODULES), d) for d in dirs]
            dirs[:] = []
    return out


class Git:
    def __init__(self, dry):
        self.dry = dry
        self._tracked = None

    def tracked(self):
        if self._tracked is None:
            out = subprocess.run(['git', 'ls-files', '-z', 'modules', 'system_models', 'reviews'], cwd=REPO, check=True,
                                 capture_output=True).stdout.decode()
            self._tracked = {os.path.join(REPO, p) for p in out.split('\0') if p}
        return self._tracked

    def mv(self, src, dst):
        if self.dry:
            return
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        subprocess.run(['git', 'mv', src, dst], cwd=REPO, check=True)
        self._tracked = None


def json_format(text):
    '''The formatter that reproduces a JSON file's text (dumps_json or json.dumps indent 2).'''
    data = json.loads(text)
    for f in (lambda d: dumps_json(d), lambda d: json.dumps(d, indent=2), lambda d: json.dumps(d, indent=2, ensure_ascii=False)):
        out = f(data)
        if out == text or out + '\n' == text:
            return lambda d, f=f, nl=text.endswith('\n') and not out.endswith('\n'): f(d) + ('\n' if nl else '')
    return lambda d: dumps_json(d) + '\n'


def rewrite_json(path, fn, dry, changed):
    with open(path) as f:
        text = f.read()
    fmt = json_format(text)
    data = json.loads(text)
    if fn(data):
        changed.append(rel(path))
        if not dry:
            with open(path, 'w') as f:
                f.write(fmt(data))


def yaml_scalar(v):
    return yaml.safe_dump(v, default_flow_style=True, width=10 ** 6).strip().removesuffix('...').strip()


def set_top_yaml_key(text, key, value):
    '''Replaces a top-level "key: scalar" line of a YAML text.'''
    new, n = re.subn(rf'(?m)^{re.escape(key)}:[^\n]*$', f'{key}: {yaml_scalar(value)}', text, count=1)
    return new if n else text


def top_yaml_block(text, key):
    '''(start, end) of a top-level key's block (its line and its indented continuation), or None.'''
    m = re.search(rf'(?m)^{re.escape(key)}:', text)
    if not m:
        return None
    nxt = re.search(r'(?m)^\S', text[m.end():])
    return m.start(), (m.end() + nxt.start()) if nxt else len(text)


def set_top_yaml_block(text, key, value, after=('reviewed', 'version')):
    dumped = yaml.safe_dump({key: value}, sort_keys=False, width=120, default_flow_style=False, allow_unicode=True)
    span = top_yaml_block(text, key)
    if span:
        return text[:span[0]] + dumped + text[span[1]:]
    for a in after:
        span = top_yaml_block(text, a)
        if span:
            return text[:span[1]] + dumped + text[span[1]:]
    return text + dumped


# ================================================================================== planning
def plan(rows):
    '''Checks the table against the tree; returns (todo rows, skipped rows, problems).'''
    have = set(all_version_dirs())
    todo, skipped, problems = [], [], []
    seen_src = set()
    for r in rows:
        ot, ov, nt, nv, _ = r
        if (ot, ov) in seen_src:
            problems.append(f'{ot}/{ov} is renamed twice')
        seen_src.add((ot, ov))
        if (ot, ov) in have:
            if os.path.exists(vdir(nt, nv)):
                problems.append(f'{nt}/{nv} exists already (target of {ot}/{ov})')
            todo.append(r)
        elif os.path.isdir(vdir(nt, nv)):
            skipped.append(r)
        else:
            problems.append(f'{ot}/{ov} not found and {nt}/{nv} does not exist')
    # the library after the move: pairs unique, type names unique, nothing left with nn
    post = {}
    moved = {(ot, ov): (nt, nv) for ot, ov, nt, nv, _ in rows}
    for t, v in sorted(have):
        nt, nv = moved.get((t, v), (t, v))
        post.setdefault((tname(nt), nv), []).append(f'{nt}/{nv}')
    for pair, where in post.items():
        if len(where) > 1:
            problems.append(f'after the move {pair} would be {where}')
    type_paths = {}
    for (t, v) in sorted(have):
        nt, _ = moved.get((t, v), (t, v))
        type_paths.setdefault(tname(nt), set()).add(nt)
    for n, ps in type_paths.items():
        if len(ps) > 1:
            problems.append(f'after the move module_type {n} would be at {sorted(ps)}')
    for (t, v), _ in [((t, v), 0) for t, v in sorted(have)]:
        nt, nv = moved.get((t, v), (t, v))
        if nv.startswith('nn') and f'{nt}/{nv}' not in KEEP_NN:
            problems.append(f'{nt}/{nv}: an nn version not in the table')
    return todo, skipped, problems


def pair_map(rows):
    return {(tname(ot), ov): (tname(nt), nv) for ot, ov, nt, nv, _ in rows}


# ================================================================================== the move
def move_version(git, row, dry, log):
    ot, ov, nt, nv, _ = row
    src, dst = vdir(ot, ov), vdir(nt, nv)
    old_stem, new_stem = f'{tname(ot)}_{ov}', f'{tname(nt)}_{nv}'
    log.append(f'git mv {rel(src)} -> {rel(dst)}')
    git.mv(src, dst)
    base = src if dry else dst
    tracked = git.tracked() if not dry else set()
    for root, dirs, files in os.walk(base):
        dirs[:] = [d for d in dirs if d not in GENERATED]
        for f in sorted(files):
            if not (f.startswith(old_stem + '_') or f.startswith(old_stem + '.')):
                continue
            p = os.path.join(root, f)
            q = os.path.join(root, new_stem + f[len(old_stem):])
            if f.endswith('.omex'):
                log.append(f'  delete stale {f}')
                if not dry:
                    os.remove(p)
                continue
            log.append(f'  {"git mv" if (dry or p in tracked) else "mv"} {f} -> {os.path.basename(q)}')
            if dry:
                continue
            if p in tracked:
                git.mv(p, q)
            else:
                os.rename(p, q)
    if dry:
        return
    # own config entry, CellML model name, spec identity
    cfg = os.path.join(dst, f'{new_stem}_modules_config.json')

    def fix_cfg(entries):
        e = entries[0]
        e['module_type'], e['module_subtype'] = tname(nt), nv
        if e.get('component_file') == f'{old_stem}_modules.cellml':
            e['component_file'] = f'{new_stem}_modules.cellml'
        return True
    rewrite_json(cfg, fix_cfg, dry, [])
    cellml = os.path.join(dst, f'{new_stem}_modules.cellml')
    if os.path.isfile(cellml):
        txt = open(cellml).read()
        new = txt.replace(f'<model name="{old_stem}"', f'<model name="{new_stem}"', 1)
        if new != txt:
            open(cellml, 'w').write(new)
    tests = os.path.join(dst, f'{new_stem}_tests.yaml')
    if os.path.isfile(tests):
        txt = open(tests).read()
        txt = set_top_yaml_key(set_top_yaml_key(txt, 'module_type', tname(nt)), 'version', nv)
        open(tests, 'w').write(txt)
    vc = os.path.join(dst, f'{new_stem}_verification_config.json')
    if os.path.isfile(vc):
        def fix_vc(d):
            d['module_type'], d['version'] = tname(nt), nv
            return True
        rewrite_json(vc, fix_vc, dry, [])


def remove_empty_types(rows, dry, log):
    '''Old type directories left with no versions (only generated files) go, then empty categories.'''
    for ot in sorted({r[0] for r in rows if r[0] != r[2]}, key=lambda p: -p.count('/')):
        d = os.path.join(MODULES, ot)
        vroot = os.path.join(d, 'versions')
        if not os.path.isdir(d):
            continue
        if os.path.isdir(vroot) and [x for x in os.listdir(vroot) if not x.startswith('.')]:
            continue
        others = [x for x in os.listdir(d) if os.path.isdir(os.path.join(d, x)) and x != 'versions']
        if others:
            continue
        log.append(f'remove emptied type {rel(d)} ({", ".join(sorted(os.listdir(d))) or "empty"})')
        if not dry:
            shutil.rmtree(d)
    for cat in ('cell/cardiomyocytes',):
        d = os.path.join(MODULES, cat)
        if os.path.isdir(d) and not [x for x in os.listdir(d) if not x.startswith('.')]:
            log.append(f'remove emptied category {rel(d)}')
            if not dry:
                os.rmdir(d)


# ================================================================================== references
def map_record(rec, pm, tkey='module_type', vkey='module_subtype'):
    k = (rec.get(tkey), rec.get(vkey))
    if k in pm:
        rec[tkey], rec[vkey] = pm[k]
        return True
    return False


def rewrite_supermodules_and_harnesses(pm, dry, changed):
    for root, dirs, files in os.walk(MODULES):
        if excluded(root):
            dirs[:] = []
            continue
        dirs[:] = [d for d in dirs if d not in GENERATED]
        for f in files:
            p = os.path.join(root, f)
            if f.endswith('_modules_config.json'):
                def fix(entries):
                    hit = False
                    for e in entries:
                        for s in e.get('submodules') or []:
                            hit |= map_record(s, pm)
                    return hit
                rewrite_json(p, fix, dry, changed)
            elif f.endswith('_verification_config.json'):
                def fix_h(d):
                    hit = False
                    for row in ((d.get('harness') or {}).get('vessel_array') or []):
                        if isinstance(row, list) and len(row) > 2 and (row[2], row[1]) in pm:
                            row[2], row[1] = pm[(row[2], row[1])]
                            hit = True
                        elif isinstance(row, dict):
                            hit |= map_record(row, pm) or map_record(row, pm, 'vessel_type', 'BC_type')
                    return hit
                rewrite_json(p, fix_h, dry, changed)


def rewrite_vessel_arrays(pm, dry, changed):
    '''One record per line: each line is edited in place, so the files keep their layout.'''
    for cat in sorted(os.listdir(SYSTEM_MODELS)):
        cdir = os.path.join(SYSTEM_MODELS, cat)
        if not os.path.isdir(cdir):
            continue
        for model in sorted(os.listdir(cdir)):
            p = os.path.join(cdir, model, f'{model}_vessel_array.json')
            if not os.path.isfile(p):
                continue
            lines = open(p).read().split('\n')
            hit = False
            for i, line in enumerate(lines):
                s = line.strip().rstrip(',')
                if not (s.startswith('{') and s.endswith('}')):
                    continue
                rec = json.loads(s)
                for tk, vk in (('module_type', 'module_subtype'), ('vessel_type', 'BC_type')):
                    k = (rec.get(tk), rec.get(vk))
                    if k in pm:
                        nt, nv = pm[k]
                        line = line.replace(f'"{tk}": "{k[0]}"', f'"{tk}": "{nt}"', 1)
                        line = line.replace(f'"{vk}": "{k[1]}"', f'"{vk}": "{nv}"', 1)
                        lines[i] = line
                        hit = True
            if hit:
                new = '\n'.join(lines)
                json.loads(new)
                changed.append(rel(p))
                if not dry:
                    open(p, 'w').write(new)
            # a multi-line layout would need the generic path: check nothing was missed
            for rec in json.load(open(p)) if not dry else []:
                k = (rec.get('module_type'), rec.get('module_subtype'))
                assert k not in pm, f'{rel(p)}: record {rec.get("name")} still names {k}'


def text_patterns(rows):
    '''[(compiled regex, replacement)] for text naming a renamed version.'''
    pats = []
    for ot, ov, nt, nv, _ in rows:
        a, b = tname(ot), tname(nt)
        W0, W1 = r'(?<![A-Za-z0-9_])', r'(?![A-Za-z0-9_])'
        pats.append((re.compile(rf'{W0}{re.escape(ot)}/versions/{re.escape(ov)}{W1}'), f'{nt}/versions/{nv}'))
        pats.append((re.compile(rf'{W0}{re.escape(a)}_{re.escape(ov)}_(?=(modules|units|tests|verification_config|references|risk)\b)'),
                     f'{b}_{nv}_'))
        pats.append((re.compile(rf'{W0}{re.escape(a)}/{re.escape(ov)}{W1}'), f'{b}/{nv}'))
        pats.append((re.compile(rf'{W0}{re.escape(a)} {re.escape(ov)}{W1}'), f'{b} {nv}'))     # prose: "inlet_flow nn_constant"
        pats.append((re.compile(rf'\({re.escape(a)}, {re.escape(ov)}\)'), f'({b}, {nv})'))
        pats.append((re.compile(rf'\({re.escape(a)}, {re.escape(ov)}, '), f'({b}, {nv}, '))
    return pats


def text_files():
    out = []
    for root, dirs, files in os.walk(MODULES):
        if excluded(root):
            dirs[:] = []
            continue
        dirs[:] = [d for d in dirs if d not in GENERATED]
        for f in files:
            if f.endswith(('_tests.yaml', '_modules_config.json', '_verification_config.json', 'SOURCES.md',
                           'source_figures.json')):
                out.append(os.path.join(root, f))
    out += [os.path.join(REVIEWS, f) for f in sorted(os.listdir(REVIEWS)) if f.endswith('.yaml')]
    for cat in sorted(os.listdir(SYSTEM_MODELS)):
        cdir = os.path.join(SYSTEM_MODELS, cat)
        if os.path.isdir(cdir):
            for model in sorted(os.listdir(cdir)):
                p = os.path.join(cdir, model, f'{model}_system.yaml')
                if os.path.isfile(p):
                    out.append(p)
    return out


def rewrite_text(rows, dry, changed):
    pats = text_patterns(rows)
    for p in text_files():
        txt = open(p).read()
        new = txt
        for rx, rep in pats:
            new = rx.sub(rep, new)
        if new != txt:
            if p.endswith('.json'):
                json.loads(new)
            changed.append(rel(p))
            if not dry:
                open(p, 'w').write(new)


# ================================================================================== port notes
def port_notes(dry, changed):
    '''Each version of a merged type gets the type's port table in its notes (decision 2).'''
    from cam_testing import library
    library._discover.cache_clear()
    for name in MERGED_TYPES:
        mt = library.load_module_type(name)
        versions = mt.versions()
        if len(versions) < 2:
            continue
        ports = {}
        for v in versions:
            p = {}
            for d in ('entrance_ports', 'exit_ports', 'general_ports'):
                for port in v.raw_entry.get(d) or []:
                    p[f'{d.split("_")[0]} {port["port_type"]}'] = ', '.join(port.get('variables') or [])
            ports[v.name] = p
        keys = sorted({k for p in ports.values() for k in p})
        same = all(set(p) == set(keys) for p in ports.values())
        lines = [f'{PORTS_NOTE_MARK}{name} ({"the same port types" if same else "they differ; allowed and documented, see modules/README.md, Versions"}):']
        for k in keys:
            lines.append(f'{k}: ' + '; '.join(f'{v} [{ports[v][k]}]' if k in ports[v] else f'{v} -' for v in ports))
        table = '\n'.join(lines)
        for v in versions:
            path = v.spec_path
            txt = open(path).read()
            tests = yaml.safe_load(txt) or {}
            notes = tests.get('notes') or ''
            if PORTS_NOTE_MARK in notes:
                notes = notes[:notes.index(PORTS_NOTE_MARK)].rstrip()
            notes = (notes + '\n\n' if notes else '') + table
            changed.append(rel(path))
            if not dry:
                open(path, 'w').write(set_top_yaml_block(txt, 'notes', notes))


# ================================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--dry-run', action='store_true', help='print the plan and the files that would change; change nothing')
    args = ap.parse_args()
    dry = args.dry_run
    rows = rename_table()
    todo, skipped, problems = plan(rows)
    print(f'{len(rows)} renames in the table: {len(todo)} to do, {len(skipped)} done already')
    if problems:
        print('PROBLEMS:\n  ' + '\n  '.join(problems))
        sys.exit(1)
    print('\nold type/version -> new type/version')
    for ot, ov, nt, nv, why in rows:
        print(f'  {ot}/{ov} -> {nt}/{nv}   [{why}]')
    print('\nkept nn:', ', '.join(f'{k} ({v})' for k, v in KEEP_NN.items()))
    git = Git(dry)
    log = []
    for r in todo:
        move_version(git, r, dry, log)
    remove_empty_types(rows, dry, log)
    pm = pair_map(rows)
    changed = []
    rewrite_supermodules_and_harnesses(pm, dry, changed)
    rewrite_vessel_arrays(pm, dry, changed)
    rewrite_text(rows, dry, changed)
    if not dry:
        port_notes(dry, changed)
        out = [{'from': [tname(ot), ov], 'to': [tname(nt), nv], 'from_path': f'{ot}/{ov}', 'to_path': f'{nt}/{nv}', 'why': why}
               for ot, ov, nt, nv, why in rows]
        with open(RENAMES_JSON, 'w') as f:
            f.write(dumps_json(out) + '\n')
        changed.append(rel(RENAMES_JSON))
    print('\nmoves:\n  ' + '\n  '.join(log))
    print(f'\nfiles rewritten ({len(set(changed))}):\n  ' + '\n  '.join(sorted(set(changed))))
    if dry:
        print('\n(dry run: nothing changed; port notes and the renames JSON are written by the real run)')


if __name__ == '__main__':
    main()
