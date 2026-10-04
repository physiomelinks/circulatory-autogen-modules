"""
Rewrite the PhLynx manifests for the versions layout
(modules/<category>/<module_type>/versions/<version>/..., nested module_types inside their parent's
directory, e.g. modules/heart/chamber/versions/vv/).

all.json lists every version available to PhLynx (a config entry with "available_to_phlynx": false
is left out of every manifest); index.json and vitalworkshop.json keep a curated set of paths under
modules/ (CURATED: categories or module_types, each with what is nested in it). Each version contributes its CellML file, config and units file, named
"<module_type> (<version>)"; a supermodule version has only its config. The PhLynx templates file
(manifests/template_modules.cellml) is listed under "templates".

    python tools/build_manifests.py
"""
import json
import os
import sys
from datetime import datetime, timezone

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)
from cam_testing.library import all_versions  # noqa: E402

MANIFESTS_DIR = os.path.join(REPO_ROOT, 'manifests')
TEMPLATE = os.path.join(MANIFESTS_DIR, 'template_modules.cellml')

# paths under modules/ (prefixes) in the curated manifests: the bond-graph vessels, the heart and the
# microvasculature: what the per-module manifests had (BG, heart, microvasculature_network), by category,
# and the lumped constitutive modules (haemodynamics/) the vessels' _lumped supermodules are built from
CURATED = ['vessels', 'heart', 'haemodynamics']
MANIFEST_SELECTION = {'index.json': CURATED, 'vitalworkshop.json': CURATED, 'all.json': None}


def _rel(path):
    return os.path.relpath(path, REPO_ROOT).replace(os.sep, '/')


def _selected(versions, prefixes):
    if prefixes is None:
        return versions
    out = []
    for v in versions:
        where = v.mtype.relpath
        if any(where == p or where.startswith(p + '/') for p in prefixes):
            out.append(v)
    return out


def build(versions):
    collections = {'modules': [], 'configs': [], 'units': []}
    for v in versions:
        name = f'{v.vessel_type} ({v.name})'
        for key, path in (('modules', v.cellml_path), ('configs', v.config_path), ('units', v.units_path)):
            if os.path.isfile(path) and not (key != 'configs' and v.is_supermodule):
                collections[key].append({'name': name, 'file': os.path.basename(path), 'path': _rel(path)})
    if os.path.isfile(TEMPLATE):
        collections['modules'].append({'name': 'templates', 'file': os.path.basename(TEMPLATE), 'path': _rel(TEMPLATE)})
    return {'generatedAt': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), 'collections': collections}


def main():
    versions = sorted((v for v in all_versions() if v.available_to_phlynx),
                      key=lambda v: (v.mtype.relpath.lower(), v.name.lower()))
    for filename, prefixes in MANIFEST_SELECTION.items():
        chosen = _selected(versions, prefixes)
        with open(os.path.join(MANIFESTS_DIR, filename), 'w') as f:
            json.dump(build(chosen), f, indent=2)
            f.write('\n')
        print(f'wrote manifests/{filename} ({len(chosen)} versions)')


if __name__ == '__main__':
    main()
