"""
Rewrite the Phlynx manifests for the per-module layout (modules/<name>/...).

index.json and vitalworkshop.json keep their curated module list (MANIFEST_MODULES);
all.json lists every module in modules/. Each listed module contributes its CellML file,
config and units file.

    python tools/build_manifests.py
"""
import json
import os
from datetime import datetime, timezone

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES_DIR = os.path.join(REPO_ROOT, 'modules')
MANIFESTS_DIR = os.path.join(REPO_ROOT, 'manifests')

# Display names Phlynx shows; modules not listed here use their directory name.
DISPLAY_NAMES = {'BG': 'Bond Graph', 'template': 'templates'}

CURATED = ['BG', 'heart', 'microvasculature_network', 'template']
MANIFEST_MODULES = {'index.json': CURATED, 'vitalworkshop.json': CURATED}


def all_module_names():
    return sorted((d for d in os.listdir(MODULES_DIR)
                   if os.path.isfile(os.path.join(MODULES_DIR, d, f'{d}_modules.cellml'))),
                  key=str.lower)


def _entry(name, filename):
    path = f'modules/{name}/{filename}'
    if not os.path.isfile(os.path.join(REPO_ROOT, path)):
        return None
    return {'name': DISPLAY_NAMES.get(name, name), 'file': filename, 'path': path}


def build(names):
    collections = {'modules': [], 'configs': [], 'units': []}
    for name in names:
        for key, filename in (('modules', f'{name}_modules.cellml'),
                              ('configs', f'{name}_modules_config.json'),
                              ('units', f'{name}_units.cellml')):
            entry = _entry(name, filename)
            if entry is not None:
                collections[key].append(entry)
    return {'generatedAt': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
            'collections': collections}


def main():
    manifests = dict(MANIFEST_MODULES)
    manifests['all.json'] = all_module_names()
    for filename, names in manifests.items():
        with open(os.path.join(MANIFESTS_DIR, filename), 'w') as f:
            json.dump(build(names), f, indent=2)
            f.write('\n')
        print(f'wrote manifests/{filename} ({len(names)} modules)')


if __name__ == '__main__':
    main()
