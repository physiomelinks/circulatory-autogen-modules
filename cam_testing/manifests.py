"""
Rewrite the PhLynx manifests (manifests/*.json) of the repo under test (cam_testing.paths) for the
versions layout (modules/<category>/<module_type>/versions/<version>/..., nested module_types inside
their parent's directory, e.g. modules/heart/chamber/versions/vv/).

all.json lists every version of the repo's modules/ (not the extra libraries': a manifest's paths
are relative to the repo). Curated manifests list the versions under some paths of modules/
(categories or module_types, each with what is nested in it), from pyproject.toml:

    [tool.cam_testing.manifests]
    "index.json" = ["vessels", "heart"]

Each version contributes its CellML file, config and units file, named "<module_type> (<version>)";
a supermodule version has only its config. The PhLynx templates file
(manifests/template_modules.cellml), where there is one, is listed under "templates".

    python -m cam_testing.manifests          # make manifests
"""
import json
import os
from datetime import datetime, timezone

from cam_testing import paths
from cam_testing.library import all_versions

TEMPLATE_NAME = 'template_modules.cellml'


def manifest_selection(roots=None):
    '''{manifest file name: path prefixes under modules/, or None for every version}: all.json, plus
    [tool.cam_testing.manifests].'''
    roots = roots or paths.roots()
    selection = {}
    for name, prefixes in (roots.settings.get('manifests') or {}).items():
        selection[name] = [prefixes] if isinstance(prefixes, str) else list(prefixes)
    selection['all.json'] = None
    return selection


def _rel(path):
    return os.path.relpath(path, paths.roots().repo_root).replace(os.sep, '/')


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
    template = os.path.join(paths.roots().manifests_dir, TEMPLATE_NAME)
    if os.path.isfile(template):
        collections['modules'].append({'name': 'templates', 'file': TEMPLATE_NAME, 'path': _rel(template)})
    return {'generatedAt': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), 'collections': collections}


def main(argv=None):
    import argparse
    argparse.ArgumentParser(description='Rewrite manifests/*.json for the repo under test').parse_args(argv)
    manifests_dir = paths.roots().manifests_dir
    os.makedirs(manifests_dir, exist_ok=True)
    versions = sorted(all_versions(), key=lambda v: (v.mtype.relpath.lower(), v.name.lower()))
    for filename, prefixes in manifest_selection().items():
        chosen = _selected(versions, prefixes)
        with open(os.path.join(manifests_dir, filename), 'w') as f:
            json.dump(build(chosen), f, indent=2)
            f.write('\n')
        print(f'wrote {_rel(os.path.join(manifests_dir, filename))} ({len(chosen)} versions)')


if __name__ == '__main__':
    main()
