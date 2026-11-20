"""
Releases of the module library: the rename map libcuflynx reads, and the release asset it fetches.

    python -m cam_testing.release renames             # writes modules/renames.json
    python -m cam_testing.release asset v2026.11.0    # writes dist/cam-modules-v2026.11.0.tar.gz and .sha256

modules/renames.json maps every old (module_type, module_subtype) -- circulatory_autogen's built-in
modules and this library's renamed versions -- to the version that replaces it
(library.legacy_renames()). It is committed, and tests/test_release.py checks it is in step.

The asset is a tar.gz of the committed files libcuflynx needs to generate models: modules/ without
the test specs, reports, figures, risk analyses and Markdown (SLIM_EXCLUDED), manifests/ and LICENSE,
under one top directory cam-modules-<tag>/, with cam_release.json: the tag, the commit, the
libcuflynx it needs ([tool.cam_testing.release] libcuflynx in pyproject.toml), and every file's
sha256. Members are sorted, with a fixed mtime, owner and mode, and gzip writes no name or time, so
building the same commit again gives the same bytes and the same sha256.
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tarfile

from cam_testing import paths

RENAMES_FILE = 'renames.json'
RELEASE_FILE = 'cam_release.json'
TAG_RE = re.compile(r'^v(\d{4})\.(\d{2})\.(\d+)$')
ASSET_CONTENT = ('modules', 'manifests', 'LICENSE')
# path parts and file name endings left out of the asset (reports, figures, analyses, test specs)
SLIM_EXCLUDED_DIRS = {'risk', 'figures', 'source_figures', 'plots', 'results'}
SLIM_EXCLUDED_SUFFIXES = ('_references_proposed.bib', '.html', '.png', '.svg', '.jpg', '.jpeg', '.gif', '.pdf', '.ttl', '.npz', '.md',
                          '_tests.yaml', '_verification_config.json')
FIXED_MTIME = 0


def repo_root():
    return paths.roots().repo_root


JUNCTION_RENAMES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'junction_renames.json')
JUNCTION_UPGRADE = 'cuflynx-upgrade junctions'


def junction_renames():
    '''The junction versions this library removed: {(module_type, version): entry}. Each entry is a
    renames.json row for libcuflynx's `cuflynx-upgrade junctions` ("by"), which also renames the
    variables ("variables", "neighbour_variables"), splits a two-inflow source ("split") or joins the
    neighbours of a massless node ("join", "to" null).'''
    with open(JUNCTION_RENAMES) as f:
        return {tuple(e['from']): e for e in json.load(f)}


def renames_document():
    '''modules/renames.json's content: every old pair and the version that replaces it, sorted. A pair
    renamed to a junction version gets that junction's entry (for cuflynx-upgrade junctions), so
    `cuflynx-upgrade module-names` never renames to a version the library no longer has.'''
    from cam_testing.library import legacy_renames
    junctions = junction_renames()
    rows = {}
    for old, new in legacy_renames().items():
        rows[old] = dict(junctions[new], **{'from': list(old)}) if new in junctions else \
            {'from': list(old), 'to': list(new)}
    for old, entry in junctions.items():
        rows.setdefault(old, entry)
    return {'schema_version': 1,
            'description': ('(module_type, module_subtype) pairs this library no longer has -> the version '
                            'that replaces each: circulatory_autogen\'s built-in modules and this library\'s '
                            'renamed versions. Entries with "by": "cuflynx-upgrade junctions" are the junction '
                            'module types, which that command converts (it renames variables too). Written by '
                            'python -m cam_testing.release renames.'),
            'renames': [rows[old] for old in sorted(rows)]}


def renames_text():
    return json.dumps(renames_document(), indent=1) + '\n'


def write_renames(root=None):
    path = os.path.join(root or repo_root(), 'modules', RENAMES_FILE)
    with open(path, 'w') as f:
        f.write(renames_text())
    return path


def check_tag(tag):
    '''A calendar tag vYYYY.MM.N (N from 0 within the month).'''
    m = TAG_RE.match(tag or '')
    if not m or not 1 <= int(m[2]) <= 12:
        raise ValueError(f'release tag {tag!r} is not vYYYY.MM.N (e.g. v2026.11.0)')
    return tag


def in_asset(path):
    '''Whether a committed file (a path relative to the repo root) goes in the release asset.'''
    parts = path.split('/')
    if parts[0] not in ASSET_CONTENT:
        return False
    if SLIM_EXCLUDED_DIRS.intersection(parts[:-1]):
        return False
    return not parts[-1].endswith(SLIM_EXCLUDED_SUFFIXES)


def _git(root, *args, **kw):
    return subprocess.run(['git', '-C', root, *args], check=True, capture_output=True, **kw).stdout


def committed_files(root, commit):
    '''{path: bytes} of the asset's files as committed at ``commit``.'''
    listing = _git(root, 'ls-tree', '-r', '-z', commit, '--', *ASSET_CONTENT).split(b'\0')
    blobs = []
    for entry in listing:
        if not entry:
            continue
        meta, path = entry.split(b'\t', 1)
        mode, kind, sha = meta.split()
        path = path.decode()
        if kind == b'blob' and mode != b'120000' and in_asset(path):
            blobs.append((path, sha.decode()))
    proc = subprocess.run(['git', '-C', root, 'cat-file', '--batch'], check=True, capture_output=True,
                          input=''.join(f'{sha}\n' for _, sha in blobs).encode())
    out, buf = {}, io.BytesIO(proc.stdout)
    for path, sha in blobs:
        header = buf.readline().split()
        size = int(header[2])
        out[path] = buf.read(size)
        buf.read(1)
    return out


def libcuflynx_requirement(root, commit='HEAD'):
    '''[tool.cam_testing.release] libcuflynx in pyproject.toml as committed at ``commit``.'''
    try:
        import tomllib
    except ModuleNotFoundError:   # Python 3.10
        import tomli as tomllib
    doc = tomllib.loads(_git(root, 'show', f'{commit}:pyproject.toml', text=True))
    req = doc.get('tool', {}).get('cam_testing', {}).get('release', {}).get('libcuflynx')
    if not req:
        raise ValueError('pyproject.toml has no [tool.cam_testing.release] libcuflynx requirement')
    return req


def _member(name, data, mode=0o644):
    info = tarfile.TarInfo(name)
    info.size, info.mtime, info.mode = len(data), FIXED_MTIME, mode
    info.uid = info.gid = 0
    info.uname = info.gname = ''
    return info


def build_asset(tag, out_dir, root=None, commit='HEAD'):
    '''Writes <out_dir>/cam-modules-<tag>.tar.gz and its .sha256 ("<hash>  <name>"); returns
    (asset path, sha256, cam_release.json content).'''
    root = root or repo_root()
    check_tag(tag)
    commit = _git(root, 'rev-parse', f'{commit}^{{commit}}', text=True).strip()
    files = committed_files(root, commit)
    if f'modules/{RENAMES_FILE}' not in files:
        raise ValueError(f'modules/{RENAMES_FILE} is not committed at {commit}')
    release = {'schema_version': 1, 'tag': tag, 'commit': commit,
               'repository': 'physiomelinks/circulatory-autogen-modules',
               'libcuflynx': libcuflynx_requirement(root, commit), 'files': len(files),
               'sha256': {p: hashlib.sha256(files[p]).hexdigest() for p in sorted(files)}}
    files[RELEASE_FILE] = (json.dumps(release, indent=1) + '\n').encode()
    top = f'cam-modules-{tag}'
    os.makedirs(out_dir, exist_ok=True)
    name = f'{top}.tar.gz'
    path = os.path.join(out_dir, name)
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode='w', format=tarfile.PAX_FORMAT) as tar:
        for p in sorted(files):
            tar.addfile(_member(f'{top}/{p}', files[p]), io.BytesIO(files[p]))
    with open(path, 'wb') as f, gzip.GzipFile(filename='', mode='wb', fileobj=f, mtime=FIXED_MTIME, compresslevel=9) as gz:
        gz.write(raw.getvalue())
    with open(path, 'rb') as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    with open(path + '.sha256', 'w') as f:
        f.write(f'{digest}  {name}\n')
    return path, digest, release


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m cam_testing.release', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('renames', help='write modules/renames.json')
    a = sub.add_parser('asset', help='build the release asset of a commit')
    a.add_argument('tag')
    a.add_argument('--commit', default='HEAD')
    a.add_argument('--out-dir', default='dist')
    args = ap.parse_args(argv)
    if args.cmd == 'renames':
        print(f'wrote {write_renames()}')
        return 0
    path, digest, release = build_asset(args.tag, args.out_dir, commit=args.commit)
    print(f'wrote {path} ({os.path.getsize(path) / 1e6:.2f} MB, {release["files"]} files)\nsha256 {digest}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
