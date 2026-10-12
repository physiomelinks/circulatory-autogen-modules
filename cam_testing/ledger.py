"""
Change-based test selection: a version's tests in a stage are skipped when nothing that could
change their outcome has changed since they last ran.

What a version's outcome in a stage depends on, hashed into its fingerprint:
  * the version's own files (results/, plots/, *.html and *.omex left out) and the files directly
    in its module_type's directory;
  * the same for every version it depends on, recursively: the neighbours of its test networks
    (harness, validation.<parameterisation>.harness) and a supermodule's submodules;
  * the system models it must reproduce (supermodule.equivalent) or is coupled in (coupled_systems);
  * shared inputs: cam_testing's code, tests/, requirements.txt, pyproject.toml, the Makefile, the
    directory schema and the installed libcuflynx commit; for the pipeline and omex stages also
    the PhLynx checkout's commit and the CUFLynx binary.

The ledger (``$CAM_LEDGER_DIR/ledger.json``, restored from the CI cache) holds, per stage and
version, the fingerprint and outcome of the last run. With ``--changed-only`` a version is skipped
when its fingerprint matches and its last outcome was a pass, or a failure of a version not
reviewed yet (those never block CI). A reviewed version that failed always runs again. The weekly
CI run tests everything and refreshes the ledger.

    python -m cam_testing.ledger fragment <stage> junit.xml fragment.json   # after a test job
    python -m cam_testing.ledger merge <ledger_dir> fragment*.json          # in the site job
    python -m cam_testing.ledger status <stage> --module cell                # what would run
"""
import argparse
import functools
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

STAGES = ('vv', 'pipeline', 'omex')
GENERATED_DIRS = {'results', 'plots', '__pycache__', 'node_modules', 'source_figures_cache'}
GENERATED_SUFFIXES = ('.html', '.omex', '.pyc')
LEDGER_NAME = 'ledger.json'
UNCHANGED_REASON = 'unchanged since'   # the start of the skip reason of a version skipped as unchanged


def ledger_dir():
    return os.environ.get('CAM_LEDGER_DIR') or os.path.join(os.path.expanduser('~'), '.cam-ledger')


# ---- hashing -----------------------------------------------------------------------------------

def _hash_file(h, path, rel):
    h.update(rel.encode() + b'\0')
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    h.update(b'\0')


def _tree_files(top):
    for d, dirs, files in os.walk(top):
        dirs[:] = sorted(x for x in dirs if x not in GENERATED_DIRS)
        for f in sorted(files):
            if not f.endswith(GENERATED_SUFFIXES):
                yield os.path.join(d, f)


@functools.lru_cache(maxsize=None)
def _dir_hash(top, recursive=True):
    '''Hash of a directory's source files (generated files left out).'''
    h = hashlib.sha256()
    if os.path.isdir(top):
        files = _tree_files(top) if recursive else sorted(
            os.path.join(top, f) for f in os.listdir(top)
            if os.path.isfile(os.path.join(top, f)) and not f.endswith(GENERATED_SUFFIXES))
        for p in files:
            _hash_file(h, p, os.path.relpath(p, top))
    return h.hexdigest()


def _libcuflynx_id():
    from importlib import metadata
    try:
        dist = metadata.distribution('libcuflynx')
    except metadata.PackageNotFoundError:
        return 'none'
    try:
        info = json.loads(dist.read_text('direct_url.json') or '{}')
    except ValueError:
        info = {}
    commit = (info.get('vcs_info') or {}).get('commit_id')
    if commit:
        return commit
    # an editable or local install: its source tree
    url = info.get('url', '')
    if url.startswith('file://') and os.path.isdir(url[7:]):
        return _dir_hash(os.path.join(url[7:], 'src'))
    return dist.version


def _git_head(path):
    try:
        return subprocess.run(['git', '-C', path, 'rev-parse', 'HEAD'], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return _dir_hash(path)


@functools.lru_cache(maxsize=None)
def shared_inputs(stage):
    '''Hash of what every version's outcome in a stage depends on.'''
    from cam_testing import paths
    repo = paths.roots().repo_root
    pkg = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    h.update(f'stage {stage}\0'.encode())
    h.update(_dir_hash(pkg).encode())
    for rel in ('tests', ):
        h.update(_dir_hash(os.path.join(repo, rel)).encode())
    for rel in ('requirements.txt', 'pyproject.toml', 'Makefile', os.path.join('modules', 'directory_schema.json')):
        p = os.path.join(repo, rel)
        if os.path.isfile(p):
            _hash_file(h, p, rel)
    h.update(_libcuflynx_id().encode())
    if stage == 'pipeline' and os.environ.get('PHLYNX_DIR'):
        h.update(_git_head(os.environ['PHLYNX_DIR']).encode())
    if stage in ('pipeline', 'omex') and os.environ.get('CUFLYNX_BIN') and os.path.isfile(os.environ['CUFLYNX_BIN']):
        hb = hashlib.sha256()
        _hash_file(hb, os.environ['CUFLYNX_BIN'], 'cuflynx')
        h.update(hb.hexdigest().encode())
    return h.hexdigest()


# ---- dependencies ------------------------------------------------------------------------------

def _network_versions(version):
    from cam_testing import harness, module_array
    from cam_testing.library import load_version
    networks = [version.spec.get('harness') or {}]
    networks += [(v or {}).get('harness') or {} for v in (version.spec.get('validation') or {}).values()
                 if isinstance(v, dict)]
    out = []
    for net in networks:
        for row in module_array.harness_rows(net) or []:
            if str(row[0]) == harness.VESSEL:
                continue
            try:
                out.append(load_version(str(row[2]), str(row[1])))
            except (OSError, ValueError, KeyError):
                pass
    return out


def _submodule_versions(version):
    from cam_testing.library import load_version
    out = []
    for s in version.submodules if version.is_supermodule else []:
        try:
            out.append(load_version(s['module_type'], s['module_subtype']))
        except (OSError, ValueError, KeyError):
            pass
    return out


def dependencies(version):
    '''Every version ``version`` depends on, recursively (itself excluded).'''
    seen, todo = {}, [version]
    while todo:
        v = todo.pop()
        for d in _network_versions(v) + _submodule_versions(v):
            if d.key != version.key and d.key not in seen:
                seen[d.key] = d
                todo.append(d)
    return [seen[k] for k in sorted(seen)]


def _system_model_dirs(version):
    from cam_testing import system
    names = [e.get('model') for e in ((version.spec.get('supermodule') or {}).get('equivalent') or [])
             if isinstance(e, dict)]
    names += list(version.spec.get('coupled_systems') or [])
    out = []
    for n in names:
        try:
            out.append(system.load_system(n).dir)
        except (KeyError, OSError, AttributeError):
            pass
    return sorted(set(out))


def _version_hash(version):
    return _dir_hash(version.dir) + _dir_hash(version.mtype.dir, recursive=False)


def fingerprint(version, stage):
    h = hashlib.sha256()
    h.update(shared_inputs(stage).encode())
    h.update(f'{version.key}\0{_version_hash(version)}\0'.encode())
    for d in dependencies(version):
        h.update(f'{d.key}\0{_version_hash(d)}\0'.encode())
    for d in _system_model_dirs(version):
        h.update(f'{os.path.basename(d)}\0{_dir_hash(d)}\0'.encode())
    return h.hexdigest()[:32]


# ---- the ledger --------------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _load(path):
    if not os.path.isfile(path):
        return {}
    with open(path) as f:
        return json.load(f)


def load():
    return _load(os.path.join(ledger_dir(), LEDGER_NAME))


def unchanged(version, stage):
    '''(True, reason) when the version's tests in ``stage`` can be skipped.'''
    entry = (load().get(stage) or {}).get(version.key)
    if not entry or entry.get('fingerprint') != fingerprint(version, stage):
        return False, None
    when = entry.get('run', 'an earlier run')
    if entry.get('outcome') == 'pass':
        return True, f'unchanged since {when}, where it passed'
    if not version.reviewed:
        return True, f'unchanged since {when}, where it failed (not reviewed yet, so not blocking)'
    return False, None


def _version_of_case(case):
    m = re.search(r'\[([^\]]+)\]$', case.get('name', ''))
    if not m:
        return None
    parts = m.group(1).split('/')
    return tuple(parts[:2]) if len(parts) >= 2 else None


def fragment(stage, junit_paths, run_label=None):
    '''The ledger entries of the versions a test job ran: fingerprint and outcome.'''
    from cam_testing.library import load_version
    outcome = {}
    for p in junit_paths:
        if not os.path.isfile(p):
            continue
        for case in ET.parse(p).getroot().iter('testcase'):
            key = _version_of_case(case)
            if key is None:
                continue
            skipped = case.find('skipped')
            if skipped is not None and UNCHANGED_REASON in (skipped.get('message') or ''):
                continue   # skipped as unchanged: its ledger entry stays as it was
            failed = case.find('failure') is not None or case.find('error') is not None
            outcome[key] = 'fail' if failed or outcome.get(key) == 'fail' else 'pass'
    run = run_label or time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())
    entries = {}
    for key, res in outcome.items():
        try:
            v = load_version(*key)
        except (OSError, ValueError, KeyError):
            continue
        entries[v.key] = {'fingerprint': fingerprint(v, stage), 'outcome': res, 'run': run}
    return {stage: entries}


def merge(directory, fragments):
    path = os.path.join(directory, LEDGER_NAME)
    data = dict(_load(path)) if os.path.isfile(path) else {}
    _load.cache_clear()
    n = 0
    for frag in fragments:
        for stage, entries in frag.items():
            data.setdefault(stage, {}).update(entries)
            n += len(entries)
    os.makedirs(directory, exist_ok=True)
    with open(path, 'w') as f:
        json.dump({s: dict(sorted(e.items())) for s, e in sorted(data.items())}, f, indent=1)
        f.write('\n')
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('fragment')
    p.add_argument('stage', choices=STAGES)
    p.add_argument('junit')
    p.add_argument('out')
    p.add_argument('--run', default=None, help='label of the run (default: now)')
    p = sub.add_parser('merge')
    p.add_argument('directory')
    p.add_argument('fragments', nargs='*')
    p = sub.add_parser('status')
    p.add_argument('stage', choices=STAGES)
    p.add_argument('--module', action='append', default=[])
    args = ap.parse_args(argv)
    if args.cmd == 'fragment':
        frag = fragment(args.stage, [args.junit], args.run)
        with open(args.out, 'w') as f:
            json.dump(frag, f, indent=1)
        print(f'{len(frag[args.stage])} versions recorded for {args.stage} in {args.out}')
    elif args.cmd == 'merge':
        frags = []
        for p in args.fragments:
            with open(p) as f:
                frags.append(json.load(f))
        print(f'{merge(args.directory, frags)} entries merged into {os.path.join(args.directory, LEDGER_NAME)}')
    elif args.cmd == 'status':
        from cam_testing.library import all_versions
        run = skip = 0
        for v in all_versions(args.module or None):
            ok, why = unchanged(v, args.stage)
            skip += ok
            run += not ok
            if not ok:
                print(f'run   {v.key}')
        print(f'{run} to run, {skip} unchanged ({args.stage}, ledger {ledger_dir()})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
