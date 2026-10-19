"""
CI shards: each category's tests split into jobs of about ten minutes or less.

A category's tests run in three stages, each its own CI job: ``vv`` (tests/test_modules.py),
``pipeline`` (PhLynx -> CUFLynx, tests/test_phlynx.py) and ``omex`` (CUFLynx archives,
tests/test_instance_omex.py). A stage with more work is split into shards (``ci/shards.json``).
A shard holds whole module_types: the pipeline and archive tests run a module_type's versions in
one batch, so splitting a module_type would repeat its batch in every shard. Module_types are
spread over the shards longest first, by their measured seconds in ``ci/test_durations.json``
(per stage; a module_type with none counts DEFAULT_SECONDS per version).

    python -m cam_testing.shards matrix                         # the CI job matrix (JSON)
    python -m cam_testing.shards list vessels pipeline 3/8      # the module_types of one shard
    python -m cam_testing.shards durations pipeline junit*.xml  # refresh ci/test_durations.json
    python -m cam_testing.shards plan                           # every shard's expected seconds

pytest and tools/build_instance_omex.py take ``--shard K/N --shard-stage STAGE`` (with --module
for the category) and keep the versions of shard K's module_types only.
"""
import argparse
import glob
import json
import os
import re
import sys
import xml.etree.ElementTree as ET

STAGES = ('vv', 'pipeline', 'omex')
DEFAULT_SECONDS = 20
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHARDS_FILE = os.path.join(REPO, 'ci', 'shards.json')
DURATIONS_FILE = os.path.join(REPO, 'ci', 'test_durations.json')


def _load(path):
    if not os.path.isfile(path):
        return {}
    with open(path) as f:
        return json.load(f)


def categories():
    '''Top-level directories of modules/ that hold module_types (as the CI's discover job).'''
    return sorted({p.split(os.sep)[1] for p in glob.glob(os.path.join('modules', '*', '**', 'versions'), recursive=True)})


def shard_count(category, stage):
    return int(_load(SHARDS_FILE).get(category, {}).get(stage, 1))


def parse_shard(text):
    k, n = (int(x) for x in str(text).split('/'))
    if not 1 <= k <= n:
        raise ValueError(f'shard {text}: need 1 <= K <= N')
    return k, n


def weights(module_types, stage):
    '''Expected seconds of each module_type's tests in a stage.'''
    from cam_testing.library import load_module_type
    measured = _load(DURATIONS_FILE).get(stage, {})
    out = {}
    for name in module_types:
        if name in measured:
            out[name] = float(measured[name])
        else:
            out[name] = DEFAULT_SECONDS * max(1, len(load_module_type(name).versions()))
    return out


def assign(module_types, stage, n):
    '''The module_types split into n shards: longest first, each to the lightest shard so far.'''
    w = weights(module_types, stage)
    bins = [[0.0, []] for _ in range(n)]
    for name in sorted(module_types, key=lambda m: (-w[m], m)):
        b = min(bins, key=lambda x: x[0])
        b[0] += w[name]
        b[1].append(name)
    return [sorted(b[1]) for b in bins], [b[0] for b in bins]


def shard_module_types(selectors, stage, shard):
    '''The module_types of shard "K/N" among those the --module selectors pick.'''
    from cam_testing.library import select_module_types
    k, n = parse_shard(shard)
    return set(assign(select_module_types(selectors), stage, n)[0][k - 1])


def matrix():
    out = []
    for c in categories():
        for stage in STAGES:
            n = shard_count(c, stage)
            for k in range(1, n + 1):
                label = f'{c} · {stage}' + (f' {k}/{n}' if n > 1 else '')
                out.append({'category': c, 'stage': stage, 'shard': f'{k}/{n}', 'label': label,
                            'id': f'{c}-{stage}-{k}of{n}'})
    return out


def _module_type_of_case(case):
    m = re.search(r'\[([^\]]+)\]$', case.get('name', ''))
    return m.group(1).split('/')[0] if m else None


def durations(stage, paths):
    '''Adds the per-module_type seconds of JUnit reports to ci/test_durations.json.'''
    data = _load(DURATIONS_FILE)
    seen = {}
    for p in paths:
        for case in ET.parse(p).getroot().iter('testcase'):
            mt = _module_type_of_case(case)
            if mt:
                seen[mt] = seen.get(mt, 0.0) + float(case.get('time') or 0)
    data.setdefault(stage, {}).update({k: round(v, 1) for k, v in seen.items()})
    data[stage] = dict(sorted(data[stage].items()))
    os.makedirs(os.path.dirname(DURATIONS_FILE), exist_ok=True)
    with open(DURATIONS_FILE, 'w') as f:
        json.dump(data, f, indent=2)
        f.write('\n')
    return seen


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('matrix')
    p = sub.add_parser('list')
    p.add_argument('category')
    p.add_argument('stage', choices=STAGES)
    p.add_argument('shard')
    p = sub.add_parser('durations')
    p.add_argument('stage', choices=STAGES)
    p.add_argument('junit', nargs='+')
    sub.add_parser('plan')
    args = ap.parse_args(argv)
    if args.cmd == 'matrix':
        print(json.dumps(matrix()))
    elif args.cmd == 'list':
        print('\n'.join(sorted(shard_module_types([args.category], args.stage, args.shard))))
    elif args.cmd == 'durations':
        seen = durations(args.stage, args.junit)
        print(f'{len(seen)} module_types updated for {args.stage} in {os.path.relpath(DURATIONS_FILE)}')
    elif args.cmd == 'plan':
        from cam_testing.library import select_module_types
        for c in categories():
            for stage in STAGES:
                n = shard_count(c, stage)
                _, secs = assign(select_module_types([c]), stage, n)
                print(f'{c:20s} {stage:8s} ' + '  '.join(f'{s / 60:5.1f}' for s in secs) + '  (min per shard)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
