"""
Builds a COMBINE archive (.omex) for every parameterisation of every version, from the library's files
(see cam_testing/omex.py): load one in CUFLynx to investigate that parameter set.

    python tools/build_parameterisation_omex.py                       # every parameterisation
    python tools/build_parameterisation_omex.py --module Lotka_Volterra --module neuron
    python tools/build_parameterisation_omex.py --version heart/Argus2026_v01
    python tools/build_parameterisation_omex.py --module vessels --shard 2/6 --shard-stage omex   # one CI shard

C++ (module_format cpp) versions are skipped: CUFLynx runs CellML models.
"""
import argparse
import os
import sys
import traceback

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from cam_testing import library, omex  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--module', action='append', default=[], help='module_type (with those nested in it) or path under modules/ (repeatable)')
    ap.add_argument('--version', action='append', default=[], help='<module_type>/<version> (repeatable)')
    ap.add_argument('--shard', default=None, metavar='K/N', help='only shard K of N of the selected module_types (cam_testing/shards.py)')
    ap.add_argument('--shard-stage', default='omex', help='whose durations balance --shard (default omex)')
    args = ap.parse_args(argv)
    versions = ([library.version_by_key(k) for k in args.version] if args.version
                else library.all_versions(args.module or None))
    if args.shard:
        from cam_testing.shards import shard_module_types
        wanted = shard_module_types(args.module, args.shard_stage, args.shard)
        versions = [v for v in versions if v.mtype.name in wanted]
    built, failed, skipped = 0, [], 0
    for v in versions:
        if v.format != 'cellml' and not v.is_supermodule:
            skipped += 1
            continue
        for inst in v.parameterisations():
            try:
                path = omex.build(v, inst)
                built += 1
                print(os.path.relpath(path, REPO))
            except Exception as e:  # noqa: BLE001 - keep going; report at the end
                failed.append((f'{v.key}/{inst.name}', f'{type(e).__name__}: {str(e).splitlines()[0][:200] if str(e) else ""}'))
                if os.environ.get('CAM_DEBUG'):
                    traceback.print_exc()
    print(f'built {built} archives; {skipped} C++ versions skipped; {len(failed)} failed')
    for key, why in failed:
        print(f'  FAILED {key}: {why}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
