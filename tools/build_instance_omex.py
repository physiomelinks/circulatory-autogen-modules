"""
Builds a COMBINE archive (.omex) for every instance of every version, from the library's files
(see cam_testing/omex.py): load one in CUFLynx to investigate that parameter set.

    python tools/build_instance_omex.py                       # every instance
    python tools/build_instance_omex.py --module Lotka_Volterra --module cell/neurons
    python tools/build_instance_omex.py --version heart/Argus2026_v01

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
    ap.add_argument('--module', action='append', default=[], help='module_type or category path (repeatable)')
    ap.add_argument('--version', action='append', default=[], help='<module_type>/<version> (repeatable)')
    args = ap.parse_args(argv)
    versions = ([library.version_by_key(k) for k in args.version] if args.version
                else library.all_versions(args.module or None))
    built, failed, skipped = 0, [], 0
    for v in versions:
        if v.format != 'cellml' and not v.is_supermodule:
            skipped += 1
            continue
        for inst in v.instances():
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
