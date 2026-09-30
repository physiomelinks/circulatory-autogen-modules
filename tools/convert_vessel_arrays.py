"""
Converts vessel arrays from CSV (name, BC_type, vessel_type, inp_vessels, out_vessels) to JSON
records in PhLynx's key names (<prefix>_vessel_array.json; see cam_testing/vessel_array.py) and
removes the CSV.

    python tools/convert_vessel_arrays.py                    # every modules/**/_vessel_array.csv
    python tools/convert_vessel_arrays.py path/to/x_vessel_array.csv ...
    python tools/convert_vessel_arrays.py --check            # exit 1 if any CSV vessel array is left
"""
import argparse
import glob
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from cam_testing import vessel_array  # noqa: E402


def csv_vessel_arrays(skip=()):
    paths = glob.glob(os.path.join(REPO, 'modules', '**', '*_vessel_array.csv'), recursive=True)
    return sorted(p for p in paths if not any(s in p for s in skip))


def convert(path):
    out = path[:-len('.csv')] + '.json'
    vessel_array.write_records(out, vessel_array.read_records(path))
    if vessel_array.read_records(out) != vessel_array.read_records(path):
        raise RuntimeError(f'{path}: JSON does not round-trip')
    os.remove(path)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('paths', nargs='*', help='CSV vessel arrays (default: every one under modules/)')
    ap.add_argument('--check', action='store_true', help='only list CSV vessel arrays still present')
    ap.add_argument('--skip', action='append', default=[], help='skip paths containing this text (repeatable)')
    args = ap.parse_args(argv)
    paths = args.paths or csv_vessel_arrays(args.skip)
    if args.check:
        print('\n'.join(os.path.relpath(p, REPO) for p in paths) or 'no CSV vessel arrays left')
        return 1 if paths else 0
    for p in paths:
        convert(p)
    print(f'converted {len(paths)} vessel arrays to JSON')
    return 0


if __name__ == '__main__':
    sys.exit(main())
