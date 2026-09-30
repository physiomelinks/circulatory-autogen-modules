"""
Rewrites the module configs (modules/<name>/<name>_modules_config.json) in PhLynx's key names:

    libcuflynx            PhLynx
    vessel_type       ->  module_type
    BC_type           ->  module_subtype
    module_file       ->  component_file
    module_type       ->  component_type     (the CellML component)

Everything else (ports, multi_port, variables_and_units, only_one_port, delay_info, ...) is
kept as it is, except that a variables_and_units kind with stray whitespace ("variable ") is
trimmed. libcuflynx reads both formats.

    python tools/convert_configs.py            # convert every module
    python tools/convert_configs.py --check    # exit 1 if any config is not in PhLynx format
"""
import argparse
import glob
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from cam_testing.library import config_format, to_phlynx_entry  # noqa: E402


def convert_entry(entry):
    e = to_phlynx_entry(entry)
    if 'variables_and_units' in e:
        e['variables_and_units'] = [[x.strip() if isinstance(x, str) else x for x in v] for v in e['variables_and_units']]
    return e


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--check', action='store_true', help='only report configs not in PhLynx format')
    args = ap.parse_args(argv)
    paths = sorted(glob.glob(os.path.join(REPO, 'modules', '**', '*_modules_config.json'), recursive=True))
    not_phlynx = []
    for path in paths:
        with open(path) as f:
            entries = json.load(f)
        if any(config_format(e) != 'phlynx' for e in entries):
            not_phlynx.append(os.path.relpath(path, REPO))
            if not args.check:
                with open(path, 'w') as f:
                    json.dump([convert_entry(e) for e in entries], f, indent=2, ensure_ascii=False)
                    f.write('\n')
    if args.check:
        print('\n'.join(not_phlynx) or 'all configs are in PhLynx format')
        return 1 if not_phlynx else 0
    print(f'converted {len(not_phlynx)} of {len(paths)} configs to PhLynx format')
    return 0


if __name__ == '__main__':
    sys.exit(main())
