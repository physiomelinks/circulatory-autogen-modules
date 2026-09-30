"""
Applying a module's confirmed review proposals.

A module prepared for review carries, in <name>_tests.yaml:
  - components[].reference_proposals: {param: {reference: "<bibkey>; note", sourced: yes|no|definitional}}
  - validation.<kind>.status: proposed      (validation data run for review, not confirmed)
  - review: {summary, questions, proposed_fixes, findings}
and <name>_references_proposed.bib with the proposed citations.

    python -m cam_testing.review apply --module heart                 # everything proposed
    python -m cam_testing.review apply --module heart --only R_ao C_ao  # just these parameters
    python -m cam_testing.review apply --module heart --no-validation   # references only

Applying writes the references into <name>_parameters.csv, moves the cited entries from the
proposed .bib into <name>_references.bib, turns proposed validation into active, and removes
what was applied from the spec. It does not set reviewed: true; do that when the module is done.
"""
import argparse
import csv
import os
import re

import yaml

from cam_testing import bib
from cam_testing.library import PARAMETER_COLUMNS, load_module


def _entry_blocks(path):
    """{key: raw text} of the entries in a .bib file."""
    if not os.path.isfile(path):
        return {}
    text = open(path).read()
    return {m.group(2): m.group(0) for m in re.finditer(r'@(\w+)\s*\{\s*([^,\s]+)\s*,.*?\n\}', text, re.S)}


def apply(name, only=None, validation=True):
    module = load_module(name)
    spec_path = os.path.join(module.dir, f'{name}_tests.yaml')
    with open(spec_path) as f:
        spec = yaml.safe_load(f)

    applied = {}      # (vessel_type, BC_type, variable) -> proposal
    for comp in spec.get('components', []):
        props = comp.get('reference_proposals') or {}
        for var in list(props):
            if only and var not in only:
                continue
            applied[(comp['vessel_type'], comp['BC_type'], var)] = props.pop(var)
        if not props:
            comp.pop('reference_proposals', None)
        if validation:
            for kind, v in (comp.get('validation') or {}).items():
                if isinstance(v, dict) and v.get('status') == 'proposed':
                    v['status'] = 'active'

    # parameters file (global constants are matched on variable name alone)
    params_path = os.path.join(module.dir, f'{name}_parameters.csv')
    rows = list(csv.DictReader(open(params_path)))
    used_keys = set()
    for r in rows:
        prop = applied.get((r['vessel_type'], r['BC_type'], r['variable_name']))
        if prop is None and r['vessel_type'] == 'global':
            prop = next((p for (vt, bc, var), p in applied.items() if var == r['variable_name']), None)
        if prop is None:
            continue
        r['data_reference'] = prop['reference']
        r['sourced'] = 'no' if str(prop.get('sourced', 'no')).strip().lower() in ('no', 'false', '0') else 'yes'
        if 'value' in prop:
            r['value'] = str(prop['value'])
        key = bib.reference_key(prop['reference'])
        if key:
            used_keys.add(key)
        used_keys.update(prop.get('also_cites') or [])
    with open(params_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=PARAMETER_COLUMNS, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

    # move cited entries from the proposed bib to the module bib
    proposed_path, bib_path = bib.proposed_bib_path(module), bib.bib_path(module)
    proposed, confirmed = _entry_blocks(proposed_path), _entry_blocks(bib_path)
    moved = [k for k in used_keys if k in proposed and k not in confirmed]
    if moved:
        with open(bib_path, 'a') as f:
            for k in moved:
                f.write('\n' + proposed[k] + '\n')
        remaining = {k: v for k, v in proposed.items() if k not in moved}
        with open(proposed_path, 'w') as f:
            f.write('\n\n'.join(remaining.values()) + ('\n' if remaining else ''))

    with open(spec_path, 'w') as f:
        yaml.safe_dump(spec, f, sort_keys=False, width=120, default_flow_style=None, allow_unicode=True)
    return len(applied), moved


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='cmd', required=True)
    ap = sub.add_parser('apply')
    ap.add_argument('--module', required=True)
    ap.add_argument('--only', nargs='*', help='only these parameters')
    ap.add_argument('--no-validation', action='store_true', help='leave proposed validation as proposed')
    args = parser.parse_args(argv)
    n, moved = apply(args.module, set(args.only) if args.only else None, not args.no_validation)
    print(f'applied {n} reference proposals; moved {len(moved)} bib entries: {", ".join(moved)}')


if __name__ == '__main__':
    main()
