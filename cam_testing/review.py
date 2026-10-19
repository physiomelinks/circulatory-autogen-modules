"""
Applying a version's confirmed review proposals.

A version prepared for review carries:
  - in <module_type>_<version>_tests.yaml: reference_proposals: {param: {reference: "<bibkey>; note",
    sourced: yes|no|definitional}} and review: {summary, questions, proposed_fixes, findings};
  - in <module_type>_<version>_verification_config.json: validation.<parameterisation>.<kind>.status: proposed
    (validation data run for review, not confirmed);
and <module_type>_<version>_references_proposed.bib with the proposed citations.

    python -m cam_testing.review apply --component heart/vp                  # everything proposed
    python -m cam_testing.review apply --module heart                        # every version of heart
    python -m cam_testing.review apply --component heart/vp --only R_ao C_ao # just these parameters
    python -m cam_testing.review apply --component heart/vp --no-validation  # references only

Applying writes the references into every parameterisation's parameters file, moves the cited entries from
the proposed .bib into the version's _references.bib, turns proposed validation into active, and
removes what was applied from the spec. It does not set reviewed: true; do that when the version is
done.
"""
import argparse
import csv
import os
import re

from cam_testing import bib
from cam_testing.library import INSTANCE_COLUMNS, all_versions, merge_spec, read_spec_files, version_by_key, write_spec


def _entry_blocks(path):
    """{key: raw text} of the entries in a .bib file."""
    if not os.path.isfile(path):
        return {}
    text = open(path).read()
    return {m.group(2): m.group(0) for m in re.finditer(r'@(\w+)\s*\{\s*([^,\s]+)\s*,.*?\n\}', text, re.S)}


def apply(version, only=None, validation=True):
    spec = merge_spec(*read_spec_files(version.dir, version.stem))

    props = spec.get('reference_proposals') or {}
    applied = {}
    for var in list(props):
        if only and var not in only:
            continue
        applied[var] = props.pop(var)
    if not props:
        spec.pop('reference_proposals', None)
    if validation:
        for inst, kinds in (spec.get('validation') or {}).items():
            for kind, v in (kinds or {}).items():
                if isinstance(v, dict) and v.get('status') == 'proposed':
                    v['status'] = 'active'

    used_keys = set()
    for prop in applied.values():
        key = bib.reference_key(prop['reference'])
        if key:
            used_keys.add(key)
        used_keys.update(prop.get('also_cites') or [])
    for inst in version.parameterisations():
        rows = list(csv.DictReader(open(inst.parameters_path)))
        for r in rows:
            prop = applied.get(r['variable_name'])
            if prop is None:
                continue
            r['data_reference'] = prop['reference']
            r['sourced'] = 'no' if str(prop.get('sourced', 'no')).strip().lower() in ('no', 'false', '0') else 'yes'
            if 'value' in prop:
                r['value'] = str(prop['value'])
        with open(inst.parameters_path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=INSTANCE_COLUMNS, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)

    # move cited entries from the proposed bib to the version's bib
    proposed_path, bib_path = bib.proposed_bib_path(version), bib.bib_path(version)
    proposed, confirmed = _entry_blocks(proposed_path), _entry_blocks(bib_path)
    moved = [k for k in used_keys if k in proposed and k not in confirmed]
    if moved:
        with open(bib_path, 'a') as f:
            for k in moved:
                f.write('\n' + proposed[k] + '\n')
        remaining = {k: v for k, v in proposed.items() if k not in moved}
        with open(proposed_path, 'w') as f:
            f.write('\n\n'.join(remaining.values()) + ('\n' if remaining else ''))

    write_spec(version.dir, version.stem, spec)
    return len(applied), moved


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='cmd', required=True)
    ap = sub.add_parser('apply')
    ap.add_argument('--module', action='append', default=[], help='every version of a module_type (or category)')
    ap.add_argument('--component', action='append', default=[], help='one version, <module_type>/<version>')
    ap.add_argument('--only', nargs='*', help='only these parameters')
    ap.add_argument('--no-validation', action='store_true', help='leave proposed validation as proposed')
    args = parser.parse_args(argv)
    if not args.module and not args.component:
        parser.error('give --module or --component')
    versions = [version_by_key(k) for k in args.component] + (all_versions(args.module) if args.module else [])
    for version in versions:
        n, moved = apply(version, set(args.only) if args.only else None, not args.no_validation)
        print(f'{version.key}: applied {n} reference proposals; moved {len(moved)} bib entries: {", ".join(moved)}')


if __name__ == '__main__':
    main()
