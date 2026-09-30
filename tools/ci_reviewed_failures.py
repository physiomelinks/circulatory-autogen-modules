"""
CI helper: which failed tests in pytest JUnit XML reports belong to reviewed versions?

A category's CI job runs every version in it; failures of versions not reviewed yet are warnings,
failures of reviewed versions fail the job. Prints each failure with its version and exits 1 if a
reviewed version failed.

    python tools/ci_reviewed_failures.py junit-tests.xml [junit-pipeline.xml ...]
"""
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing.library import load_version  # noqa: E402


def failures(paths):
    for path in paths:
        if not os.path.isfile(path):
            continue
        for case in ET.parse(path).getroot().iter('testcase'):
            if case.find('failure') is None and case.find('error') is None:
                continue
            m = re.search(r'\[([^\]]+)\]$', case.get('name', ''))
            if not m:
                yield case.get('name'), None
                continue
            parts = m.group(1).split('/')
            yield case.get('name'), tuple(parts[:2]) if len(parts) >= 2 else None


def main(argv=None):
    paths = (argv if argv is not None else sys.argv[1:])
    blocking = 0
    for name, key in failures(paths):
        reviewed = True
        if key is not None:
            try:
                reviewed = load_version(*key).reviewed
            except (OSError, ValueError):
                reviewed = True
        tag = 'error' if reviewed else 'warning'
        print(f'::{tag}::{name} failed' + ('' if reviewed else ' (version not reviewed yet, so not blocking)'))
        blocking += reviewed
    return 1 if blocking else 0


if __name__ == '__main__':
    sys.exit(main())
