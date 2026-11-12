"""
CI helper: which failed tests in pytest JUnit XML reports belong to reviewed versions?

A category's CI job runs every version in it; failures of versions not reviewed yet are warnings,
failures of reviewed versions fail the job. A system-model test (system_*_test[<category>/<model>]) is
gated the same way by the model's <model>_system.yaml reviewed key. Prints each failure with its version
or system model and exits 1 if a reviewed one failed.

    python tools/ci_reviewed_failures.py junit-tests.xml [junit-pipeline.xml ...]
"""
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing.library import load_version  # noqa: E402
from cam_testing.system import load_system  # noqa: E402


def is_reviewed(name, key):
    '''Whether the version or system model a failed test belongs to is reviewed (True when unknown).'''
    try:
        if name.startswith('system_'):
            return bool(load_system('/'.join(key)).spec.get('reviewed'))
        return load_version(*key).reviewed
    except (OSError, ValueError, KeyError):
        return True


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
        reviewed = True if key is None else is_reviewed(name, key)
        tag = 'error' if reviewed else 'warning'
        what = 'system model' if name.startswith('system_') else 'version'
        print(f'::{tag}::{name} failed' + ('' if reviewed else f' ({what} not reviewed yet, so not blocking)'))
        blocking += reviewed
    return 1 if blocking else 0


if __name__ == '__main__':
    sys.exit(main())
