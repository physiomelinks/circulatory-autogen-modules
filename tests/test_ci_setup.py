'''ci/constraints.txt: the pins of every Python package CI installs (README "CI").'''
import os
import re

from packaging.requirements import Requirement
from packaging.version import Version

try:
    import tomllib
except ModuleNotFoundError:   # Python < 3.11
    import tomli as tomllib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _pins():
    pins = {}
    with open(os.path.join(ROOT, 'ci', 'constraints.txt')) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                name, _, version = line.partition('==')
                assert version, f'not an exact pin: {line}'
                pins[re.sub(r'[-_.]+', '-', name).lower()] = Version(version)
    return pins


def test_every_pin_is_exact_and_libcuflynx_is_left_to_requirements():
    pins = _pins()
    assert 'libcuflynx' not in pins and 'cam-testing' not in pins


def test_the_pins_satisfy_the_package_dependencies():
    with open(os.path.join(ROOT, 'pyproject.toml'), 'rb') as f:
        deps = tomllib.load(f)['project']['dependencies']
    pins = _pins()
    for dep in map(Requirement, deps):
        if dep.marker and not dep.marker.evaluate({'python_version': '3.11'}):
            continue
        name = re.sub(r'[-_.]+', '-', dep.name).lower()
        assert name in pins, f'{dep.name} is not pinned in ci/constraints.txt'
        assert dep.specifier.contains(pins[name]), f'{dep.name} {pins[name]} does not satisfy {dep.specifier}'


def test_the_pins_are_part_of_the_cache_key():
    with open(os.path.join(ROOT, 'ci', 'setup_key.sh')) as f:
        assert 'ci/constraints.txt' in f.read()
