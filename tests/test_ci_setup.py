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


def test_no_workflow_uses_ubuntu_latest():
    '''The PyPI libcellml wheel needs libxml2.so.2, which ubuntu-latest (26.04) does not have.'''
    found = []
    for folder, _, names in os.walk(os.path.join(ROOT, '.github')):
        found += [os.path.join(folder, n) for n in names if n.endswith(('.yml', '.yaml'))]
    assert found
    for path in found:
        with open(path) as f:
            for number, line in enumerate(f, 1):
                assert 'ubuntu-latest' not in line.split('#')[0], f'{path}:{number}'


def test_only_the_setup_cache_job_saves_the_setup_caches():
    '''module-tests.yml: every other job that uses the setup action needs setup-cache and only restores.'''
    import yaml
    with open(os.path.join(ROOT, '.github', 'workflows', 'module-tests.yml')) as f:
        jobs = yaml.safe_load(f)['jobs']
    users = {}
    for name, job in jobs.items():
        for step in job.get('steps', []):
            if step.get('uses') == './.github/actions/setup':
                users[name] = step.get('with', {})
    assert 'setup-cache' in users and users['setup-cache'].get('save', 'true') != 'false'
    for name, with_ in users.items():
        if name != 'setup-cache':
            needs = jobs[name]['needs']
            assert 'setup-cache' in ([needs] if isinstance(needs, str) else needs), name
            assert with_.get('save') == 'false', name
            assert 'setup-cache.outputs.debs-key' in with_['debs-key'], name
            assert 'setup-cache.outputs.venv-key' in with_['venv-key'], name
