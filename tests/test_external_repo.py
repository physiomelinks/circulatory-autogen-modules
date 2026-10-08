"""
cam_testing in another repository, end to end: a tiny repo outside this one (tests/external_repo.py:
one module version whose harness uses this library's inlet_flow and arterial_simple, and a system
model whose reference is a ready CellML model) runs the suite in a subprocess, as a consumer would:
the structure checks with ``pytest --pyargs``, the module V&V tests from a one-line test file, the
system-model tests, then the report and the manifests.
"""
import json
import os
import re
import subprocess
import sys

import pytest

import external_repo
from cam_testing import paths

pytest.importorskip('libcuflynx')


def _env():
    # as a consumer's shell: no CAM_* from this session's plugin; the repo's pyproject.toml says it all
    env = {k: v for k, v in os.environ.items() if k not in (paths.ENV_ROOT, paths.ENV_LIBRARIES)}
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    return env


def _run(args, cwd, timeout=1800):
    if args[:2] == ['-m', 'pytest']:
        # cam_testing's suite files live in this checkout: pytest would load conftest.py files above it
        # (none in this repo, but a checkout nested in another, e.g. a git worktree, may have one)
        args = [*args, '--noconftest']
    p = subprocess.run([sys.executable, *args], cwd=cwd, env=_env(), capture_output=True, text=True, timeout=timeout)
    return p, p.stdout + p.stderr


def _counts(out):
    tail = [line for line in out.splitlines() if re.search(r'\d+ (passed|failed|skipped|error)', line) and ' in ' in line]
    return {k: int(n) for n, k in re.findall(r'(\d+) (passed|failed|skipped|xfailed|error|errors|deselected)', tail[-1] if tail else '')}


@pytest.fixture(scope='module')
def repo(tmp_path_factory):
    root = str(tmp_path_factory.mktemp('consumer'))
    external_repo.make(root)
    # the system model's reference: a ready CellML model (as PhLynx exports it), made here from the libraries
    from cam_testing import system as systems
    try:
        paths.configure(repo_root=root, library_dirs=[external_repo.library_modules_dir()], export=False)
        external_repo.write_reference_cellml(systems.load_system(external_repo.SYSTEM),
                                             str(tmp_path_factory.mktemp('reference_work')))
    finally:
        paths.reset()
    return root


def test_paths_resolved_from_the_consumers_pyproject(repo):
    p, out = _run(['-m', 'cam_testing.paths'], repo)
    assert p.returncode == 0, out
    assert f'repo root:          {repo}' in out
    assert external_repo.library_modules_dir() in out


def test_structure_suite(repo):
    p, out = _run(['-m', 'pytest', '--pyargs', 'cam_testing.suite.test_structure', '-q', '-p', 'no:cacheprovider',
                   '-rfE'], repo)
    assert p.returncode == 0, out[-6000:]
    c = _counts(out)
    assert c.get('passed', 0) >= 15 and not c.get('failed'), out[-3000:]
    assert 'ext_outlet_pressure/constant' in _run(
        ['-m', 'pytest', '--pyargs', 'cam_testing.suite.test_structure', '--collect-only', '-q',
         '-p', 'no:cacheprovider'], repo)[1]


def test_structure_suite_catches_a_clash_with_the_library(repo, tmp_path):
    '''Uniqueness is checked over both libraries: a component name taken from the library fails.'''
    clash = str(tmp_path / 'clash')
    external_repo.make(clash)
    cellml = os.path.join(clash, 'modules', 'boundary_conditions', external_repo.MT, 'versions', external_repo.VERSION,
                          f'{external_repo.MT}_{external_repo.VERSION}_modules.cellml')
    with open(cellml) as f:
        text = f.read()
    with open(cellml, 'w') as f:
        f.write(text.replace(external_repo.NEW_COMPONENT, external_repo.OLD_COMPONENT))
    cfg = cellml.replace('_modules.cellml', '_modules_config.json')
    with open(cfg) as f:
        entry = json.load(f)
    entry[0]['component_type'] = external_repo.OLD_COMPONENT
    with open(cfg, 'w') as f:
        json.dump(entry, f)
    p, out = _run(['-m', 'pytest', '--pyargs', 'cam_testing.suite.test_structure', '-q', '-p', 'no:cacheprovider',
                   '-k', 'uniqueness'], clash)
    assert p.returncode == 1 and f'component {external_repo.OLD_COMPONENT} is in both' in out, out[-3000:]


def test_module_suite_from_a_one_line_test_file(repo):
    tests = os.path.join(repo, 'tests')
    os.makedirs(tests, exist_ok=True)
    with open(os.path.join(tests, 'test_cam_modules.py'), 'w') as f:
        f.write('from cam_testing.suite.test_modules import *  # noqa: F401,F403\n')
    p, out = _run(['-m', 'pytest', 'tests/test_cam_modules.py', '-m', 'not slow', '--module', external_repo.MT,
                   '-p', 'no:cacheprovider', '-rfE'], repo)
    assert p.returncode == 0, out[-6000:]
    c = _counts(out)
    # run_test, invariants, BC sweep, timestep, stability, baseline (not applicable: skipped)
    assert c.get('passed', 0) >= 5 and not c.get('failed'), out[-3000:]
    results = os.path.join(repo, 'modules', 'boundary_conditions', external_repo.MT, 'versions', external_repo.VERSION,
                           'results')
    with open(os.path.join(results, 'run_test.json')) as f:
        assert json.load(f)['status'] == 'passed'


def test_system_suite_with_a_cellml_reference(repo):
    '''The system model (vessels from both libraries) generates, and reproduces its ready CellML reference.'''
    p, out = _run(['-m', 'pytest', '--pyargs', 'cam_testing.suite.test_systems', '-p', 'no:cacheprovider', '-rfEs',
                   '-v'], repo)
    assert p.returncode == 0, out[-6000:]
    assert re.search(r'system_equivalence_test\[demo/ext_chain\] PASSED', out), out[-3000:]
    assert re.search(r'system_run_test\[demo/ext_chain\] PASSED', out), out[-3000:]
    with open(os.path.join(repo, 'system_models', 'demo', 'ext_chain', 'results', 'system_equivalence_test.json')) as f:
        result = json.load(f)
    assert result['metrics']['compared'] > 0 and not result['metrics']['not_in_model']


def test_report_and_manifests(repo):
    p, out = _run(['-m', 'cam_testing.report', '--site'], repo)
    assert p.returncode == 0, out[-4000:]
    mt_dir = os.path.join(repo, 'modules', 'boundary_conditions', external_repo.MT)
    assert os.path.isfile(os.path.join(mt_dir, f'{external_repo.MT}.html'))
    assert os.path.isfile(os.path.join(repo, 'site', 'index.html'))
    p, out = _run(['-m', 'cam_testing.manifests'], repo)
    assert p.returncode == 0, out
    with open(os.path.join(repo, 'manifests', 'all.json')) as f:
        paths_listed = [e['path'] for e in json.load(f)['collections']['configs']]
    assert paths_listed == [f'modules/boundary_conditions/{external_repo.MT}/versions/{external_repo.VERSION}/'
                            f'{external_repo.MT}_{external_repo.VERSION}_modules_config.json']
