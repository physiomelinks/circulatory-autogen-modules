"""
The pytest plugin for cam_testing's test suite (cam_testing/suite/), registered through the
``pytest11`` entry point, so installing cam_testing is enough for any repo to run the suite
against its own modules/:

    pytest --pyargs cam_testing.suite                          # everything
    pytest --pyargs cam_testing.suite.test_modules --module heart
    # or a one-line test file: from cam_testing.suite.test_modules import *

It provides:
  - the repo under test (cam_testing.paths): --cam-root, --cam-library (else CAM_REPO_ROOT,
    CAM_MODULE_LIBRARY_DIRS, [tool.cam_testing] in pyproject.toml, or the working directory);
  - the selection options --module, --component, --include-unreviewed, --quick-unreviewed;
  - the parametrisation over versions, instances and supermodules (component_key, instance_key,
    supermodule_key, equivalence_key) and the component_model / instance_model fixtures;
  - the standard tests' names (run_test, verification_test_BC, ...) as test functions, and the
    markers the suite uses.
"""
import pytest

# The standard tests keep the names run_test, verification_test_BC, ... rather than test_*.
PYTHON_FUNCTIONS = ['run_test', '*_test_*', 'stability_test', 'system_*_test', 'supermodule_*_test',
                    'phlynx_*_test', 'cuflynx_*_test']
MARKERS = [
    'slow: long-running (calibration)',
    'module_vv: per-version verification and per-instance validation tests',
    'system_model: system-model tests (system_models/) and supermodule reproduces tests',
    'phlynx_pipeline: PhLynx -> .omex -> CUFLynx per-version tests (need node, PhLynx and CUFLynx)',
]
QUICK_TESTS = ('run_test', 'supermodule_structure_test', 'phlynx_export_test', 'cuflynx_simulate_test',
               'phlynx_equivalence_test', 'cuflynx_instance_omex_test')


def pytest_addoption(parser):
    group = parser.getgroup('cam_testing', 'cam_testing module library tests')
    group.addoption('--cam-root', default=None,
                    help='the repo under test (holding modules/, system_models/, manifests/); default '
                         '$CAM_REPO_ROOT, [tool.cam_testing] root in pyproject.toml, or the working directory')
    group.addoption('--cam-library', action='append', default=None,
                    help='an extra module library (a modules/ directory, or a repo holding one) its modules may '
                         'use, e.g. a circulatory-autogen-modules checkout; repeatable. Default '
                         '$CAM_MODULE_LIBRARY_DIRS or [tool.cam_testing] module_library_dirs')
    group.addoption('--module', action='append', default=[],
                    help='only test this module_type (and the module_types nested in it), or every module_type '
                         'under a path in modules/ (e.g. Lotka_Volterra, neuron, cell, cell/neuron); repeatable')
    group.addoption('--component', action='append', default=[],
                    help='only test this version, <module_type>/<version> (e.g. Lotka_Volterra/Lotka1925_v01) or '
                         '<module_type>__<version>; repeatable')
    group.addoption('--include-unreviewed', action='store_true',
                    help='also run versions whose spec has reviewed: false')
    group.addoption('--quick-unreviewed', action='store_true',
                    help='for versions not reviewed yet, run only run_test and the PhLynx pipeline (CI uses '
                         'this; the full test set runs once a version is reviewed)')
    group.addoption('--shard', default=None, metavar='K/N',
                    help='only the versions of shard K of N of the selected module_types (whole module_types, '
                         'balanced by ci/test_durations.json; see cam_testing/shards.py)')
    group.addoption('--shard-stage', default='vv', choices=('vv', 'pipeline', 'omex'),
                    help='whose durations balance --shard (default vv)')
    group.addoption('--changed-only', action='store_true',
                    help='skip versions unchanged since their last run in the ledger ($CAM_LEDGER_DIR; '
                         'cam_testing/ledger.py), stage --shard-stage; reviewed versions that failed always run')


def pytest_configure(config):
    from cam_testing import paths
    root, libs = config.getoption('--cam-root'), config.getoption('--cam-library')
    if root is not None or libs is not None:
        paths.configure(repo_root=root, library_dirs=libs)
    else:
        # resolved from the environment / pyproject.toml; exported for subprocesses (bridges, tools)
        paths.configure()
    for pattern in PYTHON_FUNCTIONS:
        if pattern not in config.getini('python_functions'):
            config.addinivalue_line('python_functions', pattern)
    for marker in MARKERS:
        config.addinivalue_line('markers', marker)


def pytest_report_header(config):
    from cam_testing import paths
    r = paths.roots()
    lines = [f'cam_testing: repo {r.repo_root} ({r.source.get("repo_root")})']
    lines += [f'cam_testing: extra library {d}' for d in r.extra_library_dirs]
    if config.getoption('--changed-only', default=False):
        from cam_testing import ledger
        lines.append(f'cam_testing: --changed-only, skipping versions unchanged since the ledger in {ledger.ledger_dir()}')
    return lines


# ---- selection ---------------------------------------------------------------------------------

def _wanted(config, version):
    wanted = set(config.getoption('--component'))
    return not wanted or version.key in wanted or version.id in wanted


def _marks(config, version):
    include_unreviewed = config.getoption('--include-unreviewed') or bool(config.getoption('--module')) \
        or bool(config.getoption('--component'))
    if version.spec.get('skip'):
        return [pytest.mark.skip(reason=version.spec['skip'])]
    if not version.reviewed and not include_unreviewed:
        return [pytest.mark.skip(reason=f'{version.key} not reviewed yet (reviewed: false in its spec)')]
    return []


_loaded = {}


def _versions(config, supermodules=None):
    '''The selected versions of the repo's modules/; supermodules=True/False: only supermodule / only
    component versions.'''
    from cam_testing.library import load_module_type, select_module_types
    shard = config.getoption('--shard')
    key = (tuple(config.getoption('--module')), _roots_key(), shard)
    if key not in _loaded:     # loaded once per session, not once per test function
        names = select_module_types(key[0])
        if shard:
            from cam_testing.shards import shard_module_types
            wanted = shard_module_types(key[0], config.getoption('--shard-stage'), shard)
            names = [n for n in names if n in wanted]
        _loaded[key] = [v for name in names for v in load_module_type(name).versions()]
    return [v for v in _loaded[key] if (supermodules is None or v.is_supermodule == supermodules) and _wanted(config, v)]


def _roots_key():
    from cam_testing import paths
    r = paths.roots()
    return r.repo_root, r.extra_library_dirs


def _selected_versions(config):
    '''(module_type, version) of every selected version, component or supermodule: the V&V and
    pipeline tests (a supermodule version is generated alone, like a component).'''
    return [pytest.param((v.vessel_type, v.name), id=v.key, marks=_marks(config, v)) for v in _versions(config)]


def _selected_instances(config):
    '''(module_type, version, instance) of every instance of a selected version.'''
    return [pytest.param((v.vessel_type, v.name, i.name), id=i.key, marks=_marks(config, v))
            for v in _versions(config) for i in v.instances()]


def _selected_supermodules(config):
    return [pytest.param((v.vessel_type, v.name), id=v.key, marks=_marks(config, v)) for v in _versions(config, True)]


def _selected_equivalences(config):
    out = []
    for v in _versions(config, True):
        for e in (v.spec.get('supermodule') or {}).get('equivalent') or []:
            out.append(pytest.param((v.vessel_type, v.name, e['model']), id=f'{v.key}/{e["model"]}',
                                    marks=_marks(config, v)))
    return out


def pytest_generate_tests(metafunc):
    if 'component_key' in metafunc.fixturenames:
        metafunc.parametrize('component_key', _selected_versions(metafunc.config), scope='module')
    if 'instance_key' in metafunc.fixturenames:
        metafunc.parametrize('instance_key', _selected_instances(metafunc.config), scope='module')
    if 'supermodule_key' in metafunc.fixturenames:
        metafunc.parametrize('supermodule_key', _selected_supermodules(metafunc.config))
    if 'equivalence_key' in metafunc.fixturenames:
        metafunc.parametrize('equivalence_key', _selected_equivalences(metafunc.config))


# ---- fixtures -----------------------------------------------------------------------------------

_models = {}


def _version_model(key, tmp_path_factory):
    # One generated model per version for the whole session: generation is the slow part.
    from cam_testing import checks
    from cam_testing.library import load_version
    if key not in _models:
        version = load_version(*key)
        _models[key] = checks.ComponentModel(version, str(tmp_path_factory.mktemp(version.id)))
    return _models[key]


@pytest.fixture
def component_model(component_key, tmp_path_factory):
    return _version_model(component_key, tmp_path_factory)


@pytest.fixture
def instance_model(instance_key, tmp_path_factory):
    '''The version's model at the instance's parameters (the default instance shares the version's).'''
    from cam_testing import checks
    from cam_testing.library import load_version
    mt, v, inst = instance_key
    version = load_version(mt, v)
    if inst == version.default_instance_name:
        return _version_model((mt, v), tmp_path_factory)
    if instance_key not in _models:
        _models[instance_key] = checks.ComponentModel(version, str(tmp_path_factory.mktemp(f'{version.id}__{inst}')),
                                                      version.instance(inst))
    return _models[instance_key]


# version keys skipped as unchanged (--changed-only): the PhLynx and CUFLynx batches leave them out
UNCHANGED = set()


def pytest_collection_modifyitems(config, items):
    quick, changed_only = config.getoption('--quick-unreviewed'), config.getoption('--changed-only')
    if not (quick or changed_only):
        return
    from cam_testing.library import load_version
    stage = config.getoption('--shard-stage')
    versions, skip_unchanged = {}, {}
    for item in items:
        callspec = getattr(item, 'callspec', None)
        if callspec is None:
            continue
        key = callspec.params.get('component_key') or callspec.params.get('instance_key') \
            or callspec.params.get('supermodule_key') or callspec.params.get('equivalence_key')
        if not isinstance(key, (tuple, list)):  # None, or NOTSET for an empty parameter set
            continue
        vkey = tuple(key[:2])
        if vkey not in versions:
            versions[vkey] = load_version(*vkey)
        version = versions[vkey]
        if changed_only:
            if vkey not in skip_unchanged:
                from cam_testing import ledger
                skip_unchanged[vkey] = ledger.unchanged(version, stage)
                if skip_unchanged[vkey][0]:
                    UNCHANGED.add(version.key)
            ok, why = skip_unchanged[vkey]
            if ok:
                item.add_marker(pytest.mark.skip(reason=why))
                continue
        # run_test and the (cheap) PhLynx -> CUFLynx pipeline run for every version
        if quick and not version.reviewed and item.originalname not in QUICK_TESTS:
            item.add_marker(pytest.mark.skip(reason=f'{"/".join(vkey)} not reviewed yet: CI runs only run_test'))
    if changed_only:
        n = sum(1 for ok, _ in skip_unchanged.values() if ok)
        config._cam_unchanged = (n, len(skip_unchanged))


def pytest_terminal_summary(terminalreporter, config):
    n = getattr(config, '_cam_unchanged', None)
    if n:
        terminalreporter.write_line(f'cam_testing --changed-only: {n[0]} of {n[1]} versions unchanged, skipped')
