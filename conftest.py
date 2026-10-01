import pytest

from cam_testing import checks
from cam_testing.library import load_version, select_module_types, load_module_type


def pytest_addoption(parser):
    parser.addoption('--module', action='append', default=[],
                     help='only test this module_type, or every module_type under a category path '
                          '(e.g. Lotka_Volterra, cell, cell/neurons); repeatable')
    parser.addoption('--component', action='append', default=[],
                     help='only test this version, <module_type>/<version> (e.g. Lotka_Volterra/nn) or '
                          '<module_type>__<version>; repeatable')
    parser.addoption('--include-unreviewed', action='store_true',
                     help='also run versions whose spec has reviewed: false')
    parser.addoption('--quick-unreviewed', action='store_true',
                     help='for versions not reviewed yet, run only run_test and the PhLynx pipeline (CI uses '
                          'this; the full test set runs once a version is reviewed)')


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
    '''The selected versions; supermodules=True/False: only supermodule / only component versions.'''
    key = tuple(config.getoption('--module'))
    if key not in _loaded:     # loaded once per session, not once per test function
        _loaded[key] = [v for name in select_module_types(key) for v in load_module_type(name).versions()]
    return [v for v in _loaded[key] if (supermodules is None or v.is_supermodule == supermodules) and _wanted(config, v)]


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


_models = {}
QUICK_TESTS = ('run_test', 'supermodule_structure_test', 'phlynx_export_test', 'cuflynx_simulate_test', 'phlynx_equivalence_test',
               'cuflynx_instance_omex_test')


def _version_model(key, tmp_path_factory):
    # One generated model per version for the whole session: generation is the slow part.
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
    mt, v, inst = instance_key
    version = load_version(mt, v)
    if inst == version.default_instance_name:
        return _version_model((mt, v), tmp_path_factory)
    if instance_key not in _models:
        _models[instance_key] = checks.ComponentModel(version, str(tmp_path_factory.mktemp(f'{version.id}__{inst}')),
                                                      version.instance(inst))
    return _models[instance_key]


def pytest_collection_modifyitems(config, items):
    if not config.getoption('--quick-unreviewed'):
        return
    reviewed = {}
    for item in items:
        callspec = getattr(item, 'callspec', None)
        if callspec is None:
            continue
        key = callspec.params.get('component_key') or callspec.params.get('instance_key') \
            or callspec.params.get('supermodule_key') or callspec.params.get('equivalence_key')
        if key is None:
            continue
        vkey = tuple(key[:2])
        if vkey not in reviewed:
            reviewed[vkey] = load_version(*vkey).reviewed
        # run_test and the (cheap) PhLynx -> CUFLynx pipeline run for every version
        if not reviewed[vkey] and item.originalname not in QUICK_TESTS:
            item.add_marker(pytest.mark.skip(reason=f'{"/".join(vkey)} not reviewed yet: CI runs only run_test'))
