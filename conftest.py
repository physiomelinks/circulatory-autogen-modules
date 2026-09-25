import pytest

from cam_testing import checks
from cam_testing.library import load_module, module_names


def pytest_addoption(parser):
    parser.addoption('--module', action='append', default=[],
                     help='only test this module (directory name under modules/); repeatable')
    parser.addoption('--component', action='append', default=[],
                     help='only test this component id, e.g. Lotka_Volterra__nn; repeatable')
    parser.addoption('--include-unreviewed', action='store_true',
                     help='also run modules whose spec has reviewed: false')
    parser.addoption('--quick-unreviewed', action='store_true',
                     help='for modules not reviewed yet, run only run_test (CI uses this; the full '
                          'test set runs once a module is reviewed)')


def _selected_components(config):
    wanted_modules = set(config.getoption('--module'))
    wanted_components = set(config.getoption('--component'))
    include_unreviewed = config.getoption('--include-unreviewed') or bool(wanted_modules)
    params = []
    for name in module_names():
        if wanted_modules and name not in wanted_modules:
            continue
        module = load_module(name)
        for component in module.components():
            if wanted_components and component.id not in wanted_components:
                continue
            marks = []
            if component.spec.get('skip'):
                marks.append(pytest.mark.skip(reason=component.spec['skip']))
            elif not module.reviewed and not include_unreviewed:
                marks.append(pytest.mark.skip(reason=f'{name} not reviewed yet (reviewed: false in its spec)'))
            params.append(pytest.param((name, component.id), id=f'{name}/{component.id}', marks=marks))
    return params


def pytest_generate_tests(metafunc):
    if 'component_key' in metafunc.fixturenames:
        metafunc.parametrize('component_key', _selected_components(metafunc.config), scope='module')


_models = {}


@pytest.fixture
def component_model(component_key, tmp_path_factory):
    # One generated model per component for the whole session: generation is the slow part.
    if component_key not in _models:
        module_name, component_id = component_key
        component = load_module(module_name).component(component_id)
        work_dir = str(tmp_path_factory.mktemp(component_id))
        _models[component_key] = checks.ComponentModel(component, work_dir)
    return _models[component_key]


def pytest_collection_modifyitems(config, items):
    if not config.getoption('--quick-unreviewed'):
        return
    reviewed = {}
    for item in items:
        callspec = getattr(item, 'callspec', None)
        if callspec is None or 'component_key' not in callspec.params:
            continue
        module_name = callspec.params['component_key'][0]
        if module_name not in reviewed:
            reviewed[module_name] = load_module(module_name).reviewed
        if not reviewed[module_name] and item.originalname != 'run_test':
            item.add_marker(pytest.mark.skip(reason=f'{module_name} not reviewed yet: CI runs only run_test'))
