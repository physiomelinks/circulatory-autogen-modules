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
