"""
cam_testing.paths: where the repo under test and the extra module libraries come from, and the
library functions across several libraries (cam_testing in another repository, README.md).
"""
import json
import os

import pytest

from cam_testing import harness, library, paths
from cam_testing.paths import ENV_LIBRARIES, ENV_ROOT, resolve


def _module_type(modules, rel, versions):
    '''modules/<rel>/versions/<v>/ with a one-entry config, for each version.'''
    name = os.path.basename(rel)
    for v in versions:
        d = os.path.join(modules, rel, 'versions', v)
        os.makedirs(os.path.join(d, 'parameterisations', 'default'), exist_ok=True)
        with open(os.path.join(d, f'{name}_{v}_modules_config.json'), 'w') as f:
            json.dump([{'module_type': name, 'module_subtype': v, 'component_file': f'{name}_{v}_modules.cellml',
                        'component_type': f'{name}_{v}', 'default_parameterisation': 'default', 'variables_and_units': []}], f)


@pytest.fixture
def two_libraries(tmp_path):
    '''A repo (its modules/: cat/own, cat/shared vA) and a library (lib/modules/: bc/neighbour,
    cat/shared vB, cat/shared/nested_in_shared).'''
    repo, lib = tmp_path / 'repo', tmp_path / 'lib'
    _module_type(str(repo / 'modules'), 'cat/own', ['v1'])
    _module_type(str(repo / 'modules'), 'cat/shared', ['vA'])
    _module_type(str(lib / 'modules'), 'bc/neighbour', ['nn'])
    _module_type(str(lib / 'modules'), 'cat/shared', ['vB'])
    _module_type(str(lib / 'modules'), 'cat/shared/nested_in_shared', ['x'])
    yield repo, lib
    paths.reset()


# ---- resolution ----------------------------------------------------------------------------------

def test_environment_gives_root_and_libraries(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    (b / 'modules').mkdir(parents=True)
    (a / 'modules').mkdir(parents=True)
    r = resolve(cwd=str(tmp_path), environ={ENV_ROOT: str(tmp_path / 'repo'),
                                             ENV_LIBRARIES: os.pathsep.join([str(a / 'modules'), str(b)])})
    assert r.repo_root == str(tmp_path / 'repo')
    assert r.modules_dir == str(tmp_path / 'repo' / 'modules')
    assert r.system_models_dir == str(tmp_path / 'repo' / 'system_models')
    assert r.site_dir == str(tmp_path / 'repo' / 'site')
    assert r.manifests_dir == str(tmp_path / 'repo' / 'manifests')
    # a library given as a repo holding modules/ is that modules/
    assert r.extra_library_dirs == (str(a / 'modules'), str(b / 'modules'))
    assert r.library_dirs == [r.modules_dir, str(a / 'modules'), str(b / 'modules')]
    assert r.source == {'repo_root': ENV_ROOT, 'extra_library_dirs': ENV_LIBRARIES}


def test_arguments_win_over_the_environment(tmp_path):
    r = resolve(repo_root=str(tmp_path / 'x'), library_dirs=[], cwd=str(tmp_path),
                environ={ENV_ROOT: str(tmp_path / 'y'), ENV_LIBRARIES: str(tmp_path / 'z')})
    assert (r.repo_root, r.extra_library_dirs) == (str(tmp_path / 'x'), ())


def test_pyproject_table(tmp_path):
    (tmp_path / 'pyproject.toml').write_text(
        '[tool.cam_testing]\nroot = "inner"\nmodule_library_dirs = ["../lib/modules"]\n'
        '[tool.cam_testing.manifests]\n"index.json" = ["vessels"]\n')
    sub = tmp_path / 'inner' / 'tests'
    sub.mkdir(parents=True)
    r = resolve(cwd=str(sub), environ={})
    assert r.repo_root == str(tmp_path / 'inner')
    assert r.extra_library_dirs == (str(tmp_path.parent / 'lib' / 'modules'),)
    assert r.settings == {'manifests': {'index.json': ['vessels']}}
    # the environment still wins
    assert resolve(cwd=str(sub), environ={ENV_ROOT: str(tmp_path)}).repo_root == str(tmp_path)


def test_pyproject_without_the_table_is_ignored(tmp_path):
    (tmp_path / 'pyproject.toml').write_text('[project]\nname = "x"\n')
    (tmp_path / 'modules').mkdir()
    (tmp_path / 'deep').mkdir()
    r = resolve(cwd=str(tmp_path / 'deep'), environ={})
    assert r.repo_root == str(tmp_path)          # the nearest directory with modules/
    assert r.extra_library_dirs == ()


def test_default_is_the_working_directory_with_modules(tmp_path):
    (tmp_path / 'modules').mkdir()
    assert resolve(cwd=str(tmp_path), environ={}).repo_root == str(tmp_path)


def test_primary_library_is_not_an_extra_one(tmp_path):
    (tmp_path / 'modules').mkdir()
    r = resolve(cwd=str(tmp_path), environ={ENV_LIBRARIES: str(tmp_path)})
    assert r.extra_library_dirs == ()


def test_this_repo_resolves_to_itself():
    '''This repo's own runs (no CAM_* settings): the checkout, as before cam_testing.paths.'''
    r = resolve(cwd=paths.PACKAGE_CHECKOUT, environ={})
    assert r.repo_root == paths.PACKAGE_CHECKOUT
    assert r.extra_library_dirs == ()
    assert r.directory_schema_path == os.path.join(paths.PACKAGE_CHECKOUT, 'modules', 'directory_schema.json')


def test_packaged_directory_schema_is_this_repos():
    with open(paths.PACKAGED_DIRECTORY_SCHEMA) as a, \
            open(os.path.join(paths.PACKAGE_CHECKOUT, 'modules', 'directory_schema.json')) as b:
        assert json.load(a) == json.load(b), 'cam_testing/data/directory_schema.json is a copy of modules/directory_schema.json'


def test_repo_without_a_schema_uses_the_packaged_one(tmp_path):
    (tmp_path / 'modules').mkdir()
    assert resolve(cwd=str(tmp_path), environ={}).directory_schema_path == paths.PACKAGED_DIRECTORY_SCHEMA


def test_configure_exports_for_subprocesses(tmp_path, monkeypatch):
    monkeypatch.delenv(ENV_ROOT, raising=False)
    monkeypatch.delenv(ENV_LIBRARIES, raising=False)
    try:
        r = paths.configure(repo_root=str(tmp_path), library_dirs=[str(tmp_path / 'lib')])
        assert os.environ[ENV_ROOT] == str(tmp_path)
        assert os.environ[ENV_LIBRARIES] == str(tmp_path / 'lib')
        assert library.MODULES_DIR == r.modules_dir == str(tmp_path / 'modules')
    finally:
        paths.reset()


# ---- the library across several libraries --------------------------------------------------------

def test_primary_and_extra_libraries(two_libraries):
    repo, lib = two_libraries
    paths.configure(repo_root=str(repo), library_dirs=[str(lib)], export=False)
    assert library.library_dirs() == [str(repo / 'modules'), str(lib / 'modules')]
    # tested / reported: the repo's own module_types
    assert library.module_type_names() == ['own', 'shared']
    assert [v.key for v in library.all_versions()] == ['own/v1', 'shared/vA']
    # loadable: every library's (a harness, supermodule or system model may use them)
    assert library.module_type_names(include_libraries=True) == ['neighbour', 'nested_in_shared', 'own', 'shared']
    assert sorted(v.key for v in library.all_versions(include_libraries=True)) == \
        ['neighbour/nn', 'nested_in_shared/x', 'own/v1', 'shared/vA', 'shared/vB']
    n = library.load_version('neighbour', 'nn')
    assert n.dir == str(lib / 'modules' / 'bc' / 'neighbour' / 'versions' / 'nn')
    assert not n.mtype.is_primary and n.mtype.relpath == 'bc/neighbour' and n.category == 'bc'
    # a module_type in both libraries: each version from the library that has it
    assert library.load_version('shared', 'vA').mtype.is_primary
    assert library.load_version('shared', 'vB').dir == str(lib / 'modules' / 'cat' / 'shared' / 'versions' / 'vB')
    assert library.module_type_dirs('shared') == [str(repo / 'modules' / 'cat' / 'shared'),
                                                  str(lib / 'modules' / 'cat' / 'shared')]
    assert library.parent_of('nested_in_shared') == 'shared'
    assert library.select_module_types(['cat']) == ['own', 'shared']
    assert set(library.version_index(include_libraries=True)) == \
        {('own', 'v1'), ('shared', 'vA'), ('shared', 'vB'), ('neighbour', 'nn'), ('nested_in_shared', 'x')}


def test_reconfigure_clears_discovery(two_libraries):
    repo, lib = two_libraries
    paths.configure(repo_root=str(repo), library_dirs=[], export=False)
    assert library.module_type_names(include_libraries=True) == ['own', 'shared']
    paths.configure(repo_root=str(repo), library_dirs=[str(lib)], export=False)
    assert 'neighbour' in library.module_type_names(include_libraries=True)


def test_former_parameterisation_names_are_still_read(tmp_path):
    '''A library on the layout before 2026-10-09: the version's parameter sets in instances/ (now
    parameterisations/), the config's "default_instance" (now "default_parameterisation") and a
    module-array record's "instance" (now "parameterisation").'''
    from cam_testing import module_array
    repo = tmp_path / 'repo'
    d = repo / 'modules' / 'cat' / 'old' / 'versions' / 'v1'
    (d / 'instances' / 'fit').mkdir(parents=True)
    (d / 'instances' / 'fit' / 'fit_parameters.csv').write_text('variable_name,units,value,data_reference\n'
                                                                'a,dimensionless,2,test\n')
    (d / 'old_v1_modules_config.json').write_text(json.dumps([{
        'module_type': 'old', 'module_subtype': 'v1', 'component_file': 'old_v1_modules.cellml',
        'component_type': 'old_v1', 'default_instance': 'fit', 'variables_and_units': []}]))
    try:
        paths.configure(repo_root=str(repo), library_dirs=[], export=False)
        v = library.load_version('old', 'v1')
        assert v.parameterisations_dir == str(d / 'instances')
        assert v.default_parameterisation_name == 'fit' and v.parameterisation_names() == ['fit']
        assert [p.value for p in v.default_parameterisation.parameters()] == ['2']
        # the former method names are aliases
        assert v.default_instance_name == 'fit' and v.instance('fit').name == 'fit'
    finally:
        paths.reset()
    record = {'name': 'a', 'module_type': 'old', 'module_subtype': 'v1', 'instance': 'fit'}
    assert module_array.parameterisation_of(record) == 'fit'
    assert module_array.normalise_record(record) == {'name': 'a', 'module_type': 'old', 'module_subtype': 'v1',
                                                     'parameterisation': 'fit', 'inp_instances': [],
                                                     'out_instances': []}


def test_libcuflynx_gets_every_library(two_libraries, monkeypatch, tmp_path):
    '''harness.generate (and the system models, the C++ run_test) pass module_library_dirs = the repo's
    modules/ then the extra libraries.'''
    repo, lib = two_libraries
    paths.configure(repo_root=str(repo), library_dirs=[str(lib)], export=False)
    seen = {}

    def fake_generate(_inp, config):
        seen.update(config)
        return True

    module = pytest.importorskip('libcuflynx.scripts.script_generate_with_new_architecture')
    monkeypatch.setattr(module, 'generate_with_new_architecture', fake_generate)
    monkeypatch.setattr(harness, '_write_resources', lambda *a, **k: None)
    version = library.load_version('own', 'v1')
    harness.generate(version, str(tmp_path / 'work'))
    assert seen['module_library_dirs'] == [str(repo / 'modules'), str(lib / 'modules')]
    assert seen['use_builtin_modules'] is False

    from cam_testing import system as systems
    s = systems.System('m', 'c', str(tmp_path / 'sys'), {})
    seen.clear()
    systems.generate(s, str(tmp_path / 'work_sys'))
    assert seen['module_library_dirs'] == [str(repo / 'modules'), str(lib / 'modules')]
    seen.clear()
    systems.generate(s, str(tmp_path / 'work_sys'), reference=True)
    assert 'module_library_dirs' not in seen      # circulatory_autogen's original: its bundled modules


def test_paths_command_prints_the_roots(two_libraries, capsys):
    repo, lib = two_libraries
    paths.configure(repo_root=str(repo), library_dirs=[str(lib)], export=False)
    paths.main([])
    out = capsys.readouterr().out
    assert str(repo) in out and str(lib / 'modules') in out


def test_plugin_registered_once(request):
    '''The pytest11 entry point and pyproject's -p name the plugin the same way (registered once).'''
    from cam_testing import pytest_plugin
    pm = request.config.pluginmanager
    assert pm.get_plugin('cam_testing.pytest_plugin') is pytest_plugin
    assert [p for p in pm.get_plugins() if p is pytest_plugin] == [pytest_plugin]
