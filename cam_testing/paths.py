"""
Where the library under test is: the repo root (holding modules/, system_models/, reviews/,
manifests/ and site/) and the extra module libraries loaded with it.

cam_testing tests, reports and builds manifests for one repo's modules/ (the *primary* library).
A repo whose modules use another library's modules (in a harness, a supermodule or a system
model; e.g. this library's boundary conditions) lists that library as an extra one: libcuflynx
gets every library in ``module_library_dirs`` (primary first) and cam_testing can load their
versions, but only the primary library is tested and reported.

Each setting is taken from the first of these that gives it:

  1. ``configure()`` (the pytest options ``--cam-root`` / ``--cam-library``);
  2. the environment: ``CAM_REPO_ROOT``, and ``CAM_MODULE_LIBRARY_DIRS`` (os.pathsep-separated);
  3. the ``[tool.cam_testing]`` table of the nearest pyproject.toml at or above the working
     directory: ``root`` and ``module_library_dirs`` (relative to that file; ``root`` defaults
     to its directory);
  4. the repo root only: the nearest directory at or above the working directory with a
     modules/ directory, else the checkout cam_testing is installed from (when it has modules/),
     else the working directory.

An extra library is a modules/ directory, or a repo holding one.

    python -m cam_testing.paths          # prints what is in use, and where it came from
"""
import functools
import os
from dataclasses import dataclass, field

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
# the checkout cam_testing is installed from (editable install or source tree); site-packages otherwise
PACKAGE_CHECKOUT = os.path.dirname(PACKAGE_DIR)
# the directory schema shipped with the package, for a repo without modules/directory_schema.json
PACKAGED_DIRECTORY_SCHEMA = os.path.join(PACKAGE_DIR, 'data', 'directory_schema.json')

ENV_ROOT = 'CAM_REPO_ROOT'
ENV_LIBRARIES = 'CAM_MODULE_LIBRARY_DIRS'
PYPROJECT_TABLE = 'cam_testing'


@dataclass(frozen=True)
class Roots:
    repo_root: str
    extra_library_dirs: tuple = ()
    source: dict = field(default_factory=dict, compare=False)   # setting -> where it came from
    settings: dict = field(default_factory=dict, compare=False)  # the rest of [tool.cam_testing]

    @property
    def modules_dir(self):
        return os.path.join(self.repo_root, 'modules')

    @property
    def system_models_dir(self):
        return os.path.join(self.repo_root, 'system_models')

    @property
    def reviews_dir(self):
        return os.path.join(self.repo_root, 'reviews')

    @property
    def manifests_dir(self):
        return os.path.join(self.repo_root, 'manifests')

    @property
    def site_dir(self):
        return os.path.join(self.repo_root, 'site')

    @property
    def directory_schema_path(self):
        '''The repo's modules/directory_schema.json, else the one shipped with cam_testing.'''
        own = os.path.join(self.modules_dir, 'directory_schema.json')
        return own if os.path.isfile(own) else PACKAGED_DIRECTORY_SCHEMA

    @property
    def library_dirs(self):
        '''Every module library, primary first: what libcuflynx gets as module_library_dirs.'''
        return [self.modules_dir, *self.extra_library_dirs]


def _library_dir(path, base):
    '''An extra library given as a modules/ directory or as a repo holding one -> the modules/ directory.'''
    p = os.path.abspath(os.path.join(base, os.path.expanduser(path)))
    inner = os.path.join(p, 'modules')
    if os.path.basename(p) != 'modules' and os.path.isdir(inner):
        return inner
    return p


def _clean_libraries(dirs, modules_dir):
    out = []
    for d in dirs:
        if os.path.realpath(d) != os.path.realpath(modules_dir) and d not in out:
            out.append(d)
    return tuple(out)


def _read_toml(path):
    try:
        import tomllib
    except ModuleNotFoundError:      # Python 3.10
        try:
            import tomli as tomllib
        except ModuleNotFoundError:
            return {}
    with open(path, 'rb') as f:
        return tomllib.load(f)


def find_pyproject_settings(start=None):
    '''(the [tool.cam_testing] table, the directory of its pyproject.toml) of the nearest
    pyproject.toml at or above ``start`` (default: the working directory); ({}, None) when the
    nearest pyproject.toml has no such table, or there is none.'''
    d = os.path.abspath(start or os.getcwd())
    while True:
        path = os.path.join(d, 'pyproject.toml')
        if os.path.isfile(path):
            table = (_read_toml(path).get('tool') or {}).get(PYPROJECT_TABLE)
            return (table, d) if isinstance(table, dict) else ({}, None)
        parent = os.path.dirname(d)
        if parent == d:
            return {}, None
        d = parent


def _nearest_with_modules(start):
    d = os.path.abspath(start)
    while True:
        if os.path.isdir(os.path.join(d, 'modules')):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def resolve(repo_root=None, library_dirs=None, cwd=None, environ=None):
    '''The Roots from the arguments, the environment, pyproject.toml and the defaults (module docstring).'''
    environ = os.environ if environ is None else environ
    cwd = os.path.abspath(cwd or os.getcwd())
    table, table_dir = find_pyproject_settings(cwd)
    source = {}

    if repo_root:
        root, source['repo_root'] = os.path.abspath(os.path.join(cwd, repo_root)), 'configure()'
    elif environ.get(ENV_ROOT):
        root, source['repo_root'] = os.path.abspath(os.path.join(cwd, environ[ENV_ROOT])), ENV_ROOT
    elif table_dir is not None:
        root = os.path.abspath(os.path.join(table_dir, table.get('root', '.')))
        source['repo_root'] = os.path.join(table_dir, 'pyproject.toml')
    elif _nearest_with_modules(cwd):
        root, source['repo_root'] = _nearest_with_modules(cwd), 'the working directory (has modules/)'
    elif os.path.isdir(os.path.join(PACKAGE_CHECKOUT, 'modules')):
        root, source['repo_root'] = PACKAGE_CHECKOUT, 'the checkout cam_testing is installed from'
    else:
        root, source['repo_root'] = cwd, 'the working directory'

    if library_dirs is not None:
        extra, source['extra_library_dirs'] = [_library_dir(p, cwd) for p in library_dirs], 'configure()'
    elif environ.get(ENV_LIBRARIES):
        extra = [_library_dir(p, cwd) for p in environ[ENV_LIBRARIES].split(os.pathsep) if p.strip()]
        source['extra_library_dirs'] = ENV_LIBRARIES
    elif table_dir is not None and table.get('module_library_dirs'):
        dirs = table['module_library_dirs']
        extra = [_library_dir(p, table_dir) for p in ([dirs] if isinstance(dirs, str) else dirs)]
        source['extra_library_dirs'] = os.path.join(table_dir, 'pyproject.toml')
    else:
        extra, source['extra_library_dirs'] = [], 'none given'

    settings = {k: v for k, v in (table or {}).items() if k not in ('root', 'module_library_dirs')}
    return Roots(root, _clean_libraries(extra, os.path.join(root, 'modules')), source, settings)


_override = {}


@functools.lru_cache(maxsize=1)
def _cached():
    return resolve(**_override)


def roots():
    '''The Roots in use (resolved once, until configure()).'''
    return _cached()


def configure(repo_root=None, library_dirs=None, export=True):
    '''
    Use this repo root and these extra libraries (None: as resolved without them). ``export``: also
    set CAM_REPO_ROOT / CAM_MODULE_LIBRARY_DIRS, so subprocesses (bridges, tools) use the same.
    Returns the new Roots. Clears cam_testing.library's discovery caches.
    '''
    _override.clear()
    if repo_root is not None:
        _override['repo_root'] = repo_root
    if library_dirs is not None:
        _override['library_dirs'] = list(library_dirs)
    _cached.cache_clear()
    r = roots()
    if export:
        os.environ[ENV_ROOT] = r.repo_root
        os.environ[ENV_LIBRARIES] = os.pathsep.join(r.extra_library_dirs)
    try:
        from cam_testing import library
        library.clear_caches()
    except ImportError:
        pass
    return r


def reset():
    '''Forget configure(): resolve from the environment and pyproject.toml again (tests use this).'''
    _override.clear()
    _cached.cache_clear()
    try:
        from cam_testing import library
        library.clear_caches()
    except ImportError:
        pass


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Print the repo root and module libraries cam_testing uses')
    parser.parse_args(argv)
    r = roots()
    print(f'repo root:          {r.repo_root}   ({r.source.get("repo_root")})')
    print(f'modules (primary):  {r.modules_dir}{"" if os.path.isdir(r.modules_dir) else "   (missing)"}')
    for d in r.extra_library_dirs:
        print(f'extra library:      {d}{"" if os.path.isdir(d) else "   (missing)"}   ({r.source.get("extra_library_dirs")})')
    if not r.extra_library_dirs:
        print('extra libraries:    none')
    print(f'system models:      {r.system_models_dir}')
    print(f'directory schema:   {r.directory_schema_path}')


if __name__ == '__main__':
    main()
