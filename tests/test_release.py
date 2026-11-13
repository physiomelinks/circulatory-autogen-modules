'''The release: modules/renames.json in step with legacy_renames(), the module pairs libcuflynx writes
itself, and a reproducible slim release asset (cam_testing.release).'''
import hashlib
import json
import os
import subprocess
import tarfile

import pytest

from cam_testing import library, release

ROOT = release.repo_root()

# (module_type, module_subtype) pairs libcuflynx names in its code (1D coupling: ModelParsers,
# convert_0d_to_1d, CVSCellMLGenerator), so every release must have them
LIBCUFLYNX_PAIRS = [('FV1D_vessel', 'vp'), ('FV1D_vessel', 'pv'), ('FV1D_vessel', 'vv'), ('FV1D_vessel', 'nn'),
                    ('FV1D_solver', 'nn'), ('FV1D_volume_sum', 'nn')]
# module_types libcuflynx's C++ generator treats specially (generators/cpp/api.py PROCESS_COORDINATORS)
LIBCUFLYNX_TYPES = ['coupler']


def committed(path):
    return subprocess.run(['git', '-C', ROOT, 'ls-files', '--error-unmatch', path], capture_output=True).returncode == 0


def test_renames_json_is_in_step():
    path = os.path.join(ROOT, 'modules', release.RENAMES_FILE)
    with open(path) as f:
        text = f.read()
    assert text == release.renames_text(), \
        'modules/renames.json is out of step with legacy_renames(): run python -m cam_testing.release renames'


def test_renames_name_versions_the_library_has():
    index = library.version_index()
    doc = release.renames_document()
    assert doc['renames'], 'no renames'
    missing = [r['to'] for r in doc['renames'] if tuple(r['to']) not in index]
    assert not missing, f'renames to versions the library does not have: {missing}'
    still = [r['from'] for r in doc['renames'] if tuple(r['from']) in index]
    assert not still, f'renamed pairs the library still has (a name would mean two versions): {still}'


def test_builtin_same_maths_versions_keep_the_builtin_component():
    '''The versions BUILTIN_SAME_MATHS maps circulatory_autogen's built-ins to are where renames.json
    sends them, and each is a version of its own (its component named <built-in>_OLD).'''
    renames = {tuple(r['from']): tuple(r['to']) for r in release.renames_document()['renames']}
    index = library.version_index()
    for old, new in library.BUILTIN_SAME_MATHS.items():
        assert renames[old] == new
        with open(index[new].path('modules_config.json')) as f:
            entry = json.load(f)[0]
        assert entry['component_type'] == f'{old[0]}_OLD'


@pytest.mark.parametrize('pair', LIBCUFLYNX_PAIRS, ids=lambda p: '/'.join(p))
def test_libcuflynx_pairs_exist(pair):
    assert pair in library.version_index(), f'libcuflynx writes {pair} itself; the library must keep it'


@pytest.mark.parametrize('module_type', LIBCUFLYNX_TYPES)
def test_libcuflynx_module_types_exist(module_type):
    assert any(t == module_type for t, _ in library.version_index()), \
        f'libcuflynx treats module_type {module_type} specially; the library must keep it'


@pytest.mark.parametrize('tag,ok', [('v2026.11.0', True), ('v2026.01.12', True), ('v2026.13.0', False),
                                    ('2026.11.0', False), ('v2026.11', False), ('v1.2.3', False)])
def test_tags_are_calendar_tags(tag, ok):
    if ok:
        assert release.check_tag(tag) == tag
    else:
        with pytest.raises(ValueError):
            release.check_tag(tag)


@pytest.mark.parametrize('path,ok', [
    ('modules/heart/versions/vp/heart_vp_modules.cellml', True),
    ('modules/heart/versions/vp/parameterisations/default/default_parameters.json', True),
    ('modules/heart/versions/vp/heart_vp_references.bib', True),
    ('modules/renames.json', True),
    ('manifests/all.json', True),
    ('LICENSE', True),
    ('modules/heart/versions/vp/heart_vp_tests.yaml', False),
    ('modules/heart/versions/vp/heart_vp_verification_config.json', False),
    ('modules/heart/versions/vp/heart_vp_references_proposed.bib', False),
    ('modules/heart/versions/vp/risk/heart_vp_risk.json', False),
    ('modules/heart/versions/vp/figures/f.png', False),
    ('modules/README.md', False),
    ('cam_testing/release.py', False),
    ('system_models/cellular/SN_simple/SN_simple_module_array.json', False),
])
def test_asset_content(path, ok):
    assert release.in_asset(path) is ok


@pytest.mark.skipif(not committed('modules/renames.json'), reason='modules/renames.json is not committed yet')
def test_asset_is_reproducible_and_complete(tmp_path):
    path, digest, rel = release.build_asset('v2026.11.0', str(tmp_path / 'a'), root=ROOT)
    path2, digest2, _ = release.build_asset('v2026.11.0', str(tmp_path / 'b'), root=ROOT)
    assert digest == digest2, 'building the same commit twice gave different bytes'
    with open(path, 'rb') as f:
        assert hashlib.sha256(f.read()).hexdigest() == digest
    with open(path + '.sha256') as f:
        assert f.read() == f'{digest}  cam-modules-v2026.11.0.tar.gz\n'
    top = 'cam-modules-v2026.11.0/'
    with tarfile.open(path) as tar:
        members = tar.getmembers()
        names = [m.name for m in members]
        assert all(n.startswith(top) and m.isfile() for n, m in zip(names, members))
        files = {n[len(top):]: tar.extractfile(m).read() for n, m in zip(names, members)}
    meta = json.loads(files.pop(release.RELEASE_FILE))
    assert meta['tag'] == 'v2026.11.0' and meta['libcuflynx'] == rel['libcuflynx']
    assert meta['files'] == len(files) == len(meta['sha256'])
    assert all(hashlib.sha256(files[p]).hexdigest() == h for p, h in meta['sha256'].items())
    assert 'modules/renames.json' in files and 'LICENSE' in files
    assert all(release.in_asset(p) for p in files)
    # every version's config, CellML and units are in it
    for v in library.version_index().values():
        rel_dir = os.path.relpath(v.dir, ROOT)
        for kind in ('modules_config.json', 'modules.cellml', 'units.cellml'):
            p = os.path.relpath(v.path(kind), ROOT)
            if committed(p):
                assert p in files, f'{p} missing from the asset ({rel_dir})'
