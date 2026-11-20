'''Versions renamed in module PRs (a library's version_renames.json) extend the legacy rename map.'''
import json

from cam_testing import library


def _with(monkeypatch, tmp_path, entries, earlier):
    (tmp_path / library.VERSION_RENAMES).write_text(json.dumps(entries))
    monkeypatch.setattr(library, 'library_dirs', lambda: [str(tmp_path)])
    monkeypatch.setattr(library, 'restructure_renames', lambda: earlier)
    monkeypatch.setattr(library.paths, 'PACKAGE_CHECKOUT', str(tmp_path))      # no restructure_map.yaml
    monkeypatch.setattr(library, 'HEART_VERSION_RENAMES', {})
    library.legacy_renames.cache_clear()


def test_an_old_name_follows_a_module_pr_rename(monkeypatch, tmp_path):
    _with(monkeypatch, tmp_path,
          [{'from': ['cardiac_clock', 'controlled'], 'to': ['cardiac_clock', 'ArgusUNPUBLISHED_v01'], 'reason': 'r'}],
          {('cardiac_clock', 'nn_controlled'): ('cardiac_clock', 'controlled')})
    try:
        renames = library.legacy_renames()
        # circulatory_autogen's name and the library's previous name both reach the new one
        assert renames[('cardiac_clock', 'nn_controlled')] == ('cardiac_clock', 'ArgusUNPUBLISHED_v01')
        assert renames[('cardiac_clock', 'controlled')] == ('cardiac_clock', 'ArgusUNPUBLISHED_v01')
    finally:
        library.legacy_renames.cache_clear()


def test_renames_chain_in_order(monkeypatch, tmp_path):
    _with(monkeypatch, tmp_path, [{'from': ['m', 'a'], 'to': ['m', 'b'], 'reason': 'r'},
                                  {'from': ['m', 'b'], 'to': ['m', 'c'], 'reason': 'r'}], {})
    try:
        assert library.legacy_renames()[('m', 'a')] == ('m', 'c')
        assert library.version_renames() == {('m', 'a'): ('m', 'b'), ('m', 'b'): ('m', 'c')}
    finally:
        library.legacy_renames.cache_clear()
