"""
A system model's reference: circulatory_autogen's original (a vessel array to generate) or a ready
CellML model simulated as it is (equivalence.reference_cellml, or reference/<model>.cellml). The
end-to-end run with a real CellML reference is in tests/test_external_repo.py.
"""
import numpy as np
import pytest

from cam_testing import system as systems


def _system(tmp_path, spec=None, files=()):
    d = tmp_path / 'm'
    (d / 'reference').mkdir(parents=True)
    for f in files:
        (d / f).write_text('x')
    return systems.System('m', 'c', str(d), spec or {})


def test_vessel_array_reference(tmp_path):
    s = _system(tmp_path, files=['reference/m_vessel_array.json', 'reference/m.cellml'])
    assert systems.reference_cellml(s) is None      # a vessel array is there: generate it


def test_implicit_cellml_reference(tmp_path):
    s = _system(tmp_path, files=['reference/m.cellml'])
    assert systems.reference_cellml(s) == str(tmp_path / 'm' / 'reference' / 'm.cellml')


def test_explicit_cellml_reference(tmp_path):
    s = _system(tmp_path, {'equivalence': {'reference_cellml': 'reference/phlynx-export.cellml'}},
                files=['reference/phlynx-export.cellml', 'reference/m_vessel_array.json'])
    assert systems.reference_cellml(s) == str(tmp_path / 'm' / 'reference' / 'phlynx-export.cellml')
    missing = _system(tmp_path / 'x', {'equivalence': {'reference_cellml': 'reference/nope.cellml'}})
    with pytest.raises(FileNotFoundError):
        systems.reference_cellml(missing)


def test_equivalence_simulates_the_cellml_reference(tmp_path, monkeypatch):
    s = _system(tmp_path, {'equivalence': {'reference_cellml': 'reference/r.cellml', 'output_map': {'a/x': 'b/x'},
                                           'ignore': {'a/t': 'time'}, 'tol': 1e-6}},
                files=['reference/r.cellml'])
    calls = []

    def fake_generate(system, work_dir, reference=False, quiet=True):
        calls.append(reference)
        return 'library.cellml'

    def fake_simulate(path, spec, solver_info, names=None, ignore_components=None):
        t = np.linspace(0, 1, 5)
        if path.endswith('r.cellml'):
            assert ignore_components == systems.REFERENCE_CELLML_IGNORED_COMPONENTS
            return t, {'a/x': t, 'a/t': t}
        assert names == ['a/t', 'b/x']
        return t, {'b/x': t + 1e-9}

    monkeypatch.setattr(systems, 'generate', fake_generate)
    monkeypatch.setattr(systems, 'simulate', fake_simulate)
    _t, ref, new, rows, missing = systems.equivalence(s, str(tmp_path / 'work'))
    assert calls == [False]                          # only the library model is generated
    assert [(r['reference'], r['model'], r['ok']) for r in rows] == [('a/x', 'b/x', True)]
    assert missing == []
