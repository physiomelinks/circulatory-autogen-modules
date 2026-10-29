'''--quick-unreviewed (CI): an unreviewed version runs only the quick tests, which include the supermodule
equivalence test, so a system model converted to an unreviewed supermodule is still checked.'''
import types

from cam_testing import library, pytest_plugin

KEY = ('heart', 'Argus2026_v01')


class Item:
    def __init__(self, name, params):
        self.originalname, self.callspec, self.markers = name, types.SimpleNamespace(params=params), []

    def add_marker(self, marker):
        self.markers.append(marker)


def test_unreviewed_supermodule_runs_its_equivalence_test(monkeypatch):
    monkeypatch.setattr(library, 'load_version', lambda *key: types.SimpleNamespace(key=key, reviewed=False))
    items = [Item('supermodule_equivalence_test', {'equivalence_key': KEY + ('closed_loop_cvs/3compartment',)}),
             Item('supermodule_structure_test', {'supermodule_key': KEY}),
             Item('run_test', {'component_key': KEY}),
             Item('verification_test_BC', {'component_key': KEY})]
    options = {'--quick-unreviewed': True, '--changed-only': False, '--shard-stage': 'vv'}
    pytest_plugin.pytest_collection_modifyitems(types.SimpleNamespace(getoption=options.__getitem__), items)
    assert {i.originalname for i in items if any(m.name == 'skip' for m in i.markers)} == {'verification_test_BC'}
