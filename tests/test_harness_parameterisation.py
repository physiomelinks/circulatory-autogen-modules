'''The version alone (no harness network) is written at the version's default parameterisation, which
need not be called "default" (e.g. valve/pp's aov).'''
import json
from types import SimpleNamespace

from cam_testing import harness


def _version(default):
    return SimpleNamespace(spec={}, BC_type='pp', vessel_type='valve', key='valve/pp', id='valve_pp',
                           default_parameterisation_name=default, parameters=lambda: [])


def _records(tmp_path, version, **kw):
    harness._write_resources(version, str(tmp_path), 'm', {}, parameters=[], **kw)
    with open(tmp_path / 'm_module_array.json') as f:
        return json.load(f)


def test_lone_version_uses_its_default_parameterisation(tmp_path):
    (rec,) = _records(tmp_path, _version('aov'))
    assert rec['name'] == harness.VESSEL and rec['parameterisation'] == 'aov'


def test_lone_version_at_the_parameterisation_under_test(tmp_path):
    (rec,) = _records(tmp_path, _version('aov'), parameterisation='miv')
    assert rec['parameterisation'] == 'miv'
    (rec,) = _records(tmp_path, _version('aov'), parameterisation=SimpleNamespace(name='trv'))
    assert rec['parameterisation'] == 'trv'
