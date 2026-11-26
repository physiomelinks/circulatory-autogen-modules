"""
A system model's parameter sets (cam_testing.system): it is generated from its default set, or the set its spec
names; a model with a calibration workflow has its default set written from the workflow's run
(write_default), and is stale (workflow_stale, system_workflow_current_test) once an input in
parameterisations/uncalibrated/ changes, e.g. when a module recalibration writes into it. No simulation.
"""
import json
import os

import pytest

from cam_testing import system as systems


def _system(tmp_path, workflow=True, spec=None):
    d = tmp_path / 'cat' / 'cell'
    for p in ('default', 'uncalibrated') if workflow else ('default',):
        (d / 'parameterisations' / p).mkdir(parents=True)
        (d / 'parameterisations' / p / f'{p}_parameters.json').write_text('{"schema_version": 1}')
    if workflow:
        u = d / 'parameterisations' / 'uncalibrated'
        (u / 'cell_calibration_workflow.json').write_text('{"name": "fit"}')
        (u / 'protocols').mkdir()
        (u / 'protocols' / 'fig1_obs_data.json').write_text('{}')
        (u / 'SOURCES.md').write_text('# sources\n')
    return systems.System('cell', 'cat', str(d), spec or {})


@pytest.fixture
def written(tmp_path, monkeypatch):
    '''A model whose default set write_default wrote from a run (libcuflynx's writer stubbed: it writes the
    path it is given).'''
    from libcuflynx.calibration_workflow import runner
    calls = []

    def write(workflow, run_dir, path, module_library_dirs=None):
        calls.append((workflow, run_dir, path))
        with open(path, 'w') as f:
            f.write('{"schema_version": 1, "x": {"value": 2, "units": "s", "data_reference": "calibrated"}}')
        return path

    monkeypatch.setattr(runner, 'write_calibrated_parameters_to', write)
    s = _system(tmp_path)
    assert systems.write_default(s, str(tmp_path / 'run')) == s.parameters_path
    return s, calls


def test_a_model_is_generated_from_its_default_set_unless_its_spec_names_another(tmp_path):
    s = _system(tmp_path)
    assert s.parameters_path == os.path.join(s.dir, 'parameterisations', 'default', 'default_parameters.json')
    other = systems.System(s.name, s.category, s.dir, {'parameterisation': 'uncalibrated'})
    assert other.parameters_path == os.path.join(s.dir, 'parameterisations', 'uncalibrated',
                                                 'uncalibrated_parameters.json')


def test_a_model_without_a_workflow_is_never_stale(tmp_path):
    s = _system(tmp_path, workflow=False)
    assert s.workflow_path is None and systems.workflow_stale(s) is None
    assert s.calibration_dir == s.dir


def test_a_default_set_not_written_by_the_workflow_is_stale(tmp_path):
    s = _system(tmp_path)
    assert s.calibration_dir == os.path.join(s.dir, 'parameterisations', 'uncalibrated')
    assert systems.workflow_stale(s) == ['default has no default_calibration.json: it was not written by the workflow']


def test_write_default_writes_the_run_and_records_its_inputs(written):
    s, calls = written
    assert calls == [(s.workflow_path, calls[0][1], s.parameters_path)]
    with open(os.path.join(os.path.dirname(s.parameters_path), systems.CALIBRATION_RECORD)) as f:
        record = json.load(f)
    # every input but notes: the run directory is not recorded (a local path)
    assert sorted(record['inputs']) == ['cell_calibration_workflow.json', 'protocols/fig1_obs_data.json',
                                        'uncalibrated_parameters.json']
    assert systems.workflow_stale(s) == []


def test_a_changed_input_makes_the_workflow_stale(written):
    s, _ = written
    u = os.path.join(s.dir, 'parameterisations', 'uncalibrated')
    with open(os.path.join(u, 'SOURCES.md'), 'a') as f:
        f.write('more notes\n')                                   # notes are not inputs
    assert systems.workflow_stale(s) == []
    with open(os.path.join(u, 'uncalibrated_parameters.json'), 'w') as f:
        f.write('{"schema_version": 1, "x": {"value": 3}}')       # e.g. a module recalibration wrote here
    with open(os.path.join(u, 'protocols', 'fig2_obs_data.json'), 'w') as f:
        f.write('{}')
    os.remove(os.path.join(u, 'protocols', 'fig1_obs_data.json'))
    assert systems.workflow_stale(s) == ['uncalibrated_parameters.json changed', 'protocols/fig2_obs_data.json added',
                                         'protocols/fig1_obs_data.json removed']


def test_write_default_needs_a_workflow(tmp_path):
    with pytest.raises(ValueError, match='has no calibration workflow'):
        systems.write_default(_system(tmp_path, workflow=False), str(tmp_path / 'run'))
