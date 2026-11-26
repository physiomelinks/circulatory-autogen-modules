# The suite lives in the package (cam_testing/suite/test_structure.py), so other repos can run it too.
from cam_testing.suite.test_structure import *  # noqa: F401,F403
from cam_testing.suite import test_structure as structure


def test_system_model_may_hold_its_calibration_workflow(tmp_path):
    """A system model's calibration workflow and its protocols/ pass the layout check; other names do not."""
    mdir = tmp_path / 'cell'
    (mdir / 'protocols').mkdir(parents=True)
    (mdir / 'reference').mkdir()
    for f in ('cell_module_array.json', 'cell_system.yaml', 'cell_calibration_workflow.json',
              'protocols/fig1_obs_data.json'):
        (mdir / f).write_text('{}')
    problems = []
    structure._check_system_model(str(mdir), problems)
    assert problems == []
    (mdir / 'workflows').mkdir()
    (mdir / 'cell_workflow.json').write_text('{}')
    structure._check_system_model(str(mdir), problems)
    assert [p.split(': ', 1)[1].split(' (')[0] for p in problems] == ['unexpected file cell_workflow.json',
                                                                       'unexpected directory workflows']
