# The suite lives in the package (cam_testing/suite/test_structure.py), so other repos can run it too.
from cam_testing.suite.test_structure import *  # noqa: F401,F403
from cam_testing.suite import test_structure as structure


def _model(tmp_path, files):
    mdir = tmp_path / 'cell'
    for f in ('cell_module_array.json', 'cell_system.yaml', *files):
        (mdir / f).parent.mkdir(parents=True, exist_ok=True)
        (mdir / f).write_text('{}')
    return mdir


def _problems(mdir):
    problems = []
    structure._check_system_model(str(mdir), problems)
    return [p.split(': ', 1)[1].split(' (')[0] for p in problems]


def test_system_model_holds_its_parameter_sets_in_parameterisations(tmp_path):
    """A system model's values are its default set, parameterisations/default/default_parameters.json; a
    parameters file beside the module array, or no default set, is refused."""
    mdir = _model(tmp_path, ['parameterisations/default/default_parameters.json', 'reference/cell.cellml'])
    assert _problems(mdir) == []
    (mdir / 'cell_parameters.json').write_text('{}')
    (mdir / 'parameterisations' / 'default' / 'default_parameters.json').unlink()
    (mdir / 'parameterisations' / 'other').mkdir()
    assert _problems(mdir) == ['unexpected file cell_parameters.json', 'missing default_parameters.json',
                               'unexpected directory other']


def test_system_model_holds_its_calibration_workflow_in_its_uncalibrated_set(tmp_path):
    """A model with a calibration workflow holds it in parameterisations/uncalibrated/, with the values it starts
    from, its obs_data, protocols/ and SOURCES.md; its default set may hold the calibration's record. Workflow
    files beside the module array, or an uncalibrated set without a workflow, are refused."""
    mdir = _model(tmp_path, ['parameterisations/default/default_parameters.json',
                             'parameterisations/default/default_calibration.json',
                             'parameterisations/uncalibrated/uncalibrated_parameters.json',
                             'parameterisations/uncalibrated/cell_calibration_workflow.json',
                             'parameterisations/uncalibrated/cell_obs_data.json',
                             'parameterisations/uncalibrated/protocols/fig1_obs_data.json',
                             'parameterisations/uncalibrated/SOURCES.md',
                             'parameterisations/uncalibrated/notes.md'])
    assert _problems(mdir) == []
    (mdir / 'cell_calibration_workflow.json').write_text('{}')
    (mdir / 'protocols').mkdir()
    (mdir / 'parameterisations' / 'uncalibrated' / 'cell_calibration_workflow.json').unlink()
    (mdir / 'parameterisations' / 'uncalibrated' / 'workflows').mkdir()
    assert _problems(mdir) == ['unexpected file cell_calibration_workflow.json', 'missing cell_calibration_workflow.json',
                               'unexpected directory workflows', 'unexpected directory protocols']
