"""
cam_testing.phlynx.component_job: the parameters PhLynx gets are named as the model names them
(<var>_<instance>, a global bare), since PhLynx matches a row to an instance by that name. The
harness writes the parameters file <vessel>/<var>. No PhLynx needed.
"""
from cam_testing import harness, phlynx
from cam_testing.library import load_version


def test_the_job_names_parameters_as_the_model_does(tmp_path):
    version = load_version('i_M', 'Argus2026_v01')
    job = phlynx.component_job(version, str(tmp_path))
    names = [p['variable_name'] for p in job['parameters']]
    assert names and not [n for n in names if '/' in n]
    own = [p.variable_name for p in version.default_parameterisation.parameters() if not p.is_global]
    assert own and f'{own[0]}_{harness.VESSEL}' in names
