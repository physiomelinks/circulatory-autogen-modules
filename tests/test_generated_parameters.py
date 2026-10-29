"""
harness.supermodule_parameters reads the parameters file libcuflynx writes beside a generated model,
<prefix>_parameters.json (libcuflynx's nested JSON), under the generated model's names.
"""
import os

import pytest

from cam_testing import harness, library

pytest.importorskip('libcuflynx')


@pytest.fixture(scope='module')
def heart(tmp_path_factory):
    version = library.load_version('heart', 'Argus2026_v01')
    return version, harness.generate(version, str(tmp_path_factory.mktemp('heart')))


def test_generated_parameters_file_is_the_nested_json(heart):
    _version, model = heart
    path = harness.generated_parameters_path(model)
    assert path == os.path.join(os.path.dirname(model), f'{os.path.basename(os.path.dirname(model))}_parameters.json')
    assert os.path.isfile(path)


def test_supermodule_parameters_are_read_from_the_generated_json(heart):
    version, model = heart
    params = {p.variable_name: p for p in harness.supermodule_parameters(version, model)}
    template = {p.variable_name: p for p in version.default_parameterisation.parameters()}
    assert params['T'].kind == 'global_constant' and params['T'].model_name == 'T'
    assert float(params['T'].value) == float(template['T'].value)
    assert params['rv/E_A'].model_name == f'E_A_{harness.VESSEL}__rv'
    assert params['rv/E_A'].kind == 'constant'
    assert float(params['rv/E_A'].value) == float(template['rv/E_A'].value)
    assert params['rv/E_A'].units == template['rv/E_A'].units
    assert params['aov/M_rg'].value != '' and float(params['aov/M_rg'].value) == 0.0     # a zero is kept
    assert params['puv/u_out'].kind == 'boundary_condition'
