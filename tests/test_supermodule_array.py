"""
A supermodule version's instances are its own module array, <module_type>_<version>_module_array.json
beside its config, and their parameter values its parameters file, <module_type>_<version>_parameters.json:
Version.submodules and Version.instance_parameters read them through libcuflynx, the structure check
reports a missing or inconsistent one, and the directory layout requires them of a supermodule version only.
"""
import json
import os

import pytest

from cam_testing import library, supermodule
from cam_testing.suite import test_structure as structure

pytest.importorskip('libcuflynx')

ENTRY = {'module_type': 'organ', 'module_subtype': 's1', 'module_format': 'supermodule',
         'default_parameterisation': 'default', 'licence': 'CC0-1.0', 'creator': [],
         'submodule_versions': [{'module_type': 'chan', 'module_subtype': 'v1'}]}
RECORDS = [{'name': 'a', 'module_type': 'chan', 'module_subtype': 'v1', 'out_instances': ['b']},
           {'name': 'b', 'module_type': 'chan', 'module_subtype': 'v1', 'inp_instances': ['a']}]


def _parameters(records, **values):
    return {'schema_version': 1, 'instances': {
        r['name']: {'parameterisation': 'default',
                    'parameters': {v: {'value': 1, 'units': 'S'} for v in values.get(r['name'], [])}}
        for r in records}}


def _version(tmp_path, records=None, parameters=None, entry=ENTRY):
    vdir = tmp_path / 'organ' / 'versions' / 's1'
    vdir.mkdir(parents=True)
    (vdir / 'organ_s1_modules_config.json').write_text(json.dumps([entry]))
    if records is not None:
        (vdir / 'organ_s1_module_array.json').write_text(json.dumps(records))
        (vdir / 'organ_s1_parameters.json').write_text(json.dumps(
            _parameters(records) if parameters is None else parameters))
    loaded = library.entry_from_file(entry, str(vdir / 'organ_s1_modules_config.json'))
    return library.Version(library.ModuleType('organ', str(tmp_path / 'organ')), 's1', loaded, {}, entry)


def test_submodules_are_read_from_the_module_array(tmp_path):
    v = _version(tmp_path, RECORDS)
    assert v.module_array_path.endswith(os.path.join('s1', 'organ_s1_module_array.json'))
    assert v.submodule_names == ['a', 'b']
    assert v.submodules[1]['inp_instances'] == ['a']


def test_structure_check_reports_a_missing_array(tmp_path):
    problems = supermodule.structure_problems(_version(tmp_path))
    assert len(problems) == 1 and 'organ_s1_module_array.json' in problems[0]


@pytest.mark.parametrize('records', [
    [dict(RECORDS[0], parameterisation='fit'), RECORDS[1]],                                    # names a parameterisation
    [dict(RECORDS[0], module_subtype='v2', out_instances=[]), dict(RECORDS[1], inp_instances=[])],   # unlisted version
], ids=['parameterisation', 'unlisted_version'])
def test_structure_check_reports_libcuflynx_refusals(tmp_path, records):
    problems = supermodule.structure_problems(_version(tmp_path, records))
    assert len(problems) == 1 and 'organ_s1_module_array.json' in problems[0]


def test_instance_parameters_are_read_from_the_parameters_file(tmp_path):
    v = _version(tmp_path, RECORDS, _parameters(RECORDS, a=['g']))
    assert v.instance_parameters_path.endswith(os.path.join('s1', 'organ_s1_parameters.json'))
    assert v.instance_parameters['a'] == {'parameterisation': 'default',
                                          'parameters': {'g': {'value': 1, 'units': 'S'}}}


def test_structure_check_reports_a_missing_parameters_file(tmp_path):
    v = _version(tmp_path, RECORDS)
    os.remove(v.instance_parameters_path)
    problems = supermodule.parameters_file_problems(v)
    assert len(problems) == 1 and 'organ_s1_parameters.json' in problems[0] and 'does not exist' in problems[0]


def test_structure_check_reports_an_instance_missing_from_the_parameters_file(tmp_path):
    problems = supermodule.parameters_file_problems(_version(tmp_path, RECORDS, _parameters(RECORDS[:1])))
    assert len(problems) == 1 and "it has no entry for ['b']" in problems[0]


def test_structure_check_reports_parameters_that_are_not_the_instance_s_variables(tmp_path):
    """A library valve (pp: rho is one of its globals) and the chamber's parameterisations."""
    entry = dict(ENTRY, submodule_versions=[{'module_type': 'valve', 'module_subtype': 'pp'}])
    records = [{'name': 'aov', 'module_type': 'valve', 'module_subtype': 'pp'}]
    parameters = _parameters(records, aov=['K_vo', 'rho', 'nope', 'x/K_vo'])
    parameters['instances']['aov']['parameterisation'] = 'missing'
    problems = supermodule.parameters_file_problems(_version(tmp_path, records, parameters, entry))
    assert [p.split(': ', 1)[1].split(' (')[0] for p in problems] == [
        'valve/pp has no parameterisation missing',
        "parameter rho is a global of valve/pp; the supermodule's parameterisation sets it",
        'parameter nope is not a variable of valve/pp',
        'parameter x/K_vo names a path, but valve/pp is not a supermodule']


@pytest.mark.parametrize('is_supermodule,files,expected', [
    (True, [], ['missing organ_s1_module_array.json', 'missing organ_s1_parameters.json']),
    (True, ['organ_s1_module_array.json', 'organ_s1_parameters.json'], []),
    (False, ['organ_s1_module_array.json', 'organ_s1_parameters.json', 'organ_s1_modules.cellml',
             'organ_s1_units.cellml'],
     ['unexpected file organ_s1_module_array.json', 'unexpected file organ_s1_parameters.json']),
])
def test_layout_requires_the_array_of_a_supermodule_only(tmp_path, is_supermodule, files, expected):
    for f in ['organ_s1_modules_config.json', 'organ_s1_tests.yaml', 'organ_s1_verification_config.json',
              'organ_s1_references.bib'] + files:
        (tmp_path / f).write_text('')
    problems = []
    structure._check_files('version', str(tmp_path), {'module_type': 'organ', 'version': 's1'}, problems,
                           is_supermodule)
    assert [p.split(': ', 1)[1].split(' (')[0] for p in problems] == expected
