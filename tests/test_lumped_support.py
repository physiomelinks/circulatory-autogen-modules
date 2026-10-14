"""
cam_testing's support for supermodule versions that stand in for a monolithic version (the lumped
vessels, modules/README.md "Lumped vessels"): available_to_phlynx in the manifests, aliased
shared_parameters, and the supermodule structure check of routes and template slots.
"""
import json
import os

import pytest

from cam_testing import harness, library, manifests, paths, supermodule


def _version(modules, rel, version, entry=None, rows=()):
    '''modules/<rel>/versions/<version>/ with a one-entry config (extra keys from ``entry``) and a
    default parameterisation with ``rows``.'''
    name = os.path.basename(rel)
    d = os.path.join(modules, rel, 'versions', version)
    os.makedirs(os.path.join(d, 'parameterisations', 'default'), exist_ok=True)
    config = {'module_type': name, 'module_subtype': version, 'component_file': f'{name}_{version}_modules.cellml',
              'component_type': f'{name}_{version}', 'default_parameterisation': 'default',
              'variables_and_units': [], **(entry or {})}
    with open(os.path.join(d, f'{name}_{version}_modules_config.json'), 'w') as f:
        json.dump([config], f)
    with open(os.path.join(d, f'{name}_{version}_modules.cellml'), 'w') as f:
        f.write('<model/>')
    with open(os.path.join(d, 'parameterisations', 'default', 'default_parameters.csv'), 'w') as f:
        f.write('variable_name,units,value,data_reference,sourced\n')
        f.writelines(','.join(map(str, r)) + '\n' for r in rows)
    return d


@pytest.fixture
def repo(tmp_path):
    '''A repo with a monolithic vessel hidden from PhLynx, its lumped twin (C -> R, sharing "l" and
    the monolithic name "C_T" for the compliance's C) and an empty template.'''
    modules = str(tmp_path / 'modules')
    sub_vars = [['C', 'm6_per_J', 'access', 'constant'], ['l', 'metre', 'access', 'constant'],
                ['u', 'J_per_m3', 'access', 'variable']]
    _version(modules, 'vessels/vessel', 'vp', {'available_to_phlynx': False})
    _version(modules, 'lumped/compliance', 'linear', {'variables_and_units': sub_vars})
    _version(modules, 'lumped/resistance', 'linear', {'variables_and_units': [['R', 'Js_per_m6', 'access', 'constant'],
                                                                             ['l', 'metre', 'access', 'constant']]})
    subs = [{'name': 'C', 'module_type': 'compliance', 'module_subtype': 'linear', 'parameterisation': 'default',
             'inp_instances': [], 'out_instances': ['R']},
            {'name': 'R', 'module_type': 'resistance', 'module_subtype': 'linear', 'parameterisation': 'default',
             'inp_instances': ['C'], 'out_instances': []}]
    _version(modules, 'vessels/vessel', 'vp_lumped',
             {'module_format': 'supermodule', 'routes': {'inputs': {'vessel_port': 'C'}, 'outputs': {'vessel_port': 'X'}},
              'shared_parameters': ['l', {'name': 'C_T', 'variable': 'C', 'submodules': ['C']}], 'submodules': subs},
             rows=[['l', 'metre', 0.1, 'test', 'no'], ['C_T', 'm6_per_J', 1e-8, 'test', 'no'],
                   ['R_R', 'Js_per_m6', 1e7, 'test', 'no']])
    _version(modules, 'lumped/lumped_vessel', 'vp_empty',
             {'module_format': 'supermodule', 'template': True,
              'submodules': [{'name': 'C', 'module_type': 'compliance', 'choices': ['linear'],
                              'inp_instances': [], 'out_instances': []}]})
    paths.configure(repo_root=str(tmp_path), library_dirs=[], export=False)
    yield tmp_path
    paths.reset()


def test_available_to_phlynx_leaves_a_version_out_of_the_manifests(repo):
    hidden = library.load_version('vessel', 'vp')
    assert not hidden.available_to_phlynx
    assert library.load_version('vessel', 'vp_lumped').available_to_phlynx
    manifests.main([])
    with open(os.path.join(str(repo), 'manifests', 'all.json')) as f:
        names = [c['name'] for c in json.load(f)['collections']['configs']]
    assert 'vessel (vp_lumped)' in names and 'vessel (vp)' not in names


def test_shared_parameters_keep_their_names(repo):
    version = library.load_version('vessel', 'vp_lumped')
    assert version.shared_parameter_names == ['l', 'C_T']
    # one model parameter per shared name, <name>_<vessel>; the others are <var>_<vessel>_<submodule>
    assert harness.supermodule_model_name(version, 'C_T', 'aorta') == 'C_T_aorta'
    assert harness.supermodule_model_name(version, 'R_R', 'aorta') == 'R_aorta_R'
    assert harness.split_model_name(version, 'C_T_aorta', 'aorta')[:3] == ('C_T', 'C_T', '')


def test_structure_check_names_a_route_to_no_submodule(repo):
    problems = supermodule.structure_problems(library.load_version('vessel', 'vp_lumped'))
    assert any('routes.outputs.vessel_port: X is not a submodule' in p for p in problems)
    # shared parameters name no submodule but are not reported as undeclared globals
    assert not [p for p in problems if 'C_T' in p or ' l ' in p]
