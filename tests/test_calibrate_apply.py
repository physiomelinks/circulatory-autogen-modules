"""
cam_testing.calibrate apply: a parameterisation's calibrated values reach the version's default
parameterisation, every instance of it in a supermodule's parameters file (whatever parameterisation the
instance started from), the system models using it (<vessel>[/<submodule>]/<var>) and the monolithic
counterparts named by a supermodule's equivalent output_map. A mock library in a temp directory, shaped like
the heart: no simulation.
"""
import json
import os
from types import SimpleNamespace

from libcuflynx.utilities.config_schemas import dump_supermodule_parameters

from cam_testing import calibrate

HEADER5 = 'variable_name,units,value,data_reference,sourced'
HEADER4 = 'variable_name,units,value,data_reference'


def _csv(path, rows, header=HEADER5, eol='\n'):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as f:
        f.write(eol.join([header] + rows) + eol)


def _version(tmp, mt, name, instances=None, spec=None):
    '''A version; a supermodule when ``instances`` ({name: (module_type, parameterisation, {var: value})}) is
    given, with its module array records and parameters file.'''
    d = os.path.join(tmp, 'modules', mt, name)

    def parameterisation(n=None):
        n = n or 'default'
        p = lambda suffix: os.path.join(d, 'parameterisations', n, f'{n}_{suffix}')
        return SimpleNamespace(name=n, parameters_path=p('parameters.csv'), params_for_id_path=p('params_for_id.csv'),
                               calibrated_parameters_path=p('calibrated_parameters.csv'),
                               calibration_path=p('calibration.json'), obs_data_path=p('obs_data.json'))
    submodules, parameters_path = [], os.path.join(d, f'{mt}_{name}_parameters.json')
    if instances is not None:
        os.makedirs(d, exist_ok=True)
        submodules = [{'name': n, 'module_type': t.split('/')[0], 'module_subtype': t.split('/')[1]}
                      for n, (t, _, _) in instances.items()]
        doc = {'schema_version': 1, 'instances': {
            n: {'parameterisation': p, 'parameters': {v: {'value': x, 'units': 'S', 'data_reference': 'its own'}
                                                      for v, x in values.items()}}
            for n, (_, p, values) in instances.items()}}
        with open(parameters_path, 'w') as f:
            f.write(dump_supermodule_parameters(doc))
    return SimpleNamespace(mtype=SimpleNamespace(name=mt), name=name, key=f'{mt}/{name}', is_supermodule=instances is not None,
                           submodules=submodules, spec=spec or {}, default_parameterisation_name='default',
                           default_parameterisation=parameterisation(), parameterisation=parameterisation, dir=d,
                           instance_parameters_path=parameters_path)


def _library(tmp):
    chamber = _version(tmp, 'chamber', 'v1')
    # ra starts from the chamber's parameterisation "ra" with its own E_A; lv from the chamber's default
    heart = _version(tmp, 'heart', 's1', instances={'ra': ('chamber/v1', 'ra', {'E_A': 3, 'k': 7}),
                                                    'lv': ('chamber/v1', 'default', {'E_A': 1, 'k': 7})},
                     spec={'supermodule': {'equivalent': [{'model': 'sys_super', 'reproduces': 'sys_mono',
                                                           'output_map': {'mono_C/x': 'C/lv/x', 'mono_C/y': 'C/other/y'}}]}})
    # a supermodule over the heart: its parameters file sets lv's E_A through the path, and nothing of ra
    body = _version(tmp, 'body', 'o1', instances={'heart': ('heart/s1', 'default', {'lv/E_A': 1})})
    single = _version(tmp, 'heart', 'single', instances={'chamber': ('chamber/v1', 'fast', {'E_A': 5})})
    mono = _version(tmp, 'heartmono', 'm1')
    fit = chamber.parameterisation('fit')
    _csv(fit.params_for_id_path, ['mod,E_A,0.1,10,E_A,uniform'], header='vessel_name,param_name,min,max,name_for_plotting,prior')
    _csv(fit.calibrated_parameters_path, ['E_A,S,2.5,calibrated,no', 'k,S,7,not fitted,no'])
    with open(fit.calibration_path, 'w') as f:
        json.dump({'calibrated_parameters': {'E_A': 2.5000001}}, f)
    with open(fit.obs_data_path, 'w') as f:
        json.dump({'obs_data_name': 'fit', 'calibration_note': 'Smith2020 Fig. 1', 'reference_key': 'Smith2020',
                   'protocol_info': {}, 'data_items': []}, f)
    _csv(chamber.default_parameterisation.parameters_path, ['E_A,S,1,old ref,no', 'k,S,7,keep,no'], eol='\r\n')
    _csv(heart.default_parameterisation.parameters_path, ['T,s,1,its global,no'])
    _csv(mono.default_parameterisation.parameters_path, ['E_A,S,1,old,no'])
    systems = []
    for key, vessels, rows in (
            ('sys_flat', [{'name': 'C_ch', 'module_type': 'chamber', 'module_subtype': 'v1', 'parameterisation': 'default'},
                          {'name': 'D_ch', 'module_type': 'chamber', 'module_subtype': 'v1', 'parameterisation': 'fit'}],
             ['C_ch/E_A,S,1,old', 'D_ch/E_A,S,1,old']),
            ('sys_super', [{'name': 'C', 'module_type': 'heart', 'module_subtype': 's1', 'parameterisation': 'default'},
                           {'name': 'E', 'module_type': 'heart', 'module_subtype': 's1', 'parameterisation': 'default'}],
             ['C/lv/E_A,S,1,old', 'E/ra/E_A,S,4,its own']),
            ('sys_mono', [{'name': 'mono_C', 'module_type': 'heartmono', 'module_subtype': 'm1', 'parameterisation': 'default'}],
             ['mono_C/E_A,S,1,old'])):
        path = os.path.join(tmp, 'system_models', key, f'{key}_parameters.csv')
        _csv(path, rows, header=HEADER4)
        systems.append({'key': key, 'dir': os.path.dirname(path), 'vessels': vessels, 'parameters_path': path})
    return chamber, [chamber, heart, body, single, mono], systems


def _instances(version):
    with open(version.instance_parameters_path) as f:
        return json.load(f)['instances']


def test_a_chamber_calibration_updates_every_chamber_instance_in_the_heart_s_parameters_file(tmp_path):
    chamber, versions, systems = _library(str(tmp_path))
    heart, body, single = versions[1:4]
    calibrate.write_plan(calibrate.apply_plan(chamber, 'fit', versions=versions, systems=systems))
    instances = _instances(heart)
    # lv, at the chamber's default
    assert instances['lv']['parameters']['E_A']['value'] == 2.5000001
    # ra too, though it started from the chamber's parameterisation "ra" with its own value; it still says so
    assert instances['ra']['parameters']['E_A']['value'] == 2.5000001
    assert instances['ra']['parameterisation'] == 'ra'
    assert instances['ra']['parameters']['k'] == {'value': 7, 'units': 'S', 'data_reference': 'its own'}
    assert instances['lv']['parameters']['E_A']['data_reference'].startswith(
        'Smith2020; Calibrated (libcuflynx) to chamber v1 parameterisations/fit/fit_obs_data.json: Smith2020 Fig. 1 (was 1;')
    assert instances['lv']['parameters']['E_A']['sourced'] == 'yes'
    # a supermodule over the heart: the path it sets; nothing added for ra
    assert _instances(body)['heart']['parameters'] == {'lv/E_A': {
        'value': 2.5000001, 'units': 'S', 'sourced': 'yes',
        'data_reference': _instances(body)['heart']['parameters']['lv/E_A']['data_reference']}}
    assert _instances(single)['chamber']['parameters']['E_A']['value'] == 2.5000001
    # the file keeps libcuflynx's layout
    with open(heart.instance_parameters_path) as f:
        assert f.read() == dump_supermodule_parameters({'schema_version': 1, 'instances': instances})


def test_apply_plan_and_write(tmp_path):
    chamber, versions, systems = _library(str(tmp_path))
    plan = calibrate.apply_plan(chamber, 'fit', versions=versions, systems=systems)
    got = {(os.path.relpath(c['path'], tmp_path), c['row']) for c in plan}
    assert got == {
        ('modules/chamber/v1/parameterisations/default/default_parameters.csv', 'E_A'),
        ('modules/heart/s1/heart_s1_parameters.json', 'ra/E_A'),
        ('modules/heart/s1/heart_s1_parameters.json', 'lv/E_A'),
        ('modules/body/o1/body_o1_parameters.json', 'heart/lv/E_A'),
        ('modules/heart/single/heart_single_parameters.json', 'chamber/E_A'),
        ('system_models/sys_flat/sys_flat_parameters.csv', 'C_ch/E_A'),
        ('system_models/sys_flat/sys_flat_parameters.csv', 'D_ch/E_A'),
        ('system_models/sys_super/sys_super_parameters.csv', 'C/lv/E_A'),
        ('system_models/sys_super/sys_super_parameters.csv', 'E/ra/E_A'),
        ('system_models/sys_mono/sys_mono_parameters.csv', 'mono_C/E_A'),
        ('modules/heartmono/m1/parameterisations/default/default_parameters.csv', 'E_A'),
    }
    assert all(c['value'] == repr(2.5000001) for c in plan)          # full precision from calibration.json
    own = next(c for c in plan if c['path'] == chamber.default_parameterisation.parameters_path)
    assert own['reference'].startswith('Smith2020; Calibrated (libcuflynx) to parameterisations/fit/fit_obs_data.json: Smith2020 Fig. 1 (was 1;')
    ra = next(c for c in plan if c['row'] == 'ra/E_A')
    assert ra['old'] == 3 and 'to chamber v1 parameterisations/fit/fit_obs_data.json' in ra['reference']

    before = {c['path']: open(c['path'], newline='').read() for c in plan}
    calibrate.write_plan([])                                             # nothing to do: nothing written
    assert all(open(p, newline='').read() == t for p, t in before.items())
    calibrate.write_plan(plan)
    text = open(chamber.default_parameterisation.parameters_path, newline='').read()
    assert '\r\n' in text and text.count('\r\n') == 3                    # line endings and other rows kept
    assert 'k,S,7,keep,no' in text
    rows = calibrate._row_values(chamber.default_parameterisation.parameters_path, {'E_A'})
    assert rows['E_A']['value'] == repr(2.5000001) and rows['E_A']['sourced'] == 'yes'
    # every vessel of the version, whatever its parameterisation, and every row through the heart
    flat = calibrate._row_values(systems[0]['parameters_path'], {'C_ch/E_A', 'D_ch/E_A'})
    assert {r['value'] for r in flat.values()} == {repr(2.5000001)}
    sup = calibrate._row_values(systems[1]['parameters_path'], {'C/lv/E_A', 'E/ra/E_A'})
    assert {r['value'] for r in sup.values()} == {repr(2.5000001)}
    # the heart's globals are not the chamber's to calibrate
    assert calibrate._row_values(versions[1].default_parameterisation.parameters_path, {'T'})['T']['value'] == '1'


def test_apply_dry_run_writes_nothing(tmp_path, monkeypatch):
    from cam_testing import library
    chamber, versions, systems = _library(str(tmp_path))
    monkeypatch.setattr(library, 'version_by_key', lambda key: chamber)
    before = {p: open(p, newline='').read() for p in (chamber.default_parameterisation.parameters_path,
                                                      versions[1].instance_parameters_path)}
    lines = []
    plan = calibrate.apply_command('chamber/v1', 'fit', dry_run=True, versions=versions, systems=systems, out=lines.append)
    assert len(plan) == 11 and lines[-1].startswith('11 row(s)') and 'dry run' in lines[-1]
    assert all(open(p, newline='').read() == t for p, t in before.items())
