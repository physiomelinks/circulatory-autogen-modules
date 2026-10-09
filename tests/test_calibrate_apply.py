"""
cam_testing.calibrate apply: a parameterisation's calibrated values reach the version's default
parameterisation, its supermodules (<var>_<submodule>), the system models using it (<var>_<vessel>[_<submodule>]) and
the monolithic counterparts named by a supermodule's equivalent output_map. A mock library in a
temp directory: no simulation.
"""
import json
import os
from types import SimpleNamespace

from cam_testing import calibrate

HEADER5 = 'variable_name,units,value,data_reference,sourced'
HEADER4 = 'variable_name,units,value,data_reference'


def _csv(path, rows, header=HEADER5, eol='\n'):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as f:
        f.write(eol.join([header] + rows) + eol)


def _version(tmp, mt, name, submodules=None, spec=None):
    d = os.path.join(tmp, 'modules', mt, name)

    def parameterisation(n=None):
        n = n or 'default'
        p = lambda suffix: os.path.join(d, 'parameterisations', n, f'{n}_{suffix}')
        return SimpleNamespace(name=n, parameters_path=p('parameters.csv'), params_for_id_path=p('params_for_id.csv'),
                               calibrated_parameters_path=p('calibrated_parameters.csv'),
                               calibration_path=p('calibration.json'), obs_data_path=p('obs_data.json'))
    return SimpleNamespace(mtype=SimpleNamespace(name=mt), name=name, key=f'{mt}/{name}', is_supermodule=submodules is not None,
                           submodules=submodules or [], spec=spec or {}, default_parameterisation_name='default',
                           default_parameterisation=parameterisation(), parameterisation=parameterisation, dir=d)


def _library(tmp):
    chan = _version(tmp, 'chan', 'v1')
    cell = _version(tmp, 'cell', 's1', submodules=[{'name': 'chan', 'module_type': 'chan', 'module_subtype': 'v1',
                                                     'parameterisation': 'default'}],
                    spec={'supermodule': {'equivalent': [{'model': 'sys_super', 'reproduces': 'sys_mono',
                                                          'output_map': {'mono_C/x': 'C_chan/x', 'mono_C/y': 'C_other/y'}}]}})
    organ = _version(tmp, 'organ', 'o1', submodules=[{'name': 'cell', 'module_type': 'cell', 'module_subtype': 's1',
                                                       'parameterisation': 'default'}])
    other = _version(tmp, 'cell', 'other', submodules=[{'name': 'chan', 'module_type': 'chan', 'module_subtype': 'v1',
                                                         'parameterisation': 'fit'}])      # not the default parameterisation: untouched
    mono = _version(tmp, 'cellmono', 'm1')
    fit = chan.parameterisation('fit')
    _csv(fit.params_for_id_path, ['mod,g,0.1,10,g,uniform'], header='vessel_name,param_name,min,max,name_for_plotting,prior')
    _csv(fit.calibrated_parameters_path, ['g,S,2.5,calibrated,no', 'k,S,7,not fitted,no'])
    with open(fit.calibration_path, 'w') as f:
        json.dump({'calibrated_parameters': {'g': 2.5000001}}, f)
    with open(fit.obs_data_path, 'w') as f:
        json.dump({'obs_data_name': 'fit', 'calibration_note': 'Smith2020 Fig. 1', 'reference_key': 'Smith2020',
                   'protocol_info': {}, 'data_items': []}, f)
    _csv(chan.default_parameterisation.parameters_path, ['g,S,1,old ref,no', 'k,S,7,keep,no'], eol='\r\n')
    _csv(cell.default_parameterisation.parameters_path, ['g_chan,S,1,old,no', 'k_chan,S,7,keep,no'])
    _csv(organ.default_parameterisation.parameters_path, ['g_cell_chan,S,1,old,no'])
    _csv(other.default_parameterisation.parameters_path, ['g_chan,S,1,old,no'])
    _csv(mono.default_parameterisation.parameters_path, ['g,S,1,old,no'])
    systems = []
    for key, vessels, rows in (
            ('sys_flat', [{'name': 'C_chan', 'module_type': 'chan', 'module_subtype': 'v1', 'parameterisation': 'default'},
                          {'name': 'D_chan', 'module_type': 'chan', 'module_subtype': 'v1', 'parameterisation': 'fit'}],
             ['g_C_chan,S,1,old', 'g_D_chan,S,1,old']),
            ('sys_super', [{'name': 'C', 'module_type': 'cell', 'module_subtype': 's1', 'parameterisation': 'default'}],
             ['g_C_chan,S,1,old']),
            ('sys_mono', [{'name': 'mono_C', 'module_type': 'cellmono', 'module_subtype': 'm1', 'parameterisation': 'default'}],
             ['g_mono_C,S,1,old'])):
        path = os.path.join(tmp, 'system_models', key, f'{key}_parameters.csv')
        _csv(path, rows, header=HEADER4)
        systems.append({'key': key, 'dir': os.path.dirname(path), 'vessels': vessels, 'parameters_path': path})
    return chan, [chan, cell, organ, other, mono], systems


def test_apply_plan_and_write(tmp_path):
    chan, versions, systems = _library(str(tmp_path))
    plan = calibrate.apply_plan(chan, 'fit', versions=versions, systems=systems)
    got = {(os.path.relpath(c['path'], tmp_path), c['row']) for c in plan}
    assert got == {
        ('modules/chan/v1/parameterisations/default/default_parameters.csv', 'g'),
        ('modules/cell/s1/parameterisations/default/default_parameters.csv', 'g_chan'),
        ('modules/organ/o1/parameterisations/default/default_parameters.csv', 'g_cell_chan'),
        ('system_models/sys_flat/sys_flat_parameters.csv', 'g_C_chan'),
        ('system_models/sys_super/sys_super_parameters.csv', 'g_C_chan'),
        ('system_models/sys_mono/sys_mono_parameters.csv', 'g_mono_C'),
        ('modules/cellmono/m1/parameterisations/default/default_parameters.csv', 'g'),
    }
    assert all(c['value'] == repr(2.5000001) for c in plan)          # full precision from calibration.json
    own = next(c for c in plan if c['path'] == chan.default_parameterisation.parameters_path)
    assert own['reference'].startswith('Smith2020; Calibrated (libcuflynx) to parameterisations/fit/fit_obs_data.json: Smith2020 Fig. 1 (was 1;')
    sup = next(c for c in plan if c['row'] == 'g_chan')
    assert 'to chan v1 parameterisations/fit/fit_obs_data.json' in sup['reference']

    before = {c['path']: open(c['path'], newline='').read() for c in plan}
    calibrate.write_plan([])                                             # nothing to do: nothing written
    assert all(open(p, newline='').read() == t for p, t in before.items())
    calibrate.write_plan(plan)
    text = open(chan.default_parameterisation.parameters_path, newline='').read()
    assert '\r\n' in text and text.count('\r\n') == 3                    # line endings and other rows kept
    assert 'k,S,7,keep,no' in text
    rows = calibrate._row_values(chan.default_parameterisation.parameters_path, {'g'})
    assert rows['g']['value'] == repr(2.5000001) and rows['g']['sourced'] == 'yes'
    flat = calibrate._row_values(systems[0]['parameters_path'], {'g_C_chan', 'g_D_chan'})
    assert flat['g_C_chan']['value'] == repr(2.5000001) and flat['g_D_chan']['value'] == '1'   # the vessel at the fit parameterisation untouched
    other = calibrate._row_values(versions[3].default_parameterisation.parameters_path, {'g_chan'})
    assert other['g_chan']['value'] == '1'


def test_apply_dry_run_writes_nothing(tmp_path, monkeypatch):
    from cam_testing import library
    chan, versions, systems = _library(str(tmp_path))
    monkeypatch.setattr(library, 'version_by_key', lambda key: chan)
    before = open(chan.default_parameterisation.parameters_path, newline='').read()
    lines = []
    plan = calibrate.apply_command('chan/v1', 'fit', dry_run=True, versions=versions, systems=systems, out=lines.append)
    assert len(plan) == 7 and lines[-1].startswith('7 row(s)') and 'dry run' in lines[-1]
    assert open(chan.default_parameterisation.parameters_path, newline='').read() == before
