"""
cam_testing.calibrate apply: a parameterisation's calibrated values reach the version's default
parameterisation, every instance of it in a supermodule version's template (its default parameterisation,
nested by instance) and in each system model's parameters file (<vessel>[/<submodule>]/<var>, whatever the
instance's values started from), and the monolithic counterparts named by a supermodule's equivalent
output_map. All of them are libcuflynx's nested JSON. A mock library in a temp directory, shaped like the
heart: no simulation.
"""
import json
import os
from types import SimpleNamespace

from libcuflynx.utilities import parameters_file as pf
from libcuflynx.utilities import params_for_id_file as pfi

from cam_testing import calibrate


def _leaf(value, reference='its own', sourced=None):
    leaf = {'value': value, 'units': 'S', 'data_reference': reference}
    if sourced is not None:
        leaf['sourced'] = sourced
    return leaf


def _write(path, body):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pf.write_parameters_file(path, dict(pf.new_document(), **body))


def _read(path, kind=pf.MODULE):
    return pf.read_parameters_file(path, kind=kind)


def _version(tmp, mt, name, instances=None, spec=None):
    '''A version; a supermodule when ``instances`` ({name: (module_type, {path: value})}) is given, with its
    module array records and its template (default parameterisation) holding those values.'''
    d = os.path.join(tmp, 'modules', mt, name)

    def parameterisation(n=None):
        n = n or 'default'
        p = lambda suffix: os.path.join(d, 'parameterisations', n, f'{n}_{suffix}')
        return SimpleNamespace(name=n, parameters_path=p('parameters.json'), params_for_id_path=p('params_for_id.json'),
                               calibrated_parameters_path=p('calibrated_parameters.json'),
                               calibration_path=p('calibration.json'), obs_data_path=p('obs_data.json'))
    default = parameterisation()
    submodules = []
    if instances is not None:
        submodules = [{'name': n, 'module_type': t.split('/')[0], 'module_subtype': t.split('/')[1]}
                      for n, (t, _) in instances.items()]
        body = {'global': {'T': _leaf(1, 'its global', 'no')}}
        for n, (_, values) in instances.items():
            node = body.setdefault(n, {})
            for path, x in values.items():
                pf.set_leaf(node, path.split('/'), _leaf(x))
        _write(default.parameters_path, body)
    return SimpleNamespace(mtype=SimpleNamespace(name=mt), name=name, key=f'{mt}/{name}', is_supermodule=instances is not None,
                           submodules=submodules, spec=spec or {}, default_parameterisation_name='default',
                           default_parameterisation=default, parameterisation=parameterisation, dir=d,
                           instance_parameters_path=default.parameters_path)


def _library(tmp):
    chamber = _version(tmp, 'chamber', 'v1')
    # ra has its own E_A (from the chamber's parameterisation "ra"); lv the chamber's default
    heart = _version(tmp, 'heart', 's1', instances={'ra': ('chamber/v1', {'E_A': 3, 'k': 7}),
                                                    'lv': ('chamber/v1', {'E_A': 1, 'k': 7})},
                     spec={'supermodule': {'equivalent': [{'model': 'sys_super', 'reproduces': 'sys_mono',
                                                           'output_map': {'mono_C/x': 'C/lv/x', 'mono_C/y': 'C/other/y'}}]}})
    # a supermodule over the heart: its template sets lv's E_A through the path, and nothing of ra
    body = _version(tmp, 'body', 'o1', instances={'heart': ('heart/s1', {'lv/E_A': 1})})
    single = _version(tmp, 'heart', 'single', instances={'chamber': ('chamber/v1', {'E_A': 5})})
    mono = _version(tmp, 'heartmono', 'm1')
    fit = chamber.parameterisation('fit')
    pfi.write_params_for_id(fit.params_for_id_path, [{'targets': ['mod/E_A'], 'min': 0.1, 'max': 10,
                                                      'name_for_plotting': 'E_A', 'prior': 'uniform'}])
    _write(fit.calibrated_parameters_path, {'E_A': _leaf(2.5, 'calibrated', 'no'), 'k': _leaf(7, 'not fitted', 'no')})
    with open(fit.calibration_path, 'w') as f:
        json.dump({'calibrated_parameters': {'E_A': 2.5000001}}, f)
    with open(fit.obs_data_path, 'w') as f:
        json.dump({'obs_data_name': 'fit', 'calibration_note': 'Smith2020 Fig. 1', 'reference_key': 'Smith2020',
                   'protocol_info': {}, 'data_items': []}, f)
    _write(chamber.default_parameterisation.parameters_path, {'E_A': _leaf(1, 'old ref', 'no'), 'k': _leaf(7, 'keep', 'no')})
    _write(mono.default_parameterisation.parameters_path, {'E_A': _leaf(1, 'old', 'no')})
    systems = []
    # each system model's parameters file holds its supermodule instances' values, copied from their versions:
    # E's ra has a value of its own; the body's B holds the heart's values nested
    for key, vessels, body_values in (
            ('sys_flat', [{'name': 'C_ch', 'module_type': 'chamber', 'module_subtype': 'v1', 'parameterisation': 'default'},
                          {'name': 'D_ch', 'module_type': 'chamber', 'module_subtype': 'v1', 'parameterisation': 'fit'}],
             {'C_ch/E_A': 1, 'D_ch/E_A': 1}),
            ('sys_super', [{'name': 'C', 'module_type': 'heart', 'module_subtype': 's1', 'parameterisation': 'default'},
                           {'name': 'E', 'module_type': 'heart', 'module_subtype': 's1', 'parameterisation': 'default'}],
             {'C/ra/E_A': 3, 'C/ra/k': 7, 'C/lv/E_A': 1, 'C/lv/k': 7,
              'E/ra/E_A': 4, 'E/ra/k': 7, 'E/lv/E_A': 1, 'E/lv/k': 7}),
            ('sys_mono', [{'name': 'mono_C', 'module_type': 'heartmono', 'module_subtype': 'm1', 'parameterisation': 'default'}],
             {'mono_C/E_A': 1}),
            ('sys_body', [{'name': 'B', 'module_type': 'body', 'module_subtype': 'o1', 'parameterisation': 'default'}],
             {'B/heart/ra/E_A': 3, 'B/heart/ra/k': 7, 'B/heart/lv/E_A': 1, 'B/heart/lv/k': 7})):
        _system_model(tmp, key, vessels, body_values)
        systems.append(key)
    # in the order above (_system_models sorts by key)
    found = {m['key']: m for m in calibrate._system_models(os.path.join(tmp, 'system_models'))}
    return chamber, [chamber, heart, body, single, mono], [found[k] for k in systems]


def _system_model(tmp, key, vessels, values, parameterisation='default'):
    '''system_models/<key>/: its module array and the parameter set ``parameterisation`` holding ``values``
    ({path: value, or a leaf}).'''
    d = os.path.join(tmp, 'system_models', key)
    path = os.path.join(d, 'parameterisations', parameterisation, f'{parameterisation}_parameters.json')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(os.path.join(d, f'{os.path.basename(key)}_module_array.json'), 'w') as f:
        json.dump(vessels, f)
    doc = pf.new_document()
    for name, x in values.items():
        pf.set_leaf(doc, name.split('/'), x if isinstance(x, dict) else _leaf(x, 'model'))
    pf.write_parameters_file(path, doc)
    return path


def test_a_chamber_calibration_updates_every_chamber_instance_in_the_heart_s_template(tmp_path):
    chamber, versions, systems = _library(str(tmp_path))
    heart, body, single = versions[1:4]
    calibrate.write_plan(calibrate.apply_plan(chamber, 'fit', versions=versions, systems=systems))
    doc = _read(heart.instance_parameters_path)
    # lv, at the chamber's default, and ra, though it has a value of its own
    assert doc['lv']['E_A']['value'] == doc['ra']['E_A']['value'] == '2.5000001'
    assert doc['ra']['k'] == {'value': '7', 'units': 'S', 'data_reference': 'its own'}
    assert doc['lv']['E_A']['data_reference'].startswith(
        'Smith2020; Calibrated (libcuflynx) to chamber v1 parameterisations/fit/fit_obs_data.json: Smith2020 Fig. 1 (was 1;')
    assert doc['global']['T']['value'] == '1'                  # the heart's globals are not the chamber's
    # a supermodule over the heart: the path its template sets; nothing added for ra
    assert set(_read(body.instance_parameters_path)['heart']) == {'lv'}
    assert _read(body.instance_parameters_path)['heart']['lv']['E_A']['value'] == '2.5000001'
    assert _read(single.instance_parameters_path)['chamber']['E_A']['value'] == '2.5000001'


def test_a_chamber_calibration_updates_the_system_models_supermodule_instances(tmp_path):
    """Each system model's parameters file holds its supermodule instances' values; every instance of the
    version in it is written, whatever its value was, through a nested supermodule's paths too."""
    chamber, versions, systems = _library(str(tmp_path))
    calibrate.write_plan(calibrate.apply_plan(chamber, 'fit', versions=versions, systems=systems))
    sup = _read(systems[1]['parameters_path'], pf.MODEL)
    for name in ('C', 'E'):
        assert sup[name]['lv']['E_A']['value'] == sup[name]['ra']['E_A']['value'] == '2.5000001'
        assert sup[name]['ra']['k'] == {'value': '7', 'units': 'S', 'data_reference': 'model'}
        assert sup[name]['lv']['E_A']['data_reference'].startswith(
            'Smith2020; Calibrated (libcuflynx) to chamber v1 parameterisations/fit/fit_obs_data.json')
    body = _read(systems[3]['parameters_path'], pf.MODEL)['B']['heart']
    assert body['lv']['E_A']['value'] == body['ra']['E_A']['value'] == '2.5000001'
    assert body['lv']['k']['value'] == '7'


def test_apply_plan_and_write(tmp_path):
    chamber, versions, systems = _library(str(tmp_path))
    plan = calibrate.apply_plan(chamber, 'fit', versions=versions, systems=systems)
    got = {(os.path.relpath(c['path'], tmp_path), c['row']) for c in plan}
    assert got == {
        ('modules/chamber/v1/parameterisations/default/default_parameters.json', 'E_A'),
        ('modules/heart/s1/parameterisations/default/default_parameters.json', 'ra/E_A'),
        ('modules/heart/s1/parameterisations/default/default_parameters.json', 'lv/E_A'),
        ('modules/body/o1/parameterisations/default/default_parameters.json', 'heart/lv/E_A'),
        ('modules/heart/single/parameterisations/default/default_parameters.json', 'chamber/E_A'),
        ('system_models/sys_flat/parameterisations/default/default_parameters.json', 'C_ch/E_A'),
        ('system_models/sys_flat/parameterisations/default/default_parameters.json', 'D_ch/E_A'),
        ('system_models/sys_super/parameterisations/default/default_parameters.json', 'C/lv/E_A'),
        ('system_models/sys_super/parameterisations/default/default_parameters.json', 'C/ra/E_A'),
        ('system_models/sys_super/parameterisations/default/default_parameters.json', 'E/lv/E_A'),
        ('system_models/sys_super/parameterisations/default/default_parameters.json', 'E/ra/E_A'),
        ('system_models/sys_mono/parameterisations/default/default_parameters.json', 'mono_C/E_A'),
        ('modules/heartmono/m1/parameterisations/default/default_parameters.json', 'E_A'),
        ('system_models/sys_body/parameterisations/default/default_parameters.json', 'B/heart/lv/E_A'),
        ('system_models/sys_body/parameterisations/default/default_parameters.json', 'B/heart/ra/E_A'),
    }
    assert all(c['value'] == repr(2.5000001) for c in plan)          # full precision from calibration.json
    own = next(c for c in plan if c['path'] == chamber.default_parameterisation.parameters_path)
    assert own['reference'].startswith('Smith2020; Calibrated (libcuflynx) to parameterisations/fit/fit_obs_data.json: Smith2020 Fig. 1 (was 1;')
    ra = next(c for c in plan if c['row'] == 'ra/E_A')
    assert ra['old'] == '3' and 'to chamber v1 parameterisations/fit/fit_obs_data.json' in ra['reference']

    before = {c['path']: open(c['path']).read() for c in plan}
    calibrate.write_plan([])                                             # nothing to do: nothing written
    assert all(open(p).read() == t for p, t in before.items())
    calibrate.write_plan(plan)
    doc = _read(chamber.default_parameterisation.parameters_path)
    assert doc['E_A']['value'] == repr(2.5000001) and doc['E_A']['sourced'] == 'yes'
    assert doc['k'] == {'value': '7', 'units': 'S', 'data_reference': 'keep', 'sourced': 'no'}
    flat = _read(systems[0]['parameters_path'], pf.MODEL)
    assert {flat['C_ch']['E_A']['value'], flat['D_ch']['E_A']['value']} == {repr(2.5000001)}
    assert _read(systems[2]['parameters_path'], pf.MODEL)['mono_C']['E_A']['value'] == repr(2.5000001)
    # a value is written as a JSON number
    with open(chamber.default_parameterisation.parameters_path) as f:
        assert json.load(f)['E_A']['value'] == 2.5000001


def test_apply_dry_run_writes_nothing(tmp_path, monkeypatch):
    from cam_testing import library
    chamber, versions, systems = _library(str(tmp_path))
    monkeypatch.setattr(library, 'version_by_key', lambda key: chamber)
    paths = [chamber.default_parameterisation.parameters_path, versions[1].instance_parameters_path] + \
        [m['parameters_path'] for m in systems]
    before = {p: open(p).read() for p in paths}
    lines = []
    plan = calibrate.apply_command('chamber/v1', 'fit', dry_run=True, versions=versions, systems=systems, out=lines.append)
    assert len(plan) == 15 and lines[-1].startswith('15 row(s)') and 'dry run' in lines[-1]
    assert all(open(p).read() == t for p, t in before.items())


def test_a_model_with_a_workflow_gets_the_values_in_its_uncalibrated_set_and_its_workflow_is_stale(tmp_path):
    """A system model with a calibration workflow gets a module recalibration in its uncalibrated set (the
    workflow's input), never in its default set (the workflow's result), and keeps a value that set gives its
    own origin; apply reports the workflow stale. A model without one gets it in its default set."""
    tmp = str(tmp_path)
    chamber, versions, _ = _library(tmp)
    vessels = [{'name': 'C_ch', 'module_type': 'chamber', 'module_subtype': 'v1', 'parameterisation': 'default'},
               {'name': 'D_ch', 'module_type': 'chamber', 'module_subtype': 'v1', 'parameterisation': 'default'}]
    default = _system_model(tmp, 'cat/sys_wf', vessels, {'C_ch/E_A': 9, 'D_ch/E_A': 9})
    own = dict(_leaf(4, 'measured by us'), origin='measured')
    uncalibrated = _system_model(tmp, 'cat/sys_wf', vessels, {'C_ch/E_A': 1, 'D_ch/E_A': own}, 'uncalibrated')
    with open(os.path.join(os.path.dirname(uncalibrated), 'sys_wf_calibration_workflow.json'), 'w') as f:
        f.write('{}')
    systems = calibrate._system_models(os.path.join(tmp, 'system_models'))
    wf = next(m for m in systems if m['key'] == 'cat/sys_wf')
    assert wf['workflow'] and wf['parameters_path'] == uncalibrated
    assert next(m for m in systems if m['key'] == 'sys_flat')['parameters_path'].endswith(
        os.path.join('sys_flat', 'parameterisations', 'default', 'default_parameters.json'))
    before = open(default).read()
    lines = []
    from cam_testing import library
    import pytest
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(library, 'version_by_key', lambda key: chamber)
        calibrate.apply_command('chamber/v1', 'fit', versions=versions, systems=systems, out=lines.append)
    doc = _read(uncalibrated, pf.MODEL)
    assert doc['C_ch']['E_A']['value'] == repr(2.5000001)
    assert doc['D_ch']['E_A'] == dict(own, value='4')                     # its own origin: kept
    assert open(default).read() == before                                # the workflow's result: untouched
    assert _read(next(m for m in systems if m['key'] == 'sys_flat')['parameters_path'], pf.MODEL)[
        'C_ch']['E_A']['value'] == repr(2.5000001)
    assert [x for x in lines if 'stale' in x] == [
        'system model cat/sys_wf: its calibration workflow is stale: rerun it to refresh its default set '
        '(python -m cam_testing.system write-default cat/sys_wf <run dir>)']


def test_an_entry_shared_by_two_variables_sets_and_applies_both(tmp_path):
    """One params_for_id entry may set several targets of different variables (SGLT1's one K_Na for K_Na_o
    and K_Na_i): the calibration's value goes to each, and apply writes each."""
    chamber, versions, systems = _library(str(tmp_path))
    fit = chamber.parameterisation('fit')
    pfi.write_params_for_id(fit.params_for_id_path, [
        {'name': 'mod/E', 'targets': ['mod/E_A', 'mod/k'], 'min': 0.1, 'max': 10, 'prior': 'uniform'},
        {'targets': ['loop/R'], 'min': 0.1, 'max': 10, 'prior': 'uniform'}])
    doc = calibrate.read_params_for_id(fit.params_for_id_path)
    assert calibrate.calibrated_values(doc, [2.5, 0.3]) == {'E_A': 2.5, 'k': 2.5, 'loop/R': 0.3}

    pfi.write_params_for_id(fit.params_for_id_path, doc['params'][:1])
    _write(fit.calibrated_parameters_path, {'E_A': _leaf(2.5, 'calibrated', 'no'), 'k': _leaf(2.5, 'calibrated', 'no')})
    with open(fit.calibration_path, 'w') as f:
        json.dump({'calibrated_parameters': {'E_A': 2.5, 'k': 2.5}}, f)
    plan = calibrate.apply_plan(chamber, 'fit', versions=versions, systems=systems)
    own = {c['row']: c['value'] for c in plan if c['path'] == chamber.default_parameterisation.parameters_path}
    assert own == {'E_A': '2.5', 'k': '2.5'}


def test_an_item_computed_from_other_items_names_no_model_output():
    """SGLT1 parent1992's K0.5 ratios divide two other items (operands [""]): only real operands are outputs."""
    items = [{'operands': ['mod/v', 'loop/u']}, {'operands': [''], 'operation': 'division',
                                                 'operation_kwargs': {'x1': 'a', 'x2': 'b'}}]
    assert calibrate.item_outputs(items) == ['loop/u', 'v']


def test_a_modifier_entry_is_refused(tmp_path):
    import pytest
    path = str(tmp_path / 'p_params_for_id.json')
    with open(path, 'w') as f:
        json.dump({'version': 1, 'params': [{'name': 'scale', 'modifies': ['mod/E_A'], 'function': 'scale',
                                             'min': 0.5, 'max': 2}]}, f)
    with pytest.raises(ValueError, match='modifiers'):
        calibrate.read_params_for_id(path)
