"""
What happens to a global constant that two modules' instances set differently, in libcuflynx
and in PhLynx (its real code, through tools/phlynx_bridge).

A global is one value for the whole model, so two instances disagreeing on one is a conflict.
The library has real ones: T is 295.15 kelvin in the Argus2026 ion channels, 310.15 kelvin in
reversal_potentials (cardiomyocyte_Paci2013_v01), and 1 second (a period, not a temperature) in
the cardiac clocks and hearts.

These tests pin down today's behaviour, so a change on either side shows up here:

* libcuflynx warns (ConflictingGlobalWarning, physiomelinks/circulatory_autogen#535), naming
  each value, its units and the record that set it, and keeps the first value.
* PhLynx has no instances: it gets the parameters file. Given two rows for one global with
  the same units, it keeps the first row without a word. Given the same name in different
  units, its export fails at flatten (libCellML's validator names the units). Given no row,
  it exports a model whose two T's have no value.

    PHLYNX_DIR=../phlynx pytest tests/test_phlynx_globals.py
"""
import copy
import json
import os
import re
import subprocess
import warnings
import zipfile

import pytest

from cam_testing import phlynx
from cam_testing.library import MODULES_DIR, load_version

pytestmark = pytest.mark.phlynx_pipeline

CHANNEL = ('i_M', 'Argus2026_v01')            # T = 295.15 kelvin
CLOCK = ('cardiac_clock', 'Liang2009_v01')    # T = 1 second, the cardiac period


def _rename(job, suffix):
    '''``job`` with every record renamed ``<name><suffix>``, its parameters to match, so two
    versions' test networks can share one model.'''
    names = [i['name'] for i in job['instances']]
    new = {n: n + suffix for n in names}
    for record in job['instances']:
        record['name'] = new[record['name']]
        for key in ('inp_instances', 'out_instances'):
            record[key] = ' '.join(new.get(n, n) for n in record[key].split())
    for row in job['parameters']:
        for name in sorted(names, key=len, reverse=True):
            if row['variable_name'].endswith('_' + name):
                row['variable_name'] = row['variable_name'][:-len(name)] + new[name]
                break
    return job


@pytest.fixture(scope='module')
def jobs(tmp_path_factory):
    work = str(tmp_path_factory.mktemp('phlynx_globals'))
    try:
        phlynx.phlynx_dir()
        phlynx._node()
    except phlynx.PipelineUnavailable as e:
        pytest.skip(f'PhLynx pipeline unavailable: {e}')
    channel = phlynx.component_job(load_version(*CHANNEL), os.path.join(work, 'channel'))
    clock = phlynx.component_job(load_version(*CLOCK), os.path.join(work, 'clock'))
    return work, channel, _rename(clock, 'B'), phlynx.library_files()


def _export(work, label, library, instances, parameters):
    '''PhLynx's report for the model, and the initial values of every variable T it exported.'''
    job_path = os.path.join(work, f'{label}.json')
    omex = os.path.join(work, f'{label}.omex')
    with open(job_path, 'w') as f:
        json.dump(dict(library, instances=instances, parameters=parameters, sim_time=1.0,
                       dt=0.01), f)
    p = subprocess.run(['node', phlynx.PHLYNX_BRIDGE, job_path, omex], capture_output=True,
                       text=True, timeout=600,
                       env=dict(os.environ, PHLYNX_DIR=phlynx.phlynx_dir()))
    if p.returncode == 3:
        pytest.skip(p.stderr.strip()[:300])
    report = json.loads(p.stdout.strip().splitlines()[-1])
    t_vars = []
    if report.get('ok'):
        with zipfile.ZipFile(omex) as z:
            for name in z.namelist():
                if name.endswith('.cellml'):
                    t_vars += re.findall(r'<variable(?=[^>]*\bname="T")[^>]*/?>',
                                         z.read(name).decode())
    return report, t_vars


def _initial_values(t_vars):
    return sorted({m for v in t_vars for m in re.findall(r'initial_value="([^"]+)"', v)})


def _t_rows(job):
    return [r for r in job['parameters'] if r['variable_name'] == 'T']


def test_libcuflynx_warns_about_a_global_its_modules_set_differently(tmp_path):
    from libcuflynx.utilities import module_instances
    if not hasattr(module_instances, 'ConflictingGlobalWarning'):
        pytest.skip('this libcuflynx has no ConflictingGlobalWarning (circulatory_autogen#535)')
    from libcuflynx.utilities.config_schemas import (load_component_registry,
                                                     load_expanded_vessel_records,
                                                     load_supermodule_registry)
    from libcuflynx.utilities.module_library import ModuleSources
    records = [{'name': n, 'module_type': mt, 'module_subtype': v, 'inp_instances': [],
                'out_instances': []} for n, (mt, v) in (('chan', CHANNEL), ('clock', CLOCK))]
    path = str(tmp_path / 'm_module_array.json')
    with open(path, 'w') as f:
        json.dump(records, f)
    files = ModuleSources({'module_library_dirs': [MODULES_DIR],
                           'use_builtin_modules': False}).config_files
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        _, rows = load_expanded_vessel_records(path, load_supermodule_registry(files),
                                               load_component_registry(files))
    conflicts = [w.message for w in caught
                 if isinstance(w.message, module_instances.ConflictingGlobalWarning)]
    assert [c.name for c in conflicts] == ['T']
    message = str(conflicts[0])
    assert '295.15 kelvin by "chan"' in message and '1 second by "clock"' in message
    assert 'units differ' in message
    # the first record's value is the one used
    assert [r['value'] for r in rows if r['variable_name'] == 'T'] == ['295.15']


def test_phlynx_keeps_the_first_of_two_rows_for_one_global_silently(jobs):
    work, channel, _, library = jobs
    rows = [r for r in channel['parameters'] if r['variable_name'] != 'T']
    (t,) = _t_rows(channel)
    for first, second in (('295.15', '310.15'), ('310.15', '295.15')):
        parameters = rows + [dict(t, value=first), dict(t, value=second)]
        report, t_vars = _export(work, f'two_rows_{first}', library, channel['instances'],
                                 parameters)
        assert report.get('ok'), report.get('error')
        assert _initial_values(t_vars) == [first]
        # no warning about the second row: PhLynx cannot tell the user T was set twice
        assert not [w for w in report.get('warnings', []) if re.search(r'\bT\b', w)]


def test_phlynx_refuses_one_global_name_in_two_units(jobs):
    # i_M's T is a temperature, the cardiac clock's a period: PhLynx joins them into one
    # global, and libCellML's validator rejects the units
    work, channel, clock, library = jobs
    instances = channel['instances'] + clock['instances']
    report, _ = _export(work, 'two_units', library, instances,
                        channel['parameters'] + clock['parameters'])
    assert not report.get('ok')
    assert report.get('stage') == 'flatten'
    assert "Variable 'T'" in report.get('error', '') and 'units' in report.get('error', '')


def test_phlynx_exports_a_global_no_row_sets_without_a_value(jobs):
    work, channel, clock, library = jobs
    instances = channel['instances'] + clock['instances']
    parameters = [r for r in channel['parameters'] + clock['parameters']
                  if r['variable_name'] != 'T']
    report, t_vars = _export(work, 'no_row', library, instances, parameters)
    # the export "succeeds", but each module's T is left unset and unconnected
    assert report.get('ok'), report.get('error')
    assert len(t_vars) == 2 and not _initial_values(t_vars)
