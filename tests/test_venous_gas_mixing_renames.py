'''circulatory_autogen's module arrays name its built-in gas-exchange modules by their nn pairs.
Importing them gives the library's versions: tissue_GE_simple/nn has the gas_flux_port, so it is
Ursino2001_v02, and venous_gas_mixing/nn is ArgusUNPUBLISHED_v01.'''
import copy

import numpy as np
import pytest

from cam_testing import harness, library, module_array
from cam_testing import system as systems

# records of circulatory_autogen's resources/control_parasymp_module_array.csv
CA_CONTROL_PARASYMP = [
    {'name': 'leg_L_GE', 'module_type': 'tissue_GE_simple', 'module_subtype': 'nn',
     'inp_instances': ['leg_L_T', 'arterial_gas_transport'], 'out_instances': ['venous_gas_mixing', 'leg_L_LC']},
    {'name': 'trunk_C_GE', 'module_type': 'tissue_GE_simple', 'module_subtype': 'nn',
     'inp_instances': ['trunk_C_T', 'arterial_gas_transport'], 'out_instances': ['venous_gas_mixing']},
    {'name': 'venous_gas_mixing', 'module_type': 'venous_gas_mixing', 'module_subtype': 'nn',
     'inp_instances': ['leg_L_GE', 'trunk_C_GE'], 'out_instances': ['venous_gas_transport']},
    {'name': 'venous_gas_transport', 'module_type': 'gas_transport_simple', 'module_subtype': 'nn',
     'inp_instances': ['venous_gas_mixing'], 'out_instances': ['pulmonary_GE']},
]


def _versions():
    return {v.key: v for v in library.all_versions()}


def test_ca_gas_exchange_records_resolve_to_the_flux_port_version_and_the_mixer():
    pairs = {rec['name']: (rec['module_type'], rec['module_subtype'])
             for rec in module_array.to_library_versions(copy.deepcopy(CA_CONTROL_PARASYMP))}
    assert pairs['leg_L_GE'] == ('tissue_GE_simple', 'Ursino2001_v02')
    assert pairs['trunk_C_GE'] == ('tissue_GE_simple', 'Ursino2001_v02')
    assert pairs['venous_gas_mixing'] == ('venous_gas_mixing', 'ArgusUNPUBLISHED_v01')
    assert pairs['venous_gas_transport'] == ('gas_transport_simple', 'Albanese2016_v01')


def test_the_resolved_versions_have_the_flux_ports():
    versions = _versions()
    if 'venous_gas_mixing/ArgusUNPUBLISHED_v01' not in versions:
        pytest.skip('venous_gas_mixing is not in this library yet (circulatory-autogen-modules#53)')
    ge = versions['tissue_GE_simple/Ursino2001_v02'].config
    mixer = versions['venous_gas_mixing/ArgusUNPUBLISHED_v01'].config
    assert 'gas_flux_port' in {p['port_type'] for p in ge['exit_ports']}
    assert 'gas_flux_port' in {p['port_type'] for p in mixer['entrance_ports']}


def test_an_unconnected_mixer_gives_zeros(tmp_path, monkeypatch):
    '''With nothing connected the sums take their defaults, 0, so the mix is an obvious 0
    rather than plausible contents.'''
    pytest.importorskip('libcuflynx')
    pytest.importorskip('myokit')
    if 'venous_gas_mixing/ArgusUNPUBLISHED_v01' not in _versions():
        pytest.skip('venous_gas_mixing is not in this library yet (circulatory-autogen-modules#53)')
    version = library.load_version('venous_gas_mixing', 'ArgusUNPUBLISHED_v01')
    monkeypatch.delitem(version.spec, 'harness')
    model = harness.generate(version, str(tmp_path))
    names = ['mod/v_mix', 'mod/C_O2_mix', 'mod/C_CO2_mix']
    t, out = systems.simulate(model, {'dt': 0.1, 'sim_time': 1.0, 'pre_time': 0.0},
                              {'rtol': 1e-8, 'atol': 1e-10}, names=names)
    for name in names:
        assert name in out, sorted(out)
        assert np.all(out[name] == 0), (name, out[name][:3])
