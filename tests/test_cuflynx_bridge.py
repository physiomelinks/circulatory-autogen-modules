# The CUFLynx bridge (cam_testing/bridges/cuflynx/simulate_omex.py), apart from a CUFLynx binary.
import importlib.util
import os
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    'simulate_omex', os.path.join(HERE, '..', 'cam_testing', 'bridges', 'cuflynx', 'simulate_omex.py'))
bridge = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bridge)

# a constants-only network as PhLynx exports it: the values in instance_parameters, mapped into the module
CELLML = """<?xml version="1.0"?>
<model xmlns="http://www.cellml.org/cellml/2.0#" name="m">
  <component name="environment">
    <variable name="t" units="second" interface="public"/>
    <variable name="clock" units="second" initial_value="0"/>
  </component>
  <component name="instance_parameters">
    <variable name="v_1" units="m3_per_s" initial_value="4.28e-06" interface="public"/>
  </component>
  <component name="mod">
    <variable name="t" units="second" interface="public"/>
    <variable name="v_1" units="m3_per_s" interface="public_and_private"/>
    <variable name="q" units="m3" interface="public_and_private"/>
  </component>
  <connection component_1="environment" component_2="mod"><map_variables variable_1="t" variable_2="t"/></connection>
  <connection component_1="instance_parameters" component_2="mod"><map_variables variable_1="v_1" variable_2="v_1"/></connection>
</model>
"""


def test_constant_values_follow_connections():
    c = bridge.constant_values(CELLML)
    assert c['mod/v_1'] == 4.28e-06 and c['instance_parameters/v_1'] == 4.28e-06
    assert 'mod/q' not in c and 'mod/t' not in c


def test_a_requested_constant_is_reported_over_the_run(tmp_path):
    omex = tmp_path / 'x.omex'
    with zipfile.ZipFile(omex, 'w') as z:
        z.writestr('manifest.xml', '<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">'
                   '<content location="./model.cellml" format="http://identifiers.org/combine.specifications/cellml" '
                   'master="true"/></omexManifest>')
        z.writestr('model.cellml', CELLML)
    res = {'time': [0.0, 0.5, 1.0], 'outputs': {'environment/clock': [0.0, 0.5, 1.0]}}
    bridge.add_constant_outputs(res, ['mod/v_1', 'mod/q'], str(omex))
    assert res['outputs']['mod/v_1'] == [4.28e-06] * 3
    assert 'mod/q' not in res['outputs'] and res['constant_outputs'] == ['mod/v_1']
