"""
libcuflynx writes a multi_port "sum" as a component multiport_sum_<module>_<variable>, whose <variable> is the total.
system.simulate() returns it as the output multiport_sum_<module>_<variable>/<variable> (the closed-loop systems'
blood-volume invariants read multiport_sum_volume_sum_q/q), and FV1D_volume_sum's exit volume_port
maps q_sum, so a volume_sum connected to it sums q_sum.
"""
import copy

import numpy as np
import pytest

from cam_testing import harness, library
from cam_testing import system as systems

MATH = 'xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/2.0#"'

# a vessel and the sum component libcuflynx writes for its multi_port "sum" variable q
CELLML = f'''<?xml version="1.0" encoding="UTF-8"?>
<model xmlns="http://www.cellml.org/cellml/2.0#" name="multiport_sum_output">
  <component name="env_module">
    <variable name="t" units="second" interface="public"/>
  </component>
  <component name="a_module">
    <variable name="t" units="second" interface="public"/>
    <variable name="x" units="dimensionless" initial_value="1" interface="public"/>
    <math {MATH}><apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>x</ci></apply>
      <cn cellml:units="per_s">1</cn></apply></math>
  </component>
  <component name="multiport_sum_s_q">
    <variable name="q" units="dimensionless"/>
    <variable name="x_a" units="dimensionless" interface="public"/>
    <math {MATH}><apply><eq/><ci>q</ci><apply><times/><cn cellml:units="dimensionless">2</cn><ci>x_a</ci></apply></apply></math>
  </component>
  <units name="per_s"><unit units="second" exponent="-1"/></units>
  <connection component_1="env_module" component_2="a_module"><map_variables variable_1="t" variable_2="t"/></connection>
  <connection component_1="a_module" component_2="multiport_sum_s_q"><map_variables variable_1="x" variable_2="x_a"/></connection>
</model>
'''

SPEC = {'dt': 0.1, 'sim_time': 1.0, 'pre_time': 0.0}
SOLVER_INFO = {'rtol': 1e-10, 'atol': 1e-12}


def test_simulate_returns_the_multiport_sum_total(tmp_path):
    pytest.importorskip('myokit')
    path = tmp_path / 'multiport_sum_output.cellml'
    path.write_text(CELLML)
    t, out = systems.simulate(str(path), SPEC, SOLVER_INFO)
    assert 'multiport_sum_s_q/q' in out, sorted(out)
    assert np.allclose(out['multiport_sum_s_q/q'], 2 * (1 + t))


def test_volume_sum_sums_the_q_sum_of_an_fv1d_volume_sum(tmp_path, monkeypatch):
    pytest.importorskip('libcuflynx')
    pytest.importorskip('myokit')
    version = library.load_version('volume_sum', 'ArgusUNPUBLISHED_v01')
    # the version's test network (two vessels into "mod"), with an FV1D_volume_sum "fv" as a third input
    network = copy.deepcopy(version.spec['harness'])
    rows = network['module_array']
    mod = next(r for r in rows if r[0] == 'mod')
    mod[3] += ' fv'
    rows.insert(0, ['fv', 'nn', 'FV1D_volume_sum', '', 'mod', 'default'])
    network['parameters'].append(['fv/q_sum', 'm3', 0.003, 'test value'])
    monkeypatch.setitem(version.spec, 'harness', network)
    model = harness.generate(version, str(tmp_path))
    t, out = systems.simulate(model, SPEC, SOLVER_INFO)
    total = out['multiport_sum_mod_q/q']
    assert np.ptp(out['va/q'] + out['vb/q']) > 0          # the vessels' volumes change
    assert np.allclose(total, out['va/q'] + out['vb/q'] + 0.003, rtol=0, atol=1e-12)
