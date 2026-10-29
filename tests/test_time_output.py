"""
system.simulate() returns time as t and never as an output. Myokit merges every component's t into one
time variable named after a component it picks per process, so a '<module>/t' output would come and go
between runs (and fail an equivalence test as "missing").
"""
import numpy as np
import pytest

from cam_testing import system as systems

MATH = 'xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/2.0#"'

# every component is a module, so whichever one Myokit names time after, the name is a module output
CELLML = f'''<?xml version="1.0" encoding="UTF-8"?>
<model xmlns="http://www.cellml.org/cellml/2.0#" name="time_named_after_a_module">
  <units name="per_s"><unit units="second" exponent="-1"/></units>
  <component name="env_module">
    <variable name="t" units="second" interface="public"/>
  </component>
  <component name="a_module">
    <variable name="t" units="second" interface="public"/>
    <variable name="x" units="dimensionless" initial_value="1"/>
    <variable name="k" units="per_s" initial_value="1"/>
    <math {MATH}><apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>x</ci></apply>
      <apply><minus/><apply><times/><ci>k</ci><ci>x</ci></apply></apply></apply></math>
  </component>
  <component name="b_module">
    <variable name="t" units="second" interface="public"/>
    <variable name="y" units="second"/>
    <math {MATH}><apply><eq/><ci>y</ci><apply><times/><cn cellml:units="dimensionless">2</cn><ci>t</ci></apply></apply></math>
  </component>
  <component name="c_module">
    <variable name="t" units="second" interface="public"/>
    <variable name="z" units="second"/>
    <math {MATH}><apply><eq/><ci>z</ci><apply><plus/><ci>t</ci><cn cellml:units="second">1</cn></apply></apply></math>
  </component>
  <connection component_1="env_module" component_2="a_module"><map_variables variable_1="t" variable_2="t"/></connection>
  <connection component_1="env_module" component_2="b_module"><map_variables variable_1="t" variable_2="t"/></connection>
  <connection component_1="env_module" component_2="c_module"><map_variables variable_1="t" variable_2="t"/></connection>
</model>
'''

SPEC = {'dt': 0.1, 'sim_time': 1.0, 'pre_time': 0.0}
SOLVER_INFO = {'rtol': 1e-10, 'atol': 1e-12}
TIME_NAMES = {'env/t', 'a/t', 'b/t', 'c/t'}


@pytest.fixture
def model(tmp_path):
    pytest.importorskip('myokit')
    path = tmp_path / 'time_named_after_a_module.cellml'
    path.write_text(CELLML)
    return str(path)


def test_time_is_not_an_output(model):
    t, out = systems.simulate(model, SPEC, SOLVER_INFO)
    assert not TIME_NAMES & set(out), sorted(out)
    assert {'a/x', 'b/y', 'c/z'} <= set(out)
    assert np.allclose(out['b/y'], 2 * t) and np.allclose(out['c/z'], t + 1)


def test_time_asked_for_by_name_is_dropped(model):
    t, out = systems.simulate(model, SPEC, SOLVER_INFO, names=sorted(TIME_NAMES | {'a/x'}))
    assert sorted(out) == ['a/x']
    assert np.allclose(out['a/x'], np.exp(-t), rtol=1e-6)
