"""
The homogenisation examples' mass-conservation invariant (tools/build_homogenisation_examples.py) on a
real run of cell_disc_lumped: it holds for the model as generated, although the cells' fill-up transient
(~5 ms) is faster than the 10 ms output step, and it fails when a term of the balance is missing or wrong.

    pytest tests/test_homogenisation_invariants.py
"""
import numpy as np
import pytest

from cam_testing import checks, system as systems

pytestmark = pytest.mark.system_model


@pytest.fixture(scope='module')
def run(tmp_path_factory):
    system = systems.load_system('homogenisation/cell_disc_lumped')
    t, out = systems.run_model(system, str(tmp_path_factory.mktemp('cell_disc_lumped')))
    env = {k.replace('/', '__'): v for k, v in out.items()}
    (mass,) = [inv for inv in system.spec['invariants'] if inv['description'].startswith('mass conservation')]
    return t, env, {'invariants': [mass]}


def _failures(run, **changed):
    t, env, spec = run
    return checks._check_invariants(spec, t, {**env, **changed}, {}, where='run')


def test_mass_invariant_holds(run):
    assert _failures(run) == []


def test_mass_invariant_catches_missing_consumption(run):
    # the balance without the consumption term: the outputs carry no M_reduced
    assert _failures(run, c0__M_reduced=np.zeros_like(run[1]['c0__M_reduced']))


def test_mass_invariant_catches_supply_error(run):
    assert _failures(run, e0__J=0.999 * run[1]['e0__J'])
