"""
The square-channel Poiseuille transport examples (modules/system/poiseuille/), checked against a direct
scipy solve of the same finite-volume equations (tools/build_poiseuille_examples.py), independent of
libcuflynx's generation and port wiring. The exact-solution and conservation checks are the systems'
invariants (tests/test_systems.py).

    pytest tests/test_poiseuille_channels.py
"""
import os
import sys

import numpy as np
import pytest

from cam_testing import system as systems

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tools'))
import build_poiseuille_examples as examples  # noqa: E402

pytestmark = pytest.mark.system_model


@pytest.mark.parametrize('model', sorted(examples.CASES))
def test_matches_direct_solve(model, tmp_path):
    system = systems.load_system(model)
    t, out = systems.run_model(system, str(tmp_path))
    P_ref, C_ref = examples.scipy_solve(examples.CASES[model], t)
    C = np.stack([out[f'{c}/C'] for c in examples.cells()])
    P = np.stack([out[f'{c}/P'] for c in examples.cells()])
    case = examples.CASES[model]
    p_scale = max(abs(case['P_in']) + abs(case['P_amp']), abs(case['P_out']), 1e-12)
    assert np.max(np.abs(C - C_ref)) <= 1e-6, f'concentration differs by {np.max(np.abs(C - C_ref)):.3g} mM'
    assert np.max(np.abs(P - P_ref)) <= 1e-6 * p_scale, f'pressure differs by {np.max(np.abs(P - P_ref)):.3g} Pa'
