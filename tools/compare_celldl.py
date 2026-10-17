"""Compares the CellML that CellDL's bg2cellml makes from a module's CellDL bond graph (tools/build_celldl.py) with the
module itself: both run from the same state with the same boundary concentrations, every internal concentration
compared over time (the bond graph's store amounts q [mol] over the volume).

    python tools/compare_celldl.py SPEC_DIR [module ...]

SPEC_DIR holds <module>.json / .cellml and <module>_run.json / .cellml (rest, and the verification config's
run_parameters). Prints, per module and run, the largest difference relative to each species' scale.
"""
import contextlib
import io
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_celldl  # noqa: E402
import build_metabolism_detailed as detailed  # noqa: E402

T_END, DT = 60.0, 1.0


def module_run(module, boundary):
    '''The module's internal concentrations [mM] at t = 0, DT, ... T_END with its boundary concentrations set.'''
    from cam_testing import checks, harness
    from cam_testing.library import load_version
    cm = checks.ComponentModel(load_version(module, detailed.VERSION), work_dir=os.path.join(
        os.environ.get('TMPDIR', '/tmp'), f'celldl_compare_{module}'))
    bg = detailed.BUILDERS[module]()[0]
    internal = [s.name for s in bg.species if not s.boundary]
    with contextlib.redirect_stdout(io.StringIO()):
        helper = harness.simulation_helper(cm.model_path, dict(cm.spec, sim_time=T_END, dt=DT))
        helper.reset_states()
        if boundary:
            helper.set_param_vals(['parameters/C_' + k + '_mod' for k in boundary], list(boundary.values()),
                                  change_states=True)
        assert helper.run() is not False
    return {s: np.asarray(helper.get_results([[harness.output_name('mod/C_' + s)]], flatten=True)[0], float).ravel()
            for s in internal}


def bond_graph_run(cellml_path, volume, species):
    '''The bond graph's store amounts over the volume [mM] at the same times.'''
    import myokit
    import myokit.formats.cellml as mc
    model = mc.CellMLImporter().model(cellml_path)
    model.validate()
    sim = myokit.Simulation(model)
    # the stores are amounts in a cell's volume (~1e-16 mol): the absolute tolerance is the module's (1e-12 mM) in mol
    sim.set_tolerance(1e-12 * volume, 1e-10)
    log = sim.run(T_END + DT / 2, log_interval=DT)
    # a store's state is q_<species>_<location> (bg2cellml names variables from symbol, species and location)
    states = {s.name().rsplit('_', 1)[0]: s.qname() for s in model.states()}
    return {s: np.asarray(log[states[f'q_{s}']], float) / volume for s in species}


def compare(spec_dir, module):
    out = {}
    with contextlib.redirect_stdout(io.StringIO()):
        _, volume = build_celldl.MODULES[module]
    for name, boundary in ((module, {}), (f'{module}_run', build_celldl.run_parameters(module))):
        mod = module_run(module, boundary)
        bg = bond_graph_run(os.path.join(spec_dir, f'{name}.cellml'), volume, list(mod))
        worst = {}
        for s in mod:
            n = min(len(mod[s]), len(bg[s]))
            scale = max(np.max(np.abs(mod[s][:n])), 1e-300)
            worst[s] = float(np.max(np.abs(bg[s][:n] - mod[s][:n])) / scale)
        s_max = max(worst, key=worst.get)
        out[name] = {'max_relative_difference': worst[s_max], 'species': s_max, 'per_species': worst,
                     'boundary': boundary}
    return out


if __name__ == '__main__':
    spec_dir = sys.argv[1]
    results = {}
    for module in sys.argv[2:] or build_celldl.MODULES:
        with contextlib.redirect_stdout(io.StringIO()):
            r = compare(spec_dir, module)
        results.update(r)
        for name, v in r.items():
            print(f'{name:22s} largest relative difference {v["max_relative_difference"]:.2e} ({v["species"]})')
    with open(os.path.join(spec_dir, 'comparison.json'), 'w') as f:
        json.dump(results, f, indent=1)
