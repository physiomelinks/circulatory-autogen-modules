'''
Check that moving a system model off the junction module types changes nothing.

For each system model whose vessel array uses a junction module_type, generate it twice from this
library: as it is, and with every junction replaced by its ordinary twin
(libcuflynx.utilities.junction_migration). Then:

* the equations: both models' state derivatives are evaluated at the same states (taken from a
  simulation of the original) and must agree to within 100 times the rounding noise of each
  derivative there (see rhs_difference). This is the pass/fail check: it is exact up to rounding,
  and independent of how the two CellML files order their equations;
* the trajectories: both are simulated (the spec's equivalence solver_info) and every variable is
  compared through the junctions' renamed port flows (v_in_sum -> v_in, v_out_sum -> v_out; the
  numbered v_out_1/v_in_1 ... have no single counterpart). This is reported, not judged: models
  with valves differ by ~1e-6 because CVODE's steps around valve events move with the order of
  the equations, and discontinuous outputs (phase wraps, floor indicators) then differ by O(1).

Run before the junction module_types are removed (both kinds of module must be in modules/):

    python tools/verify_junction_conversion.py [--write]

--write then converts the arrays in place.
'''
import argparse
import json
import os
import shutil
import sys
import tempfile

import contextlib
import io

import myokit
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing.system import compare, generate, simulate, system_models  # noqa: E402
from libcuflynx.utilities.junction_migration import RENAMED_VARIABLES, migrate_file, migrate_records  # noqa: E402

SOLVER_INFO = {'rtol': 1e-10, 'atol': 1e-12}
TOL = 1e-6


def has_parameters(system):
    return os.path.exists(os.path.join(system.dir, f'{system.name}_parameters.csv'))


def uses_junctions(system):
    path = os.path.join(system.dir, f'{system.name}_vessel_array.json')
    if not os.path.exists(path):
        return False
    with open(path) as f:
        return bool(migrate_records(json.load(f))[1])


def converted_copy(system, work_dir):
    '''A copy of the system dir whose vessel array has no junctions; returns (system, renamed).'''
    copy_dir = os.path.join(work_dir, 'converted', system.name)
    shutil.copytree(system.dir, copy_dir, ignore=shutil.ignore_patterns('results', 'plots', 'reference'))
    renamed = migrate_file(os.path.join(copy_dir, f'{system.name}_vessel_array.json'),
                           os.path.join(copy_dir, f'{system.name}_parameters.csv'))
    return type(system)(system.name, system.category, copy_dir, system.spec), renamed


# a difference of up to this many times the rounding noise is rounding
RHS_TOL = 100.0


def _myokit_model(path):
    from libcuflynx.solver_wrappers import get_simulation_helper
    with contextlib.redirect_stdout(io.StringIO()):
        helper = get_simulation_helper(model_path=path, solver='CVODE_myokit', model_type='cellml', dt=0.01,
                                       sim_time=0.1, solver_info=SOLVER_INFO, pre_time=0.0)
    return helper.model


def rhs_difference(original_path, converted_path, times=(0.05, 0.3, 0.6, 0.9)):
    '''How far the two models' state derivatives are apart at shared states, in units of the
    original's own rounding noise there.

    A derivative that is a small difference of large terms (an inertance flow, (u_in - u_out - R v)/I
    with I = 1e-6, cancels pressures of ~1e4 Pa) changes by a whole per cent when one of its terms
    moves by one ulp, which is what reordering the equations does. So each derivative is compared
    with the original evaluated again at the state perturbed by ~1e-13: the noise there. Returns
    (worst ratio of the difference to max(noise, 1e-12 of the derivative's scale), that state).
    '''
    original, converted = _myokit_model(original_path), _myokit_model(converted_path)
    key = lambda v: v.qname().replace('_module.', '/')
    index_o = {key(v): i for i, v in enumerate(original.states())}
    index_c = {key(v): i for i, v in enumerate(converted.states())}
    if set(index_o) != set(index_c):
        raise AssertionError(f'states differ: only original {sorted(set(index_o) - set(index_c))[:5]}, '
                             f'only converted {sorted(set(index_c) - set(index_o))[:5]}')
    rng = np.random.default_rng(0)
    sim = myokit.Simulation(original)
    sim.set_tolerance(1e-12, 1e-10)
    worst, worst_state = 0.0, None
    for t in times:
        sim.run(t - sim.time(), log=myokit.LOG_NONE)
        x_o = np.array(sim.state())
        x_c = np.zeros(len(x_o))
        for k, i in index_o.items():
            x_c[index_c[k]] = x_o[i]
        d_o = np.array(original.evaluate_derivatives(list(x_o)))
        d_c = np.array(converted.evaluate_derivatives(list(x_c)))
        noise = np.zeros(len(x_o))
        for _ in range(3):
            x_p = x_o * (1 + 1e-13 * rng.standard_normal(len(x_o)))
            noise = np.maximum(noise, np.abs(np.array(original.evaluate_derivatives(list(x_p))) - d_o))
        for k, i in index_o.items():
            floor = max(noise[i], 1e-12 * max(abs(d_o[i]), 1e-300))
            ratio = abs(d_o[i] - d_c[index_c[k]]) / floor
            if ratio > worst:
                worst, worst_state = ratio, k
    return worst, worst_state


def verify(system):
    work_dir = tempfile.mkdtemp(prefix=f'junction_conversion_{system.name}_')
    original_path = generate(system, os.path.join(work_dir, 'original'))
    converted, renamed = converted_copy(system, work_dir)
    converted_path = generate(converted, os.path.join(work_dir, 'converted_gen'))
    rhs = rhs_difference(original_path, converted_path)
    eq = system.spec.get('equivalence') or {}
    solver_info = eq.get('solver_info') or SOLVER_INFO
    _, original = simulate(original_path, system.spec, solver_info)
    output_map = {f'{name}/{old}': f'{name}/{new}' for name in renamed
                  for old, new in RENAMED_VARIABLES.items() if old.endswith('_sum')}
    numbered = {f'{name}/{old}': 'no single counterpart' for name in renamed
                for old in RENAMED_VARIABLES if not old.endswith('_sum')}
    wanted = sorted({output_map.get(n, n) for n in original})
    _, new = simulate(converted_path, system.spec, solver_info, names=wanted)
    ignore = dict(eq.get('ignore') or {}, **numbered)
    ignore.update({n: 'time' for n in original if n.endswith('/t')})
    rows, missing = compare(original, new, output_map, float(eq.get('tol', TOL)), ignore=ignore,
                            wrapped=eq.get('wrapped'), reciprocal=eq.get('reciprocal'))
    missing = [m for m in missing if m not in ignore]
    worst = max(rows, key=lambda r: r['difference']) if rows else None
    shutil.rmtree(work_dir, ignore_errors=True)
    return renamed, rhs, rows, missing, worst


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--write', action='store_true', help='convert the vessel arrays in place afterwards')
    parser.add_argument('systems', nargs='*', help='system names (default: every one that uses junctions)')
    args = parser.parse_args(argv)
    systems = [s for s in system_models() if uses_junctions(s) and (not args.systems or s.name in args.systems)]
    failed = False
    for system in systems:
        if not has_parameters(system):
            print(f'SKIPPED  {system.id}: no {system.name}_parameters.csv, and its junction instances have no '
                  f'values, so it does not generate before or after')
            continue
        try:
            renamed, rhs, rows, missing, worst = verify(system)
        except Exception as e:  # GenerationFailed and friends: report, don't stop the sweep
            failed = True
            print(f'ERROR    {system.id}: {str(e).strip().splitlines()[-1][:300]}')
            continue
        bad = [r for r in rows if not r['ok']]
        status = 'OK' if rhs[0] <= RHS_TOL else 'DIFFERS'
        failed |= status != 'OK'
        print(f'{status:8s} {system.id}: {len(renamed)} junctions converted; equations: worst {rhs[0]:.2g}x rounding noise '
              f'({rhs[1]}); trajectories: {len(rows)} variables, {len(bad)} beyond the spec tolerance'
              + (f', worst {worst["reference"]} {worst["difference"]:.1e}' if worst else '')
              + (f', {len(missing)} missing' if missing else ''))
    if args.write and not failed:
        for system in systems:
            migrate_file(os.path.join(system.dir, f'{system.name}_vessel_array.json'),
                         os.path.join(system.dir, f'{system.name}_parameters.csv'))
        print(f'converted {len(systems)} vessel arrays')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
