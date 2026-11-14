'''
Check that the junction removal changes no system model's equations.

Each system model is generated twice with the same libcuflynx: from a checkout of the library
before the change and from one after it. Then:

* the equations (pass/fail): both models' state derivatives are evaluated at the same states,
  taken from a simulation of the original, and must agree to within RHS_TOL times the rounding
  noise of each derivative there (see rhs_difference). This is exact up to rounding and does not
  depend on how the two CellML files order their equations;
* the trajectories (reported): both are simulated at the system's equivalence solver_info and every
  output of the original that the new model has is compared (normalised max difference, as
  cam_testing.system.compare). The junctions' port flows are mapped to their new names
  (v_in_sum -> v_in, v_out_sum -> v_out); the numbered ones (v_out_1, ...) have no single
  counterpart on the vessel and are left out.

Run it in both checkouts, then compare (libcuflynx on PYTHONPATH):

    python tools/verify_junction_conversion.py generate <out_dir> [system ...]    # in each checkout
    python tools/verify_junction_conversion.py compare <before_dir> <after_dir> [--markdown]

generate writes <out_dir>/<system>.cellml (with its parameters) and <out_dir>/status.json.
'''
import argparse
import contextlib
import io
import json
import os
import re
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing.system import compare, generate, simulate, system_models  # noqa: E402

SOLVER_INFO = {'rtol': 1e-10, 'atol': 1e-12}
RHS_TOL = 100.0     # a difference of up to this many times the rounding noise is rounding
TRAJ_TOL = 1e-6     # reported against each system's equivalence tol (default 1e-6)
RENAMED = {'v_in_sum': 'v_in', 'v_out_sum': 'v_out'}
NUMBERED = ('v_in_1', 'v_in_2', 'v_out_1', 'v_out_2', 'v_out_total', 'v_in_total')


def cmd_generate(out_dir, names):
    os.makedirs(out_dir, exist_ok=True)
    status = {}
    for system in system_models():
        if names and system.name not in names and system.id not in names:
            continue
        work = tempfile.mkdtemp(prefix=f'jr_{system.name}_')
        try:
            path, dropped = _generate(system, work)
            gen_dir = os.path.dirname(path)
            dest = os.path.join(out_dir, system.name)
            shutil.rmtree(dest, ignore_errors=True)
            shutil.copytree(gen_dir, dest)
            status[system.name] = {'ok': True, 'cellml': os.path.join(dest, os.path.basename(path)),
                                   'id': system.id, 'dropped_globals': dropped}
        except Exception as e:
            status[system.name] = {'ok': False, 'error': str(e)[-600:], 'id': system.id}
        finally:
            shutil.rmtree(work, ignore_errors=True)
        print(system.name, 'ok' if status[system.name]['ok'] else 'FAILED', flush=True)
    with open(os.path.join(out_dir, 'status.json'), 'w') as f:
        json.dump(status, f, indent=1)


UNKNOWN_GLOBAL = re.compile(r'no module in the model has a global constant "([^"]+)"')


def _generate(system, work):
    '''Generate the system; if its parameters file sets globals no module declares (stale rows the
    generator refuses), generate a copy without them instead. Returns (path, dropped names).'''
    from cam_testing.harness import GenerationFailed
    dropped, current = [], system
    for attempt in range(6):
        try:
            return generate(current, os.path.join(work, f'gen{attempt}')), dropped
        except GenerationFailed as e:
            names = sorted(set(UNKNOWN_GLOBAL.findall(str(e))) - set(dropped))
            if not names:
                raise
            dropped += names
        copy_dir = os.path.join(work, f'system_copy{attempt}', system.name)
        shutil.copytree(system.dir, copy_dir, ignore=shutil.ignore_patterns('results', 'plots', 'reference'))
        params = os.path.join(copy_dir, f'{system.name}_parameters.json')
        doc = json.load(open(params))
        for name in dropped:
            doc.get('global', {}).pop(name, None)
            doc.pop(name, None)
        with open(params, 'w') as f:
            json.dump(doc, f, indent=1)
        current = type(system)(system.name, system.category, copy_dir, system.spec)
    raise RuntimeError(f'{system.name}: still refused after dropping {dropped}')


def _myokit_model(path, solver_info):
    from libcuflynx.solver_wrappers import get_simulation_helper
    with contextlib.redirect_stdout(io.StringIO()):
        helper = get_simulation_helper(model_path=path, solver='CVODE_myokit', model_type='cellml', dt=0.01,
                                       sim_time=0.1, solver_info=solver_info, pre_time=0.0)
    return helper.model


def rhs_difference(original_path, new_path, times=(0.05, 0.3, 0.6, 0.9)):
    '''The worst ratio, over states and times, of |f_original - f_new| to the original's own rounding
    noise at that state (its derivative re-evaluated at the state perturbed by ~1e-13), floored at
    1e-12 of the derivative. Returns (ratio, state name, number of states).'''
    import myokit
    original, new = _myokit_model(original_path, SOLVER_INFO), _myokit_model(new_path, SOLVER_INFO)
    key = lambda v: v.qname().replace('_module.', '/')
    index_o = {key(v): i for i, v in enumerate(original.states())}
    index_n = {key(v): i for i, v in enumerate(new.states())}
    if set(index_o) != set(index_n):
        raise AssertionError(f'states differ: only before {sorted(set(index_o) - set(index_n))[:5]}, '
                             f'only after {sorted(set(index_n) - set(index_o))[:5]}')
    rng = np.random.default_rng(0)
    sim = myokit.Simulation(original)
    sim.set_tolerance(1e-12, 1e-10)
    worst, worst_state, evaluated = 0.0, None, 0
    for t in times:
        try:
            sim.run(t - sim.time(), log=myokit.LOG_NONE)
        except myokit.SimulationError:
            break   # the original stops here at these tolerances: compare at the states reached
        evaluated += 1
        x_o = np.array(sim.state())
        x_n = np.zeros(len(x_o))
        for k, i in index_o.items():
            x_n[index_n[k]] = x_o[i]
        d_o = np.array(original.evaluate_derivatives(list(x_o)))
        d_n = np.array(new.evaluate_derivatives(list(x_n)))
        noise = np.zeros(len(x_o))
        for _ in range(3):
            x_p = x_o * (1 + 1e-13 * rng.standard_normal(len(x_o)))
            noise = np.maximum(noise, np.abs(np.array(original.evaluate_derivatives(list(x_p))) - d_o))
        for k, i in index_o.items():
            floor = max(noise[i], 1e-12 * max(abs(d_o[i]), 1e-300))
            ratio = abs(d_o[i] - d_n[index_n[k]]) / floor
            if ratio > worst:
                worst, worst_state = ratio, k
    if not evaluated:
        raise RuntimeError(f'the original does not reach t = {times[0]} at rtol 1e-10')
    return worst, worst_state, len(index_o)


def trajectories(original_path, new_path, spec):
    eq = spec.get('equivalence') or {}
    solver_info = eq.get('solver_info') or SOLVER_INFO
    _, before = simulate(original_path, spec, solver_info)
    mapping = {}
    for name in before:
        inst, _, var = name.rpartition('/')
        if var in NUMBERED:
            continue
        mapping[name] = f'{inst}/{RENAMED.get(var, var)}'
    _, after = simulate(new_path, spec, solver_info, names=sorted(set(mapping.values())))
    # as the system equivalence test compares: the spec's ignore, wrapped and reciprocal outputs
    ignore = dict(eq.get('ignore') or {})
    ignore.update({n: 'time' for n in before if n.endswith('/t')})
    ignore.update({n: 'numbered junction port flow' for n in before if n not in mapping})
    rows, missing = compare(before, after, mapping, float(eq.get('tol', TRAJ_TOL)), ignore=ignore,
                            wrapped=eq.get('wrapped'), reciprocal=eq.get('reciprocal'))
    worst = max(rows, key=lambda r: r['difference']) if rows else None
    compared = len(rows)
    worst_name = worst['reference'] if worst else None
    worst = worst['difference'] if worst else 0.0
    return worst, worst_name, compared, missing, float(eq.get('tol', TRAJ_TOL))


def cmd_compare(before_dir, after_dir, markdown, only=()):
    before = json.load(open(os.path.join(before_dir, 'status.json')))
    after = json.load(open(os.path.join(after_dir, 'status.json')))
    specs = {s.name: s.spec for s in system_models()}
    rows = []
    for name in sorted(set(before) | set(after)):
        if only and name not in only:
            continue
        b, a = before.get(name), after.get(name)
        row = {'system': name}
        if not b or not a or not b['ok'] or not a['ok']:
            row['result'] = (f'before {"ok" if b and b["ok"] else "fails"}, after {"ok" if a and a["ok"] else "fails"}')
            row['error_before'] = (b or {}).get('error', '')[-200:]
            row['error_after'] = (a or {}).get('error', '')[-200:]
            rows.append(row)
            continue
        try:
            ratio, state, n = rhs_difference(b['cellml'], a['cellml'])
            row.update(rhs_ratio=float(ratio), rhs_state=state, states=n, rhs_ok=bool(ratio <= RHS_TOL))
        except Exception as e:
            row.update(rhs_ok=False, rhs_error=str(e)[-300:])
        try:
            d, worst, compared, missing, tol = trajectories(b['cellml'], a['cellml'], specs.get(name, {}))
            row.update(traj_max_rel=float(d), traj_worst=worst, compared=compared, missing=len(missing),
                       missing_examples=missing[:4], traj_tol=tol)
        except Exception as e:
            row.update(traj_error=str(e)[-300:])
        row['result'] = 'equivalent' if row.get('rhs_ok') else 'DIFFERS'
        rows.append(row)
        print(json.dumps(row), flush=True)
    with open(os.path.join(after_dir, 'comparison.json' if not only else 'comparison_some.json'), 'w') as f:
        json.dump(rows, f, indent=1)
    if markdown:
        print('| system | states | RHS worst / noise | trajectories max rel diff (tol) | outputs compared | result |')
        print('|---|---|---|---|---|---|')
        for r in rows:
            if 'rhs_ratio' in r or 'rhs_error' in r:
                rhs = f'{r["rhs_ratio"]:.3g} ({r["rhs_state"]})' if 'rhs_ratio' in r else r['rhs_error'][:60]
                traj = f'{r["traj_max_rel"]:.2e} ({r["traj_tol"]:g})' if 'traj_max_rel' in r else r.get('traj_error', '')[:60]
                print(f'| {r["system"]} | {r.get("states", "")} | {rhs} | {traj} | {r.get("compared", "")} | {r["result"]} |')
            else:
                print(f'| {r["system"]} | | | | | {r["result"]} |')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='cmd', required=True)
    g = sub.add_parser('generate')
    g.add_argument('out_dir')
    g.add_argument('systems', nargs='*')
    c = sub.add_parser('compare')
    c.add_argument('before_dir')
    c.add_argument('after_dir')
    c.add_argument('--markdown', action='store_true')
    c.add_argument('systems', nargs='*', help='only these systems')
    args = parser.parse_args(argv)
    if args.cmd == 'generate':
        cmd_generate(os.path.abspath(args.out_dir), args.systems)
    else:
        cmd_compare(args.before_dir, args.after_dir, args.markdown, args.systems)


if __name__ == '__main__':
    main()
