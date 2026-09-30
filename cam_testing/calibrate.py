"""
validation_test_calibrate: calibrate a version, at one of its instances, with libcuflynx parameter
identification on that instance's data, then score its prediction of held-out data.

The inputs are ordinary libcuflynx files, committed in the instance directory
(versions/<version>/instances/<instance>/), so the same calibration can be run with libcuflynx
directly:

    <instance>_obs_data.json              obs_data ("obs_data_name": "<instance>"): its data_items
                                          are fitted by CVS0DParamID; its prediction_items that
                                          carry a value are the held-out data the calibrated
                                          model is validated against, scored on
                                          prediction_window (without any, the data_items are)
    <instance>_params_for_id.csv          parameters calibrated, with bounds

and the results are written there too, and committed:

    <instance>_calibrated_parameters.csv  the instance's parameters with the calibrated values
    <instance>_calibration.json           summary: method, cost, fitted values, scores, date

The version is generated alone as a single vessel named "mod" (cam_testing.harness.VESSEL),
so operands are "mod/<variable>" and params_for_id vessel_name is "mod".

Spec (validation.<instance>.calibrate in <module_type>_<version>_verification_config.json; file
references relative to the version directory, defaulting to the instance's files):

    status: active
    source: citation / URL
    obs_data: instances/<instance>/<instance>_obs_data.json
    params_for_id: instances/<instance>/<instance>_params_for_id.csv
    prediction_window: [15, 20]           model time scored for series predictions (optional)
    initial_parameters: {a: 0.8}          optional optimiser start (default: nominal values)
    starts: [{a: 1}, {a: -1}]             optional: several starts, each calibrated and checked
    expected_parameters: {a: 0.2}         optional: known true values the fit must recover
    expected_compare: abs                 optional: compare |value| (e.g. symmetric minima)
    expected_rtol: 0.05
    z_threshold: 2.0                      constant data items: |model - value| <= z_threshold * std
    method: CMA-ES
    optimiser_options: {num_calls_to_function: 3000, seed: 1}
    metric: log_rmse                      log_rmse | nrmse
    threshold: 0.5                        on the prediction window

The obs_data files for tabular data are made with
    python -m cam_testing.calibrate from-csv --csv instances/<instance>/data.csv --time-column Year \
        --time-offset 1900 --variables x=Hare y=Lynx --window 0 15 --noise relative 0.25 \
        --name <instance> --out instances/<instance>/<instance>_obs_data.json
and held-out data is added to it as prediction_items with --validation-window (the protocol
then runs to the end of the held-out data).
"""
import contextlib
import copy
import csv
import io
import json
import os

import numpy as np

from cam_testing import harness, plots


def load_data(data_dir, v):
    with open(os.path.join(data_dir, v['data'])) as f:
        rows = list(csv.DictReader(f, skipinitialspace=True))
    rows = [{k.strip(): val for k, val in r.items()} for r in rows]
    t = np.array([float(r[v.get('time_column', 't')]) for r in rows]) - float(v.get('time_offset', 0.0))
    data = {m: np.array([float(r[col]) for r in rows]) for m, col in v['variables'].items()}
    return t, data


def score(metric, model, data):
    if metric == 'log_rmse':
        return float(np.sqrt(np.mean((np.log(model) - np.log(data)) ** 2)))
    rng = np.ptp(data) or max(np.max(np.abs(data)), 1e-300)
    return float(np.sqrt(np.mean((model - data) ** 2)) / rng)


def obs_data_from_series(t, data, window, noise, unit='dimensionless', name=None):
    """libcuflynx obs_data (series items) for data evenly spaced from model time 0."""
    mask = (t >= window[0] - 1e-9) & (t <= window[1] + 1e-9)
    tw = t[mask]
    if abs(tw[0]) > 1e-9:
        raise ValueError('series data must start at model time 0 (use --time-offset)')
    dts = np.diff(tw)
    if not np.allclose(dts, dts[0]):
        raise ValueError('series data must be evenly spaced')
    items = []
    for model_var, series in data.items():
        values = series[mask]
        if noise[0] == 'relative':
            std = np.maximum(np.abs(values) * float(noise[1]), 1e-12)
        else:
            std = np.full_like(values, float(noise[1]))
        name = harness.output_name(model_var)
        items.append({
            'data_item_name': name, 'trace_name_for_plotting': model_var, 'data_type': 'series',
            'operation': None, 'operands': [name], 'unit': unit, 'cost_type': 'gaussian_MLE',
            'weight': 1.0, 'value': [float(x) for x in values], 'std': [float(x) for x in std],
            'obs_dt': float(dts[0]), 'experiment_idx': 0, 'subexperiment_idx': 0,
        })
    out = {'obs_data_name': name} if name else {}
    out.update({'protocol_info': {'pre_times': [0.0], 'sim_times': [[float(tw[-1])]], 'params_to_change': {}},
                'data_items': items, 'prediction_items': []})
    return out


def held_out_items(obs):
    """The prediction_items of an obs_data that carry data (a value): its held-out data."""
    return [i for i in (obs.get('prediction_items') or []) if i.get('value') is not None]


def calibration_end(obs):
    """Model time the fitted data reaches: the last sample of its series data_items (the
    protocol's end when it has none). The protocol may run on, to held-out data."""
    ends = [float(i['obs_dt']) * (len(i['value']) - 1) for i in obs['data_items'] if i.get('data_type') == 'series']
    return max(ends) if ends else float(obs['protocol_info']['sim_times'][0][-1])


def calibration_obs(obs):
    """The obs_data CVS0DParamID fits: the data_items only, run to calibration_end."""
    out = copy.deepcopy(obs)
    out['prediction_items'] = []
    out['protocol_info']['sim_times'] = [[calibration_end(obs)]]
    return out


def validation_view(obs):
    """{protocol_info, data_items}: what the calibrated model is scored against: the held-out
    prediction_items when there are any, else the data_items themselves."""
    held = held_out_items(obs)
    return {'protocol_info': obs['protocol_info'], 'data_items': held or obs['data_items']}


def series_from_obs_data(obs):
    """{variable: (t, values)} from the series items of an obs_data dict."""
    out = {}
    for item in obs['data_items']:
        if item.get('data_type') != 'series':
            continue
        var = item['operands'][0].split('/', 1)[1]
        values = np.asarray(item['value'], dtype=float)
        out[var] = (np.arange(len(values)) * float(item['obs_dt']), values)
    return out


def read_params_for_id(path):
    with open(path) as f:
        rows = [{k.strip(): (v or '').strip() for k, v in r.items()} for r in csv.DictReader(f, skipinitialspace=True)]
    return [{'vessel_name': r['vessel_name'], 'param_name': r['param_name'], 'min': float(r['min']),
             'max': float(r['max']), 'name_for_plotting': r.get('name_for_plotting') or r['param_name']}
            for r in rows]


def _calibrate_once(cm, v, obs, params_for_id, start_values, tag):
    """One CVS0DParamID run from ``start_values`` ({variable: value}, or {} for nominal)."""
    from libcuflynx.param_id.paramID import CVS0DParamID

    component = cm.component
    by_var = {p.variable_name: p for p in cm.parameters()}
    cal_end = float(obs['protocol_info']['sim_times'][0][-1])
    dts = [float(i['obs_dt']) for i in obs['data_items'] if i.get('data_type') == 'series']
    dt = min(dts + [float(cm.spec['dt'])])
    model_path, work_dir = cm.model_path, cm.work_dir
    if start_values:
        # start the optimiser somewhere other than the nominal values: a model generated with them
        work_dir = os.path.join(cm.work_dir, f'calibration_{tag}')
        overrides = {harness.parameter_name(by_var[k]): float(val) for k, val in start_values.items()}
        model_path = cm.generate(work_dir, overrides=overrides)
    inp = {
        'model_path': model_path, 'model_type': 'cellml',
        'file_prefix': os.path.splitext(os.path.basename(model_path))[0],
        'param_id_method': v.get('method', 'CMA-ES'), 'sim_time': cal_end,
        'pre_time': float(obs['protocol_info'].get('pre_times', [0.0])[0]), 'dt': dt,
        'solver_info': {'solver': 'CVODE_myokit', 'rtol': 1e-8, 'atol': 1e-10},
        'optimiser_options': dict({'num_calls_to_function': 2000, 'cost_type': 'gaussian_MLE'},
                                  **(v.get('optimiser_options') or {})),
        'DEBUG': False, 'do_ad': bool(v.get('do_ad', False)),
        'param_id_output_dir': os.path.join(work_dir, 'param_id'),
        'resources_dir': os.path.join(work_dir, 'resources'), 'one_rank': True,
    }
    with contextlib.redirect_stdout(io.StringIO()):
        pid = CVS0DParamID.init_from_all_dicts(inp, obs, params_for_id)
        pid.run()
        best = np.asarray(pid.get_best_param_vals(), dtype=float).ravel()
        cost = getattr(getattr(pid, 'param_id', None), 'best_cost', None)
    cal = {p['param_name']: float(val) for p, val in zip(params_for_id, best)}
    return cal, inp['param_id_method'], (float(cost) if cost is not None and np.isfinite(cost) else None)


def _evaluate(cm, v, val_obs, calibrated, cal_end):
    """Scores of the calibrated model against the validation obs_data: series items by
    ``metric`` (calibration and prediction windows), constant items as |model - value| / std
    using libcuflynx's own operation functions."""
    from libcuflynx.param_id.operation_funcs import get_operation_funcs_dict_for_mode

    by_var = {p.variable_name: p for p in cm.parameters()}
    series = series_from_obs_data(val_obs)
    constants = [i for i in val_obs['data_items'] if i.get('data_type') == 'constant']
    t_end = float(val_obs['protocol_info']['sim_times'][0][-1])
    dts = [float(i['obs_dt']) for i in val_obs['data_items'] if i.get('data_type') == 'series']
    wanted = sorted(set(series) | {op.split('/', 1)[1] for i in constants for op in i['operands']})
    pre = float(val_obs['protocol_info'].get('pre_times', [0.0])[0])   # e.g. reach periodic steady state first
    run_spec = dict(cm.spec, sim_time=t_end, pre_time=pre, dt=min(dts + [float(cm.spec['dt'])]))
    helper = harness.simulation_helper(cm.model_path, run_spec, solver_info={'rtol': 1e-8, 'atol': 1e-10})
    overrides = {'parameters/' + harness.parameter_name(by_var[k]): val for k, val in calibrated.items()}
    tm, out = harness.run(helper, wanted, params=overrides)

    metric = v.get('metric', 'log_rmse')
    pred_window = v.get('prediction_window')
    cal_scores, pred_scores, z_scores = {}, {}, {}
    for var, (t, d) in series.items():
        m_at = np.interp(t, tm, out[var])
        in_cal = t <= cal_end + 1e-9
        cal_scores[var] = score(metric, m_at[in_cal], d[in_cal])
        if pred_window:
            in_pred = (t > pred_window[0] + 1e-9) & (t <= pred_window[1] + 1e-9)
            pred_scores[var] = score(metric, m_at[in_pred], d[in_pred])
    ops = get_operation_funcs_dict_for_mode('numpy')
    for item in constants:
        args = [out[op.split('/', 1)[1]] for op in item['operands']]
        op = item.get('operation')
        model_value = float(ops[op](*args)) if op else float(args[0][-1])
        z_scores[item['data_item_name']] = {'model': model_value, 'data': float(item['value']),
                                            'z': abs(model_value - float(item['value'])) / float(item['std'])}
    return tm, out, series, cal_scores, pred_scores, z_scores


def run(cm, v, plot_path):
    from cam_testing.checks import FAILED, PASSED, Result

    component = cm.component
    inst = cm.instance
    vdir = component.dir

    def path(key, default):
        # the spec's file reference (relative to the version directory), else the instance's convention
        return os.path.join(vdir, v[key]) if v.get(key) else default
    obs_path = path('obs_data', inst.obs_data_path)
    with open(obs_path) as f:
        obs = json.load(f)
    val_obs = validation_view(obs)
    pfi_path = path('params_for_id', inst.params_for_id_path)
    params_for_id = read_params_for_id(pfi_path)
    cal_end = calibration_end(obs)
    metric = v.get('metric', 'log_rmse')
    threshold = float(v.get('threshold', 0.5))
    z_threshold = float(v.get('z_threshold', 2.0))
    expected = v.get('expected_parameters') or {}
    expected_tol = float(v.get('expected_rtol', 0.05))
    compare_abs = v.get('expected_compare') == 'abs'

    # one calibration per start (default: from the nominal values)
    starts = v.get('starts') or [v.get('initial_parameters') or {}]
    runs, problems, figs = [], [], []
    for n, start_values in enumerate(starts):
        calibrated, method, cost = _calibrate_once(cm, v, calibration_obs(obs), params_for_id, start_values, n)
        tm, out, series, cal_s, pred_s, z_s = _evaluate(cm, v, val_obs, calibrated, cal_end)
        label = ', '.join(f'{k}={val:g}' for k, val in start_values.items()) or 'nominal'
        runs.append({'start': start_values, 'calibrated_parameters': calibrated, 'cost': cost, 'calibration_scores': cal_s,
                     'prediction_scores': pred_s, 'constant_items': z_s})
        for var, sc in pred_s.items():
            if sc > threshold:
                problems.append(f'start {label}: prediction {metric} for {var} {sc:.3g} > {threshold}')
        for name, z in z_s.items():
            if z['z'] > z_threshold:
                problems.append(f'start {label}: {name} = {z["model"]:.4g} vs {z["data"]:.4g} ({z["z"]:.2f} std)')
        for par, want in expected.items():
            got = calibrated.get(par)
            if got is None:
                continue
            g = abs(got) if compare_abs else got
            if abs(g - float(want)) > expected_tol * max(abs(float(want)), 1e-12):
                problems.append(f'start {label}: {par} = {got:.4g}, expected {"|" + par + "| = " if compare_abs else ""}{want}')
        if series:
            title = f'{component.label}: calibrated to t ≤ {cal_end:g}'
            if v.get('prediction_window'):
                title += f', predicting {v["prediction_window"][0]:g} < t ≤ {v["prediction_window"][1]:g}'
            figs.append(plots.plot_model_vs_data(
                plot_path(component, f'validation_calibrate{"" if len(starts) == 1 else f"_start{n}"}', inst), tm,
                {k: out[k] for k in series}, {k: ts for k, (ts, _) in series.items()},
                {k: d for k, (_, d) in series.items()}, cm.units(),
                title + ('' if len(starts) == 1 else f' (start {label})'), split=cal_end))

    validated = {}
    for r in runs:
        for k, val in r['calibrated_parameters'].items():
            validated.setdefault(k, []).append(val)
    files = [os.path.relpath(p, vdir) for p in (obs_path, pfi_path)]
    metrics = {'metric': metric, 'threshold': threshold, 'method': method, 'source': v.get('source', ''),
               'files': sorted(set(files), key=files.index), 'runs': runs,
               'calibrated_parameters': runs[0]['calibrated_parameters'] if len(runs) == 1 else
               {k: vals for k, vals in validated.items()},
               'validated_values': {k: vals for k, vals in validated.items()}}
    if expected:
        metrics['expected_parameters'] = expected
    summary = '; '.join(
        f"start {', '.join(f'{k}={x:g}' for k, x in r['start'].items()) or 'nominal'} -> "
        + ', '.join(f'{k}={x:.4g}' for k, x in r['calibrated_parameters'].items())
        + (f" (prediction {metric} {max(r['prediction_scores'].values()):.3g})" if r['prediction_scores'] else '')
        for r in runs)
    write_calibrated(cm, runs[0]['calibrated_parameters'], metrics, not problems, summary)
    metrics['written'] = [os.path.basename(inst.calibrated_parameters_path), os.path.basename(inst.calibration_path)]
    if problems:
        return Result('validation_test_calibrate', FAILED, '; '.join(problems[:3]), metrics, figs, problems)
    return Result('validation_test_calibrate', PASSED, summary, metrics, figs)


def write_calibrated(cm, calibrated, metrics, passed, summary):
    '''<instance>_calibrated_parameters.csv (the instance's rows, calibrated values replaced, the
    reference noting the calibration) and <instance>_calibration.json (the summary).'''
    import datetime
    from cam_testing.library import write_parameters
    inst = cm.instance
    date = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d')
    params = []
    for p in cm.parameters():
        p = copy.copy(p)
        if p.variable_name in calibrated:
            p.value = f'{calibrated[p.variable_name]:.10g}'
            p.data_reference = (f'calibrated to {os.path.basename(inst.obs_data_path)} '
                                f'({metrics.get("method")}, {date})')
            p.sourced = 'no'
        params.append(p)
    write_parameters(inst.calibrated_parameters_path, params)
    summary_json = {'instance': inst.name, 'version': cm.component.key, 'date': date, 'passed': passed,
                    'summary': summary, 'method': metrics.get('method'), 'files': metrics.get('files'),
                    'cost': (metrics.get('runs') or [{}])[0].get('cost'),
                    'calibrated_parameters': calibrated, 'runs': metrics.get('runs'),
                    'metric': metrics.get('metric'), 'threshold': metrics.get('threshold')}
    with open(inst.calibration_path, 'w') as f:
        json.dump(summary_json, f, indent=1, default=float)
        f.write('\n')


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Make libcuflynx obs_data from tabular data')
    sub = parser.add_subparsers(dest='cmd', required=True)
    fc = sub.add_parser('from-csv')
    fc.add_argument('--csv', required=True)
    fc.add_argument('--time-column', default='t')
    fc.add_argument('--time-offset', type=float, default=0.0)
    fc.add_argument('--variables', nargs='+', required=True, help='model_var=csv_column')
    fc.add_argument('--window', nargs=2, type=float, required=True)
    fc.add_argument('--noise', nargs=2, default=['relative', '0.1'], help='relative|absolute VALUE')
    fc.add_argument('--name', help='obs_data_name (the instance the file belongs to)')
    fc.add_argument('--validation-window', nargs=2, type=float,
                    help='also add this window of the data as held-out prediction_items')
    fc.add_argument('--out', required=True)
    args = parser.parse_args(argv)
    spec = {'data': os.path.abspath(args.csv), 'time_column': args.time_column, 'time_offset': args.time_offset,
            'variables': dict(x.split('=', 1) for x in args.variables)}
    t, data = load_data('/', spec)
    obs = obs_data_from_series(t, data, args.window, args.noise, name=args.name)
    if args.validation_window:
        held = obs_data_from_series(t, data, args.validation_window, args.noise)
        obs['protocol_info']['sim_times'] = held['protocol_info']['sim_times']
        obs['prediction_items'] = [
            {'data_item_name': f"{i['trace_name_for_plotting']}_validation", 'operands': i['operands'],
             'unit': i['unit'], 'trace_name_for_plotting': i['trace_name_for_plotting'], 'experiment_idx': 0,
             'data_type': 'series', 'value': i['value'], 'std': i['std'], 'obs_dt': i['obs_dt']}
            for i in held['data_items']]
    with open(args.out, 'w') as f:
        json.dump(obs, f, indent=2)
    print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
