"""
validation_test_calibrate: calibrate a component with libcuflynx parameter identification on
one data set, then score its prediction of held-out data.

The inputs are ordinary libcuflynx files, committed in the module's validation/ directory, so
the same calibration can be run with libcuflynx directly:

    validation/<component>_calibration_obs_data.json   obs_data fitted by CVS0DParamID
    validation/<component>_params_for_id.csv           parameters calibrated, with bounds
    validation/<component>_validation_obs_data.json    series over the full span, scored on
                                                       prediction_window (held-out data)

The component is generated alone as a single vessel named "mod" (cam_testing.harness.VESSEL),
so operands are "mod/<variable>" and params_for_id vessel_name is "mod".

Spec (validation.calibrate in <name>_tests.yaml):

    status: active
    source: citation / URL
    obs_data: validation/<component>_calibration_obs_data.json
    params_for_id: validation/<component>_params_for_id.csv
    validation_obs_data: validation/<component>_validation_obs_data.json
    prediction_window: [15, 20]           model time scored for the prediction
    initial_parameters: {a: 0.8}          optional optimiser start (default: nominal values)
    expected_parameters: {a: 0.2}         optional: known true values, reported alongside
    method: CMA-ES
    optimiser_options: {num_calls_to_function: 3000, seed: 1}
    metric: log_rmse                      log_rmse | nrmse
    threshold: 0.5                        on the prediction window

The obs_data files for tabular data are made with
    python -m cam_testing.calibrate from-csv --csv validation/data.csv --time-column Year \
        --time-offset 1900 --variables x=Hare y=Lynx --window 0 15 --noise relative 0.25 \
        --out validation/<component>_calibration_obs_data.json
"""
import contextlib
import csv
import io
import json
import os

import numpy as np

from cam_testing import harness, plots


def load_data(module_dir, v):
    with open(os.path.join(module_dir, v['data'])) as f:
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


def obs_data_from_series(t, data, window, noise, unit='dimensionless'):
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
    return {'protocol_info': {'pre_times': [0.0], 'sim_times': [[float(tw[-1])]], 'params_to_change': {}},
            'data_items': items, 'prediction_items': []}


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


def run(cm, v, plot_path):
    from libcuflynx.param_id.paramID import CVS0DParamID
    from cam_testing.checks import FAILED, PASSED, Result

    component = cm.component
    mdir = component.module.dir
    with open(os.path.join(mdir, v['obs_data'])) as f:
        obs = json.load(f)
    with open(os.path.join(mdir, v['validation_obs_data'])) as f:
        val_obs = json.load(f)
    params_for_id = read_params_for_id(os.path.join(mdir, v['params_for_id']))
    cal_end = float(obs['protocol_info']['sim_times'][0][-1])
    obs_dt = min(float(i['obs_dt']) for i in obs['data_items'] if i.get('data_type') == 'series')

    by_var = {p.variable_name: p for p in component.parameters()}
    model_path = cm.model_path
    work_dir = cm.work_dir
    if v.get('initial_parameters'):
        # start the optimiser somewhere other than the nominal values (e.g. a benchmark's
        # deliberately poor initial guess): a model generated with those values
        work_dir = os.path.join(cm.work_dir, 'calibration_start')
        start = {harness.parameter_name(by_var[k]): float(val) for k, val in v['initial_parameters'].items()}
        model_path = harness.generate(component, work_dir, overrides=start)
    prefix = os.path.splitext(os.path.basename(model_path))[0]
    inp = {
        'model_path': model_path, 'model_type': 'cellml', 'file_prefix': prefix,
        'param_id_method': v.get('method', 'CMA-ES'), 'sim_time': cal_end, 'pre_time': 0.0,
        'dt': min(obs_dt, float(cm.spec['dt'])),
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
    calibrated = {p['param_name']: float(val) for p, val in zip(params_for_id, best)}

    # Simulate the calibrated model over the whole validation span with this repo's harness.
    series = series_from_obs_data(val_obs)
    t_end = max(float(ts[-1]) for ts, _ in series.values())
    pred_window = v['prediction_window']
    run_spec = dict(cm.spec, sim_time=t_end, pre_time=0.0, dt=min(obs_dt, float(cm.spec['dt'])))
    helper = harness.simulation_helper(cm.model_path, run_spec, solver_info={'rtol': 1e-8, 'atol': 1e-10})
    overrides = {'parameters/' + harness.parameter_name(by_var[k]): val for k, val in calibrated.items()}
    tm, out = harness.run(helper, list(series), params=overrides)

    metric = v.get('metric', 'log_rmse')
    threshold = float(v.get('threshold', 0.5))
    cal_scores, pred_scores = {}, {}
    for var, (t, d) in series.items():
        m_at = np.interp(t, tm, out[var])
        in_cal, in_pred = t <= cal_end + 1e-9, (t > pred_window[0] + 1e-9) & (t <= pred_window[1] + 1e-9)
        cal_scores[var] = score(metric, m_at[in_cal], d[in_cal])
        pred_scores[var] = score(metric, m_at[in_pred], d[in_pred])

    fig = plots.plot_model_vs_data(plot_path(component, 'validation_calibrate'), tm, out,
                                   {k: ts for k, (ts, _) in series.items()}, {k: d for k, (_, d) in series.items()},
                                   cm.units(), f'{component.label}: calibrated to t ≤ {cal_end:g}, '
                                   f'predicting {pred_window[0]:g} < t ≤ {pred_window[1]:g}', split=cal_end)
    metrics = {'calibrated_parameters': calibrated, 'metric': metric, 'threshold': threshold,
               'calibration_scores': cal_scores, 'prediction_scores': pred_scores,
               'method': inp['param_id_method'], 'source': v.get('source', ''),
               'files': [v['obs_data'], v['params_for_id'], v['validation_obs_data']],
               'validated_values': calibrated}
    if v.get('expected_parameters'):
        metrics['expected_parameters'] = v['expected_parameters']
    worst = max(pred_scores.values())
    status = PASSED if worst <= threshold else FAILED
    return Result('validation_test_calibrate', status,
                  f'prediction {metric} {worst:.3g} {"<=" if status == PASSED else ">"} {threshold} '
                  f'(calibration {max(cal_scores.values()):.3g})', metrics, [fig])


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
    fc.add_argument('--out', required=True)
    args = parser.parse_args(argv)
    spec = {'data': os.path.abspath(args.csv), 'time_column': args.time_column, 'time_offset': args.time_offset,
            'variables': dict(x.split('=', 1) for x in args.variables)}
    t, data = load_data('/', spec)
    obs = obs_data_from_series(t, data, args.window, args.noise)
    with open(args.out, 'w') as f:
        json.dump(obs, f, indent=1)
    print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
