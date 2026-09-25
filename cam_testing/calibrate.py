"""
validation_test_calibrate: calibrate a component to part of a data set with libcuflynx
parameter identification, then check its prediction of the held-out part.

Spec (validation.calibrate in <name>_tests.yaml):

    status: active
    source: citation / URL
    data: validation/<file>.csv          time column + one column per compared output
    time_column: Year
    time_offset: 1900                     model t = data time - offset
    variables: {x: Hare, y: Lynx}         model output -> csv column
    calibration_window: [0, 15]           model time; data must be evenly spaced from 0
    prediction_window: [15, 20]
    noise: {type: relative, value: 0.25}  std used in the calibration cost
    params_for_id:                        parameters to calibrate (variable names) and bounds
      - {param: alpha, min: 0.1, max: 2.0}
    method: CMA-ES
    optimiser_options: {num_calls_to_function: 3000}
    metric: log_rmse                      log_rmse | nrmse
    threshold: 0.5                        on the prediction window
"""
import contextlib
import csv
import io
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


def _obs_data(t, data, v, window):
    mask = (t >= window[0] - 1e-9) & (t <= window[1] + 1e-9)
    tw = t[mask]
    if abs(tw[0]) > 1e-9:
        raise ValueError('calibration data must start at model time 0 (set time_offset)')
    dts = np.diff(tw)
    if not np.allclose(dts, dts[0]):
        raise ValueError('calibration data must be evenly spaced')
    noise = v.get('noise') or {'type': 'relative', 'value': 0.1}
    items = []
    for model_var, series in data.items():
        values = series[mask]
        if noise['type'] == 'relative':
            std = np.maximum(np.abs(values) * float(noise['value']), 1e-12)
        else:
            std = np.full_like(values, float(noise['value']))
        name = harness.output_name(model_var)
        items.append({
            'data_item_name': name, 'trace_name_for_plotting': model_var, 'data_type': 'series',
            'operation': None, 'operands': [name], 'unit': 'dimensionless', 'cost_type': 'gaussian_MLE',
            'weight': 1.0, 'value': values.tolist(), 'std': std.tolist(), 'obs_dt': float(dts[0]),
            'experiment_idx': 0, 'subexperiment_idx': 0,
        })
    obs = {'protocol_info': {'pre_times': [0.0], 'sim_times': [[float(tw[-1])]], 'params_to_change': {}},
           'data_items': items, 'prediction_items': []}
    return obs, float(dts[0]), float(tw[-1])


def run(cm, v, plot_path):
    from libcuflynx.param_id.paramID import CVS0DParamID
    from cam_testing.checks import FAILED, PASSED, Result

    component = cm.component
    t, data = load_data(component.module.dir, v)
    cal_window = v['calibration_window']
    pred_window = v['prediction_window']
    obs, obs_dt, cal_end = _obs_data(t, data, v, cal_window)

    by_var = {p.variable_name: p for p in component.parameters()}
    params_for_id = [{'vessel_name': harness.VESSEL, 'param_name': p['param'], 'min': float(p['min']),
                      'max': float(p['max']), 'name_for_plotting': p['param']} for p in v['params_for_id']]
    model_path = cm.model_path
    prefix = os.path.splitext(os.path.basename(model_path))[0]
    out_dir = os.path.join(cm.work_dir, 'param_id')
    inp = {
        'model_path': model_path, 'model_type': 'cellml', 'file_prefix': prefix,
        'param_id_method': v.get('method', 'CMA-ES'), 'sim_time': cal_end, 'pre_time': 0.0,
        'dt': min(obs_dt, float(cm.spec['dt'])),
        'solver_info': {'solver': 'CVODE_myokit', 'rtol': 1e-8, 'atol': 1e-10},
        'optimiser_options': dict({'num_calls_to_function': 2000, 'cost_type': 'gaussian_MLE'},
                                  **(v.get('optimiser_options') or {})),
        'DEBUG': False, 'param_id_output_dir': out_dir,
        'resources_dir': os.path.join(cm.work_dir, 'resources'), 'one_rank': True,
    }
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        pid = CVS0DParamID.init_from_all_dicts(inp, obs, params_for_id)
        pid.run()
        best = np.asarray(pid.get_best_param_vals(), dtype=float).ravel()
    calibrated = {p['param']: float(val) for p, val in zip(v['params_for_id'], best)}

    # Simulate the calibrated model over both windows with this repo's harness.
    t_end = float(pred_window[1])
    run_spec = dict(cm.spec, sim_time=t_end, pre_time=0.0, dt=min(obs_dt, float(cm.spec['dt'])))
    helper = harness.simulation_helper(model_path, run_spec, solver_info={'rtol': 1e-8, 'atol': 1e-10})
    overrides = {'parameters/' + harness.parameter_name(by_var[k]): val for k, val in calibrated.items()}
    tm, out = harness.run(helper, list(v['variables']), params=overrides)

    metric = v.get('metric', 'log_rmse')
    threshold = float(v.get('threshold', 0.5))
    in_cal = (t >= cal_window[0]) & (t <= cal_window[1])
    in_pred = (t > pred_window[0]) & (t <= pred_window[1])
    cal_scores, pred_scores = {}, {}
    for var in v['variables']:
        m_at = np.interp(t, tm, out[var])
        cal_scores[var] = score(metric, m_at[in_cal], data[var][in_cal])
        pred_scores[var] = score(metric, m_at[in_pred], data[var][in_pred])

    fig = plots.plot_model_vs_data(plot_path(component, 'validation_calibrate'), tm, out,
                                   {k: t for k in data}, data, cm.units(),
                                   f'{component.label}: calibrated to t ≤ {cal_window[1]:g}, '
                                   f'predicting {pred_window[0]:g} < t ≤ {pred_window[1]:g}',
                                   split=cal_window[1])
    metrics = {'calibrated_parameters': calibrated, 'metric': metric, 'threshold': threshold,
               'calibration_scores': cal_scores, 'prediction_scores': pred_scores,
               'method': inp['param_id_method'], 'source': v.get('source', ''),
               'validated_values': calibrated}
    worst = max(pred_scores.values())
    status = PASSED if worst <= threshold else FAILED
    return Result('validation_test_calibrate', status,
                  f'prediction {metric} {worst:.3f} {"<=" if status == PASSED else ">"} {threshold} '
                  f'(calibration {max(cal_scores.values()):.3f})', metrics, [fig])
