"""
validation_test_calibrate: calibrate a version, at one of its parameterisations, with libcuflynx
parameter identification on that parameterisation's data, then score its prediction of held-out data.

The inputs are ordinary libcuflynx files, committed in the parameterisation's
directory (versions/<version>/parameterisations/<p>/), so the same calibration can be run with libcuflynx
directly:

    <p>_obs_data.json                     obs_data ("obs_data_name": "<p>"): its data_items
                                          are fitted by CVS0DParamID; its prediction_items that
                                          carry a value are the held-out data the calibrated
                                          model is validated against, scored on
                                          prediction_window (without any, the data_items are)
    <p>_params_for_id.csv                 parameters calibrated, with bounds

and the results are written there too, and committed:

    <p>_calibrated_parameters.csv         the parameterisation's parameters with the calibrated values
    <p>_calibration.json                  summary: method, cost, fitted values, scores, date

The version is generated alone as a single vessel named "mod" (cam_testing.harness.VESSEL),
so operands are "mod/<variable>" and params_for_id vessel_name is "mod".

Spec (validation.<p>.calibrate in <module_type>_<version>_verification_config.json; file
references relative to the version directory, defaulting to the parameterisation's files):

    status: active
    source: citation / URL
    obs_data: parameterisations/<p>/<p>_obs_data.json
    params_for_id: parameterisations/<p>/<p>_params_for_id.csv
    prediction_window: [15, 20]           model time scored for series predictions (optional)
    initial_parameters: {a: 0.8}          optional optimiser start (default: nominal values)
    starts: [{a: 1}, {a: -1}]             optional: several starts, each calibrated and checked
    expected_parameters: {a: 0.2}         optional: known true values the fit must recover
    expected_compare: abs                 optional: compare |value| (e.g. symmetric minima)
    expected_rtol: 0.05
    z_threshold: 2.0                      constant data items: |model - value| <= z_threshold * std
    method: CMA-ES
    dt: 0.01                              optional: output step (default: the spec dt; series obs_dt caps it)
    optimiser_options: {num_calls_to_function: 3000, seed: 1}
    metric: log_rmse                      log_rmse | nrmse
    threshold: 0.5                        on the prediction window

The obs_data files for tabular data are made with
    python -m cam_testing.calibrate from-csv --csv parameterisations/<p>/data.csv --time-column Year \
        --time-offset 1900 --variables x=Hare y=Lynx --window 0 15 --noise relative 0.25 \
        --name <p> --out parameterisations/<p>/<p>_obs_data.json
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

from cam_testing import harness, module_array, plots


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


FEATURE_OPERATIONS = ('max', 'min', 'mean')


def feature_items(item, operations=FEATURE_OPERATIONS):
    """Scalar held-out features of a held-out series prediction_item: one constant
    prediction_item per operation (libcuflynx's max / min / mean over the experiment), its value
    the operation on the data and its std the data's at that sample (max, min) or propagated
    (mean: sqrt(sum std^2) / n). They are what SA and emulation can use as features
    (sa_options / emulator_settings include_prediction_items), and are validated like the series."""
    values = np.asarray(item['value'], dtype=float)
    std = np.broadcast_to(np.asarray(item.get('std', 0.0), dtype=float), values.shape)
    var = item['data_item_name'][:-len('_validation')] if item['data_item_name'].endswith('_validation') \
        else item['data_item_name']
    out = []
    for op in operations:
        if op == 'mean':
            value, sd = float(values.mean()), float(np.sqrt(np.sum(std ** 2)) / values.size)
        else:
            k = int(np.argmax(values) if op == 'max' else np.argmin(values))
            value, sd = float(values[k]), float(std[k])
        out.append({'data_item_name': f'{var}_{op}_validation', 'operands': list(item['operands']),
                    'operation': op, 'unit': item['unit'],
                    'trace_name_for_plotting': item.get('trace_name_for_plotting', var),
                    'item_name_for_plotting': f"{item.get('trace_name_for_plotting', var)} ({op})",
                    'experiment_idx': item.get('experiment_idx', 0),
                    'data_type': 'constant', 'value': value, 'std': sd})
    return out


def held_out_items(obs):
    """The prediction_items of an obs_data that carry data (a value): its held-out data."""
    return [i for i in (obs.get('prediction_items') or []) if i.get('value') is not None]


def calibration_end(obs):
    """Model time the fitted data reaches: the last sample of its series data_items (the
    protocol's end when it has none). The protocol may run on, to held-out data."""
    ends = [float(i['obs_dt']) * (len(i['value']) - 1) for i in obs['data_items'] if i.get('data_type') == 'series']
    return max(ends) if ends else float(obs['protocol_info']['sim_times'][0][-1])


def multi_segment(obs):
    """Whether the protocol has more than one (sub-)experiment, or changes parameters."""
    pi = obs['protocol_info']
    return sum(len(s) for s in pi['sim_times']) > 1 or bool(pi.get('params_to_change'))


def calibration_obs(obs):
    """The obs_data CVS0DParamID fits: the data_items only, run to calibration_end. A
    multi-segment protocol (sub-experiments / params_to_change, e.g. a voltage-clamp step) is
    kept as it is: its data are scalar features of its segments."""
    out = copy.deepcopy(obs)
    out['prediction_items'] = []
    if not multi_segment(obs):
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
    dt = min(dts + [float(v.get('dt', cm.spec['dt']))])   # spec dt: output step of a scalar-feature calibration
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
        constants = _constant_features(pid.param_id, best) if multi_segment(obs) else None
    cal = {p['param_name']: float(val) for p, val in zip(params_for_id, best)}
    return cal, inp['param_id_method'], (float(cost) if cost is not None and np.isfinite(cost) else None), constants


def _constant_features(param_id, best):
    """The data_items at ``best``, evaluated by libcuflynx itself over the whole protocol (each
    item on its own sub-experiment, cross-segment operation_kwargs references resolved), as the
    cost was: {data_item_name: value} for constant items, and under the key ``SERIES``
    {data_item_name: (t from the sub-experiment start, model, data, std)} for series items, the
    model aligned to the data's sample times by libcuflynx."""
    from libcuflynx.utilities.obs_data_helpers import obs_item_names

    info = param_id.obs_info
    names = list(obs_item_names(info))
    const_idx, series_idx = {}, {}
    for JJ, dtype in enumerate(info['data_types']):
        if dtype == 'constant':
            const_idx[JJ] = len(const_idx)
        elif dtype == 'series':
            series_idx[JJ] = len(series_idx)
    _, operands_list, _ = param_id.get_cost_obs_and_pred_from_params(best)
    out, k = {SERIES: {}}, 0
    with param_id.accumulating_temp_results():
        for exp_idx, n_sub in enumerate(param_id.protocol_info['num_sub_per_exp']):
            for sub_idx in range(n_sub):
                operands = operands_list[k]
                k += 1
                if operands is None:
                    continue
                with param_id.evaluating_segment(exp_idx, sub_idx):
                    obs = param_id.get_obs_output_dict(operands)
                here = lambda JJ: (int(info['experiment_idxs'][JJ]), int(info['subexperiment_idxs'][JJ])) == (exp_idx, sub_idx)
                for JJ, c in const_idx.items():
                    if here(JJ):
                        out[names[JJ]] = float(obs['const'][c])
                for JJ, c in series_idx.items():
                    if here(JJ):
                        model, data, std = param_id._align_series_to_ground_truth(np.asarray(obs['series'][c], dtype=float), c)
                        t = np.arange(len(data)) * float(info['obs_dt'][c])
                        out[SERIES][names[JJ]] = (t, np.asarray(model, dtype=float), np.asarray(data, dtype=float),
                                                  np.asarray(std, dtype=float))
    return out


SERIES = '__series__'   # the key of _constant_features' series items


def _z(model_value, item):
    return {'model': model_value, 'data': float(item['value']),
            'z': abs(model_value - float(item['value'])) / float(item['std'])}


def _evaluate(cm, v, val_obs, calibrated, cal_end, features=None):
    """Scores of the calibrated model against the validation obs_data: series items by
    ``metric`` (calibration and prediction windows), constant items as |model - value| / std
    using libcuflynx's own operation functions. ``features``: for a multi-segment protocol, the
    constant data_items' values libcuflynx evaluated over the whole protocol
    (_constant_features); they are scored directly, leaving out weight-0 helper items (items
    another one references, e.g. the two steady states of a step difference)."""
    from libcuflynx.param_id.operation_funcs import get_operation_funcs_dict_for_mode

    constants = [i for i in val_obs['data_items'] if i.get('data_type') == 'constant']
    if features is not None:
        z_scores = {i['data_item_name']: _z(features[i['data_item_name']], i) for i in constants
                    if float(i.get('weight', 1.0)) != 0 and i['data_item_name'] in features}
        # series items, by ``metric``, on the samples that carry data (leaving out samples given
        # a much larger std than the rest to exclude them, e.g. those inside a voltage jump)
        metric = v.get('metric', 'log_rmse')
        cal_scores = {}
        for name, (t, model, data, std) in features.get(SERIES, {}).items():
            keep = std < 0.5 * np.max(std) if np.ptp(std) > 0 else np.ones(len(std), bool)
            cal_scores[name] = score(metric, model[keep], data[keep])
        return None, None, features.get(SERIES, {}), cal_scores, {}, z_scores
    by_var = {p.variable_name: p for p in cm.parameters()}
    series = series_from_obs_data(val_obs)
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
        z_scores[item['data_item_name']] = _z(model_value, item)
    return tm, out, series, cal_scores, pred_scores, z_scores


def run(cm, v, plot_path):
    from cam_testing.checks import FAILED, PASSED, Result

    component = cm.component
    inst = cm.parameterisation
    vdir = component.dir

    def path(key, default):
        # the spec's file reference (relative to the version directory), else the parameterisation's convention
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
        calibrated, method, cost, features = _calibrate_once(cm, v, calibration_obs(obs), params_for_id, start_values, n)
        tm, out, series, cal_s, pred_s, z_s = _evaluate(cm, v, val_obs, calibrated, cal_end, features)
        label = ', '.join(f'{k}={val:g}' for k, val in start_values.items()) or 'nominal'
        runs.append({'start': start_values, 'calibrated_parameters': calibrated, 'cost': cost, 'calibration_scores': cal_s,
                     'prediction_scores': pred_s, 'constant_items': z_s})
        for var, sc in pred_s.items():
            if sc > threshold:
                problems.append(f'start {label}: prediction {metric} for {var} {sc:.3g} > {threshold}')
        if features is not None:
            # a multi-segment protocol has no held-out window: its series fit is checked instead
            for var, sc in cal_s.items():
                if sc > threshold:
                    problems.append(f'start {label}: fit {metric} for {var} {sc:.3g} > {threshold}')
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
        if features is not None and series:
            # one panel per series item, model and data against time from its sub-experiment's start
            figs.append(plots.plot_model_vs_data(
                plot_path(component, f'validation_calibrate{"" if len(starts) == 1 else f"_start{n}"}', inst),
                {k: t for k, (t, m, d, sd) in series.items()}, {k: m for k, (t, m, d, sd) in series.items()},
                {k: t for k, (t, m, d, sd) in series.items()},
                {k: d for k, (t, m, d, sd) in series.items()}, {},
                f'{component.label}: calibrated ({label}); time from the start of each sub-experiment'))
        elif series:
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
    '''<p>_calibrated_parameters.csv (the parameterisation's rows, calibrated values replaced,
    the reference noting the calibration) and <p>_calibration.json (the summary).'''
    import datetime
    from cam_testing.library import write_parameters
    inst = cm.parameterisation
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
    summary_json = {'parameterisation': inst.name, 'version': cm.component.key, 'date': date, 'passed': passed,
                    'summary': summary, 'method': metrics.get('method'), 'files': metrics.get('files'),
                    'cost': (metrics.get('runs') or [{}])[0].get('cost'),
                    'calibrated_parameters': calibrated, 'runs': metrics.get('runs'),
                    'metric': metrics.get('metric'), 'threshold': metrics.get('threshold')}
    with open(inst.calibration_path, 'w') as f:
        json.dump(summary_json, f, indent=1, default=float)
        f.write('\n')



# ----------------------------------------------------------------------------------------------
# calibrate apply: copy a parameterisation's calibrated values to every place that uses the version
# ----------------------------------------------------------------------------------------------
#
#     python -m cam_testing.calibrate apply <module_type>/<version> <parameterisation> [--dry-run]
#
# The calibrated parameters (those in <p>_params_for_id.csv, with their values from
# <p>_calibrated_parameters.csv) are written to:
#   1. the version's default parameterisation (row <var>);
#   2. every supermodule version using the version's default parameterisation as a submodule, recursively
#      (row <var>_<submodule path>, outer first: soma/sympathetic rho_M_i_M, a neuron using that
#      soma as "soma" rho_M_soma_i_M), where the row exists;
#   3. every system model whose module array uses the version, or a supermodule chain over it, at
#      its default parameterisation (row <var>_<vessel> or <var>_<vessel>_<submodule path>), where the row exists;
#   4. the monolithic counterparts a supermodule spec declares (supermodule.equivalent: {model,
#      reproduces, output_map}): the flattened vessel of (3) in `model` is mapped through output_map
#      to the vessel of `reproduces` ("soma_SN/w": "SN_soma_i_M/w" -> soma_SN), giving row
#      <var>_<that vessel> in that system model and row <var> in that vessel's version's default parameterisation.
# Each rewritten reference reads "Calibrated (libcuflynx) to parameterisations/<p>/<p>_obs_data.json:
# <note>" (with the version named when it is written elsewhere), and keeps the old value. <note> is the
# obs_data's "calibration_note" (else its obs_data_name); an obs_data "reference_key" (a BibTeX key of
# the version) prefixes the reference and marks the row sourced.


def _system_models(root=None):
    """[{key, dir, vessels, parameters_path}] of every system model under ``root``."""
    from cam_testing import paths
    root = root or paths.roots().system_models_dir
    out = []
    for d, _dirs, files in os.walk(root):
        for f in files:
            if f.endswith('_module_array.json'):
                name = f[:-len('_module_array.json')]
                with open(os.path.join(d, f)) as fh:
                    vessels = json.load(fh)
                out.append({'key': os.path.relpath(d, root).replace(os.sep, '/'), 'dir': d, 'vessels': vessels,
                            'parameters_path': os.path.join(d, f'{name}_parameters.csv')})
    return sorted(out, key=lambda m: m['key'])


def _csv_rows(path):
    """(lines without their line endings, line ending, header) of a parameters CSV."""
    with open(path, newline='') as f:
        raw = f.read()
    eol = '\r\n' if '\r\n' in raw else '\n'
    lines = raw.split(eol)
    return lines, eol, next(csv.reader([lines[0]]))


def _row_values(path, names):
    """{row name: {column: value}} of the rows of ``path`` named in ``names``."""
    if not os.path.isfile(path):
        return {}
    lines, _eol, header = _csv_rows(path)
    out = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        row = dict(zip(header, next(csv.reader([ln]))))
        if row.get('variable_name') in names:
            out[row['variable_name']] = row
    return out


def apply_plan(version, parameterisation_name, versions=None, systems=None):
    """The changes ``apply`` would make: [{path, row, variable, value, old, reference, sourced}].
    ``versions`` (default: the library) and ``systems`` (default: _system_models()) are what is
    searched for users of the version."""
    from cam_testing.library import all_versions, read_parameters
    inst = version.parameterisation(parameterisation_name)
    names = [p['param_name'] for p in read_params_for_id(inst.params_for_id_path)]
    if not os.path.isfile(inst.calibrated_parameters_path):
        raise FileNotFoundError(f'{inst.calibrated_parameters_path}: run the calibration first')
    calibrated = {p.variable_name: p.value for p in read_parameters(inst.calibrated_parameters_path)
                  if p.variable_name in names}
    if os.path.isfile(inst.calibration_path):
        # full precision from the summary (the CSV keeps 10 significant figures)
        with open(inst.calibration_path) as f:
            full = (json.load(f).get('calibrated_parameters') or {})
        calibrated.update({k: repr(float(v)) for k, v in full.items() if k in names and not isinstance(v, list)})
    missing = [n for n in names if n not in calibrated]
    if missing:
        raise KeyError(f'{inst.calibrated_parameters_path} has no row for {missing}')
    with open(inst.obs_data_path) as f:
        obs = json.load(f)
    note = obs.get('calibration_note') or obs.get('obs_data_name') or parameterisation_name
    key = obs.get('reference_key')
    versions = all_versions() if versions is None else versions
    systems = _system_models() if systems is None else systems
    by_key = {v.key: v for v in versions}
    by_key.setdefault(version.key, version)

    def default_of(mt, sub):
        v = by_key.get(f'{mt}/{sub}')
        return v.default_parameterisation_name if v is not None else 'default'

    # (version, submodule path) of every supermodule chain over the version's default parameterisation
    chain, todo = [], [(version.mtype.name, version.name, '')]
    while todo:
        mt, vname, sfx = todo.pop(0)
        for S in versions:
            if not S.is_supermodule:
                continue
            for sub in S.submodules:
                if (sub.get('module_type'), sub.get('module_subtype')) != (mt, vname):
                    continue
                if module_array.parameterisation_of(sub, 'default') != default_of(mt, vname):
                    continue
                path = sub['name'] + ('_' + sfx if sfx else '')
                if all((c[0].key, c[1]) != (S.key, path) for c in chain):
                    chain.append((S, path))
                    todo.append((S.mtype.name, S.name, path))

    targets = []          # (csv path, suffix, where, own)
    targets.append((version.default_parameterisation.parameters_path, '', f'{version.key} default parameterisation', True))
    for S, path in chain:
        targets.append((S.default_parameterisation.parameters_path, '_' + path, f'{S.key} default parameterisation', False))
    flattened = {}        # system key -> [flattened vessel names]
    systems_by_key = {m['key']: m for m in systems}
    chain_by_key = {}
    for S, path in chain:
        chain_by_key.setdefault(S.key, []).append(path)
    for m in systems:
        for ves in m['vessels']:
            vk = f"{ves.get('module_type')}/{ves.get('module_subtype')}"
            p_name = module_array.parameterisation_of(ves, 'default')
            if p_name != default_of(ves.get('module_type'), ves.get('module_subtype')):
                continue
            if vk == version.key:
                flattened.setdefault(m['key'], []).append(ves['name'])
            for path in chain_by_key.get(vk, []):
                flattened.setdefault(m['key'], []).append(f"{ves['name']}_{path}")
    for mkey, fl in flattened.items():
        for F in fl:
            targets.append((systems_by_key[mkey]['parameters_path'], '_' + F, f'system model {mkey}', False))
    # monolithic counterparts, through supermodule.equivalent output maps
    for S in versions:
        for e in ((S.spec.get('supermodule') or {}).get('equivalent') or []):
            rep = systems_by_key.get(e.get('reproduces'))
            if rep is None:
                continue
            omap = e.get('output_map') or {}
            for F in flattened.get(e.get('model'), []):
                mono = sorted({a.split('/', 1)[0] for a, b in omap.items() if b.split('/', 1)[0] == F})
                if len(mono) != 1:
                    continue
                targets.append((rep['parameters_path'], '_' + mono[0], f"system model {rep['key']}", False))
                ves = next((x for x in rep['vessels'] if x['name'] == mono[0]), None)
                mv = by_key.get(f"{ves.get('module_type')}/{ves.get('module_subtype')}") if ves else None
                if mv is not None:
                    targets.append((mv.default_parameterisation.parameters_path, '', f'{mv.key} default parameterisation', False))

    import datetime
    today = datetime.date.today().isoformat()
    obs_rel = os.path.relpath(inst.obs_data_path, version.dir).replace(os.sep, '/')
    plan, seen = [], set()
    for path, sfx, where, own in targets:
        rows = _row_values(path, {n + sfx for n in names})
        for n in names:
            row = rows.get(n + sfx)
            if row is None or (path, n + sfx) in seen:
                continue
            seen.add((path, n + sfx))
            to = obs_rel if own else f'{version.mtype.name} {version.name} {obs_rel}'
            ref = (f'{key}; ' if key else '') + f'Calibrated (libcuflynx) to {to}: {note} (was {row.get("value")}; applied {today})'
            plan.append({'path': path, 'where': where, 'row': n + sfx, 'variable': n, 'value': calibrated[n],
                         'old': row.get('value'), 'reference': ref, 'sourced': 'yes' if key else None})
    return plan


def write_plan(plan):
    """Rewrites the rows of ``plan`` in place, keeping each file's columns and line endings."""
    by_path = {}
    for c in plan:
        by_path.setdefault(c['path'], {})[c['row']] = c
    for path, changes in by_path.items():
        lines, eol, header = _csv_rows(path)
        for k, ln in enumerate(lines[1:], start=1):
            if not ln.strip():
                continue
            row = next(csv.reader([ln]))
            c = changes.get(row[0])
            if c is None:
                continue
            d = dict(zip(header, row))
            d['value'] = c['value']
            d['data_reference'] = c['reference']
            if c['sourced'] and 'sourced' in d:
                d['sourced'] = c['sourced']
            b = io.StringIO()
            csv.writer(b, lineterminator='').writerow([d.get(h, '') for h in header])
            lines[k] = b.getvalue()
        with open(path, 'w', newline='') as f:
            f.write(eol.join(lines))


def apply_command(version_key, parameterisation_name, dry_run=False, versions=None, systems=None, out=print):
    from cam_testing import paths
    from cam_testing.library import version_by_key
    version = version_by_key(version_key)
    plan = apply_plan(version, parameterisation_name, versions=versions, systems=systems)
    for c in plan:
        rel = os.path.relpath(c['path'], paths.roots().repo_root)
        out(f"{'would set' if dry_run else 'set'} {rel}: {c['row']} {c['old']} -> {c['value']}")
    if not dry_run:
        write_plan(plan)
    out(f'{len(plan)} row(s){" (dry run: nothing written)" if dry_run else ""}')
    return plan

def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Make libcuflynx obs_data from tabular data; apply calibrations')
    sub = parser.add_subparsers(dest='cmd', required=True)
    fc = sub.add_parser('from-csv')
    fc.add_argument('--csv', required=True)
    fc.add_argument('--time-column', default='t')
    fc.add_argument('--time-offset', type=float, default=0.0)
    fc.add_argument('--variables', nargs='+', required=True, help='model_var=csv_column')
    fc.add_argument('--window', nargs=2, type=float, required=True)
    fc.add_argument('--noise', nargs=2, default=['relative', '0.1'], help='relative|absolute VALUE')
    fc.add_argument('--name', help='obs_data_name (the parameterisation the file belongs to)')
    fc.add_argument('--validation-window', nargs=2, type=float,
                    help='also add this window of the data as held-out prediction_items')
    fc.add_argument('--validation-features', nargs='*', default=list(FEATURE_OPERATIONS),
                    help='scalar features of each held-out series added as prediction_items '
                         f'(default {" ".join(FEATURE_OPERATIONS)}; none: --validation-features)')
    fc.add_argument('--out', required=True)
    ap = sub.add_parser('apply', help="copy a parameterisation's calibrated values into the version's "
                                      'default parameterisation, its supermodules, monolithic counterparts and system models')
    ap.add_argument('version', help='<module_type>/<version>')
    ap.add_argument('parameterisation', help="the parameterisation's name (its directory under "
                    'parameterisations/)')
    ap.add_argument('--dry-run', action='store_true', help='print the changes without writing them')
    args = parser.parse_args(argv)
    if args.cmd == 'apply':
        apply_command(args.version, args.parameterisation, dry_run=args.dry_run)
        return
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
        if args.validation_features:
            obs['prediction_items'] += [f for i in list(obs['prediction_items'])
                                        for f in feature_items(i, args.validation_features)]
    with open(args.out, 'w') as f:
        json.dump(obs, f, indent=2)
    print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
