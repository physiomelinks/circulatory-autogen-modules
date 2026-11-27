"""
The calibration figures of constant data items (calibrate.constant_item_curves, plots.plot_constant_items):
an I-V experiment's items become one curve against the voltage it steps, a two-point ratio experiment's
ratio a point, its two helper currents nothing, and curves a source figure lists together share a figure.
"""
from cam_testing import calibrate, plots


def _item(name, e, s, value, std=1.0, **kw):
    return dict({'data_item_name': name, 'data_type': 'constant', 'operation': 'mean', 'operands': ['mod/I'],
                 'unit': 'A', 'value': value, 'std': std, 'weight': 1.0, 'experiment_idx': e,
                 'subexperiment_idx': s}, **kw)


V = [0.0, -0.1, 0.05]          # not in order: the curve is sorted by the input
OBS = {
    'protocol_info': {
        'pre_times': [0.0] * 3, 'sim_times': [[0.1] * 3, [0.1] * 2, [0.1] * 3],
        'experiment_labels': ['iv', 'k05', 'iv_low_na'],
        'params_to_change': {'mod/u': [V, [-0.05, -0.05], V], 'mod/c': [[1.0] * 3, [0.2, 9.0], [1.0] * 3],
                             'mod/na': [[100.0] * 3, [100.0] * 2, [10.0] * 3]}},
    'data_items': [_item(f'iv_{s}', 0, s, -s) for s in range(3)]
    + [_item('half', 1, 0, 0.0, std=1e15), _item('sat', 1, 1, 0.0, std=1e15),
       _item('ratio', 1, 1, 0.5, std=0.05, operation='division', operands=[''], unit='dimensionless',
             operation_kwargs={'x1': 'half', 'x2': 'sat'})]
    + [_item(f'low_{s}', 2, s, -s / 2) for s in range(3)],
}
Z = {i['data_item_name']: {'model': float(i['value']) + 0.1} for i in OBS['data_items']}


def test_curves_and_points():
    curves, points = calibrate.constant_item_curves(OBS, Z)
    assert [c['name'] for c in curves] == ['iv: I', 'iv_low_na: I']
    iv = curves[0]
    assert iv['input'] == 'mod/u' and iv['x'] == sorted(V)
    assert iv['items'] == ['iv_1', 'iv_0', 'iv_2'] and iv['data'] == [-1.0, 0.0, -2.0]
    assert iv['model'] == [-0.9, 0.1, -1.9] and iv['held'] == {'mod/c': 1.0, 'mod/na': 100.0}
    assert curves[1]['held']['mod/na'] == 10.0
    # the ratio is a point; the currents it divides are only its operands
    assert [p['item'] for p in points] == ['ratio']
    assert points[0]['std'] == 0.05 and points[0]['unit'] == 'dimensionless'


def test_unscored_or_trace_inputs_are_points():
    # an item left out of the scores (weight 0) is not plotted; a protocol_traces input is no axis
    z = {k: v for k, v in Z.items() if k != 'iv_2'}
    obs = dict(OBS, protocol_info=dict(OBS['protocol_info'], params_to_change={
        **OBS['protocol_info']['params_to_change'], 'mod/u': [['a', 'b', 'c'], [-0.05, -0.05], V]}))
    curves, points = calibrate.constant_item_curves(obs, z)
    assert [c['name'] for c in curves] == ['iv_low_na: I']
    assert [p['item'] for p in points] == ['iv_0', 'iv_1', 'ratio']


def test_figures(tmp_path):
    curves, points = calibrate.constant_item_curves(OBS, Z)
    path_of = lambda name: str(tmp_path / f'{plots.item_file_name(name)}.png')
    # each curve alone, and the point figure
    saved = plots.plot_constant_items(path_of, curves, points, 'title')
    assert [items for items, _ in saved] == [curves[0]['items'], curves[1]['items'], ['ratio']]
    assert [p.split('/')[-1] for _, p in saved] == ['iv__I.png', 'iv_low_na__I.png', 'constant_items.png']
    # a source figure listing both experiments' items overlays the two curves
    group = curves[0]['items'] + curves[1]['items']
    saved = plots.plot_constant_items(path_of, curves, [], 'title', groups=[group])
    assert len(saved) == 1 and saved[0][0] == group
    assert saved[0][1].endswith('iv__I_iv_low_na__I.png')
    assert all((tmp_path / p.split('/')[-1]).stat().st_size > 0 for _, p in saved)
