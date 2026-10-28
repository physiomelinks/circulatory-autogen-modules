"""
A parameterisation's validation and calibration results beside the publication figures their data came
from: one plot per data item, or one overlaid plot of the items a source figure holds; the source
figures paired with them by source_figures.json "items" (or, with no "items", the only figure pairing with
every item); source_figures.json validation and the version page's two columns.
No model is generated except where noted.
"""
import json
import os

import numpy as np

from cam_testing import checks, library, plots, report
from cam_testing.library import load_version
from cam_testing.paths import PACKAGE_CHECKOUT


# ---- one plot per item --------------------------------------------------------------------------

def test_model_vs_data_writes_one_figure_per_item(tmp_path):
    t = np.linspace(0, 1, 50)
    model = {'zeta_normal': np.sin(t), 'A_eff/normal': np.cos(t)}
    data = {k: v[::5] for k, v in model.items()}
    saved = plots.plot_model_vs_data(lambda item: str(tmp_path / f'p__{plots.item_file_name(item)}.png'),
                                     t, model, {k: t[::5] for k in model}, data, {'zeta_normal': '-'}, 'title')
    assert [n for n, _ in saved] == ['zeta_normal', 'A_eff/normal']
    assert [os.path.basename(p) for _, p in saved] == ['p__zeta_normal.png', 'p__A_eff_normal.png']
    assert all(os.path.getsize(p) > 0 for _, p in saved)
    paths, metrics = checks.item_plots(saved)
    assert paths == [p for _, p in saved] and metrics == {'plot_items': ['zeta_normal', 'A_eff/normal']}


def test_item_plot_path_names_the_parameterisation_test_and_item():
    v = load_version('chamber', 'Liang2009NoVisco_v01')
    inst = v.parameterisation('lv')
    path = checks.item_plot_path(v, 'validation_baseline', inst)('lv pressure')
    assert path == os.path.join(v.plots_dir, 'lv__validation_baseline__lv_pressure.png')


# ---- pairing ------------------------------------------------------------------------------------

FIG4 = {'file': 'f4.png', 'source': 'mynard2012simple; Fig. 4', 'items': ['A_eff_normal']}
FIG5 = {'file': 'f5.png', 'source': 'mynard2012simple; Fig. 5', 'items': ['zeta_normal']}
FIG4A = {'file': 'f4a.png', 'source': 'liang2009multi; Fig. 4a'}


def test_explicit_items_pair_figures_and_plots_unmatched_follow():
    rows = report.figure_rows(['z.png', 'a.png', 'q.png'], ['zeta_normal', 'A_eff_normal', 'aortic_flow_liang'],
                              [FIG4, FIG5, FIG4A])
    assert [([f['file'] for f in r['figures']], [p['file'] for p in r['plots']]) for r in rows] == [
        (['f5.png'], ['z.png']), (['f4.png'], ['a.png']),     # paired, in plot order
        ([], ['q.png']),                                     # an item no figure lists
        (['f4a.png'], [])]                                   # a figure without a plot here


def test_a_figure_with_several_items_takes_all_their_plots():
    fig6 = {'file': 'f6.png', 'source': 's', 'items': ['A_eff_pulmonary_normal', 'A_eff_pulmonary_hypertension']}
    rows = report.figure_rows(['n.png', 'h.png'], ['A_eff_pulmonary_normal', 'A_eff_pulmonary_hypertension'], [fig6])
    assert len(rows) == 1 and [p['item'] for p in rows[0]['plots']] == ['A_eff_pulmonary_normal',
                                                                         'A_eff_pulmonary_hypertension']
    # the overlaid plot of both
    rows = report.figure_rows(['nh.png'], [['A_eff_pulmonary_normal', 'A_eff_pulmonary_hypertension']], [fig6])
    assert rows == [{'figures': [fig6], 'plots': [{'file': 'nh.png', 'item': ['A_eff_pulmonary_normal',
                                                                              'A_eff_pulmonary_hypertension']}]}]


def test_the_only_figure_without_items_pairs_every_item():
    one = {'file': 'f7.png', 'source': 's'}
    rows = report.figure_rows(['m.png'], ['A_eff_mitral_normal'], [one])
    assert len(rows) == 1 and rows[0]['figures'] == [one] and rows[0]['plots'][0]['file'] == 'm.png'
    # two starts of a calibration: one item, two plots
    assert len(report.figure_rows(['s0.png', 's1.png'], ['x', 'x'], [one])) == 1
    # a figure overlaying two experiments: both items' plots, or their overlaid plot
    rows = report.figure_rows(['a.png', 'b.png'], ['x', 'y'], [one])
    assert len(rows) == 1 and [p['file'] for p in rows[0]['plots']] == ['a.png', 'b.png']
    assert report.figure_rows(['xy.png'], [['x', 'y']], [one]) == [
        {'figures': [one], 'plots': [{'file': 'xy.png', 'item': ['x', 'y']}]}]
    # two figures without items (nothing says which item is in which), or an unknown item (a result
    # from before plot_items): not paired
    two = [one, dict(one, file='g.png')]
    assert len(report.figure_rows(['a.png'], ['x'], two)) == 3
    assert len(report.figure_rows(['a.png', 'b.png'], ['x', 'y'], two)) == 4
    assert len(report.figure_rows(['a.png'], None, [one])) == 2
    # a figure that lists items turns the fallback off
    assert len(report.figure_rows(['a.png'], ['x'], [dict(one, items=['y'])])) == 2


def test_figure_items():
    one = {'file': 'f.png', 'source': 's'}
    assert library.figure_items([one], ['x', 'y', 'x']) == [['x', 'y']]
    assert library.figure_items([one, one], ['x', 'y']) == [[], []]
    assert library.figure_items([one], ['x', None]) == [[]]
    # an item listed by two figures goes to the first
    assert library.figure_items([dict(one, items=['y', 'x']), dict(one, items=['x', 'z'])], ['x', 'y', 'z']) == [
        ['x', 'y'], ['z']]


def _overlay(tmp_path, figures):
    t = np.linspace(0, 1, 50)
    model = {'x': np.sin(t), 'y': np.cos(t), 'z': t}
    groups = library.figure_items(figures, list(model))
    return plots.plot_model_vs_data(lambda item: str(tmp_path / f'p__{plots.item_file_name(item)}.png'),
                                    t, model, {k: t[::5] for k in model}, {k: v[::5] for k, v in model.items()},
                                    {'x': 'mmHg', 'y': 'mmHg'}, 'title', groups=groups)


def test_the_only_figure_gets_one_overlaid_plot_of_every_item(tmp_path):
    one = {'file': 'f.png', 'source': 's'}
    saved = _overlay(tmp_path, [one])
    assert [n for n, _ in saved] == [['x', 'y', 'z']] and os.path.basename(saved[0][1]) == 'p__x_y_z.png'
    assert os.path.getsize(saved[0][1]) > 0
    paths, metrics = checks.item_plots(saved)
    assert report.figure_rows(paths, metrics['plot_items'], [one]) == [
        {'figures': [one], 'plots': [{'file': saved[0][1], 'item': ['x', 'y', 'z']}]}]


def test_a_figure_listing_two_items_gets_one_overlaid_plot_of_them(tmp_path):
    fig = {'file': 'f.png', 'source': 's', 'items': ['z', 'x']}
    saved = _overlay(tmp_path, [fig])
    assert [n for n, _ in saved] == [['x', 'z'], 'y']     # at its first item; y not in the figure: alone
    paths, metrics = checks.item_plots(saved)
    rows = report.figure_rows(paths, metrics['plot_items'], [fig])
    assert [([f['file'] for f in r['figures']], [p['item'] for p in r['plots']]) for r in rows] == [
        (['f.png'], [['x', 'z']]), ([], ['y'])]


def test_two_figures_without_items_keep_one_plot_per_item(tmp_path):
    figs = [{'file': 'f.png', 'source': 's'}, {'file': 'g.png', 'source': 's'}]
    saved = _overlay(tmp_path, figs)
    assert [n for n, _ in saved] == ['x', 'y', 'z']
    paths, metrics = checks.item_plots(saved)
    rows = report.figure_rows(paths, metrics['plot_items'], figs)
    assert all(len(r['figures']) + len(r['plots']) == 1 for r in rows) and len(rows) == 5


def test_figures_alone_and_plots_alone():
    assert report.figure_rows([], None, [FIG4A]) == [{'figures': [FIG4A], 'plots': []}]
    assert report.figure_rows(['a.png'], ['x'], []) == [{'figures': [], 'plots': [{'file': 'a.png', 'item': 'x'}]}]
    assert report.figure_rows([], None, []) == []


def test_caption_split_after_its_first_sentence():
    assert report.caption_parts('LV pressure, Fig. 4a. Digitised every 0.01 s.') == {
        'caption_lead': 'LV pressure, Fig. 4a.', 'caption_rest': 'Digitised every 0.01 s.'}
    assert report.caption_parts(None) == {'caption_lead': '', 'caption_rest': ''}


# ---- source_figures.json validation -------------------------------------------------------------

def test_schema_documents_the_source_figure_keys():
    with open(os.path.join(PACKAGE_CHECKOUT, 'modules', 'directory_schema.json')) as f:
        schema = json.load(f)
    keys = schema['levels']['parameterisation']['subdirectories']['source_figures']['index_keys']
    assert {k: keys[k] for k in ('required', 'optional')} == library.SOURCE_FIGURE_KEYS


def test_source_figures_problems():
    keys = library.SOURCE_FIGURE_KEYS
    names = {'zeta_normal', 'A_eff_normal'}
    assert library.source_figures_problems([FIG4, FIG5, FIG4A], names, keys) == []
    problems = library.source_figures_problems(
        [{'file': 'a.png'}, dict(FIG4, items='A_eff_normal'), dict(FIG4, items=['nope']), dict(FIG4A, figure=1)],
        names, keys)
    assert len(problems) == 4
    assert 'missing "source"' in problems[0] and 'list of data item names' in problems[1]
    assert "['nope']" in problems[2] and 'unknown key "figure"' in problems[3]
    assert library.source_figures_problems({}, names, keys) == ['source_figures.json is not a list of figures']


def test_item_names_are_the_obs_data_items():
    inst = load_version('chamber', 'Liang2009NoVisco_v01').parameterisation('lv')
    with open(inst.obs_data_path) as f:
        obs = json.load(f)
    assert {i['data_item_name'] for i in obs['data_items'] + (obs.get('prediction_items') or [])} <= inst.item_names()


# ---- the version page ---------------------------------------------------------------------------

def test_version_page_shows_paper_and_result_columns(monkeypatch):
    v = load_version('chamber', 'Liang2009NoVisco_v01')
    inst = v.parameterisation('lv')
    figures = inst.source_figures()
    real_load = checks.load

    def load(version, test, parameterisation=None):
        if parameterisation is not None and parameterisation.name == 'lv' and test == 'validation_test_baseline':
            return checks.Result(test, checks.PASSED, 'ok', {'plot_items': ['item_a', ['item_b', 'item_c']]},
                                 ['plots/lv__validation_baseline__item_a.png', 'plots/lv__validation_baseline__item_b_item_c.png'])
        return real_load(version, test, parameterisation)
    monkeypatch.setattr(checks, 'load', load)
    monkeypatch.setattr(library.Parameterisation, 'source_figures',
                        lambda self: [dict(figures[0], items=['item_a'])] + figures[1:] if self.name == 'lv'
                        else [])
    ctx = report.version_context(v)
    lv = next(p for p in ctx['components'][0]['parameterisations'] if p['name'] == 'lv')
    rows = lv['baseline']['figure_rows']
    assert rows[0]['figures'][0]['file'] == figures[0]['file'] and rows[0]['plots'][0]['item'] == 'item_a'
    assert rows[1] == {'figures': [], 'plots': [{'file': 'plots/lv__validation_baseline__item_b_item_c.png',
                                                 'item': ['item_b', 'item_c']}]}
    assert lv['baseline']['results_heading'] == 'Model vs data' and lv['calibrate']['results_heading'] == 'Calibration results'
    html = report._env().get_template('version.html').render(**ctx)
    section = html[html.index('id="lv--validation_test_baseline"'):]
    section = section[:section.index('</article>')]
    assert 'From the paper' in section and 'Model vs data' in section
    assert section.count('class="col-paper"') == len(rows) + 1 and section.count('class="col-result"') == len(rows) + 1
    assert 'No source figure is matched to <code>item_b + item_c</code>' in section
    assert 'plot of item_b + item_c for' in section      # an overlaid plot's items
    assert section.index(figures[0]['file']) < section.index('item_a.png') < section.index('item_b_item_c.png')
    assert '<summary>Extraction notes</summary>' in section
    style = report._env().get_template('_style.html').render()
    assert style.count('--paper-col-bg') == 3 and style.count('--result-col-bg') == 3   # light, dark (system), dark
