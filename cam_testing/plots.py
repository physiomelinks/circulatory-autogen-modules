"""
Static PNG figures for the module reports (matplotlib, headless).

Colours follow the reference data-viz palette: categorical slots in fixed order for
identity, a single-hue blue ramp (ordinal steps 250 -> 700) for parameter sweeps.
"""
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE = '#fcfcfb'
TEXT_PRIMARY = '#0b0b0b'
TEXT_SECONDARY = '#52514e'
GRID = '#e1e0d9'
MUTED = '#8c8a83'
CATEGORICAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
BLUE_ORDINAL = ['#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab',
                '#184f95', '#104281', '#0d366b']
DATA_COLOUR = '#52514e'

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'axes.edgecolor': GRID, 'axes.labelcolor': TEXT_SECONDARY, 'axes.titlecolor': TEXT_PRIMARY,
    'axes.titlesize': 10, 'axes.labelsize': 9, 'axes.titleweight': 'semibold',
    'xtick.color': TEXT_SECONDARY, 'ytick.color': TEXT_SECONDARY,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'legend.frameon': False,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.6,
    'font.family': 'sans-serif',
})


def _grid(n, width=3.4, height=2.4, max_cols=3):
    cols = min(max_cols, max(n, 1))
    rows = int(np.ceil(max(n, 1) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(width * cols, height * rows), squeeze=False,
                             constrained_layout=True)
    flat = axes.ravel()
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig, flat[:n]


def _save(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def ramp(n):
    if n <= 1:
        return [BLUE_ORDINAL[4]]
    idx = np.linspace(0, len(BLUE_ORDINAL) - 1, n).round().astype(int)
    return [BLUE_ORDINAL[i] for i in idx]


def plot_outputs(path, t, outputs, units, title):
    '''One small panel per output variable: single series each, so no legend needed.'''
    names = list(outputs)
    fig, axes = _grid(len(names))
    for ax, name in zip(axes, names):
        ax.plot(t, outputs[name], color=CATEGORICAL[0])
        ax.set_title(name)
        ax.set_xlabel('time [s]')
        ax.set_ylabel(units.get(name, ''))
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11)
    return _save(fig, path)


def plot_sweep(path, sweeps, output, units, title):
    '''
    sweeps: {param_name: [(value, t, y or None)]}. One panel per swept parameter; line
    shade encodes the parameter value (light = low, dark = high); failed runs are omitted
    and noted in the panel.
    '''
    fig, axes = _grid(len(sweeps))
    for ax, (param, runs) in zip(axes, sweeps.items()):
        colours = ramp(len(runs))
        failed = 0
        for colour, (value, t, y) in zip(colours, runs):
            if y is None:
                failed += 1
                continue
            ax.plot(t, y, color=colour, label=f'{value:.3g}')
        ax.set_title(param)
        ax.set_xlabel('time [s]')
        ax.set_ylabel(f'{output} [{units.get(output, "")}]')
        ax.legend(title='value', loc='best', title_fontsize=8)
        if failed:
            ax.text(0.02, 0.02, f'{failed} run(s) failed', transform=ax.transAxes,
                    color=TEXT_PRIMARY, fontsize=8)
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11)
    return _save(fig, path)


def plot_convergence(path, dts, errors, order, observed, title):
    '''Log-log self-convergence error against step size, with the scheme's ideal slope.'''
    fig, ax = plt.subplots(figsize=(4.6, 3.4), constrained_layout=True)
    dts = np.asarray(dts, dtype=float)
    errors = np.asarray(errors, dtype=float)
    ok = errors > 0
    ax.loglog(dts[ok], errors[ok], color=CATEGORICAL[0], marker='o', markersize=5, label='measured')
    if ok.any():
        ref = errors[ok][0] * (dts[ok] / dts[ok][0]) ** order
        ax.loglog(dts[ok], ref, color=MUTED, linestyle='--', linewidth=1.2, label=f'slope {order}')
    ax.set_xlabel('step size dt [s]')
    ax.set_ylabel('max normalised difference')
    ax.set_title(f'{title}\nobserved order {observed:.2f}' if observed is not None else title)
    ax.legend(loc='lower right')
    return _save(fig, path)


def plot_timestep_traces(path, runs, output, units, title):
    '''runs: [(dt, t, y)], coarse -> fine; darker = finer step.'''
    fig, ax = plt.subplots(figsize=(4.6, 3.4), constrained_layout=True)
    for colour, (dt, t, y) in zip(ramp(len(runs)), runs):
        ax.plot(t, y, color=colour, label=f'dt={dt:g}')
    ax.set_xlabel('time [s]')
    ax.set_ylabel(f'{output} [{units.get(output, "")}]')
    ax.set_title(title)
    ax.legend(loc='best')
    return _save(fig, path)


def plot_model_vs_data(path, t_model, model, t_data, data, units, title):
    '''One panel per compared variable: model line, data as dots.'''
    names = list(model)
    fig, axes = _grid(len(names))
    for ax, name in zip(axes, names):
        ax.plot(t_model, model[name], color=CATEGORICAL[0], label='model')
        ax.plot(t_data[name], data[name], linestyle='none', marker='o', markersize=4,
                color=DATA_COLOUR, label='data')
        ax.set_title(name)
        ax.set_xlabel('time [s]')
        ax.set_ylabel(units.get(name, ''))
        ax.legend(loc='best')
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11)
    return _save(fig, path)
