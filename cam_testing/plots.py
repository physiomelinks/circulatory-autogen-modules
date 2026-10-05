"""
Static PNG figures for the module reports (matplotlib, headless).

Colours follow the reference data-viz palette: categorical slots in fixed order for
identity, a single-hue blue ramp (ordinal steps 250 -> 700) for parameter sweeps.
"""
import math
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

STYLE = {
    'text.usetex': False, 'mathtext.fontset': 'dejavusans',
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'axes.edgecolor': GRID, 'axes.labelcolor': TEXT_SECONDARY, 'axes.titlecolor': TEXT_PRIMARY,
    'axes.titlesize': 10, 'axes.labelsize': 9, 'axes.titleweight': 'semibold',
    'xtick.color': TEXT_SECONDARY, 'ytick.color': TEXT_SECONDARY,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'legend.frameon': False,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.6,
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
}
TIME_LABEL = 'time [s]'


def styled(fn):
    '''Apply this module's style per figure: libcuflynx changes the global rcParams (usetex).'''
    import functools

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with plt.rc_context(STYLE):
            return fn(*args, **kwargs)
    return wrapper


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


@styled
def plot_outputs(path, t, outputs, units, title):
    '''One small panel per output variable: single series each, so no legend needed.'''
    names = list(outputs)
    fig, axes = _grid(len(names))
    for ax, name in zip(axes, names):
        ax.plot(t, outputs[name], color=CATEGORICAL[0])
        ax.set_title(name)
        ax.set_xlabel(TIME_LABEL)
        ax.set_ylabel(units.get(name, ''))
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11)
    return _save(fig, path)


@styled
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
        ax.set_xlabel(TIME_LABEL)
        ax.set_ylabel(f'{output} [{units.get(output, "")}]')
        ax.legend(title='value', loc='best', title_fontsize=8)
        if failed:
            ax.text(0.02, 0.02, f'{failed} run(s) failed', transform=ax.transAxes,
                    color=TEXT_PRIMARY, fontsize=8)
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11)
    return _save(fig, path)


@styled
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


@styled
def plot_timestep_traces(path, runs, output, units, title):
    '''runs: [(dt, t, y)], coarse -> fine; darker = finer step.'''
    fig, ax = plt.subplots(figsize=(4.6, 3.4), constrained_layout=True)
    for colour, (dt, t, y) in zip(ramp(len(runs)), runs):
        ax.plot(t, y, color=colour, label=f'dt={dt:g}')
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel(f'{output} [{units.get(output, "")}]')
    ax.set_title(title)
    ax.legend(loc='best')
    return _save(fig, path)


@styled
def plot_model_vs_data(path, t_model, model, t_data, data, units, title, split=None):
    '''
    One panel per compared variable: model line, data as dots. ``split`` marks the end of
    the calibration window; data after it were held out. ``t_model`` may be a dict (a time base
    per compared variable).
    '''
    names = list(model)
    fig, axes = _grid(len(names))
    for ax, name in zip(axes, names):
        tm = t_model[name] if isinstance(t_model, dict) else t_model    # a dict: one time base per panel
        if split is not None:
            ax.axvspan(split, max(np.max(tm), np.max(t_data[name])), color=GRID, alpha=0.6, lw=0)
            ax.text(split, 0.98, ' held out', transform=ax.get_xaxis_transform(), va='top',
                    fontsize=8, color=TEXT_SECONDARY)
        ax.plot(tm, model[name], color=CATEGORICAL[0], label='model')
        ax.plot(t_data[name], data[name], linestyle='none', marker='o', markersize=4,
                color=DATA_COLOUR, label='data')
        ax.set_title(name)
        ax.set_xlabel(TIME_LABEL)
        ax.set_ylabel(units.get(name, ''))
        ax.legend(loc='best')
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11)
    return _save(fig, path)


def _risk_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list('risk', ['#fdf1ea', '#f4a47e', '#eb6834', '#b8461b', '#6e2308'])


def _axis_scale(ax, info, which='x'):
    if info['log']:
        (ax.set_xscale if which == 'x' else ax.set_yscale)('log')


@styled
def plot_risk_marginals(path, per_param, ranked, overall, title):
    '''Failure risk vs each parameter (binned), 95% band, most influential first.'''
    fig, axes = _grid(len(ranked), width=3.2, height=2.3, max_cols=4)
    for ax, var in zip(axes, ranked):
        info = per_param[var]
        bins = info['bins']
        mids = [math.sqrt(b['lo'] * b['hi']) if info['log'] else 0.5 * (b['lo'] + b['hi']) for b in bins]
        risk = [b['risk'] for b in bins]
        lo = [b['ci'][0] for b in bins]
        hi = [b['ci'][1] for b in bins]
        ax.fill_between(mids, lo, hi, color=CATEGORICAL[1], alpha=0.18, lw=0)
        ax.plot(mids, risk, color=CATEGORICAL[1], marker='o', markersize=4)
        ax.axhline(overall, color=MUTED, lw=1, ls='--')
        _axis_scale(ax, info)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"{var}  (η² = {info['importance']:.2f})")
        ax.set_ylabel('P(fail)')
    fig.suptitle(f'{title}\ndashed: overall P(fail) = {overall:.3f}; band: 95% interval', color=TEXT_PRIMARY, fontsize=10)
    return _save(fig, path)


@styled
def plot_risk_corner(path, x, failed, names, infos, corner_vars, title, n_bins=6):
    '''Lower triangle: failure rate for each parameter pair; diagonal: 1-D risk.'''
    idx = {n: i for i, n in enumerate(names)}
    k = len(corner_vars)
    fig, axes = plt.subplots(k, k, figsize=(2.3 * k + 1.2, 2.1 * k), squeeze=False, constrained_layout=True)
    cmap = _risk_cmap()
    mappable = None

    def edges(info):
        lo, hi = info['range']
        return np.geomspace(lo, hi, n_bins + 1) if info['log'] else np.linspace(lo, hi, n_bins + 1)
    for r, vy in enumerate(corner_vars):
        for c, vx in enumerate(corner_vars):
            ax = axes[r][c]
            if c > r:
                ax.set_visible(False)
                continue
            ix = infos[vx]
            if r == c:
                bins = ix['bins']
                mids = [math.sqrt(b['lo'] * b['hi']) if ix['log'] else 0.5 * (b['lo'] + b['hi']) for b in bins]
                ax.plot(mids, [b['risk'] for b in bins], color=CATEGORICAL[1], marker='o', markersize=3)
                ax.set_ylim(-0.02, 1.02)
                _axis_scale(ax, ix)
                ax.set_ylabel('P(fail)' if c == 0 else '')
            else:
                iy = infos[vy]
                ex, ey = edges(ix), edges(iy)
                tot, _, _ = np.histogram2d(x[:, idx[vx]], x[:, idx[vy]], bins=[ex, ey])
                bad, _, _ = np.histogram2d(x[failed, idx[vx]], x[failed, idx[vy]], bins=[ex, ey])
                with np.errstate(invalid='ignore', divide='ignore'):
                    rate = np.where(tot > 0, bad / tot, np.nan)
                mappable = ax.pcolormesh(ex, ey, rate.T, cmap=cmap, vmin=0, vmax=1, shading='flat')
                _axis_scale(ax, ix, 'x')
                _axis_scale(ax, iy, 'y')
                ax.grid(False)
                if c == 0:
                    ax.set_ylabel(vy)
            if r == k - 1:
                ax.set_xlabel(vx)
    if mappable is not None:
        fig.colorbar(mappable, ax=[a for row in axes for a in row if a.get_visible()], shrink=0.6, label='P(fail)')
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11)
    return _save(fig, path)
