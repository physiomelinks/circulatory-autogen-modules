"""Digitise Liang 2009 Fig. 4a (normal LV pressure) and Fig. 4e (normal LV P-V loop) from the 300 dpi page render.

The "Normal" curves are the thick black lines (the thin black femoral-stenosis curve lies on them).
Axes are calibrated from the tick marks found in the image."""
import json, sys
import numpy as np
from PIL import Image

D = '/home/farg967/.claude/jobs/3ade3515/tmp/fig'
page = np.asarray(Image.open(f'{D}/p751-09.png').convert('RGB')).astype(int)


def black(a):
    return (a[..., 0] < 90) & (a[..., 1] < 90) & (a[..., 2] < 90)


def axes(box):
    '''Left axis column and bottom axis row of a plot in page box (x0, y0, x1, y1): the longest black runs.'''
    x0, y0, x1, y1 = box
    b = black(page[y0:y1, x0:x1])
    col = int(np.argmax(b.sum(axis=0)))
    row = int(np.argmax(b.sum(axis=1)))
    return x0 + col, y0 + row, b


def ticks_x(axis_row, x0, x1, below=12):
    '''Columns of the tick marks hanging below the bottom axis.'''
    strip = black(page[axis_row + 3:axis_row + below, x0:x1]).sum(axis=0)
    cols = np.where(strip >= below - 6)[0] + x0
    groups = np.split(cols, np.where(np.diff(cols) > 2)[0] + 1)
    return [float(g.mean()) for g in groups if len(g)]


def ticks_y(axis_col, y0, y1, left=12):
    strip = black(page[y0:y1, axis_col - left:axis_col - 3]).sum(axis=1)
    rows = np.where(strip >= left - 6)[0] + y0
    groups = np.split(rows, np.where(np.diff(rows) > 2)[0] + 1)
    return [float(g.mean()) for g in groups if len(g)]


def curve(xa, ya, x0, x1, y0, y1):
    '''Mean row of the black pixels in each column of the plot area (the thick curve).'''
    out = []
    for x in range(x0, x1):
        rows = np.where(black(page[y0:y1, x:x + 1])[:, 0])[0] + y0
        if len(rows):
            out.append((x, float(np.median(rows))))
    return out


# ---- Fig. 4a: LV pressure, x 0-1 s, y 0-200 mmHg ------------------------------------------
box_a = (150 + 40, 240 + 20, 150 + 1000, 240 + 640)
xa, ya, _ = axes(box_a)
tx = ticks_x(ya, xa - 5, 150 + 1000)
ty = ticks_y(xa, 240 + 20, ya + 5)
print('4a axis col/row', xa, ya, 'x ticks', [round(t) for t in tx], 'y ticks', [round(t) for t in ty])
# major ticks: x at 0.0, 0.5, 1.0 (minor at 0.25, 0.75); y every 20 mmHg (labels 0..200 by 40)
x_of = lambda px: np.interp(px, [tx[0], tx[-1]], [0.0, 1.0])
y_of = lambda py: np.interp(py, [ty[0], ty[-1]], [200.0, 0.0])
row_of = lambda mmhg: ty[-1] - mmhg*(ty[-1]-ty[0])/200.0
pts = curve(xa, ya, int(tx[0]) + 3, int(tx[-1]) - 1, int(row_of(118)), ya - 2)   # below the legend
t = np.array([x_of(x) for x, _ in pts]); p = np.array([y_of(y) for _, y in pts])
# resample every 0.01 s
tg = np.round(np.arange(0.0, 1.0001, 0.01), 3)
pg = np.interp(tg, t, p)
print('4a: peak %.1f mmHg at %.3f s; min %.1f at %.3f s; p(0)=%.1f p(1)=%.1f' % (pg.max(), tg[pg.argmax()], pg.min(), tg[pg.argmin()], pg[0], pg[-1]))

# ---- Fig. 4e: LV P-V loop, x 30-140 ml, y 0-200 mmHg --------------------------------------
box_e = (100 + 150, 1660 + 120, 100 + 1080, 1660 + 740)
xe, ye, _ = axes(box_e)
txe = ticks_x(ye, xe - 5, 100 + 1080)
tye = ticks_y(xe, 1660 + 120, ye + 5)
print('4e axis col/row', xe, ye, 'x ticks', [round(t) for t in txe], 'y ticks', [round(t) for t in tye])
json.dump({'4a': {'t': tg.tolist(), 'p': pg.round(2).tolist(), 'xticks': tx, 'yticks': ty},
           '4e_axes': {'xticks': txe, 'yticks': tye, 'axis': [xe, ye]}}, open(f'{D}/digitised.json', 'w'))

# P-V loop: the normal loop's vertical edges (ESV left, EDV right), below the legend (< 160 mmHg)
v_of = lambda px: np.interp(px, [txe[0], txe[-1]], [40.0, 140.0])
pe_of = lambda py: np.interp(py, [tye[0], tye[-1]], [200.0, 0.0])
r0 = int(tye[-1] - 150*(tye[-1]-tye[0])/200)       # 150 mmHg
r1 = int(tye[-1]) - 3
b = black(page[r0:r1, xe + 3:int(txe[-1])])
counts = b.sum(axis=0)
cols = np.where(counts > 0.25*(r1 - r0))[0] + xe + 3
groups = np.split(cols, np.where(np.diff(cols) > 3)[0] + 1)
edges = [(float(g.mean()), int(counts[g - xe - 3].max())) for g in groups if len(g)]
print('4e vertical edges (ml, height px):', [(round(float(v_of(c)), 1), h) for c, h in edges])
