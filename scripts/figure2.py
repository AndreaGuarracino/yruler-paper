#!/usr/bin/env python3
"""Fig. 2 in Alejandro's v3 layout (manuscript/figures/main/fig2/Fig2_Yruler_v3.svg), every panel drawn from the tables.
usage: python3 scripts/figure2.py [out.pdf]   (default manuscript/figures/main/fig2_biology.pdf)
Layout units are the SVG units (540 wide, y down); 1 unit = 0.75 pt in the PDF, as Fig. 1.
Nothing is cut from the SVG: its six panels all carry data. Positions, sizes and colours follow it; numbers come from the tables. Axis titles italic, as every figure (palette.save)."""
import re, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon, Rectangle
import figures_1kgp as F
from palette import PAD, crop_ink, CLASS, RATIO, INK, MUTED, LINE, GREY, ORANGE, GRID, SUPERPOP, LINEAGE, REF

ROOT = Path(__file__).resolve().parents[1]
M = ROOT / 'results/model'
W, Y0, Y1, S = 540, 4, 582, .75   # page width, SVG y band, pt per unit
AT, AB, AK = 45.2, 220.9, 1.15   # panel a rows stretched 15%
BT, BB, BK = 434.9, 544.4, 1.2   # d-f row stretched 20%
ga = lambda y: y if y <= AT else (AT + (y - AT) * AK if y <= AB else y + (AB - AT) * (AK - 1))
g = lambda y: ga(y) if y <= BT else (ga(BT) + (y - BT) * BK if y <= BB else ga(BB) + (BB - BT) * (BK - 1) + (y - BB))   # SVG y to page y
FX = 523.3 - (544.2 - 434.9) * BK   # f left edge, square axes
Y1 = g(Y1)
BAND_A, BAND, DIA, FRAME = '#EDEEEE', '#F2F4F5', '#7D8A91', '#231F20'   # v3 row bands, diamonds, frames
HI = {'DYZ3': '#CFE0EE', 'DYZ19': '#F1D2D9', 'c': '#ECF3F8'}   # v3 highlighted rows
UNIT = 2404.4 / 3478.9   # bp per Hallast 2026 unit, as figures_1kgp
fig = plt.figure(figsize=(W * S / 72, (Y1 - Y0) * S / 72))
fx, fy = (lambda x: x / W), (lambda y: 1 - (g(y) - Y0) / (Y1 - Y0))
rd = F.read_tsv


def axes(x0, y0, x1, y1):
    """Axes at SVG coordinates."""
    return fig.add_axes([fx(x0), fy(y1), (x1 - x0) / W, fy(y0) - fy(y1)])


def put(x, y, s, u=8, **k):
    """Text at SVG coordinates, y the baseline, u units."""
    k.setdefault('color', INK); k.setdefault('va', 'baseline')
    return fig.text(fx(x), fy(y), s, fontsize=max(u * S, 5), **k)


def line(xs, ys, **k):
    fig.add_artist(plt.Line2D([fx(x) for x in xs], [fy(y) for y in ys], transform=fig.transFigure, **k))


def rect(x0, y0, x1, y1, c, z=0):
    fig.add_artist(Rectangle((fx(x0), fy(y1)), (x1 - x0) / W, fy(y0) - fy(y1), transform=fig.transFigure, color=c, lw=0, zorder=z))


def style(ax, sides=('bottom',), u=8):
    for k, v in ax.spines.items(): v.set_visible(k in sides); v.set_linewidth(.5 * S); v.set_color(FRAME)
    ax.tick_params(which='both', width=.5 * S, length=2, labelsize=max(u * S, 5), pad=1.5, colors=INK); ax.tick_params(which='minor', length=1.2)


def power(x, y, base, exp, u=8):
    """'base × 10^exp' with a 5 pt exponent."""
    t = put(x, y, f'{base} × 10', u); fig.canvas.draw(); b = t.get_window_extent().transformed(fig.transFigure.inverted())
    put(b.x1 * W + .3, y - .4 * u, f'−{exp}', 5 / S)


def tree(ax, root, leaves, h, tip, end):
    """Collapsed SNP tree; tips at x tip, leaders to end."""
    m = max(n['x'] + n['extent'] for n in leaves); ax.set_xlim(0, m)
    def draw(n):
        if n['children']:
            ys = [c['y'] for c in n['children']]; ax.plot([n['x']] * 2, [min(ys), max(ys)], color=LINE, lw=.7 * S)
            for c in n['children']: ax.plot([n['x'], c['x']], [c['y']] * 2, color=LINE, lw=.7 * S); draw(c)
            return
        c = LINEAGE.get('NO' if n['group'].startswith('NO') else n['group'][0], LINE); e = n['x'] + n['extent']
        if n['extent'] > 1e-8: ax.add_patch(Polygon([[n['x'], n['y']], [e, n['y'] - h], [e, n['y'] + h]], fc=to_rgba(c, .75), ec=c, lw=.3 * S))
        else: ax.plot(e, n['y'], 'o', color=c, ms=2.4 * S, mew=0)
        ax.plot([e, m * end / tip], [n['y']] * 2, color='#939598', lw=.5 * S, ls=(0, (1, 1.6)), clip_on=False)
    draw(root); ax.axis('off'); return m


def scale_bar(m, x0, x1, bx, by, val, exp, base):
    """Scale bar of val substitutions per site."""
    L = val * (x1 - x0) / m; line([bx, bx + L], [by, by], color=LINE, lw=.75 * S)
    power(bx + L + 4, by + 2.8, base, exp); put(bx - 6, by + 12.5, 'substitutions per site')
    print(f'  tree scale: {(x1 - x0) / m:.4g} units per substitution, bar {val:g} = {L:.1f} units', flush=True)


def panel_a(valid):
    rows = {r['sample']: r for r in valid}
    root, leaves = F.collapse_tree(F.prune_tree(F.parse_newick(F.TREE_SOURCE), set(rows)), rows); n = len(leaves)
    if n != 33: sys.exit(f'ERROR a: {n} tree rows, the v3 layout has 33')
    y0, p, T, B = 48.3, (219.5 - 50.6) / 32, 45.2, 220.9; lim = ((B - y0) / p, (T - y0) / p)   # row 0 centre, pitch
    tx = axes(18.3, T, 66, B); tx.set_ylim(*lim); m = tree(tx, root, leaves, .29, 66 - 18.3, 68.6 - 18.3)
    scale_bar(m, 18.3, 66, 55.8, 229.2, 5e-5, 5, 5)
    for i, nd in enumerate(leaves):
        if i % 2 == 0: rect(68.6, y0 + p * (i - .5), 127.2, y0 + p * (i + .5), BAND_A)
        g = nd['group']; put(69.1, y0 + p * i + 2.4, f'{g} ({len(nd["members"])})', 6.67, fontweight='bold' if g.startswith('J2a1a1') or g in ('R1b', 'O2a2b1a1') else 'normal')
    print(f'a: {n} rows, {sum(len(nd["members"]) for nd in leaves)} males: ' + ', '.join(f'{nd["group"]} {len(nd["members"])}' for nd in leaves), flush=True)
    sp = axes(139.1, T, 169.8, B); sp.set_ylim(*lim); sp.set_xlim(0, 100)
    for i, nd in enumerate(leaves):
        left = 0
        for k, col in SUPERPOP.items():
            w = 100 * sum(r['superpopulation'] == k for r in nd['members']) / len(nd['members']); sp.barh(i, w, left=left, height=.6, color=col, lw=0); left += w
    style(sp); sp.set_yticks([]); sp.set_xticks([0, 100], ['0', '100']); put(154.5, 240.3, 'Males (%)', 10, ha='center', style='italic')
    for k, x, y in (('AFR', 148.5, 13.7), ('EAS', 134.7, 23.9), ('AMR', 165.9, 23.3), ('SAS', 134.7, 32.9), ('EUR', 165.7, 33.0)):
        rect(x - 10, y - 4.8, x - 2.7, y - .2, SUPERPOP[k], 2); put(x, y, k)
    put(69.1, 41.9, 'Y lineage (n)', fontweight='bold'); put(135.8, 42.4, 'Superpop.', fontweight='bold')
    for r in valid: r['unit_ratio'] = r['ratio_bp'] * UNIT
    cols = [('DYZ3', 'DYZ3_bp', 1e6, 181.2, 'log', (.0735, 3.21), [.1, .2, .5, 1, 2], 'Length (Mb, log)', lambda g: g in ('R1b', 'O2a2b1a1')),
            ('DYZ19', 'DYZ19_bp', 1e3, 268.8, 'linear', (-1.5, 557), [0, 200, 400], 'Length (kb)', lambda g: g.startswith('J2a1a1')),
            ('Yq12', 'Yq12_bp', 1e6, 356.4, 'linear', (-.12, 61.7), [0, 20, 40, 60], 'Length (Mb)', lambda g: False),
            ('DYZ1:DYZ2', 'unit_ratio', 1, 444.1, 'log2', (.4457, 2.396), [.5, 1, 2], 'Unit ratio (log)', lambda g: False)]
    rng = np.random.default_rng(2710)
    for name, key, sc, x0, scale, xl, ticks, lab, hi in cols:
        ax = axes(x0, T, x0 + 75.8, B); ax.set_ylim(*lim); c = CLASS.get(name, RATIO); nref = nmed = 0
        allv = np.array([r[key] / sc for r in valid]); out = int(((allv < xl[0]) | (allv > xl[1])).sum())
        if out: v3, xl = xl, (min(xl[0], allv.min() * .95), max(xl[1], allv.max() * 1.05)); print(f'  {name}: {out} males outside the v3 axis {v3[0]:.3g}-{v3[1]:.3g} (data {allv.min():.3g}-{allv.max():.3g}), widened to {xl[0]:.3g}-{xl[1]:.3g}', flush=True)
        for i, nd in enumerate(leaves):
            if hi(nd['group']): ax.axhspan(i - .5, i + .5, color=HI[name], lw=0, zorder=0)
            elif i % 2 == 0: ax.axhspan(i - .5, i + .5, color=BAND_A, lw=0, zorder=0)
            v = np.array([r[key] / sc for r in nd['members']]); y = i + rng.uniform(-.3, .3, len(v))
            t = np.array([r['sample'] in F.TRUTH[F.TRUTH_OF[key]] for r in nd['members']]); inc = np.array([F.in_summary(r, key) for r in nd['members']])
            ax.scatter(v[~t], y[~t], s=1.6, color=c, alpha=.55, lw=0, zorder=3); ax.scatter(v[t], y[t], zorder=4, **{**REF, 'edgecolors': DIA, 'linewidths': .45 * S, 'alpha': .9}); nref += int(t.sum())
            if inc.sum() >= 5:
                q1, md, q3 = np.percentile(v[inc], [25, 50, 75]); nmed += 1
                ax.plot([q1, q3], [i, i], color=INK, alpha=.55, lw=2 * S, solid_capstyle='butt', zorder=5); ax.plot([md, md], [i - .33, i + .33], color=INK, lw=1.1 * S, zorder=6)
        if scale == 'log': ax.set_xscale('log'); ax.xaxis.set_minor_locator(mt.LogLocator(subs=range(2, 10))); ax.xaxis.set_minor_formatter(mt.NullFormatter())
        if scale == 'log2': ax.set_xscale('log', base=2); ax.xaxis.set_minor_locator(mt.NullLocator()); ax.axvline(1, color=RATIO, lw=.6 * S, ls=(0, (3.7, 1.6)), zorder=2)
        ax.set_xlim(*xl); ax.set_xticks(ticks, [f'{t:g}' for t in ticks]); ax.set_yticks([]); style(ax); ax.patch.set_alpha(0)
        for t in ticks: ax.axvline(t, color='white', lw=.6 * S, zorder=1)
        put(x0 + 37.9, 42.2, name, 10, ha='center', fontweight='bold', color=c); put(x0 + 37.9, 240.3, lab, 10, ha='center', style='italic')
        print(f'  a {name}: {len(allv)} males, {nref} reference diamonds, medians in {nmed} of {n} rows (n >= 5)', flush=True)
    g = {nd['group']: i for i, nd in enumerate(leaves)}; a1, a2 = g['O2a2b1a1'], g['O2a2b1a2']
    med = lambda i: np.median([r['DYZ3_bp'] for r in leaves[i]['members'] if F.in_summary(r, 'DYZ3_bp')])
    fold = med(a2) / med(a1); ya, yb = y0 + p * a1, y0 + p * a2; line([243.5, 245, 245, 243.5], [ya, ya, yb, yb], color=INK, lw=.5 * S)
    put(247.3, (ya + yb) / 2 + 2.8, f'{fold:.1f}×', fontweight='bold'); print(f'  a fold O2a2b1a2 / O2a2b1a1 DYZ3 (1kGP): {med(a2) / 1e3:.0f} / {med(a1) / 1e3:.0f} kb = {fold:.2f}', flush=True)
    put(17.3, 25, 'a', 20, fontweight='bold'); put(350.5, 18.9, 'Satellite Length Differences Across Lineages', 12, ha='center', fontweight='bold')
    keys = [(6, lambda x: line([x + 2], [27.2], marker='D', ms=3.6 * S, mfc='none', mec=DIA, mew=.45 * S), 'Reference assembly'),
            (17, lambda x: line([x, x + 14], [27.2] * 2, color=INK, alpha=.55, lw=2 * S, solid_capstyle='butt'), 'Interquartile range'),
            (5, lambda x: line([x + 1] * 2, [25, 29.4], color=INK, lw=1.1 * S), 'Median (n ≥ 5)')]
    T = [put(0, 29.8, s) for _, _, s in keys]; fig.canvas.draw(); w = [t.get_window_extent().width / fig.bbox.width * W for t in T]
    x = 350.5 - (sum(w) + sum(k for k, _, _ in keys) + 2 * 10) / 2   # centred under the title, 10-unit gaps
    for (k, draw, _), t, wi in zip(keys, T, w): draw(x); t.set_x(fx(x + k)); x += k + wi + 10


def bracket(ax, xa, xb, ta, tb, label, hu, w='bold'):
    """Fold bracket above two groups; hu axes height."""
    f = lambda t: ax.transAxes.inverted().transform(ax.transData.transform((0, t)))[1]
    ya, yb = f(ta) + 3.5 / hu, f(tb) + 3.5 / hu; h = max(max(ya, yb) + 1.5 / hu, min(max(ya, yb) + 5.5 / hu, 1 - 10.4 / hu)); tr = ax.get_xaxis_transform()   # label 2 units under the frame
    ax.plot([xa, xa, xb, xb], [ya, h, h, yb], transform=tr, color=MUTED, lw=.5 * S, clip_on=False)
    ax.text((xa + xb) / 2, h + 2.6 / hu, label, transform=tr, ha='center', va='baseline', fontsize=8 * S, color=INK, fontweight=w)


def pairs(ax, G, xs, cols, hu, ycount, ycoh, s, jit, half):
    """Two strips per cohort, medians, fold brackets."""
    rng = np.random.default_rng(5)
    for i, (coh, A, B) in enumerate(G):
        for (v, t), x, (c, a) in zip((A, B), xs[2 * i: 2 * i + 2], cols):
            v, t = np.array(v), np.array(t, bool); j = x + rng.uniform(-jit, jit, len(v))
            ax.scatter(j[~t], v[~t], s=s, color=c, alpha=a, lw=0, zorder=3); ax.scatter(j[t], v[t], zorder=4, **{**REF, 's': s * 1.6, 'edgecolors': DIA, 'linewidths': .5 * S, 'alpha': .9})
            ax.plot([x - half, x + half], [np.median(v)] * 2, color=INK, lw=1.3 * S, zorder=5); put(x, ycount, f'{len(v):,}', ha='center')
        f = np.median(B[0]) / np.median(A[0]); bracket(ax, xs[2 * i], xs[2 * i + 1], max(A[0]), max(B[0]), f'{f:.1f}×', hu)
        put((xs[2 * i] + xs[2 * i + 1]) / 2, ycoh, coh, ha='center', fontweight='bold')
        print(f'  {coh}: {len(A[0])} vs {len(B[0])}, median {np.median(A[0]):.4g} vs {np.median(B[0]):.4g}, fold {f:.3f}, diamonds {int(sum(A[1]) + sum(B[1]))}', flush=True)


def panel_b(valid, Q):
    one = [r for r in valid if F.in_summary(r, 'DYZ3_bp')]; tr = F.TRUTH['dyz3']
    G = [('1kGP', *[([r['DYZ3_bp'] / 1e6 for r in one if (r['display_group'] == 'R1b') == y], [r['sample'] in tr for r in one if (r['display_group'] == 'R1b') == y]) for y in (True, False)])]
    G += [(c, *[([float(r['DYZ3_bp']) / 1e6 for r in Q if r['cohort'] == c and r['call'].startswith('R1b') == y], [False] * sum(r['cohort'] == c and r['call'].startswith('R1b') == y for r in Q)) for y in (True, False)]) for c in ('HGDP', 'SGDP')]
    ax = axes(47, 262.5, 246.5, 384.1); ax.set_yscale('log'); ax.set_ylim(.0887, 5.16); ax.set_xlim(47, 246.5)
    print('b: DYZ3, R1b vs all other males', flush=True)
    pairs(ax, G, [61.8, 89.9, 134.1, 160.1, 202.5, 229.1], ((CLASS['DYZ3'], .75), (GREY, .6)), 384.1 - 262.5, 392.5, 401.7, 5.5, 8.5, 9)
    ax.set_xticks([]); ax.set_yticks([.1, .2, .5, 1, 2], ['0.1', '0.2', '0.5', '1', '2']); ax.yaxis.set_minor_formatter(mt.NullFormatter()); style(ax, ('left', 'bottom', 'top', 'right'))
    put(27.2, 372.1, 'DYZ3', 10, rotation=90, rotation_mode='anchor', fontweight='bold', color=CLASS['DYZ3'], style='italic'); put(27.2, 372.1 - 25.6, ' length (Mb, log)', 10, rotation=90, rotation_mode='anchor', style='italic')
    put(21.6, 256.8, 'b', 20, fontweight='bold'); put(146.75, 259.4, 'R1b vs. all other males', 12, ha='center', fontweight='bold')


def panel_c(valid, Q):
    rows = F.subtree_rows(valid, 'O2a'); root, leaves = F.collapse_tree(F.prune_tree(F.parse_newick(F.TREE_SOURCE), set(rows)), rows)
    D = {nd['group']: {'1kGP': [(r['DYZ3_bp'] / 1e6, r['sample'] in F.TRUTH['dyz3']) for r in nd['members'] if F.in_summary(r, 'DYZ3_bp')],
                       **{c: [(float(r['DYZ3_bp']) / 1e6, False) for r in Q if r['cohort'] == c and r['call'].startswith(nd['group'])] for c in ('HGDP', 'SGDP')}} for nd in leaves}
    n, T, B, p = len(leaves), 267.6, 389.2, (389.2 - 267.6) / 6; cy = lambda i: T + p * (i + .5)
    if n != 6: sys.exit(f'ERROR c: {n} clades, the v3 layout has 6')
    tx = axes(264, T, 312.6, B); tx.set_ylim(n - .5, -.5); m = tree(tx, root, leaves, .3, 312.6 - 264, 314 - 264)
    scale_bar(m, 264, 312.6, 284.8, 397.5, 1e-5, 5, 1)
    ax = axes(352.7, T, 508, B); ax.set_ylim(n - .5, -.5); ax.set_xlim(-.18, 2.081); ax.patch.set_alpha(0); c = CLASS['DYZ3']; rng = np.random.default_rng(2710)
    for i, nd in enumerate(leaves):
        g = nd['group']; put(315.1, cy(i) - 1.6, g, 7, fontweight='bold' if g == 'O2a2b1a1' else 'normal'); put(315.1, cy(i) + 5.8, '/'.join(str(len(D[g][k])) for k in ('1kGP', 'HGDP', 'SGDP')), 7, color=MUTED)
        if g == 'O2a2b1a1': rect(314, cy(i) - p / 2, 508, cy(i) + p / 2, HI['c'], -1)
        if i: ax.axhline(i - .5, color='#6D6E71', lw=.5 * S, zorder=2)
        for k, coh in enumerate(('1kGP', 'HGDP', 'SGDP')):
            y = i + (k - 1) / 3
            if k != 1 and g != 'O2a2b1a1': ax.axhspan(y - 1 / 6, y + 1 / 6, color=BAND, lw=0, zorder=0)
            if not D[g][coh]: continue
            v = np.array([a for a, _ in D[g][coh]]); t = np.array([b for _, b in D[g][coh]], bool); q = np.percentile(v, [10, 25, 50, 75, 90])
            ax.bxp([dict(whislo=q[0], q1=q[1], med=q[2], q3=q[3], whishi=q[4], fliers=[])], positions=[y], widths=.2, vert=False, patch_artist=True, manage_ticks=False, showcaps=False,
                   boxprops=dict(facecolor=to_rgba(c, .25), edgecolor=c, lw=.5 * S), whiskerprops=dict(color=c, lw=.5 * S), medianprops=dict(color=INK, lw=1 * S), zorder=3)
            j = y + rng.uniform(-.07, .07, len(v)); ax.scatter(v[~t], j[~t], s=1.2, color=c, alpha=.5, lw=0, zorder=4); ax.scatter(v[t], j[t], zorder=5, **{**REF, 'edgecolors': c, 'linewidths': .8 * S})
    print('c: DYZ3 per clade (1kGP/HGDP/SGDP), 1kGP diamonds: ' + ', '.join(f'{g} {"/".join(str(len(d[k])) for k in d)} ({sum(b for _, b in d["1kGP"])})' for g, d in D.items()), flush=True)
    ax.set_xticks([0, .5, 1, 1.5, 2], ['0.0', '0.5', '1.0', '1.5', '2.0']); ax.set_xticks(np.arange(0, 2.01, .25), minor=True); ax.set_yticks([]); style(ax)
    for t in (0, .5, 1, 1.5, 2): ax.axvline(t, color=GRID, lw=.4 * S, zorder=1)
    i1, i2 = [i for i, nd in enumerate(leaves) if nd['group'] in ('O2a2b1a1', 'O2a2b1a2')]; line([510.8, 512.8, 512.8, 510.8], [cy(i1), cy(i1), cy(i2), cy(i2)], color=INK, lw=.5 * S)
    for k, coh in enumerate(('1kGP', 'HGDP', 'SGDP')):
        a, b = (np.median([v for v, _ in D[leaves[i]['group']][coh]]) for i in (i1, i2)); put(514.4, (cy(i1) + cy(i2)) / 2 - 5 + 6.9 * k, f'{b / a:.1f}×', fontweight='bold')
        put(351.4, cy(0) - p / 3 + p / 3 * k + 3.0, coh, 7); print(f'  c fold O2a2b1a2 / O2a2b1a1 {coh}: {b * 1e3:.0f} / {a * 1e3:.0f} kb = {b / a:.3f}', flush=True)
    put(404.1, 409.6, 'DYZ3', 10, fontweight='bold', color=c, style='italic'); put(404.1 + 25.6, 409.6, ' length (Mb)', 10, style='italic')
    put(260.7, 256.8, 'c', 20, fontweight='bold'); put(430.35, 261.6, 'O2a2b1a1 vs. O2a2b1a2', 12, ha='center', fontweight='bold')
    put(293.2, 259.3, 'Clade', fontweight='bold'); put(264.9, 269.3, 'n: 1kGP/HGDP/SGDP', fontweight='bold')


def panel_d(Q):
    C = rd(M / 'cohort_1kgp.tsv'); keep = ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain', 'par2_gain')   # as figures_1kgp.biology_inputs
    J = [r for r in C if r['sex'] == 'male' and r['estimate_basis'] == 'ruler' and r['y_state'] in keep and r['lineage'] == 'J']; Jq = [r for r in Q if r['call'].startswith('J')]
    kb = lambda R: ([float(r['DYZ19_bp']) / 1e3 for r in R], [r['sample'] in F.TRUTH['dyz19'] for r in R])
    G = [('1kGP', kb([r for r in J if r['haplogroup'].startswith('J2a1a1')]), kb([r for r in J if not r['haplogroup'].startswith('J2a1a1')]))]
    G += [(c, kb([r for r in Jq if r['cohort'] == c and r['call'].startswith('J2a1a1')]), kb([r for r in Jq if r['cohort'] == c and not r['call'].startswith('J2a1a1')])) for c in ('HGDP', 'SGDP')]
    DR = 152; k = (DR - 38.7) / (167.3 - 38.7); dx = lambda x: 38.7 + (x - 38.7) * k   # d narrowed for e
    ax = axes(38.7, 435.3, DR, 544.3); ax.set_ylim(0, 405); ax.set_xlim(38.7, DR); print('d: DYZ19, J2a1a1 vs other J', flush=True)
    pairs(ax, G, [dx(x) for x in (48.3, 67.1, 93.7, 112.4, 137.9, 156.6)], ((CLASS['DYZ19'], .8), (GREY, .65)), 544.3 - 435.3, 552.9, 561.7, 4.5, 5 * k, 6.5 * k)
    ax.set_xticks([]); ax.set_yticks([0, 100, 200, 300]); style(ax, ('left', 'bottom', 'top', 'right'))
    from matplotlib.textpath import TextPath; from matplotlib.font_manager import FontProperties
    wt = max(TextPath((0, 0), t, size=8 * S, prop=FontProperties(family=plt.rcParams['font.sans-serif'])).get_extents().width for t in ('J2a1a1', 'Other J')) / S
    xm = (dx(93.7) + dx(137.9)) / 2 - (4.3 + wt - 1.6) / 2; print(f'  d legend centred between legs {dx(93.7):.1f} and {dx(137.9):.1f}: marker x {xm:.1f}, text width {wt:.1f}', flush=True)   # equidistant
    for y, c, s in ((499.4, CLASS['DYZ19'], 'J2a1a1'), (509.2, GREY, 'Other J')): line([xm], [y], marker='o', ms=3.2 * S, color=c, mew=0); put(xm + 4.3, y + 2.3, s, 8)
    put(16.6, 526.6, 'DYZ19', 10, rotation=90, rotation_mode='anchor', fontweight='bold', color=CLASS['DYZ19'], style='italic'); put(16.6, 526.6 - 27.3, ' length (kb)', 10, rotation=90, rotation_mode='anchor', style='italic')
    put(22.8, 429.7, 'd', 20, fontweight='bold'); put((38.7 + DR) / 2, 430.8, 'J2a1a1 vs. other J', 12, ha='center', fontweight='bold')


def panel_e():
    U = rd(M / 'dyz1_dyz2_lineages.tsv'); A = {r['measure'].split('_')[1]: float(r['median']) for r in rd(M / 'dyz1_dyz2_truth_vs_yruler.tsv') if r['measure'].startswith('lineage_') and '_truth_' in r['measure']}
    C = rd(M / 'cohort_1kgp.tsv'); r2 = {r['haplogroup'][:2] for r in C if r['sex'] == 'male' and r['lineage'] == 'R'}   # table lineage R: all R2?
    name = lambda l: 'R2' if l == 'R' and r2 == {'R2'} else l; print(f'e: {len(U)} lineages; 1kGP lineage R haplogroups start with {sorted(r2)}', flush=True)
    U.sort(key=lambda r: float(r['unit_ratio_median'])); lo, hi = min(float(r['ci95_low']) for r in U), max(float(r['ci95_high']) for r in U)
    ax = axes(195, 435.2, 346, 544.4); xl = (min(.7768, lo - .01), max(1.1685, hi + .01)); ax.set_xlim(*xl); ax.set_ylim(len(U) - 1 + .7, -.7); c = CLASS['DYZ3']
    for i, r in enumerate(U):
        if i % 2 == 0: ax.axhspan(i - .5, i + .5, color=BAND, lw=0, zorder=0)
        ax.plot([float(r['ci95_low']), float(r['ci95_high'])], [i, i], color=c, lw=1.1 * S, zorder=3); ax.plot(float(r['unit_ratio_median']), i, 'o', ms=2.6, color=c, mew=0, zorder=4)
        if r['lineage'] in A: ax.plot(A[r['lineage']], i, 'x', ms=3.4, color=ORANGE, mew=.8 * S, zorder=5)
        print(f'  e {name(r["lineage"])} ({r["men"]}): median {r["unit_ratio_median"]} ({r["ci95_low"]}-{r["ci95_high"]})' + (f', assemblies {A[r["lineage"]]}' if r['lineage'] in A else ''), flush=True)
    ax.axvline(1, color=RATIO, lw=.6 * S, ls=(0, (3.7, 1.6)), zorder=2)
    for t in (.8, .9, 1, 1.1): ax.axvline(t, color=GRID, lw=.4 * S, zorder=1)
    ax.set_yticks(range(len(U)), [f'{name(r["lineage"])} ({r["men"]})' for r in U]); ax.set_xticks([.8, .9, 1, 1.1], ['0.8', '0.9', '1.0', '1.1']); style(ax, ('left', 'bottom', 'top', 'right')); ax.tick_params(axis='y', length=0)
    line([201.3], [531.6], marker='x', ms=3.4, color=ORANGE, mew=.8 * S); put(204.9, 534.5, 'Assembly median'); line([201.3], [539.6], marker='o', ms=2.6, color=c, mew=0); put(205, 542.5, 'Read median')
    put(270.5, 568.5, 'DYZ1:DYZ2 unit ratio (median, 95% interval)', 10, ha='center', style='italic')
    put(166, 429.7, 'e', 20, fontweight='bold'); put(270.5, 429.6, 'Unit ratio by lineage', 12, ha='center', fontweight='bold')


def panel_f():
    """DYZ1 against DYZ2, lineage-centred log2 (old f data)."""
    from matplotlib.legend_handler import HandlerTuple
    _, C, _, N, inc = F.biology_inputs(); lin = {r['sample']: r['lineage'] for r in C}
    kin = {r['sample']: r['kinship_component_id'] for r in rd(ROOT / 'results/estimates_1kgp.relatedness_samples.tsv')}; rep = {}
    for r in sorted((r for r in inc if r['y_state'] != 'par2_gain'), key=lambda r: r['sample']): rep.setdefault(kin.get(r['sample'], r['sample']), r)   # one male per kinship component
    def centred(P):
        a, b, l = np.log2([p[0] for p in P]), np.log2([p[1] for p in P]), np.array([p[2] for p in P])
        for k in set(l): a[l == k] -= a[l == k].mean(); b[l == k] -= b[l == k].mean()
        return b, a
    P = [r for r in rep.values() if r['reference_male'] != 'yes']; tri = np.array([r['sample'] in F.TRUTH['yq12'] for r in P])
    T = {r['assembly']: r for r in rd(ROOT / 'data/method/yq12_truth_hallast.tsv')}
    men = centred([(float(r['DYZ1_bp']), float(r['DYZ2_bp']), r['lineage']) for r in P]); asm = centred([(float(r['dyz1_units']), float(r['dyz2_units']), lin[s]) for s, r in T.items() if s in lin and s != 'HG01261'])
    sd = {k: re.search(r'^(\d+): .*ratio ([0-9.]+) \(95% ([0-9.]+) to ([0-9.]+)\)', N[f'composition_sd_{k}']).groups() for k in ('assemblies', '1kgp', 'hgdp', 'sgdp')}
    if (len(men[0]), len(asm[0])) != (int(sd['1kgp'][0]), int(sd['assemblies'][0])): sys.exit(f'ERROR f: {len(men[0])} males, {len(asm[0])} assemblies, table {sd["1kgp"][0]}, {sd["assemblies"][0]}')
    r, a = float(sd['1kgp'][1]), float(sd['assemblies'][1]); big = max(np.abs(np.concatenate([*men, *asm]))); L = max(3, .5 * np.ceil(big / .5))
    print(f'f: {len(men[0])} 1kGP males ({int(tri.sum())} Yq12 reference diamonds), {len(asm[0])} assemblies; SD ratio from points {np.std(men[1]) / np.std(men[0]):.3f} and {np.std(asm[1]) / np.std(asm[0]):.3f}, table {r} and {a}; largest |log2| {big:.2f}, axis +-{L:g}', flush=True)
    ax = axes(FX, 434.9, 523.3, 544.2); ax.set_xlim(-L, L); ax.set_ylim(-L, L)
    ax.scatter(men[0][~tri], men[1][~tri], s=2, c=GREY, lw=0, alpha=.6, zorder=3, rasterized=True); ax.scatter(men[0][tri], men[1][tri], zorder=3, rasterized=True, **{**REF, 'edgecolors': DIA, 'linewidths': .45 * S})
    ax.scatter(*asm, s=6, marker='x', c=ORANGE, linewidths=.8 * S, alpha=.9, zorder=4, rasterized=True)
    ax.plot([-L, L], [-L, L], color=MUTED, lw=.6 * S, ls=(0, (3.7, 1.6)), zorder=2); ax.plot([-L, L], [-L * r, L * r], color=INK, lw=1.1 * S, zorder=5); ax.plot([-L, L], [-L * a, L * a], color=ORANGE, lw=1.1 * S, zorder=5)   # slopes: SD ratios through 0
    lab = {-3: '1/8×', -2: '1/4×', -1: '1/2×', 0: '1×', 1: '2×', 2: '4×', 3: '8×'}; t = [k for k in lab if abs(k) <= L]; tl = [lab[k] for k in t]; ax.set_xticks(t, tl); ax.set_yticks(t, tl); style(ax, ('left', 'bottom', 'top', 'right'))   # b-e sizes
    for k in t: ax.axvline(k, color=GRID, lw=.4 * S, zorder=1); ax.axhline(k, color=GRID, lw=.4 * S, zorder=1)
    h = lambda c, m, ms, mew=0, ls='-': (plt.Line2D([], [], color=c, lw=1.1 * S, ls=ls), plt.Line2D([], [], ls='none', marker=m, ms=ms, color=c, mew=mew))
    lg = ax.legend([h(ORANGE, 'x', 2.5, .8 * S), h(INK, 'o', 2.2), plt.Line2D([], [], color=MUTED, lw=.6 * S, ls=(0, (3.7, 1.6)))], [f'Assemblies ({len(asm[0])}), {a:g}', f'1kGP reads ({len(men[0]):,}), {r:g}', 'Equal variation, 1'],
              handler_map={tuple: HandlerTuple(ndivide=1)}, loc='upper left', bbox_to_anchor=(-.015, 1.012), fontsize=8 * S, frameon=False, borderpad=.2, borderaxespad=.3, handlelength=1.3, handletextpad=.5, labelspacing=.25)   # one line each, no backing
    fig.canvas.draw(); (bx0, by0), (bx1, by1) = ax.transData.inverted().transform(lg.get_window_extent())   # legend box in data units
    print(f'  f legend box x {bx0:.2f} to {bx1:.2f}, y {by0:.2f} to {by1:.2f}; lines inside: {[k for k in (1, r, a) if k * bx1 >= by0]}', flush=True)
    B = [.719, .16, .25, .19]; ins = ax.inset_axes(B); x0, y0 = -L + 2 * L * B[0], -L + 2 * L * B[1]
    print(f'  f points under the inset: {int(sum(((x >= x0) & (x <= x0 + 2 * L * B[2]) & (y >= y0) & (y <= y0 + 2 * L * B[3])).sum() for x, y in (men, asm)))}', flush=True)
    for i, (n, k) in enumerate((('Assemblies', 'assemblies'), ('1kGP', '1kgp'), ('HGDP', 'hgdp'), ('SGDP', 'sgdp'))):
        v, lo, hi = map(float, sd[k][1:]); c = ORANGE if i == 0 else INK; ins.plot([lo, hi], [3 - i] * 2, color=c, lw=.75 * S * 1.3); ins.plot(v, 3 - i, 'x' if i == 0 else 'o', ms=2.5 if i == 0 else 2, mew=.8 * S, color=c)
        print(f'  f inset {n}: n {sd[k][0]}, SD ratio {v} (95% {lo}-{hi})', flush=True)
    ins.axvline(1, color=MUTED, lw=.6 * S, ls=(0, (3.7, 1.6))); ins.set_yticks(range(4), ['SGDP', 'HGDP', '1kGP', 'Assemblies']); ins.set_xlim(.55, 1.05); ins.set_ylim(-.5, 3.5)
    ins.set_xticks([.7, 1], ['0.7', '1.0']); ins.set_xticks([.8, .9], minor=True); style(ins, ('left', 'bottom')); ins.tick_params(length=1.5, pad=1); ins.patch.set_alpha(.85)
    ins.set_xlabel('SD ratio', fontsize=8 * S, labelpad=1.5, style='italic', color=INK, x=.5 + 2.5 / (.25 * (523.3 - FX)))   # 2.5 units right
    put(FX - 21, 489.55, 'DYZ1 / lineage geometric mean', 10, rotation=90, rotation_mode='anchor', ha='center', style='italic'); put((FX + 523.3) / 2, 568.5, 'DYZ2 / lineage geometric mean', 10, ha='center', style='italic')
    put(FX - 39.5, 429.7, 'f', 20, fontweight='bold'); put((FX + 523.3) / 2, 430.2, 'DYZ2 varies more than DYZ1', 12, ha='center', fontweight='bold')


def check():
    """Every text at least 5 pt, inside the page."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); T = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    small = [(t.get_text(), t.get_fontsize()) for t in T if t.get_fontsize() < 5 - 1e-9]
    out = [t.get_text() for t in T if (lambda b: b.x0 < 0 or b.y0 < 0 or b.x1 > fig.bbox.width or b.y1 > fig.bbox.height)(t.get_window_extent(r))]
    b = fig.get_tightbbox(r); w, h = fig.get_size_inches(); over = {k: round(float(v) * 25.4, 2) for k, v in (('left', -b.x0), ('bottom', -b.y0), ('right', b.x1 - w), ('top', b.y1 - h)) if v > 0}
    print(f'check: {len(T)} texts, smallest {min(t.get_fontsize() for t in T):.2f} pt, below 5 pt: {small}, texts outside the page: {out}, drawing past the page (mm): {over}', flush=True)
    if small or out or over: sys.exit('ERROR text below 5 pt or content outside the page')


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'manuscript/figures/main/fig2_biology.pdf'
    _, valid, _ = F.load_snapshot(); Q = rd(M / 'hgsg_calls.tsv'); print(f'inputs: {len(valid)} 1kGP males on the ruler, {len(Q)} HGDP/SGDP calls', flush=True)
    panel_a(valid); panel_b(valid, Q); panel_c(valid, Q); panel_d(Q); panel_e(); panel_f(); check()
    fig.savefig(out, transparent=False, bbox_inches='tight', pad_inches=PAD); crop_ink(out); w, h = fig.get_size_inches() * 25.4   # 1.5 mm ink margin, as palette.save
    print(f'written {out}: {w * 72 / 25.4:.1f} x {h * 72 / 25.4:.1f} pt ({w:.1f} x {h:.1f} mm)', flush=True)


if __name__ == '__main__':
    main()
