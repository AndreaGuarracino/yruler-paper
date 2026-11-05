#!/usr/bin/env python3
"""Fig. 1: panels a-b cut from Alejandro's Fig1_Yruler_v4.svg, panels c-f drawn from the model tables in its style.
usage: python3 scripts/figure1.py [out dir] | python3 scripts/figure1.py ed4   (ed4: Extended Data Fig. 4 in the Fig. 1c style)   (default manuscript/figures/main; TOP: c drops points above TOP %, counted as triangles; needs inkscape and pypdf)
writes <out>/fig1/fig1{a..f}_*.pdf (one per panel) and <out>/fig1_method.pdf/.svg (the panels composed).
Layout units are the SVG units (540 wide, y down); 1 unit = 0.75 pt in the PDF, as in the Inkscape export of v4."""
import collections as co, csv, os, re, subprocess, sys, tempfile, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.colors import to_rgba
from pypdf import PdfReader, PdfWriter, Transformation
from palette import CLASS, RATIO, DOSE

ROOT = Path(__file__).resolve().parents[1]
SVG = ROOT / 'manuscript/figures/main/fig1/Fig1_Yruler_v4.svg'   # Alejandro's v4, source of b
SVG_A = ROOT / 'manuscript/figures/main/fig1/Fig1A_Yruler_v3 (1).svg'   # Alejandro's compact a (2026-11-04)
A_BOX = (10.83, 103)   # a viewBox y range: 5.5 units above the title, 1 below the lowest label
M = ROOT / 'results/model'
W, H, S = 540, 870, 0.75   # page width and height (units), pt per unit
PANELS = {'a': ('chrY', 0, 0, 540, 114), 'b': ('workflow', 0, 96, 540, 342), 'c': ('cv', 0, 342, 540, 532),   # name, region x0 y0 x1 y1
          'd': ('coverage', 0, 532, 204, 720), 'e': ('pedigrees', 204, 532, 540, 720), 'f': ('specificity', 0, 720, 540, 870)}
LOGO = (234, 99, 322, 136)   # Yruler logo box, part of b
A_UP = 96 - 68.47   # b top on the page: logo letters 5 units under the DYZ18S bar of a (bar bottom y 87.3 in SVG_A, page 76.47)
NAMES = ['DYZ3', 'DYZ17', 'DYZ19', 'Yq12', 'DYZ18S', 'DYZ1', 'DYZ2', 'DYZ1:DYZ2']
COL = {**CLASS, 'DYZ1:DYZ2': RATIO}
KEY = {n: n.lower() for n in NAMES}
FONTS = [("'ArialMT';", "'Arial';"), ("'Arial-BoldMT';", "'Arial';font-weight:bold;"), ("'Arial-ItalicMT';", "'Arial';font-style:italic;"),
         ("'Arial-BoldItalicMT';", "'Arial';font-weight:bold;font-style:italic;"), ("'HelveticaNeue-Light';", "'Arial';font-weight:300;"),
         ("'HelveticaNeue-Italic';", "'Arial';font-style:italic;"), ("'HelveticaNeue-MediumItalic';", "'Arial';font-weight:500;font-style:italic;")]
DRAWN = {'path', 'text', 'line', 'rect', 'circle', 'ellipse', 'polygon', 'polyline', 'g', 'image', 'use'}
GREY, SHADE, GRIDC = '#55666F', '#F4F6F7', '#E8EDEF'
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
FIG, R = None, None   # current panel figure and region


def svg_panels(tmp, out):
    """Panels a and b from the SVG, as PDF."""
    ns = 'http://www.w3.org/2000/svg'; ET.register_namespace('', ns); ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
    tree = ET.parse(SVG); root = tree.getroot()
    for i, e in enumerate(root): e.set('id', e.get('id') or f'top{i}')   # keep gradient ids
    tree.write(tmp / 'ids.svg')
    q = subprocess.run(['inkscape', str(tmp / 'ids.svg'), '--query-all'], capture_output=True, text=True, check=True).stdout
    ids = {e.get('id') for e in root}; box = {l.split(',')[0]: [float(v) for v in l.split(',')[1:]] for l in q.splitlines() if l.split(',')[0] in ids}
    def part(e):
        if e.tag.split('}')[1] not in DRAWN or e.get('id') not in box: return 'keep'
        x, y, w, h = box[e.get('id')]
        logo = x >= LOGO[0] and y >= LOGO[1] and x + w <= LOGO[2] and y + h <= LOGO[3]
        return 'c-f' if y >= PANELS['b'][4] else 'b' if logo or y >= PANELS['a'][4] - 1 else 'a'
    parts = {e.get('id'): part(e) for e in root}
    s = SVG_A.read_text(); s = re.sub(r'viewBox="[^"]*"', f'viewBox="0 {A_BOX[0]} 540 {A_BOX[1] - A_BOX[0]}"', s, 1); s = re.sub(r'enable-background:new [^;]*;', f'enable-background:new 0 {A_BOX[0]} 540 {A_BOX[1] - A_BOX[0]};', s, 1)
    for a, b in FONTS: s = s.replace(a, b)
    f = out / f'fig1a_{PANELS["a"][0]}.pdf'; (tmp / 'a.svg').write_text(s)
    subprocess.run(['inkscape', str(tmp / 'a.svg'), '--export-type=pdf', f'--export-filename={f}'], check=True, capture_output=True); print(f'panel a: {SVG_A.name}, y {A_BOX[0]}-{A_BOX[1]} -> {f}', flush=True)
    for p in 'b':
        _, x0, y0, x1, y1 = PANELS[p]; t = ET.parse(tmp / 'ids.svg'); r = t.getroot()
        for e in list(r):
            if parts[e.get('id')] not in ('keep', p): r.remove(e)
        low = max((box[e.get('id')][1] + box[e.get('id')][3] for e in r if parts[e.get('id')] == p), default=0)
        high = min((box[e.get('id')][1] for e in r if parts[e.get('id')] == p), default=0)
        if low > y1 or high < y0: sys.exit(f'ERROR panel {p}: elements span y {high:.1f}-{low:.1f}, outside its region {y0}-{y1}')
        r.set('viewBox', f'{x0} {y0} {x1 - x0} {y1 - y0}'); r.set('style', f'enable-background:new {x0} {y0} {x1 - x0} {y1 - y0};')
        s = ET.tostring(r, encoding='unicode')
        for a, b in FONTS: s = s.replace(a, b)
        f = out / f'fig1{p}_{PANELS[p][0]}.pdf'; (tmp / f'{p}.svg').write_text(s)
        subprocess.run(['inkscape', str(tmp / f'{p}.svg'), '--export-type=pdf', f'--export-filename={f}'], check=True, capture_output=True)
        print(f'panel {p}: {sum(v == p for v in parts.values())} SVG elements, y {high:.1f}-{low:.1f} -> {f}', flush=True)
    print(f'SVG: {len(parts)} top-level elements, {sum(v == "c-f" for v in parts.values())} of panels c-e dropped (redrawn here)', flush=True)


SHRINK = {'c': (362.05, 501.05, 14), 'd': (571, 695.5, 14), 'e': (571, 695.5, 14), 'f': (762, 844, 20)}   # plot band (design y) and units removed from it
ROWS = {'c': 0, 'd': 1, 'e': 1, 'f': 2}   # rows below b; shrinks of rows above shift a row up
P_ = None   # current panel


def g(y):
    """Design y to drawn y in the current panel."""
    a, b, sh = SHRINK.get(P_, (0, 1, 0)); return y if y <= a else y - sh * min(y - a, b - a) / (b - a)


def panel(p):
    """New figure for panel p."""
    global FIG, R, P_
    _, x0, y0, x1, y1 = PANELS[p]; R = (x0, y0, x1 - x0, y1 - y0); P_ = p
    FIG = plt.figure(figsize=(R[2] * S / 72, (g(y1) - g(y0)) * S / 72))


XS = 1   # x stretch: ED figures 180 mm wide with Fig. 1 font sizes
fx = lambda x: (x * XS - R[0]) / R[2]
fy = lambda y: 1 - (g(y) - g(R[1])) / (g(R[1] + R[3]) - g(R[1]))   # bands shrunk by SHRINK


def axes(x0, y0, x1, y1):
    """Axes at SVG coordinates."""
    return FIG.add_axes([fx(x0), fy(y1), (x1 - x0) * XS / R[2], fy(y0) - fy(y1)])   # drawn height (SHRINK)


def put(x, y, s, size, **k):
    """Text at SVG coordinates, y is the baseline."""
    FIG.text(fx(x), fy(y), s, fontsize=size * S, va='baseline', **k)


def mark(xs, ys, **k):
    """Line or marker at SVG coordinates."""
    FIG.add_artist(plt.Line2D([fx(x) for x in xs], [fy(y) for y in ys], transform=FIG.transFigure, **k))


def frame(ax, sides='all'):
    for k, v in ax.spines.items(): v.set_visible(sides == 'all' or k in sides); v.set_linewidth(.5)
    ax.tick_params(width=.5, length=2.5, labelsize=8 * S, pad=1.5)


def save(p, out):
    f = out / f'fig1{p}_{PANELS[p][0]}.pdf'; FIG.savefig(f, transparent=True); plt.close(FIG)
    print(f'panel {p}: {R[2] * S:.1f} x {(g(R[1] + R[3]) - g(R[1])) * S:.1f} pt -> {f}', flush=True)


def cv_errors():
    """Per male mean absolute CV error."""
    U = co.defaultdict(dict)
    for r in rd(M / 'cv_units.tsv'): U[r['class']][(r['repeat'], r['male'])] = float(r['error_pct'])
    U['dyz1:dyz2'] = {k: 100 * ((1 + U['dyz1'][k] / 100) / (1 + U['dyz2'][k] / 100) - 1) for k in U['dyz1'] if k in U['dyz2']}
    out = {}
    for n in NAMES:
        m = co.defaultdict(list)
        for (rep, male), v in U[KEY[n]].items(): m[male].append(abs(v))
        if {len(v) for v in m.values()} != {10}: sys.exit(f'ERROR {n}: not 10 repeats per male')
        out[n] = ({k: float(np.mean(v)) for k, v in m.items()}, float(np.median(np.abs(list(U[KEY[n]].values())))))
    return out


def lolo_errors():
    """Per male excluded-lineage absolute error."""
    U = co.defaultdict(dict)
    for r in rd(M / 'lolo_units.tsv'): U[r['class']][r['male']] = float(r['error_pct'])
    U['dyz1:dyz2'] = {m: 100 * ((1 + U['dyz1'][m] / 100) / (1 + U['dyz2'][m] / 100) - 1) for m in U['dyz1'] if m in U['dyz2']}
    return {n: ({m: abs(v) for m, v in U[KEY[n]].items()}, float(np.median(np.abs(list(U[KEY[n]].values()))))) for n in NAMES}


def boxes(E, ylabel, title, letter, cvmed=None, ybot=501.05, x1=521.95, tx=None):
    """Fig. 1c box plot on the current panel figure, log y, every male."""
    ax = axes(37.15, 362.05, x1, ybot); rng = np.random.default_rng(1); J2 = ('HG02492', 'HG02735', 'NA20809')   # short J2a1a1 DYZ19 arrays
    allv = np.concatenate([np.array(list(per.values())) for per, _ in E.values()]); lo, hi, yl = 10 ** np.floor(np.log10(allv.min())), 50, 33   # log axis from the decade below the smallest error; median labels at yl
    assert allv.max() < yl, allv.max(); print(f'{title}: errors {allv.min():.4f}-{allv.max():.2f}%, y axis {lo:g}-{hi}', flush=True)
    for i, n in enumerate(NAMES):
        per, med = E[n]; males = sorted(per); v = np.array([per[m] for m in males]); x = i + rng.uniform(-.24, .24, len(v)); c = COL[n]
        ax.scatter(x, v, s=3.2, color=c, alpha=.55, lw=0, zorder=2)
        q = np.percentile(v, [10, 25, 50, 75, 90])
        ax.bxp([dict(whislo=q[0], q1=q[1], med=q[2], q3=q[3], whishi=q[4], fliers=[])], positions=[i], widths=.55, patch_artist=True,
               boxprops=dict(facecolor=to_rgba(c, .28), edgecolor=c, lw=.75), medianprops=dict(color='k', lw=1.2),
               whiskerprops=dict(color='k', lw=.6), capprops=dict(color='k', lw=.6), capwidths=.22, manage_ticks=False, zorder=3)
        if cvmed: ax.scatter(i + .36, cvmed[n], marker='s', s=10, facecolor='white', edgecolor='k', lw=.6, zorder=6)
        ax.text(i, yl, f'{med:.2f}%', color=c, ha='center', va='center', fontsize=8 * S, zorder=6, bbox=dict(facecolor='white', alpha=.75, edgecolor='none', pad=.6))
        ax.annotate(f'{n}\n(n = {len(v)})', (i, 0), xycoords=('data', 'axes fraction'), xytext=(0, -3), textcoords='offset points',
                    ha='center', va='top', color=c, fontsize=8 * S, fontweight='bold', linespacing=1.15)
        if n == 'DYZ19':
            for m in J2:
                if m in per: ax.scatter(x[males.index(m)], per[m], marker='D', s=12, facecolor='white', edgecolor='k', lw=.6, zorder=5)
        print(f'{title} {n}: {len(v)} males, median of predictions {med:.2f}%, p10-p90 {q[0]:.2f}-{q[4]:.2f}, min {v.min():.3f}, max {v.max():.2f}' + (f', CV median {cvmed[n]:.2f}%' if cvmed else ''), flush=True)
    ax.set_xlim(-.576, 7.53); ax.set_yscale('log'); ax.set_ylim(lo, hi); ax.set_xticks([]); T = [t for t in (.01, .02, .05, .1, .2, .5, 1, 2, 5, 10, 20) if t >= lo]; ax.set_yticks(T, [f'{t:g}' for t in T])   # 1-2-5 per decade
    ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=range(2, 10))); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.tick_params(which='minor', width=.4, length=1.5)
    for t in ax.get_yticklabels(): t.set_fontstyle('italic')
    frame(ax)
    ax.set_ylabel(ylabel, fontsize=10 * S, fontstyle='italic', labelpad=2)
    if letter: put(21.8, 356.9, letter, 20, fontweight='bold')
    put(tx or (37.15 + x1) / 2, 357, title, 12, fontweight='bold', ha='center')   # tx: clear of the diamond key
    yd = 343 if cvmed else 354   # legends stacked when the CV square is shown
    mark([42], [yd], marker='D', ms=4.5, mfc='white', mec='k', mew=.6); put(48.5, yd + 3, '= short J2a1a1 DYZ19 arrays', 8)
    if cvmed: mark([42], [354], marker='s', ms=3.6, mfc='white', mec='k', mew=.6); put(48.5, 357, '= cross-validation median (Fig. 1c)', 8)


def cost(x0, x1, uc):   # uc: design units per box step in c, so boxes match
    """Runtime and peak memory of the released model, 130 local CRAMs."""
    B = rd(ROOT / 'results/checks/runtime_memory.tsv'); rng = np.random.default_rng(4)
    t, m = np.array([float(r['wall_s']) / 60 for r in B]), np.array([float(r['peak_rss_kb']) * 1024 / 1e9 for r in B])   # GNU time %M is KiB
    a = axes(x0, 362.05, x1, 501.05); b = a.twinx(); k = uc / ((x1 - x0) / 1.72)   # c scale over d scale (xlim span 1.72)
    for ax, v, i, c in ((a, t, 0, DOSE), (b, m, 1, GREY)):
        ax.scatter(i + rng.uniform(-.24 * k, .24 * k, len(v)), v, s=3.2, color=c, alpha=.55, lw=0, zorder=2); q = np.percentile(v, [10, 25, 50, 75, 90])
        ax.bxp([dict(whislo=q[0], q1=q[1], med=q[2], q3=q[3], whishi=q[4], fliers=[])], positions=[i], widths=.55 * k, patch_artist=True, manage_ticks=False, zorder=3,
               boxprops=dict(facecolor=to_rgba(c, .28), edgecolor=c, lw=.75), medianprops=dict(color='k', lw=1.2), whiskerprops=dict(color='k', lw=.6), capprops=dict(color='k', lw=.6), capwidths=.22 * k)
        print(f'cost {["runtime min", "peak memory GB"][i]}: {len(v)} CRAMs, median {q[2]:.2f}, p10-p90 {q[0]:.2f}-{q[4]:.2f}, min {v.min():.2f}, max {v.max():.2f}', flush=True)
    a.set_xlim(-.36, 1.36); a.set_xticks([0, 1], ['Runtime', 'Memory']); a.set_ylim(1.1, 2.3); b.set_ylim(3.45, 3.53); a.set_yticks([1.2, 1.4, 1.6, 1.8, 2, 2.2], ['1.2', '1.4', '1.6', '1.8', '2.0', '2.2']); b.set_yticks([3.46, 3.48, 3.5, 3.52], ['3.46', '3.48', '3.50', '3.52'])   # zoomed on the data (axes do not start at 0)   # axes end above the top tick, clear of the title
    for ax, side in ((a, 'left'), (b, 'right')): frame(ax, ('bottom', side))
    for tk in a.get_xticklabels(): tk.set_fontweight('bold')   # black
    a.set_ylabel('Runtime (min)', fontsize=10 * S, fontstyle='italic', labelpad=2); b.set_ylabel('Memory (GB)', fontsize=10 * S, fontstyle='italic', labelpad=2)
    put((x0 + x1) / 2, 357, 'Performance', 12, fontweight='bold', ha='center'); put(x0 - 40, 357, 'd', 20, fontweight='bold')   # letter on the title baseline, clear of the tick labels


def panel_c(out):
    panel('c'); boxes(cv_errors(), 'Mean absolute error (%, log scale)', 'Cross-validation error', 'c', x1=382, tx=232); cost(426, 496, (382 - 37.15) / (7.53 + .576)); save('c', out)   # c xlim span 8.106


def ed4(path):
    """Extended Data Fig. 4: excluded-lineage errors, Fig. 1c style."""
    global R, XS
    global P_
    XS = 180 / 25.4 * 72 / S / 540; panel('c'); P_ = 'ed4'; R = (0, 342, 540 * XS, 270); FIG.set_size_inches(R[2] * S / 72, R[3] * S / 72)   # ED4: no shrink
    E = lolo_errors(); mx = max(max(per.values()) for per, _ in E.values()); print(f'ED4 largest error {mx:.2f}%', flush=True)
    boxes(E, 'Absolute error (%, log scale)', 'Excluded-lineage error', '', ybot=581)   # every point shown, log y as Fig. 1c
    from palette import PAD, crop_ink
    FIG.savefig(path, bbox_inches='tight', pad_inches=PAD); crop_ink(path); plt.close(FIG); XS = 1; print(f'written {path}', flush=True)


def panel_d(out):
    panel('d'); D = {(r['depth'], r['class']): float(r['median_abs']) for r in rd(M / 'depth.tsv')}
    X = {r['fraction']: float(r['median_autosomal_depth']) for r in rd(ROOT / 'results/figures/figure1_depth_summary.tsv')}
    levels = [('d1x', X['d1x']), ('d2x', X['d2x']), ('d5x', X['d5x']), ('d10x', X['d10x']), ('d15x', X['d15x']), ('d20x', X['d20x']), ('full', X['1'])]
    ax = axes(32, 571, 198, 695.5)
    for n in NAMES[:-1]:
        y = [D[(d, KEY[n])] for d, _ in levels]
        ax.plot([x for _, x in levels], y, '-o', color=COL[n], lw=.9 * S, ms=2.6, zorder=3)
        print(f'd {n}: ' + ', '.join(f'{x:.1f}x {v:.2f}' for (_, x), v in zip(levels, y)), flush=True)
    ax.set_xscale('log'); ax.set_xlim(.804, 39.4); ax.set_ylim(0, 5)
    ax.set_xticks([1, 2, 5, 10, 20, X['1']], ['1', '2', '5', '10', '20', f'{X["1"]:.0f}'])
    frame(ax); ax.tick_params(which='minor', width=.4, length=1.5)
    ax.set_xlabel('Median autosomal coverage (×)', fontsize=10 * S, fontstyle='italic', labelpad=2)
    ax.set_ylabel('Median absolute error (%)', fontsize=10 * S, fontstyle='italic', labelpad=4)
    put(22.3, 550.3, 'e', 20, fontweight='bold'); put(54.1, 546, 'Sequencing coverage', 12, fontweight='bold')   # baseline as the e title
    for i, n in enumerate(NAMES[:-1]):
        x, y = [48.7, 89.6, 130.5, 176.7][i // 2], [558, 567][i % 2]
        mark([x - 11, x - 3], [y - 2.6] * 2, color=COL[n], lw=.9 * S); mark([x - 7], [y - 2.6], marker='o', ms=2.6, color=COL[n]); put(x, y, n, 8)
    save('d', out)


def pairs():
    """Parent-offspring percent differences, batch group."""
    E = {r['sample']: r for r in rd(M / 'cohort_1kgp.tsv')}; out = []
    for p, o in (l.split()[:2] for l in open(ROOT / 'data/catalog_d376/fs_pairs.tsv') if l.strip()):
        if p not in E or o not in E or any(E[s]['sex'] != 'male' or E[s]['estimate_basis'] != 'ruler' for s in (p, o)): continue   # both males with a Y
        b = (E[p]['batch'], E[o]['batch']); g = 'same' if b[0] == b[1] else 'cross' if b == ('2504', '698') else 'reverse'
        val = lambda s, n: float(E[s]['DYZ1_bp']) / float(E[s]['DYZ2_bp']) if n == 'DYZ1:DYZ2' else float(E[s][f'{n}_bp'])
        out.append(dict(parent=p, child=o, group=g, **{n: 100 * (val(o, n) / val(p, n) - 1) for n in NAMES}))
    return out


def panel_e(out):
    panel('e'); P = pairs(); g = co.Counter(r['group'] for r in P); rng = np.random.default_rng(2)
    print(f'e: {len(P)} pairs, groups {dict(g)}', flush=True)
    ax = axes(259.3, 564.9, 485.1, 695.5);   # x axis level with d
    ax.set_xscale('symlog', linthresh=5, linscale=1); ax.set_xlim(-220, 220); ax.set_ylim(7.5, -.5)
    special = {'HG02966': '^', 'HG02773': 'o'}
    for i, n in enumerate(NAMES):
        c = COL[n]
        if i % 2 == 0: ax.axhspan(i - .5, i + .5, color=SHADE, lw=0, zorder=0)
        meds = []
        for grp, dy, mk, s, a in (('cross', -.2, 'o', 3, .4), ('same', .2, '^', 7, .9)):
            Rw = [r for r in P if r['group'] == grp or (grp == 'cross' and r['group'] == 'reverse')]   # reverse pairs drawn, not summarised
            v = np.array([r[n] for r in Rw]); plain = np.array([r['parent'] not in special and r['child'] not in special for r in Rw])
            ax.scatter(v[plain], i + dy + rng.uniform(-.09, .09, plain.sum()), s=s, marker=mk, color=c, alpha=a, lw=0, zorder=2)
            q = np.percentile([r[n] for r in P if r['group'] == grp], [25, 50, 75]); meds.append(q[1])
            ax.plot([q[0], q[2]], [i + dy] * 2, color=c, lw=2.6 * S, solid_capstyle='butt', zorder=3)
            ax.plot([q[1]] * 2, [i + dy - .2, i + dy + .2], color='k', lw=1.2 * S, zorder=4)
        for r in P:
            for m, mk in special.items():
                if m in (r['parent'], r['child']): ax.scatter(r[n], i + (-.2 if r['group'] != 'same' else .2), marker=mk, s=16, facecolor='white', edgecolor='k', lw=.75, zorder=5)
        ax.annotate(n, (0, i), xycoords=('axes fraction', 'data'), xytext=(-3, 0), textcoords='offset points', ha='right', va='center', color=c, fontsize=8 * S, fontweight='bold')
        ax.annotate(f'{meds[0]:+.2f}/{meds[1]:+.2f}', (1, i), xycoords=('axes fraction', 'data'), xytext=(3, 0), textcoords='offset points',
                    ha='left', va='center', color=c, fontsize=8 * S, fontweight='bold')
        print(f'e {n}: median cross {meds[0]:+.2f}, same {meds[1]:+.2f}', flush=True)
    ticks = [-100, -20, -5, 0, 5, 20, 100]; ax.set_xticks(ticks, [f'{t}'.replace('-', '−') for t in ticks]); ax.set_yticks([])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for t in ticks: ax.axvline(t, color=GRIDC, lw=.6 * S, zorder=1)
    ax.axvline(0, color=GREY, lw=.6 * S, ls=(0, (3, 2)), zorder=1)
    frame(ax, ('left', 'bottom'))
    ax.set_xlabel('Relative offspring–parent difference (%)', fontsize=10 * S, fontstyle='italic', labelpad=1)
    put(224.2, 551.1, 'f', 20, fontweight='bold'); put(372, 546, f'Parent–offspring agreement ({len(P)} pairs)', 12, fontweight='bold', ha='center')   # title size as c, d, f
    for x, mk, s in ((323.5, '^', 'HG02966'), (378, 'o', 'HG02773')):
        mark([x], [555.2], marker=mk, ms=4, mfc='white', mec='k', mew=.75); put(x + 5.5, 557.9, f'= {s}', 8)
    put(511, 553.2, 'Median (%)', 8, color=GREY, fontweight='bold', ha='center'); put(511, 562.7, 'Cross/Same', 8, color=GREY, fontweight='bold', ha='center')
    save('e', out)


def panel_f(out):
    panel('f')
    F = [r for r in rd(M / 'cohort_fem.tsv') if r['cohort'] == '1kGP']; males = [r for r in rd(M / 'cohort_1kgp.tsv') if r['sex'] == 'male']
    fy_, my = np.array([float(r['y_dose']) for r in F]), np.array([float(r['y_dose']) for r in males]); cut, top = .02, .03   # no-Y cut (user 2026-11-04: above every female, max 0.0191 in HGDP), log part end
    print(f'f: 1kGP females {len(fy_)} (max Y dosage {fy_.max():.4f}, above the cut {cut}: {(fy_ > cut).sum()}, at or below 0.0001: {(fy_ <= 1e-4).sum()}), males {len(my)} (min {my.min():.4f}, at or below the cut: {(my <= cut).sum()})', flush=True)
    gap = ((fy_ > top) & (fy_ < .2)).sum() + ((my > top) & (my < .2)).sum(); lo, hi = (my >= .2).sum(), (my > 2.2).sum()
    print(f'f: broken x axis, samples in the hidden {top}-0.2 gap: {gap}, males from 0.2: {lo}, above 2.2 (clipped into the last bin): {hi}', flush=True)
    if gap: sys.exit('ERROR samples in the hidden part of the dosage axis')
    T1, B, X1, X2 = 751, 844, 135, 141; HB = g(B) - g(T1)   # top, bottom, axis break; HB drawn axes height
    a1, a2 = axes(52, T1, X1, B), axes(X2, T1, 270, B)   # log part (no Y), linear part (Y)
    for ax, bins, clip in ((a1, np.logspace(-4, np.log10(top), 26), (1.01e-4, top)), (a2, np.arange(.2, 2.225, .025), (.2, 2.2))):
        for v, c, lab in ((fy_, '#B07AA1', f'Recorded females (n = {len(fy_)})'), (my, DOSE, f'Recorded males (n = {len(my):,})')):
            w = v[(v >= clip[0] * .99) & (v <= top)] if ax is a1 else v[v >= .2]
            ax.hist(np.clip(w, *clip), bins, color=c, label=lab, zorder=2)
        ax.set_yscale('log'); ax.set_ylim(.7, 5e3); frame(ax, ('left', 'bottom') if ax is a1 else ('bottom',))
    a1.set_xscale('log'); a1.set_xlim(1e-4, top); a1.set_xticks([1e-4, 1e-3, 1e-2], ['0.0001', '0.001', '0.01']); a1.minorticks_off()
    a1.set_yticks([1, 10, 100, 1000], ['1', '10', '100', '1000']); a1.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    a2.set_xlim(.2, 2.2); a2.set_xticks([.5, 1, 1.5, 2], ['0.5', '1', '1.5', '2']); a2.set_xticks([x / 10 for x in range(2, 23) if x % 5], minor=True)   # minor ticks every 0.1
    a2.tick_params(which='minor', width=.4, length=1.5); a2.set_yticks([]); a2.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for x in (X1, X2): mark([x - 1.5, x + 1.5], [B + 2.5, B - 2.5], color='k', lw=.5)   # axis break
    a1.axvline(cut, color=GREY, lw=.6 * S, ls=(0, (3, 2)), zorder=1); a1.text(cut / 1.15, 1.6e3 * 10 ** (3 / HB * np.log10(5e3 / .7)), f'{cut:g}', fontsize=8 * S, color='k', va='center', ha='right')
    for r in males:
        if float(r['y_dose']) <= cut:   # recorded male without Y
            e = np.logspace(-4, np.log10(top), 26); j = np.searchsorted(e, float(r['y_dose'])) - 1   # its bin
            a1.annotate(r['sample'], (np.sqrt(e[j] * e[j + 1]), np.sqrt(.7)), xytext=(.0043, (g(B) - g(865)) / HB), textcoords=('data', 'axes fraction'), annotation_clip=False, fontsize=8 * S, color='k', va='baseline', ha='center',   # on the 'Y dosage' baseline
                        arrowprops=dict(arrowstyle='-|>', color='k', lw=.5, mutation_scale=4, shrinkA=1, shrinkB=0))
            print(f"f: labelled {r['sample']} (recorded male, Y dosage {r['y_dose']}, X dosage {r['x_dose']}, {r['sex_chromosomes']})", flush=True)
    a2.legend(loc='upper left', fontsize=8 * S, handlelength=.8, frameon=False, bbox_to_anchor=(.0, 1 + 5 / HB), borderaxespad=0, borderpad=0, labelspacing=.3)   # 5 units above the axes top, clear of the bars
    put(161, 865, 'Y dosage', 10, fontstyle='italic', ha='center'); a1.set_ylabel('Samples', fontsize=10 * S, fontstyle='italic', labelpad=2)
    put(161, 743, '1000 Genomes samples', 10, fontweight='bold', ha='center')
    T = {r['class']: r for r in rd(M / 'cv_females.tsv')}; CV = {r['class']: float(r['median_abs']) for r in rd(M / 'cv.tsv')}
    ax = axes(312, 762, 530, B); n7 = NAMES[:-1]
    for i, n in enumerate(n7):
        r = T[KEY[n]]; md, mx = float(r['median_pct_of_median_truth']), float(r['max_pct_of_median_truth']); c = COL[n]
        ax.bar(i - .19, mx, .36, facecolor='white', edgecolor=c, hatch='//////', lw=.75, zorder=2)   # females hatched, males as in c
        ax.bar(i + .19, CV[KEY[n]], .36, facecolor=to_rgba(c, .28), edgecolor=c, lw=.75, zorder=2)
        print(f'f {n}: female median {md}% max {mx}% of the median length ({r["n_predictions"]} predictions); CV median error {CV[KEY[n]]}%', flush=True)
    ax.set_ylim(0, 4.2); ax.set_xlim(-.6, len(n7) - .4); ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xticks(range(len(n7)), n7)
    for t, n in zip(ax.get_xticklabels(), n7): t.set_color(COL[n]); t.set_fontweight('bold')
    frame(ax, ('left', 'bottom')); ax.tick_params(axis='x', length=0)
    ax.set_ylabel('% of median\nsatellite length', fontsize=10 * S, fontstyle='italic', labelpad=2, linespacing=1.1)
    put(421, 743, f'Signal in {int(T["dyz3"]["n_predictions"]) // 10} held-out females', 10, fontweight='bold', ha='center')
    for x, txt, h in ((330, 'highest female signal', '//////'), (423, 'error in males (CV median)', None)):
        FIG.add_artist(matplotlib.patches.Rectangle((fx(x - 3), fy(754.8)), 6 / R[2], fy(748.8) - fy(754.8), transform=FIG.transFigure, facecolor='white' if h else to_rgba(GREY, .28), edgecolor=GREY, hatch=h, lw=.75))
        put(x + 5.5, 754.4, txt, 8)
    put(22.3, 735, 'g', 20, fontweight='bold'); put(270, 735, 'Specificity', 12, fontweight='bold', ha='center')
    save('f', out)


def main():
    if sys.argv[1:2] == ['ed4']: return ed4(ROOT / 'manuscript/figures/ed/ed4_lineage_validation.pdf')
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'manuscript/figures/main'; pd = out / 'fig1'; pd.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as t: svg_panels(Path(t), pd)
    panel_c(pd); panel_d(pd); panel_e(pd); panel_f(pd)
    up = lambda p: (A_UP if p != 'a' else 0) + sum(SHRINK[q][2] for q in ('c', 'd', 'f') if ROWS.get(q, -1) < ROWS.get(p, -1))   # shift from a and rows above
    Hn = H - A_UP - sum(SHRINK[q][2] for q in ('c', 'd', 'f')) + 3; w = PdfWriter(); page = w.add_blank_page(W * S, Hn * S)   # +3 units: 1.5 mm bottom margin
    for p, (name, x0, y0, x1, y1) in PANELS.items():
        sub = PdfReader(pd / f'fig1{p}_{name}.pdf').pages[0]; sw, sh = float(sub.mediabox.width), float(sub.mediabox.height); hn = A_BOX[1] - A_BOX[0] if p == 'a' else y1 - y0 - SHRINK.get(p, (0, 0, 0))[2]
        if abs(sw - (x1 - x0) * S) > .5 or abs(sh - hn * S) > .5: sys.exit(f'ERROR panel {p}: {sw:.1f} x {sh:.1f} pt, region {(x1 - x0) * S:.1f} x {hn * S:.1f}')
        page.merge_transformed_page(sub, Transformation().translate(x0 * S, (Hn - (y0 - up(p)) - hn) * S))
    f = out / 'fig1_method.pdf'; w.write(f); from palette import crop_ink; crop_ink(f)   # 1.5 mm ink margin
    subprocess.run(['inkscape', str(f), '--export-type=svg', f'--export-filename={f.with_suffix(".svg")}'], check=True, capture_output=True)
    print(f'written {f} and .svg: {W * S:.1f} x {Hn * S:.1f} pt ({W * S / 72 * 25.4:.1f} x {Hn * S / 72 * 25.4:.1f} mm), {len(PANELS)} panels', flush=True)


if __name__ == '__main__':
    main()
