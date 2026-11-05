#!/usr/bin/env python3
"""Supplementary Fig. 1: aligned read depth, which uses no Yruler marker, against the Yruler estimates of DYZ3 and DYZ19 (Fig. 1 style).
usage: python3 scripts/suppfig1_markerfree.py [out.pdf]   (default manuscript/figures/supp/suppfig1_markerfree.pdf)
inputs: results/checks/dyz3_check.tsv (aligned_ratio: GRCh38 DYZ3 block over chrY:14.5-15.5 Mb),
        results/checks/dyz19_depth.tsv (depth_ratio_core: chrY:20,054,500-20,351,000 over chrY:14.5-15.5 Mb, um_dyz19depth.py with v9bBx),
        results/checks/alignchk.tsv (same reads aligned by DRAGEN 3.7.6 and NYGC BWA: aligned depth per class, um_alignchk.py, 7 classes),
        Yruler estimates and population set from figures_1kgp.load_snapshot (final model v9bBx)."""
import csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figures_1kgp as f
from palette import CLASS, LINEAGE, MUTED, RATIO
from matplotlib.colors import to_rgba

ROOT = Path(__file__).resolve().parents[1]
S, H = 0.75, 589   # pt per font unit, page height in Fig. 1 units
XS = 180 / 25.4 * 72 / S / 540; W = 540 * XS   # 180 mm wide, Fig. 1 font sizes
GREY = '#B8C0C6'
rd = lambda p: {r['sample']: r for r in csv.DictReader(open(p), delimiter='\t')}
fig = plt.figure(figsize=(W * S / 72, H * S / 72))


def axes(x0, y0, x1, y1):
    return fig.add_axes([x0 * XS / W, 1 - y1 / H, (x1 - x0) * XS / W, (y1 - y0) / H])


def put(x, y, s, size, **k):
    fig.text(x * XS / W, 1 - y / H, s, fontsize=size * S, va='baseline', **k)


def style(ax, xl, yl):
    for v in ax.spines.values(): v.set_linewidth(.5); v.set_visible(True)   # full frame as Fig. 1d
    ax.tick_params(width=.5, length=2.5, labelsize=8 * S, pad=1.5)
    ax.set_xlabel(xl, fontsize=10 * S, fontstyle='italic', labelpad=2); ax.set_ylabel(yl, fontsize=10 * S, fontstyle='italic', labelpad=2)


def panel(ax, groups, xk, yk, scale, title, letter, x0, x1, xl, yl):
    allx, ally = [], []
    for lab, rows, c, a, s, z in groups:
        x = np.array([float(r[xk]) / scale for r in rows]); y = np.array([float(r[yk]) for r in rows]); allx += list(x); ally += list(y)
        ax.scatter(x, y, s=s, color=c, alpha=a, lw=0, zorder=z, label=f'{lab} ({len(rows):,})')
        print(f'{title} {lab}: {len(rows)} males, estimate {x.min():.3g}-{x.max():.3g}, read depth ratio {y.min():.3f}-{y.max():.3f}', flush=True)
    rho = spearmanr(allx, ally).statistic; n = len(allx)
    ax.text(.03, .96, f'Spearman ρ = {rho:.4f} ({n:,} males)', transform=ax.transAxes, va='top', fontsize=8 * S)
    ax.legend(loc='lower right', fontsize=8 * S, frameon=False, handletextpad=.2, borderaxespad=.3, markerscale=1.5)
    style(ax, xl, yl); put((x0 + x1) / 2, 22, title, 12, fontweight='bold', ha='center'); put(x0 - 44, 22, letter, 20, fontweight='bold')
    print(f'{title}: Spearman {rho:.4f}, {n} males', flush=True)
    return rho, n


def aligner(ax):
    """Aligned depth, DRAGEN against BWA on the same reads."""
    A = list(csv.DictReader(open(ROOT / 'results/checks/alignchk.tsv'), delimiter='\t')); cl = ['DYZ3', 'DYZ17', 'DYZ19', 'Yq12', 'DYZ18S', 'DYZ1', 'DYZ2', 'DYZ1:DYZ2']; rng = np.random.default_rng(3)   # Fig. 1c order
    V = {(r['sample'], r['class']): r for r in A}; males = sorted({r['sample'] for r in A})
    ratio = np.array([100 * ((float(V[m, 'dyz1']['dragen']) / float(V[m, 'dyz2']['dragen'])) / (float(V[m, 'dyz1']['bwa']) / float(V[m, 'dyz2']['bwa'])) - 1) for m in males])   # DYZ1:DYZ2 depth ratio
    for i, n in enumerate(cl):
        v = ratio if n == 'DYZ1:DYZ2' else np.array([float(r['dragen_vs_bwa_pct']) for r in A if r['class'] == n.lower()]); c = RATIO if n == 'DYZ1:DYZ2' else CLASS[n]
        ax.scatter(i + rng.uniform(-.22, .22, len(v)), v, s=4, color=c, alpha=.55, lw=0, zorder=2)
        q = np.percentile(v, [10, 25, 50, 75, 90])
        ax.bxp([dict(whislo=q[0], q1=q[1], med=q[2], q3=q[3], whishi=q[4], fliers=[])], positions=[i], widths=.45, patch_artist=True, manage_ticks=False, zorder=3,
               boxprops=dict(facecolor=to_rgba(c, .28), edgecolor=c, lw=.75), medianprops=dict(color='k', lw=1.2), whiskerprops=dict(color='k', lw=.6), capprops=dict(color='k', lw=.6), capwidths=.22)
        print(f'aligner {n}: {len(v)} males, DRAGEN vs BWA aligned depth median {q[2]:.2f}%, range {v.min():.2f} to {v.max():.2f}%', flush=True)
    ax.set_xticks(range(len(cl)), cl); ax.set_xlim(-.6, len(cl) - .4)
    for t, n in zip(ax.get_xticklabels(), cl): t.set_color(RATIO if n == 'DYZ1:DYZ2' else CLASS[n]); t.set_fontweight('bold')
    style(ax, '', 'DRAGEN − BWA aligned depth (%)'); put(291, 347, f'Aligner dependence ({len(v)} males)', 12, fontweight='bold', ha='center'); put(8, 347, 'c', 20, fontweight='bold')


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'manuscript/figures/supp/suppfig1_markerfree.pdf'
    _, valid, _ = f.load_snapshot(); E = {r['sample']: r for r in valid if r.get('population_included') in (True, 'True', 'yes')}
    print(f'population males {len(E)}', flush=True)
    D3 = rd(ROOT / 'results/checks/dyz3_check.tsv'); D19 = rd(ROOT / 'results/checks/dyz19_depth.tsv')
    r3 = [dict(D3[s], est=E[s]['DYZ3_bp'], hg=E[s]['haplogroup'], lin=E[s]['lineage']) for s in E if s in D3]
    if len(r3) != len(E): print(f'WARNING {len(E) - len(r3)} population males without a DYZ3 read depth', flush=True)
    o1 = [r for r in r3 if r['hg'].startswith('O2a2b1a1')]; o2 = [r for r in r3 if r['hg'].startswith('O2a2b1a2')]; r1b = [r for r in r3 if r['lin'] == 'R1b']
    rest = [r for r in r3 if r not in o1 and r not in o2 and r not in r1b]
    a = axes(52, 34, 262, 299)   # square: 210 x units = 265 y units
    panel(a, [('other lineages', rest, GREY, .5, 5, 1), ('R1b', r1b, LINEAGE['R'], .7, 7, 2), ('O2a2b1a2', o2, LINEAGE['O'], .9, 12, 3), ('O2a2b1a1', o1, CLASS['DYZ3'], .9, 12, 4)],
          'est', 'aligned_ratio', 1e6, 'DYZ3', 'a', 52, 262, 'Yruler DYZ3 estimate (Mb)', 'Aligned read depth, DYZ3 / control')
    r19 = [dict(r, est=r['DYZ19_bp']) for s, r in D19.items() if s in E]
    for nm, rr, k in (('DYZ3', r3, 'aligned_ratio'), ('DYZ19', r19, 'depth_ratio_core')):   # least-squares line for the caption
        x, y = np.array([float(r['est']) for r in rr]), np.array([float(r[k]) for r in rr]); sl, ic = np.polyfit(x, y, 1); r2 = 1 - np.var(y - sl * x - ic) / np.var(y)
        print(f'{nm} fit: {k} = estimate / {1 / sl / 1e3:.1f} kb + {ic:.3f}, R2 {r2:.4f}, {len(x)} males, median |residual| {100 * np.median(np.abs(y - sl * x - ic) / y):.2f}%', flush=True)
    j1 = [r for r in r19 if r['haplogroup'].startswith('J2a1a1')]; jo = [r for r in r19 if r['lineage'] == 'J' and r not in j1]; ot = [r for r in r19 if r['lineage'] != 'J']
    b = axes(320, 34, 530, 299)
    panel(b, [('other lineages', ot, GREY, .5, 5, 1), ('other J', jo, LINEAGE['J'], .6, 8, 2), ('J2a1a1', j1, CLASS['DYZ19'], .9, 10, 3)],
          'est', 'depth_ratio_core', 1e3, 'DYZ19', 'b', 320, 530, 'Yruler DYZ19 estimate (kb)', 'Aligned read depth, DYZ19 core / control')
    aligner(axes(52, 359, 530, 569))
    from palette import PAD, crop_ink
    fig.savefig(out, bbox_inches='tight', pad_inches=PAD); crop_ink(out); print(f'written {out}: {W * S / 72 * 25.4:.1f} x {H * S / 72 * 25.4:.1f} mm', flush=True)


if __name__ == '__main__':
    main()
