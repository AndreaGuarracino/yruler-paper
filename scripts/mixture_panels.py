#!/usr/bin/env python3
"""Extended Data Fig. 7 b-d: read-level male-female mixtures (v2 release).

draw(axes, mix, coh, tab, band): b relative difference in measured Y dosage; c, d length change against measured Y dosage at 5x and 1x,
grey band below the no-Y cut (Yruler reports no lengths there). Called by figures_1kgp.specificity_supplement.
mixtures: 6 male-female pairs, male read fractions 0, 0.01, 0.02, 0.1, 0.5, 1, mixture depths 5, 1 and 0.2x, 3 seeds (um_mixreads.sh, um_mixreads_sum.py).
expected Y dosage = male fraction x Y dosage of the male in the cohort table; change = 100 x (mixture estimate / whole-genome estimate of the male - 1).
"""
import csv, numpy as np
from palette import CLASS, FS, GRID, head
C = ['DYZ3', 'DYZ17', 'DYZ19', 'Yq12', 'DYZ18S', 'DYZ1', 'DYZ2']; DC = {10.0: '#4a9a6e', 5.0: '#1f4e79', 1.0: '#d36a27'}


def draw(axes, mix, coh, tab, band=.025, xlim=(5e-3, 1.2)):
    """Three mixture panels; summary rows to tab."""
    R = list(csv.DictReader(open(mix), delimiter='\t')); X = {r['sample']: r for r in csv.DictReader(open(coh), delimiter='\t')}
    runs = {(r['male'], r['female'], r['fraction'], r['depth'], r['seed']): r for r in R if r['class'] == 'DYZ3'}   # one row per run for the dosage
    D = [d for d in (10.0, 5.0, 1.0) if any(float(k[3]) == d for k in runs)]; FR = sorted({float(k[2]) for k in runs if float(k[2]) > 0}); rows = []; a = axes[0]
    for d in D:   # measured Y dosage relative to expected
        sel = [(k, r) for k, r in runs.items() if float(k[3]) == d and float(k[2]) > 0]
        e = np.array([float(k[2]) * float(X[k[0]]['y_dose']) for k, r in sel]); m = np.array([float(r['y_dose']) for k, r in sel])
        style = {10.0: dict(marker='^', s=14, color=DC[d], linewidths=0), 5.0: dict(marker='o', s=14, color=DC[d], linewidths=0),
                 1.0: dict(marker='x', s=20, color=DC[d], linewidths=.7)}[d]
        a.scatter(e, 100 * (m / e - 1), alpha=.85, label=f'{d:g}x mixture depth', **style)
        for f in FR:
            s = np.array([abs(float(r['y_dose']) - float(k[2]) * float(X[k[0]]['y_dose'])) for k, r in sel if float(k[2]) == f])
            rows.append(('dosage', d, f, len(s), f'{np.median(s):.4f}', f'{s.max():.4f}'))
            print(f'mixtures {d:g}x fraction {f:g}: {len(s)} runs, median |measured - expected| {np.median(s):.4f}, max {s.max():.4f}', flush=True)
        f0 = np.array([float(r['y_dose']) for k, r in runs.items() if float(k[3]) == d and float(k[2]) == 0]); rows.append(('female_only', d, 0, len(f0), f'{np.median(f0):.4f}', f'{f0.max():.4f}'))
        print(f'mixtures {d:g}x female only: {len(f0)} runs, Y dosage median {np.median(f0):.4f}, max {f0.max():.4f}', flush=True)
    a.axvspan(xlim[0], band, color='#e6e6e6', lw=0, zorder=0); a.axhline(0, color='grey', lw=.5, zorder=0)
    a.axvline(band, color='grey', lw=.6, ls=(0, (3, 2)), zorder=1)
    a.annotate(f'{band:g}', (band, 1), xycoords=('data', 'axes fraction'), xytext=(-3, -4), textcoords='offset points', ha='right', va='top', fontsize=FS['point'])
    a.set_xscale('log'); a.set_xlim(*xlim); a.set_ylim(-30, 30); a.set_yticks(range(-30, 31, 10)); a.set_xticks([.001, .01, .1, 1], ['0.001', '0.01', '0.1', '1']); a.minorticks_off()
    a.set_xlabel('Expected Y dosage (male fraction x male Y dosage)', fontstyle='italic'); a.set_ylabel('Relative Y dosage\ndifference (%)', fontstyle='italic'); a.legend(fontsize=FS['legend'], frameon=False, loc='upper right')
    head(a, 'Y dosage of male-female read mixtures')
    for ax, d in zip(axes[1:], D):   # length change against Y dosage, one panel per depth (10x, 5x, 1x)
        for n in C:
            x, med, lo, hi = [], [], [], []
            for f in FR:
                rr = [r for r in R if r['class'] == n and float(r['depth']) == d and float(r['fraction']) == f and float(r['y_dose']) > band]   # Yruler reports lengths only above the no-Y cut
                if len(rr) < 9: continue   # fraction mostly at or below the cut
                v = np.array([float(r['change_pct']) for r in rr])
                x.append(np.median([float(r['y_dose']) for r in rr])); med.append(np.median(v)); lo.append(np.percentile(v, 10)); hi.append(np.percentile(v, 90))
                rows.append(('change', d, f, len(v), n, f'{np.median(v):+.2f}', f'{np.percentile(v, 10):+.2f}', f'{np.percentile(v, 90):+.2f}'))
            ax.plot(x, med, color=CLASS[n], lw=1, marker='o', ms=2, label=n); ax.fill_between(x, lo, hi, color=CLASS[n], alpha=.1, lw=0)
            print(f'mixtures {d:g}x {n}: ' + ', '.join(f'dosage {x_:.3f} {m:+.1f}%' for x_, m in zip(x, med)), flush=True)
        ax.axvspan(xlim[0], band, color='#e6e6e6', lw=0, zorder=0)
        ax.axvline(band, color='grey', lw=.6, ls=(0, (3, 2)), zorder=1)
        ax.annotate(f'{band:g}', (band, 1), xycoords=('data', 'axes fraction'), xytext=(-3, -4), textcoords='offset points', ha='right', va='top', fontsize=FS['point'])
        ax.axhline(0, color='grey', lw=.5); ax.set_xscale('log'); ax.set_xlim(*xlim); ax.set_xticks([.001, .01, .1, 1], ['0.001', '0.01', '0.1', '1']); ax.minorticks_off(); ax.set_ylim(-30, 40); ax.set_yticks(range(-30, 41, 10))
        ax.set_ylabel('Change from the\nwhole genome (%)', fontstyle='italic'); ax.set_xlabel('Measured Y dosage', fontstyle='italic'); head(ax, f'Satellite estimates at {d:g}x mixture depth')
    for ax in axes[1:]:
        ax.legend(fontsize=FS['legend'], frameon=False, loc='upper right', ncol=4)
    for ax in axes: ax.spines[['top', 'right']].set_visible(False); ax.grid(color=GRID, lw=.5); ax.set_axisbelow(True)
    with open(tab, 'w') as o: o.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
    print('written', tab, flush=True)
