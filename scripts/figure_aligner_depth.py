#!/usr/bin/env python3
"""Supplementary figure: aligned read density under two aligners, same reads.

usage: python3 scripts/figure_aligner_depth.py <alignchk.tsv> <out.pdf>
"""
import sys, csv, collections as co
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
IN, OUT = sys.argv[1:3]; G = co.defaultdict(list); N = [('dyz3', 'DYZ3', '#277da8'), ('dyz2', 'DYZ2', '#8c7139'), ('dyz19', 'DYZ19', '#9663a0'), ('dyz1', 'DYZ1', '#338565')]
for r in csv.DictReader(open(IN), delimiter='\t'): G[r['class']].append(float(r['dragen_vs_bwa_pct']))
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'pdf.fonttype': 42})
fig, ax = plt.subplots(figsize=(4.6, 3.1))
for i, (k, n, c) in enumerate(N):
    v = np.array(G[k]); ax.scatter(i + np.random.default_rng(i).uniform(-.28, .28, len(v)), v, s=8, c=c, lw=0, alpha=.65); m = np.median(v); ax.plot([i - .38, i + .38], [m, m], color='#22313a', lw=1.3); print(n, len(v), f'median {m:+.2f}', flush=True)
ax.axhline(0, color='#55666f', lw=.9, ls='--'); ax.text(3.45, .35, 'yruler counts: identical', ha='right', va='bottom', fontsize=7, color='#22313a')
ax.set_xticks(range(4), [n for _, n, _ in N]); ax.set_ylabel('Aligned read density, DRAGEN vs BWA (%)'); ax.spines[['top', 'right']].set_visible(False); ax.tick_params(length=2.5)
fig.tight_layout(); fig.savefig(OUT); print('written', OUT, flush=True)
