#!/usr/bin/env python3
"""Supplementary figure: aligned read density under two aligners, same reads.

usage: python3 scripts/figure_aligner_depth.py <alignchk.tsv> <out.pdf>
"""
from palette import CLASS, MM, FS, save, ORDER, INK, MUTED
import sys, csv, collections as co
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
IN, OUT = sys.argv[1:3]; G = co.defaultdict(list); N = [(k.lower(), k, CLASS[k]) for k in ORDER if k in ('DYZ3', 'DYZ2', 'DYZ19', 'DYZ1')]
for r in csv.DictReader(open(IN), delimiter='\t'): G[r['class']].append(float(r['dragen_vs_bwa_pct']))
fig, ax = plt.subplots(figsize=(88 * MM, 60 * MM)); fig.subplots_adjust(left=.14, right=.98, top=.97, bottom=.1)
for i, (k, n, c) in enumerate(N):
    v = np.array(G[k]); ax.scatter(i + np.random.default_rng(i).uniform(-.28, .28, len(v)), v, s=4, c=c, lw=0, alpha=.65); m = np.median(v); ax.plot([i - .38, i + .38], [m, m], color=INK, lw=1); print(n, len(v), f'median {m:+.2f}', flush=True)
ax.axhline(0, color=MUTED, lw=.6, ls='--'); ax.text(3.45, .35, 'Yruler counts: identical', ha='right', va='bottom', fontsize=FS['point'], color=INK)
ax.set_xticks(range(4), [n for _, n, _ in N]); ax.set_ylabel('Aligned read density, DRAGEN vs BWA (%)'); save(fig, OUT)
