#!/usr/bin/env python3
"""Fig. 1b: CV errors per class and of the DYZ1:DYZ2 ratio; ratio errors and estimates by lineage.

usage: python3 scripts/figure1b_ratio.py <cv_units.tsv> <cohort_1kgp.tsv> <out.pdf> <lineage table out.tsv>
DYZ1:DYZ2 ratio error per prediction = (1 + e_DYZ1) / (1 + e_DYZ2) - 1 from the same male, repeat and fold.
Unit ratio = base-pair ratio x U2 / U1, Hallast 2026 units (one HMMER hit): median bp per hit DYZ1 3,478.9, DYZ2 2,404.4 in 83 assemblies (their filter).
"""
import sys, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CVU, COH, OUT, TAB = sys.argv[1:5]
GRAY, ORANGE, INK, MUTED = '#a8a7a2', '#eb6834', '#0b0b0b', '#52514e'
CL = [('dyz3', 'DYZ3'), ('yq12', 'Yq12'), ('dyz1', 'DYZ1'), ('dyz2', 'DYZ2'), ('dyz18s', 'DYZ18S'), ('dyz19', 'DYZ19'), ('dyz17', 'DYZ17')]
E = collections.defaultdict(dict)
for l in list(open(CVU))[1:]:
    c, r, m, e = l.rstrip('\n').split('\t'); E[c][(r, m)] = float(e)
pm = lambda d: {m: np.mean([v for (r, mm), v in d.items() if mm == m]) for m in {k[1] for k in d}}   # per-male mean over repeats
R = {k: 100 * ((1 + E['dyz1'][k] / 100) / (1 + E['dyz2'][k] / 100) - 1) for k in E['dyz1'] if k in E['dyz2']}
L = [l.rstrip('\n').split('\t') for l in open(COH)]; H = L[0]; X = [dict(zip(H, x)) for x in L[1:]]
lin = {x['sample']: x['lineage'] for x in X}
inc = [x for x in X if x['estimate_basis'] == 'ruler' and x['sex'] != 'female' and x['y_state'] in ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain')]
print(f'CV predictions: ratio {len(R)}; cohort men included {len(inc)}', flush=True)
plt.rcParams.update({'font.size': 7.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED})
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7), constrained_layout=True, gridspec_kw={'width_ratios': [1.5, 1, 1.15]})
pma = lambda d: {m: np.mean([abs(v) for (r, mm), v in d.items() if mm == m]) for m in {k[1] for k in d}}   # per-male mean |error|
a = ax[0]; data = [np.array(list(pma(E[c]).values())) for c, _ in CL] + [np.array(list(pma(R).values()))]; names = [n for _, n in CL] + ['DYZ1:DYZ2']
for i, v in enumerate(data):   # per-male mean |error| over 10 repeats
    col = ORANGE if i == len(data) - 1 else GRAY
    a.scatter(i + np.random.default_rng(i).uniform(-0.28, 0.28, len(v)), v, s=3, c=col, lw=0, alpha=0.8)
    a.plot([i - 0.32, i + 0.32], [np.median(v)] * 2, color=INK, lw=1)
    print(f'{names[i]}: males {len(v)}, median of per-male mean |error| {np.median(v):.2f}%, p90 {np.percentile(v, 90):.2f}%, max {v.max():.2f}%', flush=True)
a.set_xticks(range(len(names)), names, rotation=45, ha='right', fontsize=6.5); a.set_ylabel('|CV error| per male (%; mean of 10 repeats)'); a.set_ylim(0, None)
g = collections.defaultdict(list)
for m, v in pm(R).items(): g[lin.get(m, '?')].append(v)
big = [k for k in sorted(g, key=lambda k: -len(g[k])) if len(g[k]) >= 4 and k != '?']; oth = [v for k in g if k not in big for v in g[k]]
a = ax[1]; groups = [(k, g[k]) for k in big] + [('other', oth)]
for i, (k, v) in enumerate(groups):
    a.scatter(i + np.random.default_rng(10 + i).uniform(-0.25, 0.25, len(v)), v, s=4, c=ORANGE, lw=0, alpha=0.85); a.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=INK, lw=1)
    print(f'ratio signed CV error, lineage {k}: males {len(v)}, median {np.median(v):+.2f}%', flush=True)
a.axhline(0, color=MUTED, lw=0.6, ls=':'); a.set_xticks(range(len(groups)), [f'{k}\n{len(v)}' for k, v in groups], fontsize=6.5); a.set_ylabel('DYZ1:DYZ2 signed CV error (%)')
U1, U2 = 3478.9, 2404.4   # bp per Hallast unit (hit), DYZ1 and DYZ2
u = {x['sample']: (float(x['DYZ1_bp']) / U1) / (float(x['DYZ2_bp']) / U2) for x in inc}
gl = collections.defaultdict(list)
for x in inc: gl[x['lineage']].append(x['sample'])
import csv
KIN = {r['sample']: r['kinship_component_id'] for r in csv.DictReader(open('results/estimates_1kgp.relatedness_samples.tsv'), delimiter='\t')}
rng = np.random.default_rng(20261101); rows = []
for k, v in gl.items():
    if len(v) < 20 or k == '.': continue
    cp = collections.defaultdict(list)
    for s in v: cp[KIN.get(s, s)].append(u[s])   # resample kinship components
    cl = list(cp.values()); b = [np.median(np.concatenate([cl[i] for i in rng.integers(0, len(cl), len(cl))])) for _ in range(2000)]; rows.append((np.median([u[s] for s in v]), k, len(v), np.percentile(b, 2.5), np.percentile(b, 97.5)))
rows.sort()
with open(TAB, 'w') as o:
    o.write('lineage\tmen\tunit_ratio_median\tci95_low\tci95_high\n'); o.writelines(f'{k}\t{n}\t{m:.4f}\t{lo:.4f}\t{hi:.4f}\n' for m, k, n, lo, hi in rows)
a = ax[2]
for i, (m, k, n, lo, hi) in enumerate(rows):
    a.plot([lo, hi], [i, i], color=MUTED, lw=1.2); a.scatter([m], [i], s=14, c=ORANGE, ec='white', lw=0.5, zorder=3)
a.axvline(1, color=MUTED, lw=0.6, ls=':'); a.set_yticks(range(len(rows)), [f'{k} ({n})' for m, k, n, lo, hi in rows], fontsize=6.5); a.set_xlabel('DYZ1:DYZ2 unit ratio, 1kGP males\n(median, 95% bootstrap interval)')
for i, t in enumerate('abc'): ax[i].text(-0.02, 1.04, t, transform=ax[i].transAxes, fontweight='bold', fontsize=9, ha='right')
fig.savefig(OUT); print('written', OUT, TAB, flush=True)
