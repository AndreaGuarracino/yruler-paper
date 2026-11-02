#!/usr/bin/env python3
"""CV error shared across classes against marker dropout, library G+C and x dosage; class trim effect.

usage: python3 scripts/shared_error_check.py <results dir> <fzero.tsv> <trimchk.tsv> <gc_curves.tsv>
factor = mean log CV error of DYZ1, DYZ17, DYZ18S per man; dropout estimator: class mean over present markers / ruler mean over present markers, applied to CV errors without a new fit
"""
import sys, csv, math, collections as co
import numpy as np
from scipy.stats import spearmanr
D, FZ, TR, GC = sys.argv[1:5]; rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
E = co.defaultdict(lambda: co.defaultdict(list))
for r in rd(f'{D}/cv_units.tsv'): E[r['class']][r['male']].append(math.log1p(float(r['error_pct']) / 100))
M = {c: {s: float(np.mean(v)) for s, v in d.items()} for c, d in E.items()}; men = sorted(set(M['dyz1']) & set(M['dyz17']) & set(M['dyz18s'])); F = np.array([np.mean([M[c][s] for c in ('dyz1', 'dyz17', 'dyz18s')]) for s in men])
print(f'men {len(men)}; SD of the shared factor {100 * F.std():.2f}%')
Z = {r['sample']: r for r in rd(FZ)}; X = {r['sample']: r for r in rd(f'{D}/cohort_1kgp.tsv')}; G = {r['sample']: r for r in rd(GC)}; I = {r['sample']: r for r in rd('results/model_checks_2026-10-29/insert_1kgp.tsv')}
men2 = [i for i, s in enumerate(men) if s in Z and s in G and s in X and s in I]; f = F[men2]; S = [men[i] for i in men2]; print('men with all covariates', len(S))
V = {'ruler frac_zero': [float(Z[s]['ruler_fz']) for s in S], 'dyz1 frac_zero': [float(Z[s]['dyz1_fz']) for s in S], 'dyz17 frac_zero': [float(Z[s]['dyz17_fz']) for s in S], 'dyz18s frac_zero': [float(Z[s]['dyz18s_fz']) for s in S], 'low G+C curve (gc8)': [float(G[s]['curve_gc8']) for s in S], 'x_dose': [float(X[s]['x_dose']) for s in S], 'overlap share': [float(I[s]['overlap_share']) for s in S], 'ruler mean count': [float(Z[s]['ruler_mean']) for s in S]}
for k, v in V.items(): r = spearmanr(f, v); print(f'{k}: Spearman {r[0]:+.3f}, P {r[1]:.2g}')
b = np.array([x['batch'] == '2504' for x in (X[s] for s in S)]); print(f'batch 2504: n {int(b.sum())}, factor median {100 * np.median(f[b]):+.2f}%; batch 698: n {int((~b).sum())}, {100 * np.median(f[~b]):+.2f}%')
st = lambda a: f'{np.mean(a):.2f}/{np.median(a):.2f}/{np.percentile(a, 90):.2f}/{np.max(a):.2f}'
for c in ('dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17'):
    R = [(r['male'], float(r['error_pct'])) for r in rd(f'{D}/cv_units.tsv') if r['class'] == c and r['male'] in Z]; s = {m: (1 - float(Z[m]['ruler_fz'])) / (1 - float(Z[m][f'{c}_fz'])) for m, _ in R}; k = np.median(list(s.values()))
    a = np.abs([e for _, e in R]); n = np.abs([100 * ((1 + e / 100) * s[m] / k - 1) for m, e in R]); print(f'present-only {c}: mean/median/p90/max {st(a)} -> {st(n)}')
T = {r['sample']: r for r in rd(TR)}
for c in ('dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17'):
    R = [(r['male'], float(r['error_pct'])) for r in rd(f'{D}/cv_units.tsv') if r['class'] == c and r['male'] in T]; s = {m: 1 - float(T[m][f'{c}_share']) for m, _ in R}; k = np.median(list(s.values())); v = np.array(list(s.values()))
    a = np.abs([e for _, e in R]); n = np.abs([100 * ((1 + e / 100) * s[m] / k - 1) for m, e in R]); print(f'trim {c}: share of the sum above the cap median {100 * (1 - np.median(v)):.3f}%, max {100 * (1 - v.min()):.3f}%; {st(a)} -> {st(n)}')
