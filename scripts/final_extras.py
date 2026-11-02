#!/usr/bin/env python3
"""Extra tables of a final model: CV median bootstrap, prediction intervals, lineage medians, DYZ1:DYZ2 truth vs yruler.

usage: python3 scripts/final_extras.py <results dir> <indep_setB.tsv> <truth_regions.tsv>
reads cv_units.tsv, holdout_21_units.tsv, cohort_1kgp.tsv in the results dir; interval half-width w: |error| as % of the estimate (truth within estimate x (1 +- w)); writes cv_median_bootstrap.tsv, prediction_intervals.tsv,
lineage_medians_bootstrap.tsv, dyz1_dyz2_truth_vs_yruler.tsv there. Bootstrap: 2,000 resamples, seed 20261101.
DYZ1:DYZ2 unit ratio in Hallast 2026 units (bp per hit): DYZ1 3,478.9, DYZ2 2,404.4.
"""
import sys, csv, collections
import numpy as np
from scipy.stats import spearmanr

D, IND, REG = sys.argv[1:4]; NB, SEED, U1, U2 = 2000, 20261101, 3478.9, 2404.4
CL = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17']
def log(*a): print(*a, flush=True)
def rows(f): return list(csv.DictReader(open(f), delimiter='\t'))
def put(f, H, R):
    with open(f'{D}/{f}', 'w') as o: o.write('\t'.join(H) + '\n'); o.writelines('\t'.join(map(str, r)) + '\n' for r in R)
    log('written', f'{D}/{f}', len(R), 'rows')
def boot(v, f, rng): v = np.asarray(v); b = [f(v[rng.integers(0, len(v), len(v))]) for _ in range(NB)]; return np.percentile(b, 2.5), np.percentile(b, 97.5)
E = collections.defaultdict(dict)
for r in rows(f'{D}/cv_units.tsv'): E[r['class']][(r['repeat'], r['male'])] = float(r['error_pct'])
H21 = collections.defaultdict(list)
for r in rows(f'{D}/holdout_21_units.tsv'): H21[r['class']].append(float(r['error_pct']))
pm = lambda d: {m: np.mean([abs(v) for (r, mm), v in d.items() if mm == m]) for m in sorted({k[1] for k in d})}   # per-male mean |error|
rng = np.random.default_rng(SEED); B, PI = [], []
for c in CL:
    v = list(pm(E[c]).values()); lo, hi = boot(v, np.median, rng); B.append((c, len(v), f'{np.median(v):.2f}', f'{lo:.2f}', f'{hi:.2f}'))
    ae = lambda v: np.abs(v) / (1 + np.asarray(v) / 100)   # |error| as % of the estimate
    a = ae(list(E[c].values())); w = np.percentile(a, [80, 90, 95]); h = ae(H21[c])
    PI.append((c, len(a), *(f'{x:.2f}' for x in w), len(h), int((h <= w[1]).sum())))
    log(f'{c}: men {len(v)}, median per-male mean |error| {np.median(v):.2f}% (95% {lo:.2f}-{hi:.2f}); 90% half-width {w[1]:.2f}%, holdout_21 within {int((h <= w[1]).sum())} of {len(h)}')
put('cv_median_bootstrap.tsv', ['class', 'n_males', 'median_of_male_mean_abs_pct', 'ci95_low', 'ci95_high'], B)
put('prediction_intervals.tsv', ['class', 'cv_predictions', 'half_width_80_pct', 'half_width_90_pct', 'half_width_95_pct', 'holdout21_n', 'holdout21_within_90'], PI)
X = rows(f'{D}/cohort_1kgp.tsv'); lin = {x['sample']: x['lineage'] for x in X}
inc = [x for x in X if x['estimate_basis'] == 'ruler' and x['sex'] != 'female' and x['y_state'] in ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain', 'par2_gain')]
G = collections.defaultdict(list)
for x in inc: G[x['lineage']].append(x)
KIN = {r['sample']: r['kinship_component_id'] for r in rows('results/estimates_1kgp.relatedness_samples.tsv')}
LN = sorted(k for k in G if len(G[k]) >= 10 and k not in ('.', '?', '')); L = []
log(f'cohort men included {len(inc)}; lineages with >= 10 men {len(LN)} ({sum(len(G[k]) for k in LN)} men)')
for k in LN:
    for col in ('DYZ3', 'Yq12', 'DYZ1', 'DYZ2', 'DYZ19', 'DYZ17'):
        w = [x for x in G[k] if col in ('DYZ3', 'DYZ19', 'DYZ17') or x['y_state'] != 'par2_gain']; v = [float(x[f'{col}_bp']) for x in w]; cp = collections.defaultdict(list)
        for x in w: cp[KIN.get(x['sample'], x['sample'])].append(float(x[f'{col}_bp']))   # resample kinship components
        cl = [np.array(c) for c in cp.values()]; bs = [np.median(np.concatenate([cl[i] for i in rng.integers(0, len(cl), len(cl))])) for _ in range(NB)]; lo, hi = np.percentile(bs, 2.5), np.percentile(bs, 97.5); L.append((k, len(v), col, f'{np.median(v):.0f}', f'{lo:.0f}', f'{hi:.0f}'))
put('lineage_medians_bootstrap.tsv', ['lineage', 'n', 'class', 'median_bp', 'ci95_low', 'ci95_high'], L)
A12 = {x[0]: x[3] for x in (l.rstrip('\n').split('\t') for l in open(REG)) if x[1] == 'yq12'}; T = {r['assembly']: r for r in rows(IND)}
M = sorted({m for (r, m) in E['dyz1']} & {m for (r, m) in E['dyz2']}); miss = [m for m in M if A12.get(m) not in T]
if miss: sys.exit(f'ERROR no truth row for {miss}')
tr = {m: (float(T[A12[m]]['dyz1_units']) / U1) / (float(T[A12[m]]['dyz2_units']) / U2) for m in M}
cv = {m: np.mean([tr[m] * (1 + E['dyz1'][k] / 100) / (1 + E['dyz2'][k] / 100) for k in E['dyz1'] if k[1] == m and k in E['dyz2']]) for m in M}   # mean over repeats
dv = [100 * (cv[m] / tr[m] - 1) for m in M]; q = lambda v: [len(v), f'{np.mean(v):.3f}', f'{np.median(v):.3f}', f'{np.percentile(v, 10):.3f}', f'{np.percentile(v, 90):.3f}', f'{np.min(v):.3f}', f'{np.max(v):.3f}']
coh = [(float(x['DYZ1_bp']) / U1) / (float(x['DYZ2_bp']) / U2) for x in inc if x['y_state'] != 'par2_gain']; sd = lambda v: np.std(np.log(v), ddof=1)
R = [['truth_hallast_unit_ratio'] + q(list(tr.values())), ['yruler_cv_prediction_unit_ratio'] + q(list(cv.values())), ['yruler_over_truth_minus_1_pct'] + q(dv),
     ['spearman_yruler_cv_vs_truth', len(M), f'{spearmanr([cv[m] for m in M], [tr[m] for m in M])[0]:.4f}', '', '', '', '', ''],
     ['sd_log_ratio_cv_error', len(M), f'{sd([cv[m] / tr[m] for m in M]):.4f}', '', '', '', '', ''], ['sd_log_ratio_1kgp_included', len(coh), f'{sd(coh):.4f}', '', '', '', '', ''],
     ['sd_log_ratio_truth_assemblies', len(M), f'{sd(list(tr.values())):.4f}', '', '', '', '', ''], ['error_variance_share_of_1kgp', len(M), f'{sd([cv[m] / tr[m] for m in M]) ** 2 / sd(coh) ** 2:.4f}', '', '', '', '', '']]
for k in sorted({lin.get(m, '?') for m in M}):
    g = [m for m in M if lin.get(m, '?') == k]
    if len(g) >= 10 and k != '?': R += [[f'lineage_{k}_truth_unit_ratio', len(g), '', f'{np.median([tr[m] for m in g]):.3f}', '', '', '', ''], [f'lineage_{k}_yruler_cv_unit_ratio', len(g), '', f'{np.median([cv[m] for m in g]):.3f}', '', '', '', '']]
for r in R: log('  ', *r)
put('dyz1_dyz2_truth_vs_yruler.tsv', ['measure', 'men', 'mean', 'median', 'p10', 'p90', 'min', 'max'], R)
