#!/usr/bin/env python3
"""San Roman 2024 full method (Teitz 2018 GC correction + Get_DYZ1_Depth.R, their code unchanged) on our 1kGP 30x CRAMs against Yruler, same males, assembly truth.

usage: python3 scripts/sanroman_full.py <their table> <yruler cv_units.tsv> <truth tsv> <ped> <out.tsv> <cohort_1kgp.tsv>
their table: results/checks/sanroman_30x_teitz.tsv (DYZ1_norm_coverage = DYZ1 depth / single-copy depth, per male).
San Roman: one K per fold, 10 repeats of 5-fold CV grouped by family (rng seed 7 + repeat), as sanroman_compare.py. Yruler: its own cross-validation errors (cv_units.tsv).
Their DYZ1 signal is also scored against Yq12 and DYZ2 truth, as a proxy. Errors in % of the assembly length: mean, median, p90, max of the absolute error.
Paired: per-male mean absolute error (males with the smaller one, two-sided Wilcoxon signed-rank P), Spearman of the per-male mean signed errors of the two methods.
cv_truth_signal: CV of truth / signal, no scale factor (Yruler signal = released-model estimate, cohort_1kgp.tsv, one constant per class cancels).
"""
import sys, collections
import numpy as np
from scipy.stats import spearmanr, wilcoxon

SR, CVU, TRUTH, PED, OUT, COH = sys.argv[1:7]
def log(*a): print(*a, flush=True)
def rows(p):
    L = [l.rstrip('\n').split('\t') for l in open(p)]
    return [dict(zip(L[0], x)) for x in L[1:]]
sr = {r['sample']: float(r['dyz1_norm_coverage']) for r in rows(SR)}
T = collections.defaultdict(dict)
for r in rows(TRUTH): T[r['class']][r['sample']] = float(r['truth_bp'])
Y = collections.defaultdict(lambda: collections.defaultdict(list))
for r in rows(CVU): Y[r['class']][r['male']].append(float(r['error_pct']))
fam = {x.split()[1]: x.split()[0] for x in open(PED) if x.strip() and not x.startswith('FamilyID')}
E = {r['sample']: r for r in rows(COH)}; NM = {'dyz1': 'DYZ1', 'yq12': 'Yq12', 'dyz2': 'DYZ2'}
cv = lambda v: np.std(v, ddof=1) / np.mean(v) * 100
st = lambda v: (np.mean(np.abs(v)), np.median(np.abs(v)), np.percentile(np.abs(v), 90), np.max(np.abs(v)))
out = []
for c in ('dyz1', 'yq12', 'dyz2'):
    M = sorted(m for m in T[c] if m in sr and m in Y[c])
    fams = sorted({fam.get(m, m) for m in M}); P = collections.defaultdict(list)
    for rep in range(10):
        g = list(fams); np.random.default_rng(7 + rep).shuffle(g); fold = {f: i % 5 for i, f in enumerate(g)}
        for k in range(5):
            tr = [m for m in M if fold[fam.get(m, m)] != k]; K = sum(T[c][m] for m in tr) / sum(sr[m] for m in tr)
            for m in M:
                if fold[fam.get(m, m)] == k: P[m].append(100 * (K * sr[m] / T[c][m] - 1))
    a = st([e for m in M for e in P[m]]); b = st([e for m in M for e in Y[c][m]])
    rho = spearmanr([sr[m] for m in M], [T[c][m] for m in M]).correlation
    cvr = cv([T[c][m] / sr[m] for m in M]); cvy = cv([T[c][m] / float(E[m][f'{NM[c]}_bp']) for m in M])
    ma = np.array([np.mean(np.abs(P[m])) for m in M]); ya = np.array([np.mean(np.abs(Y[c][m])) for m in M]); wp = wilcoxon(ma, ya).pvalue
    es = spearmanr([np.mean(P[m]) for m in M], [np.mean(Y[c][m]) for m in M]).correlation
    log(f'{c}: {len(M)} males, {len(fams)} families | San Roman CV |error| mean {a[0]:.2f} median {a[1]:.2f} p90 {a[2]:.2f} max {a[3]:.2f}, Spearman {rho:.3f}, CV of truth/signal {cvr:.2f}% | Yruler CV |error| mean {b[0]:.2f} median {b[1]:.2f} p90 {b[2]:.2f} max {b[3]:.2f}, CV of truth/signal {cvy:.2f}%')
    log(f'  per-male mean |error| smaller: San Roman {int((ma < ya).sum())}, Yruler {int((ya < ma).sum())}, Wilcoxon P {wp:.2g}; Spearman of per-male signed errors {es:.4f}')
    out.append([c, len(M), len(fams), 'San Roman 2024 (Teitz GC), their code', *[f'{x:.2f}' for x in a], f'{rho:.4f}', f'{cvr:.2f}', int((ma < ya).sum()), f'{wp:.2g}', f'{es:.4f}'])
    out.append([c, len(M), len(fams), 'Yruler, cross-validation', *[f'{x:.2f}' for x in b], '', f'{cvy:.2f}', int((ya < ma).sum()), f'{wp:.2g}', f'{es:.4f}'])
with open(OUT, 'w') as o:
    o.write('class\tmales\tfamilies\tmethod\tabs_error_mean\tabs_error_median\tabs_error_p90\tabs_error_max\tspearman_signal_truth\tcv_truth_signal\tmales_smaller_mean_abs_error\twilcoxon_p\tspearman_signed_errors\n')
    for r in out: o.write('\t'.join(str(x) for x in r) + '\n')
log(f'wrote {OUT}')
