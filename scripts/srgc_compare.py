#!/usr/bin/env python3
"""San Roman depth with Teitz GC correction against yruler, same males and CV folds.

usage: python3 scripts/srgc_compare.py <srgc out dir> <depth_f1.tsv> <out.tsv>
depth measure of a class = mean GC-weighted coverage over covered bases of its blocks / that of chrY:14.5-15.5 Mb;
blocks on chrY only (as San Roman) or on every contig; Yq12 also from the San Roman DYZ1 region and as K-scaled DYZ1 + DYZ2 + rest;
K per fold = sum truth / sum measure on the training males of the yruler fold; errors of yruler from the same folds and males
"""
import sys, csv, glob, collections as co, numpy as np
D, F1, OUT = sys.argv[1:4]; rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
M = co.defaultdict(dict)   # sample -> measure name -> value
for f in sorted(glob.glob(f'{D}/*.tsv')):
    R = rd(f); s = R[0]['sample']; ctl = next(float(r['mean_covered']) for r in R if r['region'] == 'control'); A = co.defaultdict(lambda: [0.0, 0])
    for r in R:
        if r['region'] == 'control': continue
        for scope in ('chrY', 'all'):
            if scope == 'chrY' and r['chrom'] != 'chrY': continue
            k = (r['region'], scope); A[k][0] += float(r['mean_covered']) * int(r['covered_bases']); A[k][1] += int(r['covered_bases'])
    for (reg, scope), (w, n) in A.items(): M[s][f'{reg}|{scope}'] = w / n / ctl if n else np.nan
print(len(M), 'males with depth', flush=True)
T = co.defaultdict(dict); E = co.defaultdict(dict); FO = {}
for r in rd(F1):
    T[r['class']][r['sample']] = float(r['reference_bp']); E[r['class']][(r['repeat'], r['sample'])] = abs(float(r['error_pct'])); FO[r['class'], r['repeat'], r['sample']] = r['fold']
T['yq12rest'] = {s: T['yq12'][s] - T['dyz1'][s] - T['dyz2'][s] for s in T['yq12'] if s in T['dyz1'] and s in T['dyz2']}
st = lambda a: [len(a), *(f'{v:.2f}' for v in (np.mean(a), np.median(a), np.percentile(a, 90), np.max(a)))] if len(a) else [0, '', '', '', '']
def cv(c, comps):   # comps: [(truth class, measure name)], estimate = sum of K x measure
    men = sorted(s for s in T[c] if s in M and all(s in T[t] and np.isfinite(M[s].get(m, np.nan)) and M[s].get(m, 0) > 0 for t, m in comps))
    reps = sorted({k[0] for k in E[c]}); err, yerr = [], []
    for rp in reps:
        ms = [s for s in men if (c, rp, s) in FO]
        for f in sorted({FO[c, rp, s] for s in ms}):
            tr = [s for s in ms if FO[c, rp, s] != f]
            K = [sum(T[t][s] for s in tr) / sum(M[s][m] for s in tr) for t, m in comps]
            for s in (s for s in ms if FO[c, rp, s] == f):
                est = sum(k * M[s][m] for k, (t, m) in zip(K, comps)); err.append(100 * abs(est / T[c][s] - 1)); yerr.append(E[c][rp, s])
    return men, err, yerr
rows = []
for c, scope in [(c, sc) for c in ('dyz3', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17', 'yq12') for sc in ('chrY', 'all')]:
    variants = [(f'{c} blocks', [(c, f'{c}|{scope}')])]
    if c == 'yq12': variants += [('San Roman DYZ1 region', [('yq12', 'sanroman_dyz1|chrY')])] if scope == 'chrY' else []; variants += [('DYZ1 + DYZ2 + rest blocks', [('dyz1', f'dyz1|{scope}'), ('dyz2', f'dyz2|{scope}'), ('yq12rest', f'yq12rest|{scope}')])]
    if c == 'dyz1' and scope == 'chrY': variants += [('San Roman DYZ1 region', [('dyz1', 'sanroman_dyz1|chrY')])]
    for name, comps in variants:
        men, err, yerr = cv(c, comps)
        rows.append([c, scope, name, len(men), *st(err)[1:], *st(yerr)[1:]]); print(*rows[-1], flush=True)
with open(OUT, 'w') as o:
    o.write('class\tblocks\tdepth_measure\tmales\tdepth_mean\tdepth_median\tdepth_p90\tdepth_max\tyruler_mean\tyruler_median\tyruler_p90\tyruler_max\n')
    for r in rows: o.write('\t'.join(str(v) for v in r) + '\n')
print('written', OUT)
