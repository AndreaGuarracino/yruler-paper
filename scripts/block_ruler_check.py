#!/usr/bin/env python3
"""Ruler blocks in normal men: length shift of each man and CV errors with a median-of-blocks ruler.

usage: python3 scripts/block_ruler_check.py <results dir> <ruler_ypq label counts, 8 comma-separated> <out.tsv>
shift = own median block dose / key-weighted mean block dose - 1 = length bias of the summed ruler; CV estimates are divided by 1 + shift, without a new fit
"""
import sys, csv, collections as co
import numpy as np
D, W, OUT = sys.argv[1:4]; rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
W = np.array(W.split(','), float); W /= W.sum(); B = [f'xdr{i}_dose' for i in range(1, 9)]
X = {r['sample']: r for r in rd(f'{D}/cohort_1kgp.tsv')}; P = {r['sample']: r for r in rd(f'{D}/yprofile_1kgp.tsv')}
def shift(s): v = np.array([float(P[s][b]) for b in B]); return np.median(v) / (W * v).sum()
A = {s: shift(s) for s in P if X[s]['sex'] == 'male' and X[s]['y_state'] == 'normal'}; v = 100 * (np.array(list(A.values())) - 1)
print(f'normal men {len(v)}; |shift| median {np.median(np.abs(v)):.3f}%, p99 {np.percentile(np.abs(v), 99):.3f}%, above 2% {int((np.abs(v) > 2).sum())}, min {v.min():.2f}%, max {v.max():.2f}%', flush=True)
E = co.defaultdict(list)
for r in rd(f'{D}/cv_units.tsv'):
    if r['class'] != 'yq12rest': E[r['class']].append((r['male'], float(r['error_pct'])))
st = lambda a: (np.mean(a), np.median(a), np.percentile(a, 90), np.max(a))
with open(OUT, 'w') as o:
    o.write('class\tpredictions\tmen_with_profile\tmean_old\tmean_new\tmedian_old\tmedian_new\tp90_old\tp90_new\tmax_old\tmax_new\n')
    for c, R in E.items():
        R = [(m, e) for m, e in R if m in P]; s = {m: shift(m) for m, _ in R}; k = np.median(list(s.values())); a = np.abs([e for _, e in R]); b = np.abs([100 * ((1 + e / 100) / s[m] * k - 1) for m, e in R]); q = [x for p in zip(st(a), st(b)) for x in p]
        o.write(f'{c}\t{len(R)}\t{len(s)}\t' + '\t'.join(f'{x:.2f}' for x in q) + '\n'); print(c, ' '.join(f'{x:.2f}' for x in q), flush=True)
J = [m for m in {m for m, _ in E['yq12']} if X[m]['lineage'] == 'J']
for m in sorted(J): e = [e for q, e in E['yq12'] if q == m]; print(m, X[m]['haplogroup'][:9], f'Yq12 median error {np.median(e):+.2f}%, shift {100 * (shift(m) - 1):+.2f}%, after {np.median(e) - 100 * (shift(m) - 1):+.2f}%')
