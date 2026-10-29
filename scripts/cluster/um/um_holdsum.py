# Holdout of the model: per-class source, Yq12 composite
# usage: python3 um_holdsum.py <base tsv> <floor tsv> <floor classes,comma> <out tsv>
import sys, numpy as np
P = '/path/to/project'; O = f'{P}/markers_u25'; X = f'{P}/truth_qc/x'; sys.path.insert(0, f'{P}/scripts'); import truth_cl
B, F, FC, OUT = sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4]
T = truth_cl.load(O)
for l in open(f'{X}/rule_w4.txt'): c, s = l.split()[:2]; T[c].pop(s, None)
for l in open(f'{X}/truth_all.tsv'):
    c, s, v = l.split()[:3]
    if s in T[c]: T[c][s] = float(v)
rd = lambda f: {(c, s): float(e) for c, s, e in (l.split() for l in open(f))}
b, f = rd(B), rd(F); E = {k: (f if k[0] in FC else b)[k] for k in set(b) | set(f) if k in (f if k[0] in FC else b)}
for (c, s) in [k for k in E if k[0] == 'yq12']:
    if ('dyz2', s) in b and ('dyz2', s) in f: E[c, s] = 100 * ((T['yq12'][s] * (1 + b[c, s] / 100) - T['dyz2'][s] * (1 + b['dyz2', s] / 100) + T['dyz2'][s] * (1 + f['dyz2', s] / 100)) / T['yq12'][s] - 1)
    else: print(f'WARNING yq12 {s}: no DYZ2 holdout error, left out', flush=True); del E[c, s]
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'yq2k7': 'YqHet2.7', 'dyz18sub': 'DYZ18sub', 'dyz19': 'DYZ19'}
with open(OUT, 'w') as o:
    o.write('class\tmales\tmean_abs_pct\tmedian_abs_pct\tp90_abs_pct\tmax_abs_pct\tsigned_median_pct\tworst_male\n')
    for c in NM:
        v = {s: e for (k, s), e in E.items() if k == c}; a = np.abs(list(v.values())); w = max(v, key=lambda s: abs(v[s]))
        o.write(f'{NM[c]}\t{len(a)}\t{a.mean():.2f}\t{np.median(a):.2f}\t{np.percentile(a, 90):.2f}\t{a.max():.2f}\t{np.median(list(v.values())):+.2f}\t{w}\n')
print(open(OUT).read(), flush=True)
