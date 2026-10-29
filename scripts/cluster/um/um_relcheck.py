# Release check: yruler output vs Python model
# usage: python3 um_relcheck.py <sample> <yruler out dir> [model dir]; needs markers_u25/cnpy/<sample>.npy
import sys, numpy as np
P = '/path/to/project'; O = f'{P}/markers_u25'; sys.path.insert(0, f'{P}/scripts'); import um_model
s, YD = sys.argv[1:3]; MD = sys.argv[3] if len(sys.argv) > 3 else f'{P}/model'
U = np.load(f'{O}/u_keys.npy'); fm, M = um_model.load(MD); c = np.load(f'{O}/cnpy/{s}.npy').astype(np.float64)
L, S, Y = um_model.lengths(c, fm, M, um_model.gcount(U), fm[f'ysc_{M[0]["class"]}'])
rows = lambda f: [l.rstrip('\n').split('\t') for l in open(f'{YD}/{f}') if not l.startswith(('#', 'sample'))]
yl = {x[1]: float(x[2]) for x in rows('lengths.tsv')}; ys = {x[2]: (float(x[7]), float(x[4])) for x in rows('summary.tsv')}
bad = int(set(yl) != set(L)); print(f'{s}: reported classes yruler {sorted(yl)} python {sorted(L)}', flush=True)
for k in S: ok = abs(ys[k][0] - S[k]) <= 0.05 + 1e-6; bad += not ok; print(f'{s} {k}: sum_gc yruler {ys[k][0]:.1f} python {S[k]:.3f} {"OK" if ok else "DIFF"}', flush=True)   # print rounding
ok = ys['ysinglecopy'][1] == Y; bad += not ok; print(f'{s} ruler sum: yruler {ys["ysinglecopy"][1]:.0f} python {Y:.0f} {"OK" if ok else "DIFF"}', flush=True)
for k in L: ok = abs(yl.get(k, np.nan) - L[k]) <= 0.5 + 1e-6; bad += not ok; print(f'{s} {k}: length yruler {yl.get(k, np.nan):.0f} python {L[k]:.3f} bp {"OK" if ok else "DIFF"}', flush=True)
print('PASS' if bad == 0 else f'FAIL {bad}'); sys.exit(1 if bad else 0)
