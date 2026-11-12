# Release check: yruler output vs Python model
# usage: [UM=markers_v6 COHORT=<um_apply table>] python3 um_relcheck.py <sample> <yruler out dir> [model dir]; counts from <out dir>/counts.tsv, else markers_u25/cnpy/<sample>.npy; COHORT: y_dose of the pipeline (slice counts) vs yruler (|diff| <= 0.005)
import os, sys, numpy as np
P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_u25")}'; sys.path.insert(0, f'{P}/scripts'); import um_model
s, YD = sys.argv[1:3]; MD = sys.argv[3] if len(sys.argv) > 3 else f'{P}/model'
U = np.load(f'{O}/u_keys.npy'); fm, M = um_model.load(MD)
if os.path.exists(f'{YD}/counts.tsv'):   # yruler per-marker counts of this run (thinned runs have no cnpy)
    km = np.loadtxt(f'{YD}/counts.tsv', dtype='U31', usecols=0); n = np.loadtxt(f'{YD}/counts.tsv', dtype=np.float64, usecols=1)
    LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]; v = LUT[np.frombuffer(''.join(km).encode(), np.uint8).reshape(-1, 31)]
    k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
    for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
    kc = np.minimum(k, rc); ix = np.searchsorted(U, kc); ok = (ix < len(U)) & (U[np.minimum(ix, len(U) - 1)] == kc)
    c = np.zeros(len(U)); c[ix[ok]] = n[ok]; print(f'{s}: counts from {YD}/counts.tsv, {len(n)} markers, sum {n.sum():.0f}; not in u_keys (chr2depth) {(~ok).sum()}', flush=True)
else: c = np.load(f'{O}/cnpy/{s}.npy').astype(np.float64)
e = um_model.lengths(c, fm, M, um_model.gcount(U), fm[f'ysc_{M[0]["class"]}'])
if e is None: sys.exit(f'{s}: no G+C curve in Python')
L, S, Y = e
rows = lambda f: [l.rstrip('\n').split('\t') for l in open(f'{YD}/{f}') if not l.startswith(('#', 'sample'))]
LF = 'estimates.tsv' if os.path.exists(f'{YD}/estimates.tsv') else 'lengths.tsv'; LH = open(f'{YD}/{LF}').readline().rstrip('\n').split('\t'); LT = rows(LF); col = lambda x, n: x[LH.index(n)] if n in LH else '.'   # estimates.tsv, else lengths.tsv, by header (relative column since 2026-11-01)
yl = {x[1]: float(x[2]) for x in LT}; bs = {{'ruler': 'y_found', 'haploid_depth': 'no_y'}.get(v, v) for v in (col(x, 'y_status' if 'y_status' in LH else 'estimate_basis') for x in LT)}; ys = {x[2]: (float(x[7]), float(x[4])) for x in rows('summary.tsv')}   # y_status, old values mapped to new
bad = int(set(yl) != set(L)); print(f'{s}: reported classes yruler {sorted(yl)} python {sorted(L)}', flush=True)
for k in S: ok = abs(ys[k][0] - S[k]) <= 0.05 + 1e-6; bad += not ok; print(f'{s} {k}: sum_gc yruler {ys[k][0]:.1f} python {S[k]:.3f} {"OK" if ok else "DIFF"}', flush=True)   # print rounding
ok = ys['ysinglecopy'][1] == Y; bad += not ok; print(f'{s} ruler sum: yruler {ys["ysinglecopy"][1]:.0f} python {Y:.0f} {"OK" if ok else "DIFF"}', flush=True)
for k in L: ok = abs(yl.get(k, np.nan) - L[k]) <= 0.5 + 1e-6; bad += not ok; print(f'{s} {k}: length yruler {yl.get(k, np.nan):.0f} python {L[k]:.3f} bp {"OK" if ok else "DIFF"}', flush=True)
if 'relative' in LH:   # relative = sum_gc(class) / ruler sum, the row class itself
    for x in LT:
        if col(x, 'relative') == '.': print(f'{s} {x[1]}: composed class, no relative value (terms {col(x, "terms") if "terms" in LH else "?"})', flush=True); continue
        k = x[1]; r = ys[k][0] / ys['ysinglecopy'][1] if k in ys else np.nan; ok = abs(float(col(x, 'relative')) / r - 1) <= 1e-4; bad += not ok; print(f'{s} {k}: relative yruler {col(x, "relative")} from summary {r:.6g} {"OK" if ok else "DIFF"}', flush=True)
ok = bs == {'y_found'}; bad += not ok; print(f'{s} y_status yruler {sorted(bs)} ({LF}) python y_found {"OK" if ok else "DIFF"}', flush=True)
if os.environ.get('COHORT'):   # y_dose: pipeline (window keys in region slices, YDN scale) vs yruler (whole input, model y_dose_scale)
    R = [l.rstrip('\n').split('\t') for l in open(os.environ['COHORT'])]; py = {r[0]: r[R[0].index('y_dose')] for r in R[1:]}.get(s, 'nan'); yd = float(col(LT[0], 'y_dose'))
    ok = abs(float(py) - yd) <= 0.005; bad += not ok; print(f'{s} y_dose yruler {yd:.4f} pipeline {py} diff {yd - float(py):+.4f} {"OK" if ok else "DIFF"}', flush=True)
print('PASS' if bad == 0 else f'FAIL {bad}'); sys.exit(1 if bad else 0)
