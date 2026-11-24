# Release: write = markers and model.tsv for yruler (was um_release.py); check = yruler output vs the Python model (was um_relcheck.py)
# usage: python3 release.py write <model dir> | check <sample> <yruler out dir> <model dir> | ontscale <model dir> <male list>
import sys
MODE = sys.argv.pop(1) if len(sys.argv) > 1 else ''
if MODE == 'write':
    # Release markers and model.tsv for yruler
    # usage: UM=markers_v1 DENOMK=<window keys npy>[,<npy>...] YDN=<reference raw y_dose> NOY=<no-Y cut, scaled y_dose>] python3 release.py write [model dir]; writes <model dir>/release/markers.fa and model.tsv; denominator classes depth1..depthN = the window keys (default: chr2depth = chr2 50-60 markers); YDN -> #y_dose_scale in model.tsv
    import os, sys, numpy as np
    P = os.environ['P']; O = f'{P}/{os.environ["UM"]}'; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import model
    MD = sys.argv[1]; R = f'{MD}/release'; os.makedirs(R, exist_ok=True)
    def log(*a): print(*a, flush=True)
    U = np.load(f'{O}/u_keys.npy'); fm, M = model.load(MD)
    if any('_' in r['class'] for r in M): raise SystemExit('ERROR class names with _ (yruler reads the header name before _)')
    LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]
    def enc(lines):   # canonical keys
        v = LUT[np.frombuffer(''.join(lines).encode(), np.uint8).reshape(-1, 31)]; k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
        for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
        for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
        return np.minimum(k, rc)
    def dec(keys):   # uint64 keys -> 31-mers
        b = np.zeros((len(keys), 31), np.uint8)
        for i in range(31): b[:, 30 - i] = ((keys >> np.uint64(2 * i)) & np.uint64(3)).astype(np.uint8)
        return np.frombuffer(b'ACGT', np.uint8)[b]
    if not os.environ.get('DENOMK'): sys.exit('ERROR DENOMK (chr8 window keys) is required')
    DW = [(f'depth{i + 1}', np.unique(np.load(f).astype(np.uint64))) for i, f in enumerate(os.environ['DENOMK'].split(','))]
    y0 = fm[f'ysc_{M[0]["class"]}']
    for r in M:
        if not np.array_equal(fm[f'ysc_{r["class"]}'], y0): raise SystemExit(f'ERROR ruler of {r["class"]} differs')
        if len(fm[f'm_{r["class"]}']) != int(r['markers']) or len(np.unique(fm[f'm_{r["class"]}'])) != int(r['markers']): raise SystemExit(f'ERROR marker count {r["class"]}')
    S = [(r['class'], U[np.sort(fm[f'm_{r["class"]}'])]) for r in M] + [('ysinglecopy', U[np.sort(y0)])] + DW   # sorted canonical keys
    with open(f'{R}/markers.fa', 'wb') as o:
        for c, ix in S:
            s = dec(ix); o.write(f'>{c}\n'.encode()); o.write(np.hstack([s, np.full((len(s), 1), 10, np.uint8)]).tobytes())
            if not np.array_equal(np.sort(enc([s.tobytes().decode()])), ix): raise SystemExit(f'ERROR round trip {c}')   # all k-mers
            log(f'{c}: {len(ix)} markers')
    for (a, x), (b, y) in [(p, q) for i, p in enumerate(S) for q in S[i + 1:]]:
        n = len(np.intersect1d(x, y))
        if n: log(f'shared {a} {b}: {n}')
    with open(f'{R}/model.tsv', 'w') as o:
        o.write('# length_bp = K * sum_gc(class) / sum(ysinglecopy); K per technology (-x); length: reported length as a signed sum of class lengths, none = not reported\n' + (f'#y_dose_scale={float(os.environ["YDN"]):.6f}\n' if os.environ.get('YDN') else '') + (f'#no_y_dose={float(os.environ["NOY"]):.6f}\n' if os.environ.get('NOY') else '') + 'class\tK\ttech\tlength\n')
        for r in M: o.write(f'{r["class"]}\t{r["K"]}\tillumina\t{r["length"]}\n')
    log('written', f'{R}/markers.fa', f'{R}/model.tsv')
elif MODE == 'check':
    # Release check: yruler output vs Python model
    # usage: UM=markers_v1 [COHORT=<cohort.py table>] python3 release.py check <sample> <yruler out dir> [model dir]; counts from <out dir>/counts.tsv, else <UM>/counts/<sample>.npy; COHORT: y_dose of the pipeline (slice counts) vs yruler (|diff| <= 0.005)
    import os, sys, numpy as np
    P = os.environ['P']; O = f'{P}/{os.environ["UM"]}'; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import model
    s, YD = sys.argv[1:3]; MD = sys.argv[3]
    U = np.load(f'{O}/u_keys.npy'); fm, M = model.load(MD)
    if os.path.exists(f'{YD}/counts.tsv'):   # yruler per-marker counts of this run (thinned runs have no cnpy)
        km = np.loadtxt(f'{YD}/counts.tsv', dtype='U31', usecols=0); n = np.loadtxt(f'{YD}/counts.tsv', dtype=np.float64, usecols=1)
        LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]; v = LUT[np.frombuffer(''.join(km).encode(), np.uint8).reshape(-1, 31)]
        k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
        for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
        for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
        kc = np.minimum(k, rc); ix = np.searchsorted(U, kc); ok = (ix < len(U)) & (U[np.minimum(ix, len(U) - 1)] == kc)
        c = np.zeros(len(U)); c[ix[ok]] = n[ok]; print(f'{s}: counts from {YD}/counts.tsv, {len(n)} markers, sum {n.sum():.0f}; not in u_keys (chr2depth) {(~ok).sum()}', flush=True)
    else: c = np.load(f'{os.environ["IN"]}/counts/males/{s}.npy').astype(np.float64)
    e = model.lengths(c, fm, M, model.gcount(U), fm[f'ysc_{M[0]["class"]}'])
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
elif MODE == 'ontscale':
    # Nanopore y_dose_scale: median raw y_dose of the listed R10 males above half their median; raw = 2 x ONT ruler mean (counts_nanopore) / chr8 window depth_mean of the same reads; check: counts_nanopore equals the ont_countsB counts
    # usage: python3 release.py ontscale <model dir> <male list> <depth dir: <male>/summary.tsv with depth_mean of the chr8 window> (writes <model dir>/ydn.txt, ydn.tsv)
    import os, numpy as np
    P = os.environ['P']; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import model
    M, L, DD = sys.argv[1:4]; fm, MM = model.load(M); y0 = fm[f'ysc_{MM[0]["class"]}']; fb, MB = model.load(f'{os.environ["IN"]}/counts/nanopore_release_b_model'); yb = fb[f'ysc_{MB[0]["class"]}']
    print(f'depth dir {DD}; ONT ruler {len(y0)} keys, B ruler (ont_countsB) {len(yb)} keys, shared {len(np.intersect1d(y0, yb))}', flush=True); rows = []
    for l in open(L):
        s = l.split('\t')[0]; c = np.load(f'{os.environ["IN"]}/counts/nanopore/{s}.npy').astype(np.float64)
        H = dict(x.split('=', 1) for x in open(f'{os.environ["IN"]}/counts/nanopore_release_b/{s}/summary.tsv').readline().lstrip('#').rstrip('\n').split('\t')); dB = float(H['depth_mean'])
        d = float(dict(x.split('=', 1) for x in open(f'{DD}/{s}/summary.tsv').readline().lstrip('#').rstrip('\n').split('\t'))['depth_mean'])
        rows.append((s, 2 * c[y0].sum() / len(y0) / d, 2 * c[yb].sum() / len(yb) / dB, float(H['y_dose']), d))
    df = max(abs(x[2] - x[3]) for x in rows); print(f'{len(rows)} males; B-ruler raw y_dose from counts_nanopore vs ont_countsB summaries: max |diff| {df:.2e}', flush=True)
    if df > 1e-6: sys.exit('ERROR counts_nanopore and ont_countsB counts differ')
    r = np.array([x[1] for x in rows]); k = r[r > 0.5 * np.median(r)]; ydn = float(np.median(k))
    print(f'raw y_dose, ONT ruler: median {np.median(r):.6f}, range {r.min():.4f}-{r.max():.4f}, {len(k)} above half the median -> y_dose_scale {ydn:.6f}; B ruler median {np.median([x[2] for x in rows]):.6f}', flush=True)
    open(f'{M}/ydn.txt', 'w').write(f'{ydn:.6f}\n')
    with open(f'{M}/ydn.tsv', 'w') as o: o.write('sample\traw_y_dose_ont_ruler\traw_y_dose_b_ruler\tont_countsB_y_dose\tdepth_mean\n'); o.writelines('\t'.join([x[0]] + [f'{v:.6f}' for v in x[1:]]) + '\n' for x in rows)
else: sys.exit(f'ERROR mode {MODE}: write, check or ontscale')
