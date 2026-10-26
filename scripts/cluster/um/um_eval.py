# Uniform-rule markers: grouped repeated 5-fold CV, fixed held-out, production set
import os, sys, math, glob, subprocess as sp, numpy as np
P = '/path/to/project'; O = f'{P}/markers_u'; Z = '/path/to/conda/zstd/1.5.7/bin/zstd'
THR = (0.5, 0.9); REP, FOLD = int(os.environ.get('REP', 10)), 5; FULL = 0.4   # superset presence
CL = ['dyz3', 'yq12', 'ysinglecopy', 'dyz1', 'dyz2', 'dyz18s', 'yq2k7', 'dyz19']
REG = {'dyz18s': ('dyz18', 'yq3k1')}   # merged class: union of regions
def log(*a): print(*a, flush=True)
# superset keys, G+C, chr2depth markers
U = np.load(f'{O}/u_keys.npy'); NU = len(U)
gc = np.zeros(NU, np.uint8)
for i in range(31): b = (U >> np.uint64(2 * i)) & np.uint64(3); gc += ((b == 1) | (b == 2)).astype(np.uint8)
LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]
def enc_lines(lines):
    v = LUT[np.frombuffer(''.join(lines).encode(), np.uint8).reshape(-1, 31)]; k = np.zeros(len(v), np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
    rc = np.zeros(len(v), np.uint64)
    for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
    return np.minimum(k, rc)
fa = open(f'{O}/markers_u.fa').read().split('\n'); c2 = [fa[i + 1] for i in range(0, len(fa) - 1, 2) if fa[i].startswith('>chr2depth_')]
CHR2 = np.zeros(NU, bool); CHR2[np.searchsorted(U, enc_lines(c2))] = True; log('superset', NU, 'chr2depth', CHR2.sum())
# truth per male and class
sys.path.insert(0, f'{P}/scripts'); import truth_cl   # dyz3, yq12 from truth.load; units from indep_setB, dyz19_clean
T = truth_cl.load(O)
# regions: male -> assembly per class
A = {}
for l in open(f'{O}/regions.tsv'): a, c = l.split('\t')[:2]; A.setdefault(a.replace('_R2', ''), {})[c] = a
# per-marker counts per male (cached matrix)
order = [l.strip() for l in open(f'{O}/counts_order.txt')]; oi = np.searchsorted(U, enc_lines(order))
if not (U[oi] == enc_lines(order)).all(): raise SystemExit('ERROR counts order not in superset')
males = sorted(d.split('/')[-1] for d in glob.glob(f'{P}/yruler_u/*') if os.path.exists(f'{d}/counts.col.zst'))
cache = f'{O}/counts_{len(males)}.npy'
if os.path.exists(cache): C = np.load(cache, mmap_mode='r')
else:
    C = np.zeros((len(males), NU), np.uint32)
    for j, s in enumerate(males):
        v = np.array(sp.run([Z, '-dcq', f'{P}/yruler_u/{s}/counts.col.zst'], check=True, capture_output=True).stdout.split(), np.uint32)
        if len(v) != len(order): raise SystemExit(f'ERROR {s}: {len(v)} counts vs {len(order)} markers')
        C[j, oi] = v
    np.save(cache, C)
MI = {s: j for j, s in enumerate(males)}; log('males with counts', len(males))
# bits per assembly, unpacked lazily
B = {}
def bits(a, c):
    if a not in B: z = np.load(f'{O}/bits/{a}.npz'); B[a] = {k: z[k] for k in z.files}
    r = c if c in B[a] else None
    if c in REG: return np.logical_or.reduce([np.unpackbits(B[a][x])[:NU].astype(bool) for x in REG[c] if x in B[a]])
    return np.unpackbits(B[a][r])[:NU].astype(bool) if r else np.zeros(NU, bool)
def asm_of(s, c): return A.get(s, {}).get(REG.get(c, (c,))[0])
# candidates per class: superset membership on all assemblies (presence >= 40%)
CAND = {}
for c in CL:
    al = [asm_of(s, c) for s in A if asm_of(s, c)]; pres = np.zeros(NU, np.uint16)
    for a in al: pres += bits(a, c)
    CAND[c] = np.flatnonzero(pres >= math.ceil(FULL * len(al))); log(f'{c}: assemblies {len(al)}, candidates {len(CAND[c])}')
def build(train, c, thr):   # one rule: presence >= thr in training regions, absent from training Y outside the class
    al = [asm_of(s, c) for s in train if asm_of(s, c)]; idx = CAND[c]; pres = np.zeros(len(idx), np.uint16); out = np.zeros(len(idx), bool)
    for a in al: r = bits(a, c)[idx]; pres += r; out |= bits(a, 'whole')[idx] & ~r
    return idx[(pres >= math.ceil(thr * len(al))) & ~out]
def gcsum(j, m, ysc):   # yruler gc_sums: G+C curve from single-copy Y markers
    y = C[j, ysc].astype(np.float64); g = gc[ysc]; med = np.median(y)
    tm = lambda a, lo, hi: a[(a >= lo * med) & (a <= hi * med)].mean() if ((a >= lo * med) & (a <= hi * med)).any() else np.nan
    al = tm(y, 0.5, 1.5)
    if med < 10 or not al > 0: return C[j, m].sum(dtype=np.float64)
    d = [(b, tm(y[g == b], 0.25, 2.5) / al) for b in range(32) if (g == b).sum() >= 1000]
    r = np.array([min(d, key=lambda e: abs(e[0] - b))[1] for b in range(32)]); f = (np.bincount(g, minlength=32) * r).sum() / len(g)
    return (C[j, m] / r[gc[m]]).sum() * f
def zval(s, mk):
    j = MI[s]; ysc = mk['ysinglecopy']; Y = C[j, ysc].sum(dtype=np.float64)
    return {c: gcsum(j, mk[c], ysc) / Y for c in mk if c != 'ysinglecopy'}
def fit(train, test, thr_fixed=None):   # markers + K on train; threshold by K-only LOO on train
    out = {}
    for c in [x for x in CL if x != 'ysinglecopy']:
        tr = [s for s in train if s in T[c] and s in MI]; te = [s for s in test if s in T[c] and s in MI]
        if len(tr) < 5 or not te: continue
        best = None
        for thr in ([thr_fixed] if thr_fixed else THR):
            mk = {'ysinglecopy': build(train, 'ysinglecopy', thr), c: build(train, c, thr)}
            z = {s: zval(s, mk)[c] for s in tr + te}; St, Sz = sum(T[c][s] for s in tr), sum(z[s] for s in tr)
            loo = np.median([abs((St - T[c][s]) / (Sz - z[s]) * z[s] / T[c][s] - 1) for s in tr])
            if best is None or loo < best[0]: best = (loo, thr, St / Sz, z, len(mk[c]))
        loo, thr, K, z, n = best; out[c] = (thr, n, {s: 100 * (K * z[s] / T[c][s] - 1) for s in te})
    return out
def summ(e):
    a = np.abs(np.array(list(e.values()))); return f'n {len(a)} median |err| {np.median(a):.2f}% p90 {np.percentile(a, 90):.2f}% bias {np.median(list(e.values())):+.2f}%'
fam = truth_cl.families(f'{O}/ped.txt'); everyone = sorted(set(MI) & set().union(*[set(T[c]) for c in T]))
grp = lambda s: fam.get(s, s)
mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
if mode in ('all', 'oldsplit'):   # old split 22/15 with the new rule
    for c, f in (('dyz3', 'holdout_pub.dyz3.txt'), ('yq12', 'holdout_pub.yq12.txt')):
        h = {l.strip() for l in open(f'{O}/{f}') if l.strip()}; r = fit([s for s in everyone if s not in h], [s for s in everyone if s in h])
        if c in r: log(f'oldsplit {c}: threshold {r[c][0]}, markers {r[c][1]}, {summ(r[c][2])}')
if mode in ('all', 'holdout'):   # fixed held-out: Miga/Altemose males, relatives out of training
    h = {l.strip() for l in open(f'{O}/holdout_lit.txt') if l.strip()}; hg = {grp(s) for s in h}
    r = fit([s for s in everyone if grp(s) not in hg], [s for s in everyone if s in h])
    for c in r: log(f'holdout {c}: threshold {r[c][0]}, markers {r[c][1]}, {summ(r[c][2])}')
if mode in ('all', 'cv'):   # repeated grouped 5-fold
    rng = np.random.default_rng(7); G = sorted({grp(s) for s in everyone}); E = {}
    for rep in range(REP):
        g = rng.permutation(G); f = {x: i % FOLD for i, x in enumerate(g)}
        for k in range(FOLD):
            r = fit([s for s in everyone if f[grp(s)] != k], [s for s in everyone if f[grp(s)] == k])
            for c in r:
                E.setdefault(c, {}).setdefault(rep, {}).update(r[c][2]); E[c].setdefault('thr', []).append(r[c][0])
        log(f'cv repeat {rep + 1} done')
    for c in E:
        med = [np.median(np.abs(list(E[c][rep].values()))) for rep in range(REP)]
        log(f'cv {c}: males {len(E[c][0])}, median |err| per repeat mean {np.mean(med):.2f}% (min {min(med):.2f}, max {max(med):.2f}), threshold 0.5 in {E[c]["thr"].count(0.5)} of {len(E[c]["thr"])} folds')
if mode in ('all', 'production'):   # all truth males: final markers and K
    r = fit(everyone, everyone)
    for c in r: log(f'production {c}: threshold {r[c][0]}, markers {r[c][1]}')
