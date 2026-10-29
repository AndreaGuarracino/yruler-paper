# Uniform-rule markers: grouped repeated 5-fold CV, fixed held-out, production set
import os, sys, math, glob, collections, subprocess as sp, numpy as np
P = '/path/to/project'; O = f'{P}/{os.environ.get("UM", "markers_u")}'; Z = '/path/to/conda/zstd/1.5.7/bin/zstd'
THR = tuple(map(float, os.environ.get('THRS', '0.5,0.9').split(','))); REP, FOLD = int(os.environ.get('REP', 10)), 5; FULL = float(os.environ.get('FULL', 0.4))   # superset presence, as in um_super.py
CL = [c for c in os.environ.get('CLASSES', 'dyz3,yq12,ysinglecopy,dyz1,dyz2,dyz18s,yq2k7,dyz18sub,dyz19').split(',')]
REG = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7')}   # merged classes: union of regions
def log(*a): print(*a, flush=True)
if len(sys.argv) > 1 and sys.argv[1] == 'lolosum':   # pool lineage files
    E = collections.defaultdict(dict)
    for fn in sorted(glob.glob(f'{O}/{os.environ.get("LOLODIR", "lolo")}/*.tsv')):
        for l in open(fn): x, c, s_, th, n, e = l.rstrip('\n').split('\t'); E[c][s_] = float(e)
    for c in E: a = np.abs(list(E[c].values())); log(f'lolo pooled {c}: n {len(a)} median |err| {np.median(a):.2f}% p90 {np.percentile(a, 90):.2f}% bias {np.median(list(E[c].values())):+.2f}%')
    sys.exit(0)
if len(sys.argv) > 1 and sys.argv[1] == 'cvsum':   # summary over repeat files
    E = collections.defaultdict(lambda: collections.defaultdict(list)); TH = collections.defaultdict(list)
    for fn in sorted(glob.glob(f'{O}/{os.environ.get("CVDIR", "cv")}/rep*.tsv')):
        for l in open(fn):
            rep, k, c, s, thr, n, e = l.rstrip('\n').split('\t'); E[c][rep].append(abs(float(e))); TH[c].append(thr)
    for c in E:
        med = [np.median(v) for v in E[c].values()]
        log(f'cv {c}: repeats {len(med)}, males per repeat {len(next(iter(E[c].values())))}, median |err| mean {np.mean(med):.2f}% (min {min(med):.2f}, max {max(med):.2f}), choices (presence, outside) {dict(collections.Counter(TH[c]))}')
    sys.exit(0)
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
CN = os.environ.get('CNPY', 'cnpy')   # counts dir: cnpy (short reads), cnpy_hifi, cnpy_ont
males = sorted(f[:-4] for f in os.listdir(f'{O}/{CN}') if f.endswith('.npy'))   # per-male counts from um_cnpy.py
C = np.stack([np.load(f'{O}/{CN}/{s}.npy', mmap_mode='r') for s in males])
MI = {s: j for j, s in enumerate(males)}; log('males with counts', len(males), 'from', CN)
# bits per assembly, unpacked lazily
B = {}
def bits(a, c):
    if a not in B: z = np.load(f'{O}/bits/{a}.npz'); B[a] = {k: z[k] for k in z.files}
    r = c if c in B[a] else None
    if c in REG: return np.logical_or.reduce([np.unpackbits(B[a][x])[:NU].astype(bool) for x in REG[c] if x in B[a]])
    return np.unpackbits(B[a][r])[:NU].astype(bool) if r else np.zeros(NU, bool)
def asm_of(s, c): return A.get(s, {}).get(REG.get(c, (c,))[0])
# candidates per class: superset membership on all assemblies (presence >= 40%)
CAND, NALL, BIND = {}, {}, set()
for c in CL:
    al = [asm_of(s, c) for s in A if asm_of(s, c)]; pres = np.zeros(NU, np.uint16); NALL[c] = len(al)
    for a in al: pres += bits(a, c)
    CAND[c] = np.flatnonzero(pres >= math.ceil(FULL * len(al))); log(f'{c}: assemblies {len(al)}, candidates {len(CAND[c])}')
RB, WB = {}, {}   # class bits and outside-class Y bits at the class candidates, per assembly
for c in CL:
    al = [asm_of(s, c) for s in A if asm_of(s, c)]
    RB[c] = {a: bits(a, c)[CAND[c]] for a in al}; WB[c] = {a: bits(a, 'out_' + c)[CAND[c]] for a in al}   # Y sequence outside the class regions
B.clear(); log('bits cached')
OUTF = tuple(map(float, os.environ.get('OUTFS', '0,0.05,0.1,0.25').split(',')))   # outside-class presence allowed: any assembly (0) or below this fraction
TRIM = float(os.environ.get('TRIM', 0))   # per sample: drop this fraction of markers at each tail of count/depth
NYS = os.environ.get('NYSCORE')   # unified outside: non-Y genomes count as outside occurrences
if NYS: NY, NYN, NYR = np.load(NYS).astype(np.uint16), int(os.environ['NYN']), np.load(os.environ['NYRULER']) > 0; log(f'unified outside: {NYN} non-Y genomes from {NYS}; ruler drops {int(NYR.sum())} keys')
REFA = os.environ.get('REFASM')   # single-reference baseline: markers from this assembly only
if REFA: log('single-reference baseline: markers and ruler from', REFA)
def build(train, c, thr, of=0.0):   # one rule: present in >= thr of training class regions, outside the class in < of of them (of 0: in none)
    al = [asm_of(s, c) for s in ([REFA] if REFA else train) if asm_of(s, c)]; pres = np.zeros(len(CAND[c]), np.uint16); out = np.zeros(len(CAND[c]), np.uint16)
    if math.ceil(thr * len(al)) < math.ceil(FULL * NALL[c]) and (c, thr, len(al)) not in BIND:   # 40% prefilter on all assemblies can exclude markers
        BIND.add((c, thr, len(al))); log(f'WARNING prefilter can bind: {c} thr {thr}, training assemblies {len(al)} of {NALL[c]}')
    for a in al: pres += RB[c][a]; out += WB[c][a]
    ok = pres >= math.ceil(thr * len(al))
    if not NYS: return CAND[c][ok & (out < max(1, math.ceil(of * len(al))))]
    if c == 'ysinglecopy': return CAND[c][ok & (out < 1) & ~NYR[CAND[c]]]   # ruler: CHM13 rule only
    return CAND[c][ok & (out + NY[CAND[c]] < max(1, math.ceil(of * (len(al) + NYN))))]   # Y outside the class + non-Y genomes
NOGC = os.environ.get('NOGC') == '1'   # no G+C correction
if NOGC: log('G+C correction off')
def gcbins(y, g):   # G+C curve: ruler counts pooled per G+C bin, counts above 10 x (mean + 1) dropped; needs 1e6 ruler counts, bins 1,000 markers and 1e4 counts
    m = y.mean(); k = y <= 10 * (m + 1)
    if y.sum() < 1e6 or not m > 0: return None
    al = y[k].mean(); d = [(b, y[k & (g == b)].mean() / al) for b in range(32) if (g == b).sum() >= 1000 and y[k & (g == b)].sum() >= 1e4]
    return d or None
def gccurve(j, ysc):   # yruler gc_sums: G+C curve from single-copy Y markers -> (r per G+C, f, Y sum)
    y = C[j, ysc].astype(np.float64); g = gc[ysc]; Y = y.sum(); d = None if NOGC else gcbins(y, g)
    if not d: return np.ones(32), 1.0, Y
    r = np.array([min(d, key=lambda e: abs(e[0] - b))[1] for b in range(32)]); return r, (np.bincount(g, minlength=32) * r).sum() / len(g), Y
def tsum(v):   # trimmed sum, rescaled to all markers (TRIM 0: plain sum)
    if not TRIM or len(v) < 100: return v.sum()
    lo, hi = np.quantile(v, [TRIM, 1 - TRIM]); k = (v >= lo) & (v <= hi); return v[k].sum() * len(v) / k.sum()
RUL = os.environ.get('RULERS', 'match').split(',')   # denominators: match (Y presence = class p), yQQ (Y at QQ%), chr2, yQQchr2 (geometric mean)
RSEL = os.environ.get('RULERS_SEL', ','.join(RUL)).split(',')   # rulers eligible for the nested choice; others reported only
rq = lambda R: None if R in ('match', 'chr2') else int(R[1:].replace('chr2', '')) / 100
C2I = np.flatnonzero(CHR2)
# k-mer filters: drop markers whose score in a source exceeds t
FILTS = os.environ.get('FILTS')   # 'score.npy=t1,t2;score2.npy=t3': several sources, labels s<k>:<t>
if FILTS:
    SRC = [x.split('=') for x in FILTS.split(';')]; FRS = [np.load(p) for p, _ in SRC]
    FT = ['inf'] + [f's{k}:{t}' for k, (_, ts) in enumerate(SRC) for t in ts.split(',')]
else:
    FT = os.environ.get('FEMT', 'inf').split(','); FRS = [np.load(os.environ.get('FILT', f'{O}/fem_ratio_filter.npy'))] if any(t != 'inf' for t in FT) else []
    FT = [t if t == 'inf' else f's0:{t}' for t in FT]
MINGC = {k: int(v) for k, v in (x.split(':') for x in os.environ.get('MINGC', '').split(',') if x)}   # class: min G+C of 31
if MINGC: log('minimum marker G+C', MINGC)
MGS = [int(x) for x in os.environ.get('MINGCS', '').split(',') if x]   # nested G+C floor grid
if MGS: log('nested minimum marker G+C grid', MGS)
def keep(idx, t):
    if t == 'inf': return idx
    k, v = t[1:].split(':'); return idx[FRS[int(k)][idx] <= float(v)]
def fit(train, test, thr_fixed=None):   # markers + K on train; thresholds and ruler by K-only LOO on train
    if TRIM: raise SystemExit('ERROR TRIM is not supported with binned sums')
    out = {}; Z = {}; YS = {}; who = sorted(s for s in set(train) | set(test) if s in MI); th = [thr_fixed] if thr_fixed else list(THR)
    for q in sorted(set(th) | {rq(R) for R in RUL if rq(R) is not None}):
        y = build(train, 'ysinglecopy', q)
        for t in FT: YS[(q, t)] = keep(y, t)
    GY = {k: {s: gccurve(MI[s], YS[k]) for s in who} for k in YS}; G2 = {s: gccurve(MI[s], C2I) for s in who} if any('chr2' in R for R in RUL) else {}
    def den(R, p, t, s):   # (G+C curve, f, denominator) of ruler R
        if R == 'chr2': return G2[s]
        r, f, Y = GY[(p if R == 'match' else rq(R), t)][s]; return (r, f, math.sqrt(Y * G2[s][2])) if R.endswith('chr2') else (r, f, Y)
    def arrays(key): R, t = key[2], key[3]; return {'ysc': YS[(key[0] if R == 'match' else rq(R), t)] if R != 'chr2' else np.zeros(0, np.int64), 'chr2': R.endswith('chr2')}
    for p in th:
        for c in [x for x in CL if x != 'ysinglecopy']:
            for of in OUTF:
                m0 = build(train, c, p, of)
                for t in FT:
                    m1 = keep(m0, t)
                    if c in MINGC: m1 = m1[gc[m1] >= MINGC[c]]   # AT-robust markers
                    for g in MGS or [None]:
                        m = m1 if g is None else m1[gc[m1] >= g]
                        if not len(m): continue
                        Bm = {s: np.bincount(gc[m], weights=C[MI[s], m], minlength=32) for s in who}   # counts per G+C bin
                        for R in RUL:
                            z = {}
                            for s in who: r, f, D = den(R, p, t, s); z[s] = (Bm[s] / r).sum() * f / D
                            Z[((p, of, R, t) + (() if g is None else (g,)), c)] = (z, len(m), m)
    ALL = {}
    for c in [x for x in CL if x != 'ysinglecopy']:
        tr = [s for s in train if s in T[c] and s in MI]; te = [s for s in test if s in T[c] and s in MI]
        if len(tr) < 5 or not te: continue
        cd = {}
        for key in [k for k, cc in Z if cc == c]:
            z, n, mm = Z[(key, c)]; St, Sz = sum(T[c][s] for s in tr), sum(z[s] for s in tr)
            if any(Sz == z[s] for s in tr): continue
            cd[key] = (np.median([abs((St - T[c][s]) / (Sz - z[s]) * z[s] / T[c][s] - 1) for s in tr]), St / Sz, z, n, mm)
        ALL[c] = (cd, te)
    res = lambda c, key: (key, ALL[c][0][key][3], {s: 100 * (ALL[c][0][key][1] * ALL[c][0][key][2][s] / T[c][s] - 1) for s in ALL[c][1]}, ALL[c][0][key][1], ALL[c][0][key][4], arrays(key))
    pick = lambda c, keys: min(keys, key=lambda k: ALL[c][0][k][0])   # first minimum, insertion order
    for c, (cd, te) in ALL.items():
        out[c] = res(c, pick(c, [k for k in cd if k[2] in RSEL]))   # nested: thresholds and ruler per class
        if len(RUL) > 1:
            for R in RUL:
                ks = [k for k in cd if k[2] == R]
                if ks: out[f'{c}|{R}'] = res(c, pick(c, ks))   # fixed ruler
        if len(MGS) > 1:
            for g in MGS:
                ks = [k for k in cd if k[4] == g and k[2] in RSEL]
                if ks: out[f'{c}|g{g}'] = res(c, pick(c, ks))   # fixed G+C floor
        if len(FT) > 1:
            for t in FT:
                ks = [k for k in cd if k[3] == t and k[2] in RSEL]
                if ks: out[f'{c}|t{t}'] = res(c, pick(c, ks))   # fixed filter
    if len(RUL) > 1 and ALL:   # one ruler for all classes, chosen on train: lowest mean of per-class best LOO
        sc = {R: np.mean([min(cd[k][0] for k in cd if k[2] == R) for cd, _ in ALL.values()]) for R in RSEL if all(any(k[2] == R for k in cd) for cd, _ in ALL.values())}
        Rg = min(sc, key=sc.get)
        for c, (cd, _) in ALL.items(): out[f'{c}|global'] = res(c, pick(c, [k for k in cd if k[2] == Rg]))
    if len(FT) > 1 and ALL:   # one female threshold for all classes, chosen on train
        sc = {t: np.mean([min(cd[k][0] for k in cd if k[3] == t and k[2] in RSEL) for cd, _ in ALL.values()]) for t in FT if all(any(k[3] == t and k[2] in RSEL for k in cd) for cd, _ in ALL.values())}
        tg = min(sc, key=sc.get)
        for c, (cd, _) in ALL.items(): out[f'{c}|tglobal'] = res(c, pick(c, [k for k in cd if k[3] == tg and k[2] in RSEL]))
    return out
def summ(e):
    a = np.abs(np.array(list(e.values()))); return f'n {len(a)} median |err| {np.median(a):.2f}% p90 {np.percentile(a, 90):.2f}% bias {np.median(list(e.values())):+.2f}%'
fam = truth_cl.families(f'{O}/ped.txt'); everyone = sorted(set(MI) & set().union(*[set(T[c]) for c in T]))
grp = lambda s: fam.get(s, s)
if os.environ.get('TRUTHX'):   # drop truth pairs, folds unchanged
    for l in open(os.environ['TRUTHX']):
        c, s = l.split()[:2]; log(f'TRUTHX drop {c} {s}: {"dropped" if T[c].pop(s, None) is not None else "not in truth"}')
if os.environ.get('TRUTHSET'):   # replace truth values: class, male, bp
    for l in open(os.environ['TRUTHSET']):
        c, s, v = l.split()[:3]
        if s in T[c]: log(f'TRUTHSET {c} {s}: {T[c][s]:.0f} -> {float(v):.0f}'); T[c][s] = float(v)
mode =sys.argv[1] if len(sys.argv) > 1 else 'all'
if mode in ('all', 'oldsplit'):   # old split 22/15 with the new rule
    for c, f in (('dyz3', 'holdout_pub.dyz3.txt'), ('yq12', 'holdout_pub.yq12.txt')):
        h = {l.strip() for l in open(f'{O}/{f}') if l.strip()}; r = fit([s for s in everyone if s not in h], [s for s in everyone if s in h])
        if c in r: log(f'oldsplit {c}: threshold {r[c][0]}, markers {r[c][1]}, {summ(r[c][2])}')
if mode in ('all', 'holdout'):   # fixed held-out: Miga/Altemose males, relatives out of training
    h = {l.strip() for l in open(f'{O}/holdout_lit.txt') if l.strip()}; hg = {grp(s) for s in h}
    r = fit([s for s in everyone if grp(s) not in hg], [s for s in everyone if s in h])
    for c in r: log(f'holdout {c}: threshold {r[c][0]}, markers {r[c][1]}, {summ(r[c][2])}')
    if os.environ.get('HOLDTSV'): open(os.environ['HOLDTSV'], 'w').write(''.join(f'{c}\t{s}\t{e}\n' for c in r for s, e in sorted(r[c][2].items())))   # per-male holdout errors
CA = {}   # alternative test counts (ALTDIR/<male>/counts.col.zst, same marker order)
if os.environ.get('ALTDIR'):
    oi = np.load(f'{O}/counts_order_idx.npy')
    for s_ in everyone:
        fn = f"{P}/{os.environ['ALTDIR']}/{s_}/counts.col.zst"
        if not os.path.exists(fn): log('ALTDIR missing', s_); continue
        v = np.array(sp.run([Z, '-dcq', fn], check=True, capture_output=True).stdout.split(), np.uint32)
        if len(v) != len(oi): raise SystemExit(f'ERROR {fn}: {len(v)} counts vs {len(oi)} markers')
        CA[s_] = np.zeros(NU, np.uint32); CA[s_][oi] = v
    log('ALTDIR', os.environ['ALTDIR'], 'males', len(CA))
if mode in ('all', 'cv'):   # repeated grouped 5-fold; one repeat per REPSTART, errors per male saved
    G = sorted({grp(s) for s in everyone}); r0 = int(os.environ.get('REPSTART', 0)); os.makedirs(f'{O}/{os.environ.get("CVDIR", "cv")}', exist_ok=True)
    for rep in range(r0, r0 + REP):
        g = np.random.default_rng(7 + rep).permutation(G); f = {x: i % FOLD for i, x in enumerate(g)}
        KS = [int(x) for x in os.environ.get('FOLDS', ','.join(map(str, range(FOLD)))).split(',')]   # folds of this task
        with open(f'{O}/{os.environ.get("CVDIR", "cv")}/rep{rep}' + ('' if len(KS) == FOLD else '_k' + '_'.join(map(str, KS))) + ('_c' + '-'.join(x for x in CL if x != 'ysinglecopy') if 'CLASSES' in os.environ else '') + '.tsv', 'w') as o:   # one file per class split
            for k in KS:
                te = [s for s in everyone if f[grp(s)] == k]; keep_rows = {}
                for s_ in te:   # test males read from ALTDIR (e.g. a depth fraction), training males at full depth
                    if s_ in CA: keep_rows[s_] = C[MI[s_]].copy(); C[MI[s_]] = CA[s_]
                r = fit([s for s in everyone if f[grp(s)] != k], te)
                for s_, v in keep_rows.items(): C[MI[s_]] = v
                for c in r:
                    for s, e in r[c][2].items(): o.write(f'{rep}\t{k}\t{c}\t{s}\t{r[c][0]}\t{r[c][1]}\t{e:.4f}\n')
        log(f'cv repeat {rep} done')
if mode == 'lolo':   # leave one lineage out: major clade, R1a/R1b split
    L = dict(l.rstrip('\n').split('\t')[:2] for l in open(f'{O}/lineage.tsv') if l.strip())
    lab = lambda s: None if s not in L else L[s][:3] if L[s].startswith('R1') else L[s][0]
    log('lineages', dict(collections.Counter(lab(s) for s in everyone))); E = collections.defaultdict(dict)
    LD = f'{O}/{os.environ.get("LOLODIR", "lolo")}'; os.makedirs(LD, exist_ok=True); XS = sorted({lab(s) for s in everyone} - {None})
    for x in [v for v in os.environ.get('LIN', ','.join(XS)).split(',') if v in XS]:   # lineages of this task
        te = [s for s in everyone if lab(s) == x]; tg = {grp(s) for s in te}
        r = fit([s for s in everyone if lab(s) not in (x, None) and grp(s) not in tg], te)
        with open(f'{LD}/{x}.tsv', 'w') as o:
            for c in r:
                E[c].update(r[c][2]); log(f'lolo {x} {c}: threshold {r[c][0]}, markers {r[c][1]}, {summ(r[c][2])}')
                for s_, e in r[c][2].items(): o.write(f'{x}\t{c}\t{s_}\t{r[c][0]}\t{r[c][1]}\t{e:.4f}\n')
    for c in E: log(f'lolo pooled {c} (this task): {summ(E[c])}')
if mode in ('all', 'production'):   # all truth males: final markers and K
    r = fit(everyone, everyone)
    for c in r: log(f'production {c}: threshold {r[c][0]}, markers {r[c][1]}, K {r[c][3]:.6g}; in-sample (train = test) {summ(r[c][2])}')
    # frozen model: superset indices per class, K, G+C markers
    r = {c: v for c, v in r.items() if '|' not in c}
    FM = os.environ.get('FINAL', 'final_model')   # output name
    np.savez(f'{O}/{FM}.npz', **{f'm_{c}': r[c][4] for c in r}, **{f'ysc_{c}': r[c][5]['ysc'] for c in r})
    with open(f'{O}/{FM}.tsv', 'w') as o:
        o.write('class\tpresence\toutside\truler\tfemale_t\tmarkers\tK\tysinglecopy_markers\tchr2' + ('\tmin_gc' if MGS else '') + '\n')
        for c in r: o.write(f'{c}\t{r[c][0][0]}\t{r[c][0][1]}\t{r[c][0][2]}\t{r[c][0][3]}\t{r[c][1]}\t{float(r[c][3])!r}\t{len(r[c][5]["ysc"])}\t{r[c][5]["chr2"]}' + (f'\t{r[c][0][4]}' if MGS else '') + '\n')
    log('model written', f'{O}/{FM}.npz', f'{O}/{FM}.tsv')
