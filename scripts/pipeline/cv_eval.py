# Marker choice and K per class: grouped repeated 5-fold CV, lineage-out, fixed held-out, length range, production model
# usage: (environment of v1.sh) python3 cv_eval.py cv|lolo|holdout|range|production
import os, sys, math, collections, subprocess as sp, numpy as np
P = os.environ['P']; O = f'{P}/{os.environ["UM"]}'; IN = os.environ['IN']; Z = 'zstd'
THR = tuple(map(float, os.environ.get('THRS', '0.5,0.9').split(','))); REP, FOLD = int(os.environ.get('REP', 10)), 5; FULL = float(os.environ.get('FULL', 0.4))   # superset presence, as in um_super.py
CL = [c for c in os.environ.get('CLASSES', 'dyz3,yq12,ysinglecopy,dyz1,dyz2,dyz18s,yq2k7,dyz18sub,dyz19,dyz17').split(',')]
REG = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7'), 'dyz18suba': ('dyz18', 'yq3k1', 'yq2k7a')}   # merged classes: union of regions
def log(*a): print(*a, flush=True)
mode = sys.argv[1] if len(sys.argv) > 1 else ''
if mode not in ('cv', 'lolo', 'holdout', 'range', 'production'): raise SystemExit('ERROR usage: cv_eval.py cv|lolo|holdout|range|production')
# superset keys, G+C
U = np.load(f'{O}/u_keys.npy'); NU = len(U)
gc = np.zeros(NU, np.uint8)
for i in range(31): b = (U >> np.uint64(2 * i)) & np.uint64(3); gc += ((b == 1) | (b == 2)).astype(np.uint8)
log('superset', NU)
# truth per male and class
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import truth   # dyz3, yq12 from truth.load; units from indep_setB, dyz19_clean
T = truth.load(f'{IN}/truth')
# regions: male -> assembly per class
A = {}
for l in open(f'{O}/regions.tsv'): a, c = l.split('\t')[:2]; A.setdefault(a.replace('_R2', ''), {})[c] = a
# per-marker counts per male (cached matrix)
CN = os.environ.get('CNPY', 'cnpy')   # counts dir: short reads or nanopore
males = sorted(f[:-4] for f in os.listdir(f'{O}/{CN}') if f.endswith('.npy'))   # per-male counts from um_cnpy.py
C = np.stack([np.load(f'{O}/{CN}/{s}.npy', mmap_mode='r') for s in males])
MI = {s: j for j, s in enumerate(males)}; log('males with counts', len(males), 'from', CN)
FEMNPZ = os.environ.get('FEMNPZ')   # female counts at class keys
if FEMNPZ:
    z = np.load(FEMNPZ); FX, FIDX, FN, FD2 = z['F'], z['FIDX'], [str(x) for x in z['names']], z['d2'].astype(np.float64); FPOS = np.full(NU, -1, np.int64); FPOS[FIDX] = np.arange(len(FIDX)); log('females', len(FN), 'class keys', len(FIDX))
# bits per assembly, unpacked lazily
B = {}
def bits(a, c):
    if a not in B: z = np.load(f'{O}/bits/{a}.npz'); B[a] = {k: z[k] for k in z.files}
    r = c if c in B[a] else None
    if c in ('yq12rest2', 'out_yq12rest2'):   # rest of Yq12 without DYZ18S: Yq12 outside DYZ1, DYZ2 and DYZ18S; outside it: Y outside Yq12, DYZ1, DYZ2, DYZ18S
        y, d1, d2, d18 = (bits(a, x) for x in ('yq12', 'dyz1', 'dyz2', 'dyz18s')); return y & ~d1 & ~d2 & ~d18 if c == 'yq12rest2' else bits(a, 'out_yq12') | d1 | d2 | d18
    if c in ('yq12rest', 'out_yq12rest'):   # rest of Yq12: Yq12 sequence outside DYZ1 and DYZ2; outside it: Y outside Yq12, DYZ1, DYZ2
        y, d1, d2 = (bits(a, x) for x in ('yq12', 'dyz1', 'dyz2')); return y & ~d1 & ~d2 if c == 'yq12rest' else bits(a, 'out_yq12') | d1 | d2
    if c in REG: return np.logical_or.reduce([np.unpackbits(B[a][x])[:NU].astype(bool) for x in REG[c] if x in B[a]])
    return np.unpackbits(B[a][r])[:NU].astype(bool) if r else np.zeros(NU, bool)
def asm_of(s, c): return A.get(s, {}).get(REG.get(c, ('yq12',) if c in ('yq12rest', 'yq12rest2') else (c,))[0])
# candidates per class: superset membership on all assemblies (presence >= FULL)
CAND, NALL = {}, {}
for c in CL:
    al = [asm_of(s, c) for s in A if asm_of(s, c)]; pres = np.zeros(NU, np.uint16); NALL[c] = len(al)
    for a in al: pres += bits(a, c)
    CAND[c] = np.flatnonzero(pres >= math.ceil(FULL * len(al))); log(f'{c}: assemblies {len(al)}, candidates {len(CAND[c])}')
OUTLAB = os.environ.get('OUTLAB') == '1'   # outside = Y outside Yq12 or another annotated Yq12 class: unannotated Yq12 sequence is neither inside nor outside
OLC = {'dyz1': ('dyz2', 'dyz18', 'yq3k1', 'yq2k7a'), 'dyz2': ('dyz1', 'dyz18', 'yq3k1', 'yq2k7a'), 'dyz18s': ('dyz1', 'dyz2', 'yq2k7a')}
if OUTLAB: log('OUTLAB: outside of dyz1, dyz2, dyz18s = Y outside Yq12 or another annotated Yq12 class')
def outb(a, c): return (bits(a, 'out_yq12') | np.logical_or.reduce([bits(a, k) for k in OLC[c]])) if OUTLAB and c in OLC else bits(a, 'out_' + c)
RB, WB = {}, {}   # class bits and outside-class Y bits at the class candidates, per assembly
for c in CL:
    al = [asm_of(s, c) for s in A if asm_of(s, c)]
    RB[c] = {a: bits(a, c)[CAND[c]] for a in al}; WB[c] = {a: outb(a, c)[CAND[c]] for a in al}   # Y sequence outside the class regions
B.clear(); log('bits cached')
OUTF = tuple(map(float, os.environ.get('OUTFS', '0,0.05,0.1,0.25').split(',')))   # outside-class presence allowed: any assembly (0) or below this fraction
NYS = os.environ.get('NYSCORE')   # non-Y genome counts: any k-mer seen in a non-Y genome is dropped from markers and ruler (NYSTRICT=1)
if NYS and os.environ.get('NYSTRICT') != '1': raise SystemExit('ERROR NYSCORE needs NYSTRICT=1')
RX = np.load(os.environ['RULEX']) > 0 if os.environ.get('RULEX') else None   # RULEX: keys dropped from the ruler
if RX is not None: log('ruler exclusion', os.environ['RULEX'], int(RX.sum()), 'keys')
CX = np.load(os.environ['CLSX']) > 0 if os.environ.get('CLSX') else None   # CLSX: keys dropped from every class except the ruler
if CX is not None: log('class exclusion', os.environ['CLSX'], int(CX.sum()), 'keys')
if NYS: NY, NYN = np.load(NYS).astype(np.uint16), int(os.environ['NYN']); log(f'non-Y genomes: {NYN} from {NYS}; strict: any k-mer seen in one of them is dropped from markers and ruler')
def build(train, c, thr, of=0.0):   # one rule: present in >= thr of training class regions, outside the class in < of of them (of 0: in none)
    al = [asm_of(s, c) for s in train if asm_of(s, c) and (c, s) not in XM]; pres = np.zeros(len(CAND[c]), np.uint16); out = np.zeros(len(CAND[c]), np.uint16)
    if math.ceil(thr * len(al)) < math.ceil(FULL * NALL[c]):   # the prefilter on all assemblies could exclude a marker: not allowed
        raise SystemExit(f'ERROR prefilter can bind: {c} thr {thr}, training assemblies {len(al)} of {NALL[c]}')
    for a in al: pres += RB[c][a]; out += WB[c][a]
    ok = pres >= math.ceil(thr * len(al))
    if not NYS: return CAND[c][ok & (out < max(1, math.ceil(of * len(al))))]
    r_ = CAND[c][ok & (out < (1 if c == 'ysinglecopy' else max(1, math.ceil(of * len(al))))) & (NY[CAND[c]] == 0)]; return r_[~RX[r_]] if RX is not None and c == 'ysinglecopy' else r_[~CX[r_]] if CX is not None and c != 'ysinglecopy' else r_   # absent from every non-Y genome; of for Y outside the class only
def gcbins(y, g):   # G+C curve: ruler counts pooled per G+C bin, counts above 10 x (mean + 1) dropped; needs 1e6 ruler counts, bins 1,000 markers and 1e4 counts
    m = y.mean(); k = y <= 10 * (m + 1)
    if y.sum() < 1e6 or not m > 0: return None
    al = y[k].mean(); d = [(b, y[k & (g == b)].mean() / al) for b in range(32) if (g == b).sum() >= 1000 and y[k & (g == b)].sum() >= 1e4]
    return d or None
def gccurve(j, ysc):   # yruler gc_sums: G+C curve from single-copy Y markers -> (r per G+C, f, Y sum)
    y = C[j, ysc].astype(np.float64); g = gc[ysc]; Y = y.sum(); d = gcbins(y, g)
    if not d: return np.ones(32), 1.0, Y
    r = np.array([min(d, key=lambda e: abs(e[0] - b))[1] for b in range(32)]); return r, (np.bincount(g, minlength=32) * r).sum() / len(g), Y
RUL = os.environ.get('RULERS', 'match').split(',')   # denominators: match (Y presence = class p), yQQ (Y at QQ%)
if any(R != 'match' and not R[1:].isdigit() for R in RUL): raise SystemExit(f'ERROR rulers {RUL}: only match and yQQ')
rq = lambda R: None if R == 'match' else int(R[1:]) / 100
FT = 'inf'   # no k-mer filter; the rule key keeps this field (female_t in model.tsv)
MINGC = {k: int(v) for k, v in (x.split(':') for x in os.environ.get('MINGC', '').split(',') if x)}   # class: min G+C of 31
if MINGC: log('minimum marker G+C', MINGC)
def fit(train, test, thr_fixed=None, allrules=False):   # markers + K on train; thresholds and ruler by K-only LOO on train; allrules: also every fixed rule as <class>|r<p>_<o>
    out = {}; Z = {}; YS = {}; who = sorted(s for s in set(train) | set(test) if s in MI); th = [thr_fixed] if thr_fixed else list(THR)
    for q in sorted((set(th) if 'match' in RUL else set()) | {rq(R) for R in RUL if rq(R) is not None}):   # only the ruler thresholds in use
        YS[(q, FT)] = build(train, 'ysinglecopy', q)
    GY = {k: {s: gccurve(MI[s], YS[k]) for s in who} for k in YS}
    def arrays(key): R, t = key[2], key[3]; return {'ysc': YS[(key[0] if R == 'match' else rq(R), t)]}
    for p in th:
        for c in [x for x in CL if x != 'ysinglecopy']:
            na = sum(1 for s in train if asm_of(s, c))
            if math.ceil(p * na) < math.ceil(FULL * NALL[c]): log(f'{c}: presence {p} skipped, {na} training assemblies of {NALL[c]} (candidate prefilter could bind)'); continue
            for of in OUTF:
                m = build(train, c, p, of)
                if c in MINGC: m = m[gc[m] >= MINGC[c]]   # G+C floor
                if not len(m): continue
                Bm = {s: np.bincount(gc[m], weights=C[MI[s], m], minlength=32) for s in who}   # counts per G+C bin
                for R in RUL:
                    z = {}
                    for s in who: r, f, D = GY[(p if R == 'match' else rq(R), FT)][s]; z[s] = (Bm[s] / r).sum() * f / D
                    Z[((p, of, R, FT), c)] = (z, len(m), m)
    ALL = {}
    for c in [x for x in CL if x != 'ysinglecopy']:
        tr = [s for s in train if s in T[c] and s in MI]; te = [s for s in test if s in T[c] and s in MI]
        if len(tr) < 5 or not te: continue
        cd = {}
        for key in [k for k, cc in Z if cc == c]:
            z, n, mm = Z[(key, c)]; Tk = {s: T[c][s] for s in tr + te}; St, Sz = sum(Tk[s] for s in tr), sum(z[s] for s in tr)
            if any(Sz == z[s] for s in tr): continue
            cd[key] = (np.median([abs((St - Tk[s]) / (Sz - z[s]) * z[s] / Tk[s] - 1) for s in tr]), St / Sz, z, n, mm, Tk)
        ALL[c] = (cd, te)
    res = lambda c, key: (key, ALL[c][0][key][3], {s: 100 * (ALL[c][0][key][1] * ALL[c][0][key][2][s] / ALL[c][0][key][5][s] - 1) for s in ALL[c][1]}, ALL[c][0][key][1], ALL[c][0][key][4], arrays(key), ALL[c][0][key][5])   # [6]: truth used per male
    pick = lambda c, keys: min(keys, key=lambda k: ALL[c][0][k][0])   # first minimum, insertion order
    for c, (cd, te) in ALL.items():
        if not cd: log(f'ERROR {c}: no marker set left (all rules give no markers)'); continue
        out[c] = res(c, pick(c, list(cd)))   # thresholds and ruler per class, by train LOO
        if allrules: out.update({f'{c}|r' + '_'.join(map(str, k)): res(c, k) for k in cd})   # rule named by the full key (presence, outside, ruler, filter)
    return out
INNER = int(os.environ.get('INNER', 0)); NEST = None   # inner k-fold CV inside each training set (markers, ruler, K and truth rebuilt per inner fold)
def inner_fold(j):   # |errors| per fixed rule on inner fold j of NEST = (train, fold map)
    try: return inner_rule(j)
    except SystemExit as e: raise RuntimeError(f'inner fold {j}: {e}')   # a pool worker that exits hangs the parent
def inner_rule(j):
    tr, f = NEST; return {k: [abs(e) for e in v[2].values()] for k, v in fit([s for s in tr if f[grp(s)] != j], [s for s in tr if f[grp(s)] == j], allrules=True).items() if '|r' in k}
SCORE = os.environ.get('INNERSCORE', 'median'); SCOREF = {'median': np.median, 'mean': np.mean, 'p90': lambda v: np.percentile(v, 90)}[SCORE]   # inner-CV rule score
def nested(train, test):   # class rows: rule per class chosen by the inner CV; <class>|grule: one rule for all classes chosen by the inner CV
    global NEST; gi = sorted({grp(s) for s in train}); f = {x: i % INNER for i, x in enumerate(np.random.default_rng(11).permutation(gi))}; S = collections.defaultdict(list); NF = collections.Counter(); NEST = (train, f)
    np_ = min(INNER, int(os.environ.get('SLURM_CPUS_PER_TASK', 1)))
    if np_ > 1:   # inner folds in forked processes: counts shared, only the errors come back
        from multiprocessing import get_context
        with get_context('fork').Pool(np_) as pool: RS = pool.map(inner_fold, range(INNER))
    else: RS = [inner_fold(j) for j in range(INNER)]
    for r in RS:
        for k, v in r.items(): S[k] += v; NF[k] += 1
    sc = {k: SCOREF(v) for k, v in S.items() if NF[k] == INNER}; out = fit(train, test, allrules=True); res = {}   # only rules scored in every inner fold; INNERSCORE: median (default), mean or p90 of the inner |errors|
    for c in sorted({k.split('|')[0] for k in out if '|r' in k}):
        ks = [k for k in out if k.startswith(c + '|r') and k in sc]
        if ks: kb = min(ks, key=sc.get); res[c] = out[kb]; log(f'inner CV {c}: rule {kb.split("|")[1]}, inner {SCORE} |err| {sc[kb]:.2f}%')
    miss = [c for c in out if '|' not in c and c not in res]
    if miss: log(f'ERROR inner CV: no rule scored in all {INNER} inner folds for {miss}')
    ok = sorted(r for r in {k.split('|')[1] for k in sc} if all(f'{c}|{r}' in sc for c in res))
    if ok:
        g = min(ok, key=lambda r: np.mean([sc[f'{c}|{r}'] for c in res])); log(f'inner CV one rule for all classes: {g}, mean of inner class medians {np.mean([sc[f"{c}|{g}"] for c in res]):.2f}%')
        for c in list(res): res[f'{c}|grule'] = out[f'{c}|{g}']
    return res
fitc = lambda train, test: nested(train, test) if INNER else fit(train, test)
def summ(e):
    a = np.abs(np.array(list(e.values()))); return f'n {len(a)} median |err| {np.median(a):.2f}% p90 {np.percentile(a, 90):.2f}% bias {np.median(list(e.values())):+.2f}%'
fam = truth.families(os.environ.get('PED', f'{O}/ped.txt')); log('CV groups from', os.environ.get('PED', f'{O}/ped.txt')); everyone = sorted(set(MI) & set().union(*[set(T[c]) for c in T]))
grp = lambda s: fam.get(s, s)
if os.environ.get('EXCL'):   # external test males, not trained
    ex = set(open(os.environ['EXCL']).read().split()); everyone = [s for s in everyone if s not in ex]; log('EXCL', sorted(ex), 'males left', len(everyone))
if os.environ.get('TRUTHX'):   # drop truth pairs, folds unchanged
    for l in open(os.environ['TRUTHX']):
        c, s = l.split()[:2]; log(f'TRUTHX drop {c} {s}: {"dropped" if T[c].pop(s, None) is not None else "not in truth"}')
XM = set()   # (class, male) left out of marker selection too: MARKX=1 applies the TRUTHX exclusions to the marker regions
if os.environ.get('MARKX') == '1' and os.environ.get('TRUTHX'):
    for l in open(os.environ['TRUTHX']):
        c, s = l.split()[:2]; XM.add((c, s))
        if c == 'yq12': XM |= {('yq12rest', s), ('yq12rest2', s)}   # the rest classes sample the Yq12 region
    log('MARKX: marker regions left out', sorted(XM))
if 'yq12rest' in CL: T['yq12rest'] = {s: T['yq12'][s] - T['dyz1'][s] - T['dyz2'][s] for s in T['yq12'] if s in T['dyz1'] and s in T['dyz2']}   # truth of the rest of Yq12
if 'yq12rest2' in CL: T['yq12rest2'] = {s: T['yq12'][s] - T['dyz1'][s] - T['dyz2'][s] - T['dyz18s'][s] for s in T['yq12'] if all(s in T[k] for k in ('dyz1', 'dyz2', 'dyz18s'))}; log('yq12rest2 truth males', len(T['yq12rest2']))
if mode == 'holdout':   # fixed held-out list (HOLDLIST), relatives out of training
    h = {l.strip() for l in open(os.environ.get('HOLDLIST', f'{O}/holdout_lit.txt')) if l.strip()}; hg = {grp(s) for s in h}
    r = fitc([s for s in everyone if grp(s) not in hg], sorted(s for s in h if s in MI))
    for c in r: log(f'holdout {c}: threshold {r[c][0]}, markers {r[c][1]}, {summ(r[c][2])}')
    if os.environ.get('HOLDTSV'): open(os.environ['HOLDTSV'], 'w').write(''.join(f'{c}\t{s}\t{e}\t{r[c][6][s]:.0f}\n' for c in r for s, e in sorted(r[c][2].items())))   # per-male holdout errors
CA = {}   # alternative test counts (ALTDIR/<male>/counts.col.zst, same marker order; counts.col.npy: the same values in binary)
if os.environ.get('ALTDIR'):
    oi = np.load(f'{O}/counts_order_idx.npy')
    for s_ in everyone:
        fn = f"{P}/{os.environ['ALTDIR']}/{s_}/counts.col.zst"
        if not os.path.exists(fn): log('ALTDIR missing', s_); continue
        fz = fn[:-4] + '.npy'
        v = np.load(fz) if os.path.exists(fz) else np.array(sp.run([Z, '-dcq', fn], check=True, capture_output=True).stdout.split(), np.uint32)
        if len(v) != len(oi): raise SystemExit(f'ERROR {fn}: {len(v)} counts vs {len(oi)} markers')
        CA[s_] = np.zeros(NU, np.uint32); CA[s_][oi] = v
    log('ALTDIR', os.environ['ALTDIR'], 'males', len(CA))
if mode == 'cv':   # repeated grouped 5-fold; one repeat per REPSTART, errors per male saved
    G = sorted({grp(s) for s in everyone} | ({grp(s) for s in FN} if FEMNPZ else set())); r0 = int(os.environ.get('REPSTART', 0)); os.makedirs(f'{O}/{os.environ.get("CVDIR", "cv")}', exist_ok=True)
    for rep in range(r0, r0 + REP):
        g = np.random.default_rng(7 + rep).permutation(G); f = {x: i % FOLD for i, x in enumerate(g)}
        KS = [int(x) for x in os.environ.get('FOLDS', ','.join(map(str, range(FOLD)))).split(',')]   # folds of this task
        with open(f'{O}/{os.environ.get("CVDIR", "cv")}/rep{rep}' + ('' if len(KS) == FOLD else '_k' + '_'.join(map(str, KS))) + ('_c' + '-'.join(x for x in CL if x != 'ysinglecopy') if 'CLASSES' in os.environ else '') + '.tsv', 'w') as o, open(os.devnull if not FEMNPZ else o.name[:-4] + '_fem.tsv', 'w') as fo:   # one file per class split
            for k in KS:
                te = [s for s in everyone if f[grp(s)] == k]; keep_rows = {}; FTE = [i for i, s in enumerate(FN) if f[grp(s)] == k] if FEMNPZ else []
                for s_ in te:   # test males read from ALTDIR (e.g. a depth fraction), training males at full depth
                    if s_ in CA: keep_rows[s_] = C[MI[s_]].copy(); C[MI[s_]] = CA[s_]
                r = fitc([s for s in everyone if f[grp(s)] != k], te)
                for s_, v in keep_rows.items(): C[MI[s_]] = v
                for c in r:
                    for s, e in r[c][2].items(): o.write(f'{rep}\t{k}\t{c}\t{s}\t{r[c][0]}\t{r[c][1]}\t{e:.4f}\t{r[c][6][s]:.0f}\n')
                    m, ysn = r[c][4], len(r[c][5]['ysc'])   # test females: bp at haploid depth
                    bad = int((FPOS[m] < 0).sum()) if FTE else 0
                    if bad: log(f'ERROR {c}: {bad} markers without female counts, test females not written')
                    for i in ([] if bad else FTE): fo.write(f'{rep}\t{k}\t{c}\t{FN[i]}\t{r[c][0]}\t{r[c][1]}\t{float(r[c][3]) * FX[i, FPOS[m]].sum(dtype=np.float64) / (ysn * FD2[i] / 2):.1f}\n')
        log(f'cv repeat {rep} done')
if mode == 'range':   # extrapolation: train on one length tertile of the class, test on others (relatives of test males out of training)
    cf = [x for x in CL if x != 'ysinglecopy'][0]; c0 = os.environ.get('RANGEBY', cf); xs = sorted((T[c0][s], s) for s in everyone if s in T[c0]); n = len(xs)   # RANGEBY: class whose true length sets the tertiles (yq12 for the composed Yq12)
    q = lambda a, b: {s for _, s in xs[round(a * n):round(b * n)]}   # males between two quantiles of the true length
    log(f'range {c0}: {n} males with truth; lengths {xs[0][0]:.0f}-{xs[-1][0]:.0f}; thirds at {xs[n // 3][0]:.0f} and {xs[2 * n // 3][0]:.0f} bp, sixths at {xs[n // 6][0]:.0f} and {xs[5 * n // 6][0]:.0f} bp')
    RD = f'{O}/{os.environ.get("RANGEDIR", "range")}'; os.makedirs(RD, exist_ok=True)
    with open(f'{RD}/{cf}.tsv', 'w') as o:
        for name, (trs, tes) in (('short_to_long', (q(0, 2 / 3), q(2 / 3, 1))), ('long_to_short', (q(1 / 3, 1), q(0, 1 / 3))), ('middle_to_extremes', (q(1 / 6, 5 / 6), q(0, 1 / 6) | q(5 / 6, 1)))):   # training two thirds, test outside its range
            tg = {grp(s) for s in tes}; tr = [s for s in everyone if s in trs and grp(s) not in tg]; te = sorted(tes)
            r = fitc(tr, te); log(f'range {c0} {name}: training {len(tr)} males ({len(trs) - len(tr)} relatives of test males out), test {len(te)}')
            for c in r:
                for s, e in r[c][2].items(): o.write(f'{name}\t{c}\t{s}\t{T[c0][s]:.0f}\t{r[c][6][s]:.0f}\t{e:.4f}\t{c0}\n')   # truth of the sorting class
    log('range done', RD)
if mode == 'lolo':   # leave one lineage out: major clade, R1a/R1b split
    L = dict(l.rstrip('\n').split('\t')[:2] for l in open(f'{IN}/lists/lineage.tsv') if l.strip())
    lab = lambda s: None if s not in L else L[s][:3] if L[s].startswith('R1') else L[s][0]
    log('lineages', dict(collections.Counter(lab(s) for s in everyone))); E = collections.defaultdict(dict)
    LD = f'{O}/{os.environ.get("LOLODIR", "lolo")}'; os.makedirs(LD, exist_ok=True); XS = sorted({lab(s) for s in everyone} - {None})
    for x in [v for v in os.environ.get('LIN', ','.join(XS)).split(',') if v in XS]:   # lineages of this task
        te = [s for s in everyone if lab(s) == x]; tg = {grp(s) for s in te}
        r = fitc([s for s in everyone if lab(s) not in (x, None) and grp(s) not in tg], te)
        with open(f'{LD}/{x}' + ('_c' + '-'.join(c for c in CL if c != 'ysinglecopy') if 'CLASSES' in os.environ else '') + '.tsv', 'w') as o:   # one file per lineage and class split
            for c in r:
                E[c].update(r[c][2]); log(f'lolo {x} {c}: threshold {r[c][0]}, markers {r[c][1]}, {summ(r[c][2])}')
                for s_, e in r[c][2].items(): o.write(f'{x}\t{c}\t{s_}\t{r[c][0]}\t{r[c][1]}\t{e:.4f}\t{r[c][6][s_]:.0f}\n')
    for c in E: log(f'lolo pooled {c} (this task): {summ(E[c])}')
if mode == 'production':   # all truth males: final markers and K
    r = fitc(everyone, everyone)
    for c in r: log(f'production {c}: threshold {r[c][0]}, markers {r[c][1]}, K {r[c][3]:.6g}; in-sample (train = test) {summ(r[c][2])}')
    # frozen model: superset indices per class, K, G+C markers
    r = {c: v for c, v in r.items() if '|' not in c}
    FM = os.environ.get('FINAL', 'final_model')   # output name
    np.savez(f'{O}/{FM}.npz', **{f'm_{c}': r[c][4] for c in r}, **{f'ysc_{c}': r[c][5]['ysc'] for c in r})
    with open(f'{O}/{FM}.tsv', 'w') as o:
        o.write('class\tpresence\toutside\truler\tfemale_t\tmarkers\tK\tysinglecopy_markers\n')
        for c in r: o.write(f'{c}\t{r[c][0][0]}\t{r[c][0][1]}\t{r[c][0][2]}\t{r[c][0][3]}\t{r[c][1]}\t{float(r[c][3])!r}\t{len(r[c][5]["ysc"])}\n')
    log('model written', f'{O}/{FM}.npz', f'{O}/{FM}.tsv')
