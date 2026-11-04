# Where class marker copies lie on Y
# usage: python3 um_v5pos.py <model npz or per-class prefix> <u_keys npy> <out dir> [workers] | sup <out npz> (class candidates) | sum <scan dir> <out tsv> | ref <out dir> (depth denominators) | denfa <dir> <base> <fa> <out fa> | den <counts dir> <marker fa> <out tsv> | pgprep/pgc/pgscan/pgsum (pangenome copies of the depth sets) | cand <blocks> <dir> [step] (region keys) | copies <cand npz> <dir> (per-key copies in / outside the class) | lift <blocks.tsv> <bed> <paf list> <index> (CHM13 -> GRCh38, IMPG) | blocks <windows.tsv> <n> <out> | tile <GRCh38 chrom> <paf list> <index> <out> | tilesum <dir> <out> | proj <GRCh38 bed> <paf list> <agc> <out dir> [ROW=n] (IMPG) | projsum <out dir>; every markers_v5 assembly -> <out>/<asm>.tsv (copies in region, other Yq12, rest of Y) and <asm>.bins.tsv (100 kb bins outside the region)
import os, sys, gzip, numpy as np
from multiprocessing import Pool
LUT = np.full(256, 4, np.uint8); LUT[np.frombuffer(b'ACGTacgt', np.uint8)] = [0, 1, 2, 3, 0, 1, 2, 3]
def contigs(f):   # (name, sequence bytes)
    n, s = None, []
    for l in gzip.open(f, 'rb') if f.endswith('.gz') else open(f, 'rb'):
        if l.startswith(b'>'):
            if n: yield n, b''.join(s)
            n, s = l[1:].split()[0].decode(), []
        else: s.append(l.strip())
    if n: yield n, b''.join(s)
def contigs_bytes(b):   # (name, sequence) from FASTA bytes
    for r in b.split(b'>')[1:]:
        h, _, q = r.partition(b'\n'); yield h.split()[0].decode(), q.replace(b'\n', b'')
def canon(seq):   # canonical 31-mer codes from 1, 2, 4, 8, 16-base blocks (9 passes), valid mask
    b = LUT[np.frombuffer(seq, np.uint8)]; n = len(b) - 30
    if n <= 0: return np.zeros(0, np.uint64), np.zeros(0, bool)
    x = np.minimum(b, 3).astype(np.uint64); F, R = {1: x}, {1: np.uint64(3) - x}
    for w in (2, 4, 8, 16): h = w // 2; F[w] = (F[h][:-h] << np.uint64(2 * h)) | F[h][h:]; R[w] = R[h][:-h] | (R[h][h:] << np.uint64(2 * h))
    f, r, p = F[16], R[16], 16
    for w in (8, 4, 2, 1): f = (f[:len(F[w]) - p] << np.uint64(2 * w)) | F[w][p:]; r = r[:len(R[w]) - p] | (R[w][p:] << np.uint64(2 * p)); p += w
    bad = np.concatenate([[0], np.cumsum(b == 4)]); return np.minimum(f[:n], r[:n]), bad[31:] - bad[:n] == 0
def rowkeys(a):   # canonical codes of fixed-width 31-mer rows (uint8 ASCII, n x >= 31)
    b = LUT[a[:, :31]].astype(np.uint64); f = np.zeros(len(a), np.uint64); r = np.zeros(len(a), np.uint64)
    for q in range(31): f = (f << np.uint64(2)) | b[:, q]; r |= (np.uint64(3) - b[:, q]) << np.uint64(2 * q)
    return np.minimum(f, r)
def orderkeys(f):   # canonical code per count row of a yruler .order file (one row per distinct k-mer)
    OR = np.memmap(f, np.uint8, 'r').reshape(-1, 32); return np.concatenate([rowkeys(np.asarray(OR[x:x + 5000000])) for x in range(0, len(OR), 5000000)])
def inside(p, iv):   # positions inside half-open intervals
    if not iv: return np.zeros(len(p), bool)
    m = []
    for a, b in sorted(iv):   # merge overlapping and nested intervals
        if m and a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
        else: m.append([a, b])
    s = np.array([x[0] for x in m]); e = np.array([x[1] for x in m]); i = np.searchsorted(s, p, side='right') - 1
    return (i >= 0) & (p < e[np.maximum(i, 0)])
def best(WT, L, n):   # top n non-overlapping autosome blocks of L windows without N or censat, by mean unique share
    b = []
    for g in [f'chr{i}' for i in range(1, 23)]:
        sh, nb, cen = WT[g]; ok = (nb == 0) & (cen == 0); c = np.convolve(np.where(ok, sh, 0), np.ones(L), 'valid') / L; bad = np.convolve(~ok, np.ones(L), 'valid'); b += [(c[i], g, i) for i in range(len(c)) if bad[i] == 0]   # every window valid
    b.sort(reverse=True); top = []
    for m, g, i in b:
        if all(g != h or abs(i - j) >= L for _, h, j in top): top.append((m, g, i))
        if len(top) == n: break
    return top
if sys.argv[1] == 'blocks':   # [SD=<bed> [SDX=1]] blocks <windows.tsv> <n> <out>: chr2, chr17 and the top n blocks of 4.5 and 10 Mb from a windows table; SD bp per block, SDX=1 also excludes SD windows
    WT = {}; n = int(sys.argv[3])
    for l in list(open(sys.argv[2]))[1:]: g, b0, sh, nb, cen = l.split(); WT.setdefault(g, ([], [], [])); [a.append(float(v)) for a, v in zip(WT[g], (sh, nb, cen))]
    WT = {g: tuple(np.array(a) for a in v) for g, v in WT.items()}; SD = {}
    for l in open(os.environ['SD']) if os.environ.get('SD') else ():   # SD=<bed>: segmental duplication bp per window
        x = l.split('\t'); SD.setdefault(x[0], []).append((int(x[1]), int(x[2])))
    for g, iv in SD.items():
        m = []
        for a, b in sorted(iv):
            if m and a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
            else: m.append([a, b])
        SD[g] = m; sb = np.zeros(len(WT[g][0])) if g in WT else None
        for a, b in m if sb is not None else ():
            for k in range(a // 100000, min(len(sb) - 1, (b - 1) // 100000) + 1): sb[k] += min(b, (k + 1) * 100000) - max(a, k * 100000)
        if sb is not None and os.environ.get('SDX'): WT[g] = (WT[g][0], WT[g][1], WT[g][2] + sb)   # SDX=1: exclude SD windows
    sdbp = lambda g, s0, e0: sum(max(0, min(b, e0) - max(a, s0)) for a, b in SD.get(g, []))
    with open(sys.argv[4], 'w') as w:
        w.write('set\tchrom\tstart\tend\tmean_unique_share\tsd_bp\n')
        for nm, (g, s0, e0) in {'chr2': ('chr2', 50000000, 60000000), 'chr17': ('chr17', 72950000, 77450000)}.items(): w.write(f'{nm}\t{g}\t{s0}\t{e0}\t{WT[g][0][s0 // 100000:e0 // 100000].mean():.4f}\t{sdbp(g, s0, e0)}\n')
        for L in (45, 100):
            for r, (m, g, i) in enumerate(best(WT, L, n)): w.write(f'best{L / 10:g}Mb_{r + 1}\t{g}\t{i * 100000}\t{(i + L) * 100000}\t{m:.4f}\t{sdbp(g, i * 100000, (i + L) * 100000)}\n')
    print(open(sys.argv[4]).read(), flush=True); sys.exit(0)
def chm13_sorted(P):   # CHM13 without chrY: contigs and all 31-mers sorted
    C = [(n, q) for n, q in contigs(f'{P}/ext/ref/chm13v2.0.fa') if n != 'chrY']; A = np.empty(sum(max(0, len(q) - 30) for n, q in C), np.uint64); o = 0
    for n, q in C: k, ok = canon(q); A[o:o + len(k)] = np.where(ok, k, np.uint64(2 ** 64 - 1)); o += len(k)
    A.sort(); print(f'CHM13 without chrY: {len(C)} contigs, {len(A)} 31-mers sorted', flush=True); return C, A
if sys.argv[1] == 'cand':   # cand <blocks.tsv> <out dir> [step]: 31-mers once in CHM13 without chrY at every step-th position of each region (1: all) -> <base>_keys.npy, _set.npy, _sets.txt, base cand or cand<step>
    P = os.environ['P']; OD = sys.argv[3]; st = int(sys.argv[4]) if len(sys.argv) > 4 else 1; BS = 'cand' if st == 1 else f'cand{st}'; os.makedirs(OD, exist_ok=True); KK, SS, ST = [], [], []
    NOC = os.environ.get('NOCHM') == '1'   # NOCHM=1: no once-in-CHM13 filter (CHM13 counted as one haplotype later)
    if NOC: SQ = {n: q for n, q in contigs(os.environ.get('REFFA', f'{P}/ext/ref/chm13v2.0.fa')) if n != 'chrY'}   # REFFA: other reference (NOCHM=1 only)
    else: C, A = chm13_sorted(P); SQ = dict(C)
    for b in [l.split() for l in list(open(sys.argv[2]))[1:]]:
        k, ok = canon(SQ[b[1]][int(b[2]):int(b[3]) + 30]); k = np.unique(k[::st][ok[::st]])
        if not NOC: k = k[np.searchsorted(A, k, 'right') - np.searchsorted(A, k, 'left') == 1]
        ST.append(b[0]); KK.append(k); SS.append(np.full(len(k), len(ST) - 1, np.uint16)); print(f'{b[0]}: {len(k)} k-mers{"" if NOC else " once in CHM13 without chrY"} (every {st} positions)', flush=True)
    k = np.concatenate(KK); sid = np.concatenate(SS); o = np.argsort(k, kind='stable'); k, sid = k[o], sid[o]; d = np.concatenate([[True], k[1:] != k[:-1]])
    print(f'keys {len(k)}, in two regions {(~d).sum()} (kept in the first)', flush=True); np.save(f'{OD}/{BS}_keys.npy', k[d]); np.save(f'{OD}/{BS}_set.npy', sid[d]); open(f'{OD}/{BS}_sets.txt', 'w').write('\n'.join(ST) + '\n')
    os.system(f'cp {sys.argv[2]} {OD}/blocks.tsv'); sys.exit(0)
if sys.argv[1] == 'pgscan':   # pgscan <agc> <haplotype: set name or .fa.gz> <keys npy> <out npz>: genome-wide copies of each key in one haplotype (uint8, capped at 255)
    import subprocess; K = np.load(sys.argv[4]); ix = []; L = 0
    fa = gzip.open(sys.argv[3]).read() if sys.argv[3].endswith('.gz') else subprocess.run(['agc', 'getset', sys.argv[2], sys.argv[3]], check=True, capture_output=True).stdout   # haplotype: AGC set name or FASTA path
    for n, q in contigs_bytes(fa):
        for c0 in range(0, max(1, len(q) - 30), 50000000):   # 50 Mb chunks overlapping by 30 bp: memory bound
            k, ok = canon(q[c0:c0 + 50000030]); k = k[ok]; i = np.minimum(np.searchsorted(K, k), len(K) - 1); ix.append(i[K[i] == k])
        L += len(q)
    c = np.minimum(np.bincount(np.concatenate(ix), minlength=len(K)), 255).astype(np.uint8); np.savez_compressed(sys.argv[5] + '.tmp.npz', c=c); os.replace(sys.argv[5] + '.tmp.npz', sys.argv[5])
    print(f'{sys.argv[3]}: {L} bp, keys {len(K)}: one copy {(c == 1).mean():.4f}, absent {(c == 0).mean():.4f}, more {(c >= 2).mean():.4f}', flush=True); sys.exit(0)
if sys.argv[1] == 'blockproj':   # [ROWS=a-b] blockproj <blocks.tsv> <paf list> <impg index> <out tsv>: IMPG projection of each GRCh38 block in every haplotype (pieces merged within 10 kb)
    import subprocess; B = [l.split() for l in list(open(sys.argv[2]))[1:]]; IM = 'impg'
    if os.environ.get('ROWS'): a, b = map(int, os.environ['ROWS'].split('-')); B = B[a:b]
    CL = {l.split()[0]: int(l.split()[1]) for l in open(os.environ['P'] + '/ext/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai')}
    with open(sys.argv[5] + '.tmp', 'w') as o:
        for nm, g, s0, e0 in B:
            e0 = str(min(int(e0), CL[g]))   # blocks are cut at the chromosome end
            if int(e0) <= int(s0): print(f'{nm}: beyond {g} end {CL[g]}, skipped', flush=True); continue
            r = subprocess.run([IM, 'query', '--alignment-list', sys.argv[3], '-i', sys.argv[4], '-r', f'GRCh38#0#{g}:{int(s0) + 1}-{e0}', '-d', '10000', '-o', 'bed', '-t', os.environ.get('SLURM_CPUS_PER_TASK', '4')], check=True, capture_output=True).stdout.decode().splitlines()
            r = [x.split('\t')[:3] for x in r if not x.startswith('GRCh38#')]; o.writelines(f'{nm}\t{c}\t{a0}\t{b0}\n' for c, a0, b0 in r)
            print(f'{nm} {g}:{s0}-{e0}: pieces {len(r)}, haplotypes {len({c.rsplit("#", 1)[0] for c, _, _ in r})}', flush=True)
    os.replace(sys.argv[5] + '.tmp', sys.argv[5]); sys.exit(0)
if sys.argv[1] == 'pgloc':   # pgloc <agc> <haplotype: set name or .fa.gz> <dir with cand20_keys/set/sets> <projection tsv> <out npz>: copies of each key in the haplotype (c) and whether the key has one copy lying in the IMPG projection of its block (g, packed bits)
    import subprocess; D = sys.argv[4]; K = np.load(f'{D}/cand20_keys.npy'); sid = np.load(f'{D}/cand20_set.npy'); ST = open(f'{D}/cand20_sets.txt').read().split(); SI = {t: i for i, t in enumerate(ST)}
    fa = gzip.open(sys.argv[3]).read() if sys.argv[3].endswith('.gz') else subprocess.run(['agc', 'getset', sys.argv[2], sys.argv[3]], check=True, capture_output=True).stdout
    lp = np.zeros(len(K), np.int64); lc = np.full(len(K), -1, np.int32); CT = []; L = 0; ix = []
    for n, q in contigs_bytes(fa):
        CT.append(n)
        for c0 in range(0, max(1, len(q) - 30), 50000000):   # 50 Mb chunks overlapping by 30 bp
            k, ok = canon(q[c0:c0 + 50000030]); p = np.nonzero(ok)[0]; k = k[ok]; i = np.minimum(np.searchsorted(K, k), len(K) - 1); h = K[i] == k
            ix.append(i[h]); lp[i[h]] = c0 + p[h]; lc[i[h]] = len(CT) - 1   # position: the last hit (the only one when one copy)
        L += len(q)
    del fa; c = np.bincount(np.concatenate(ix), minlength=len(K)) if ix else np.zeros(len(K), np.int64); del ix; CI = {n: j for j, n in enumerate(CT)}; good = np.zeros(len(K), bool); good0 = np.zeros(len(K), bool); bo = np.argsort(sid, kind='stable'); bnd = np.searchsorted(sid[bo], np.arange(len(ST) + 1)); npj = 0
    for l in open(sys.argv[5]):
        t, ctg, a0, b0 = l.split('\t'); j = CI.get(ctg)
        if j is None or t not in SI: continue
        x = bo[bnd[SI[t]]:bnd[SI[t] + 1]]; m = (c[x] == 1) & (lc[x] == j); good[x[m & (lp[x] >= int(a0) - 10000) & (lp[x] < int(b0) + 10000)]] = True; good0[x[m & (lp[x] >= int(a0)) & (lp[x] < int(b0))]] = True; npj += 1   # g: within 10 kb of the projected ends, g0: strictly inside
    np.savez_compressed(sys.argv[6] + '.tmp.npz', c=np.minimum(c, 255).astype(np.uint8), g=np.packbits(good), g0=np.packbits(good0)); os.replace(sys.argv[6] + '.tmp.npz', sys.argv[6])
    print(f'{sys.argv[3]}: {L} bp, contigs {len(CT)}, projection pieces {npj}; keys {len(K)}: one copy {(c == 1).mean():.4f}, one copy in the projection of its block {good.mean():.4f} (strictly inside {good0.mean():.4f}), absent {(c == 0).mean():.4f}, more {(c >= 2).mean():.4f}', flush=True); sys.exit(0)
if sys.argv[1] == 'copies':   # [UM=markers_v5] copies <cand npz> <out dir> [workers]: per assembly and class, copies of each candidate key in the class region and elsewhere on Y
    P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_v5")}'; CF, OUT = sys.argv[2:4]; NW = int(sys.argv[4]) if len(sys.argv) > 4 else 8; U = np.load(f'{O}/u_keys.npy'); CN = np.load(CF)
    CS = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7'), 'dyz18suba': ('dyz18', 'yq3k1', 'yq2k7a')}; CK = {c: U[CN[c]] for c in CN.files if c != 'ysinglecopy'}   # u_keys sorted: keys in candidate order
    REG = {}
    for l in open(f'{O}/regions.tsv'): a, c, g, s0, e0 = l.rstrip('\n').split('\t')[:5]; REG.setdefault(a, []).append((c, g, int(s0), int(e0)))
    def cp(x):
        a, f = x
        if os.path.exists(f'{OUT}/{a}.npz'): return f'{a} done'
        R = REG.get(a, []); Z = {f'{c}_{w}': np.zeros(len(k), np.uint16) for c, k in CK.items() for w in ('in', 'out')}
        for ctg, q in contigs(f):
            k, ok = canon(q); pos = np.nonzero(ok)[0]; k = k[ok]
            for c, K in CK.items():
                i = np.minimum(np.searchsorted(K, k), len(K) - 1); h = K[i] == k
                if not h.any(): continue
                ir = inside(pos[h], [(s0, e0) for cc, g, s0, e0 in R if g == ctg and cc in CS.get(c, (c,))])
                Z[f'{c}_in'] += np.bincount(i[h][ir], minlength=len(K)).astype(np.uint16); Z[f'{c}_out'] += np.bincount(i[h][~ir], minlength=len(K)).astype(np.uint16)
        np.savez_compressed(f'{OUT}/{a}.tmp.npz', **Z); os.replace(f'{OUT}/{a}.tmp.npz', f'{OUT}/{a}.npz')
        return f'{a} ' + ' '.join(f'{c}:{int(Z[c + "_in"].sum())}/{int(Z[c + "_out"].sum())}' for c in CK)
    os.makedirs(OUT, exist_ok=True); A = [l.split()[:2] for l in open(f'{O}/assemblies.tsv')]; print(f'candidates ' + ' '.join(f'{c} {len(k)}' for c, k in CK.items()), flush=True)
    if not all(np.all(np.diff(k.astype(np.int64)) > 0) for k in CK.values()): sys.exit('ERROR candidate keys not sorted')
    with Pool(NW) as p:
        for r in p.imap_unordered(cp, A): print(r, flush=True)
    sys.exit(0)
if sys.argv[1] == 'ref':   # ref <out dir>: CHM13 31-mers once without chrY per 100 kb, best autosome blocks, depth marker FASTA
    P = os.environ['P']; OD = sys.argv[2]; os.makedirs(OD, exist_ok=True); W = 100000; log = lambda *a: print(*a, flush=True)
    CEN = {}
    for l in open(f'{P}/ext/ref/chm13v2.0_censat_v2.1.bed'):
        if not l.startswith('track'): x = l.split('\t'); CEN.setdefault(x[0], []).append((int(x[1]), int(x[2])))
    C, A = chm13_sorted(P)
    once = lambda k: np.searchsorted(A, k, 'right') - np.searchsorted(A, k, 'left') == 1
    WT, KS = {}, {}
    with open(f'{OD}/windows.tsv', 'w') as w:
        w.write('chrom\tstart\tunique_share\tn_bp\tcensat_bp\n')
        for n, q in C:
            k, ok = canon(q); p = np.arange(0, len(k), 10); u = ok[p] & once(k[p]); nw = len(q) // W + 1; cen = np.zeros(nw)
            for s0, e0 in CEN.get(n, []):
                for b in range(s0 // W, min(nw - 1, (e0 - 1) // W) + 1): cen[b] += min(e0, (b + 1) * W) - max(s0, b * W)
            sh = np.bincount(p // W, u, nw) / np.maximum(np.bincount(p // W, None, nw), 1); nb = np.bincount(np.nonzero(LUT[np.frombuffer(q, np.uint8)] == 4)[0] // W, None, nw)
            WT[n] = (sh, nb, cen); w.writelines(f'{n}\t{b * W}\t{sh[b]:.4f}\t{nb[b]:.0f}\t{cen[b]:.0f}\n' for b in range(nw))
    R = {'chr2depth': ('chr2', 50000000, 60000000), 'chr17': ('chr17', 72950000, 77450000)}   # yruler needs the classes chr2depth and ysinglecopy
    with open(f'{OD}/blocks.tsv', 'w') as w:
        w.write('set\tchrom\tstart\tend\tmean_unique_share\n')
        for nm, (g, s0, e0) in R.items(): w.write(f'{nm}\t{g}\t{s0}\t{e0}\t{WT[g][0][s0 // W:e0 // W].mean():.4f}\n')
        for L in (45, 100):
            top = best(WT, L, int(os.environ.get('TOPN', 5)))
            for r, (m, g, i) in enumerate(top): w.write(f'best{L * W // 100000 / 10:g}Mb_{r + 1}\t{g}\t{i * W}\t{(i + L) * W}\t{m:.4f}\n')
            R[f'best{L * W // 100000 / 10:g}Mb'] = (top[0][1], top[0][2] * W, (top[0][2] + L) * W)
    log(open(f'{OD}/blocks.tsv').read()); SQ = dict(C); seen = np.zeros(0, np.uint64); U = np.load(f'{P}/markers_v5/u_keys.npy'); fm = np.load(f'{P}/model_v5/model.npz')
    for nm, (g, s0, e0) in list(R.items()) + [('auto', (None, None, None))]:
        if g: k, ok = canon(SQ[g][s0:e0 + 30]); k = k[ok]
        else: k = np.concatenate([(lambda k, ok: k[::300][ok[::300]])(*canon(SQ[h])) for h in [f'chr{i}' for i in range(1, 23)]])   # every 300th position
        k = np.unique(k[once(k)]); KS[nm] = np.setdiff1d(k, seen, assume_unique=True); seen = np.union1d(seen, KS[nm]); log(f'{nm}: {len(KS[nm])} k-mers once in CHM13 without chrY')
    KS['ysinglecopy'] = np.unique(U[fm['ysc_dyz3']]); log(f'ysinglecopy (Y ruler): {len(KS["ysinglecopy"])} k-mers')
    with open(f'{OD}/denom.fa', 'w') as w:
        for nm, k in KS.items():
            for c0 in range(0, len(k), 1000000):   # decode in chunks
                b = np.array([(k[c0:c0 + 1000000] >> np.uint64(2 * (30 - i))) & np.uint64(3) for i in range(31)]).T.astype(np.uint8); S = np.frombuffer(b'ACGT', np.uint8)[b]
                w.writelines(f'>{nm}_{c0 + i}\n{S[i].tobytes().decode()}\n' for i in range(len(S)))
    log('written', f'{OD}/denom.fa'); sys.exit(0)
if sys.argv[1] == 'lift':   # lift <blocks.tsv> <out bed> <paf list> <impg index>: CHM13 regions -> GRCh38 with IMPG (CHM13 vs GRCh38 PAF), largest merged piece on the same chromosome
    import subprocess; B = [l.split() for l in list(open(sys.argv[2]))[1:]]; IM = 'impg'
    with open(sys.argv[3], 'w') as o:
        for nm, g, s0, e0 in ((b[0], b[1], int(b[2]), int(b[3])) for b in B):
            SR, DS = os.environ.get('SRC', 'CHM13'), os.environ.get('DST', 'GRCh38')   # SRC=GRCh38 DST=CHM13: the other way
            r = subprocess.run([IM, 'query', '--alignment-list', sys.argv[4], '-i', sys.argv[5], '-r', f'{SR}#0#{g}:{s0 + 1}-{e0}', '-d', '100000', '-o', 'bed', '-t', '4'], check=True, capture_output=True).stdout.decode().splitlines()
            iv = sorted((int(x.split()[1]), int(x.split()[2])) for x in r if x.split()[0] == f'{DS}#0#{g}')
            if not iv: print(f'ERROR {nm}: no GRCh38 projection', flush=True); continue
            a0, z0 = max(iv, key=lambda x: x[1] - x[0]); o.write(f'{g}\t{a0}\t{z0}\t{nm}\n')
            print(f'{nm} {SR} {g}:{s0}-{e0} -> {DS} {g}:{a0}-{z0} (length {z0 - a0} vs {e0 - s0}; pieces {len(iv)})', flush=True)
    sys.exit(0)
if sys.argv[1] == 'proj':   # proj <GRCh38 bed> <paf list> <agc[,fasta,..]> <out dir>: IMPG projection per haplotype: coverage, pieces, length, identity, indels, k-mers
    import subprocess, re, collections; BED, PL, AGC, D = sys.argv[2:6]; IM = 'impg'; os.makedirs(D, exist_ok=True)
    q = lambda r, a: subprocess.run([IM, 'query', '--alignment-list', PL, '-i', f'{D}/hprc.impg', '-r', r, '--no-merge', '-t', os.environ.get('SLURM_CPUS_PER_TASK', '8')] + a, check=True, stdout=subprocess.PIPE, stderr=open(f'{D}/impg.log', 'a')).stdout.decode()
    BL = [l for l in open(BED) if l.strip()]; BL = [BL[int(os.environ['ROW']) - 1]] if os.environ.get('ROW') else BL   # ROW: one region per array task
    o = open(f'{D}/proj_{BL[0].split()[3] if len(BL) == 1 else "all"}.tsv', 'w'); o.write('region\thaplotype\tpieces\tcontigs\ttarget_cover\tlength_ratio\tidentity\tsnv_per_kb\tindels_lt50\tindels_ge50\tindel_bp_ge50\tlargest_indel\tkmers_once\tkmers_absent\tkmers_multi\n')
    for l in BL:
        g, s0, e0, nm = l.split()[:4]; s0, e0 = int(s0), int(e0); L = e0 - s0; r = f'GRCh38#0#{g}:{s0 + 1}-{e0}'; H = collections.defaultdict(list)
        for x in q(r, ['-o', 'paf']).splitlines():
            x = x.split('\t')
            if not x[0].startswith('GRCh38'): H['#'.join(x[0].split('#')[:2])].append(x)
        SQ = collections.defaultdict(list); h = None
        for x in q(r, ['-o', 'fasta', '--sequence-files'] + AGC.split(',')).split('\n'):
            if x.startswith('>'): h = 'GRCh38' if x[1:].startswith('GRCh38') else '#'.join(x[1:].split('#')[:2]); SQ[h].append([])   # one list of lines per piece
            elif h: SQ[h][-1].append(x.encode())
        if 'GRCh38' not in SQ: sys.exit(f'ERROR {nm}: no GRCh38 sequence in the IMPG fasta')
        seq = lambda h: b'N'.join(b''.join(p) for p in SQ.get(h, []))   # pieces joined by N, lines joined directly
        k, ok = canon(seq('GRCh38')); u, n = np.unique(k[ok], return_counts=True); KS = u[n == 1]   # k-mers once in the GRCh38 region
        for h, R in sorted(H.items()):
            iv = sorted((int(x[7]), int(x[8])) for x in R); cov = 0; p = s0
            for a, b in iv: cov += max(0, min(b, e0) - max(a, p)); p = max(p, b)
            ops = [(int(n), c) for x in R for n, c in re.findall(r'(\d+)([=XID])', next((t[5:] for t in x[12:] if t.startswith('cg:Z:')), ''))]
            M = sum(n for n, c in ops if c == '='); X = sum(n for n, c in ops if c == 'X'); ind = [n for n, c in ops if c in 'ID']; big = [n for n in ind if n >= 50]
            k, ok = canon(seq(h)); i = np.minimum(np.searchsorted(KS, k[ok]), len(KS) - 1); c = np.bincount(i[KS[i] == k[ok]], minlength=len(KS))
            km = f'{(c == 1).mean():.4f}\t{(c == 0).mean():.4f}\t{(c >= 2).mean():.4f}'
            o.write(f'{nm}\t{h}\t{len(R)}\t{len({x[0] for x in R})}\t{cov / L:.4f}\t{sum(int(x[3]) - int(x[2]) for x in R) / L:.4f}\t{M / max(M + X + len(ind), 1):.5f}\t{1000 * X / max(M + X, 1):.3f}\t{len(ind) - len(big)}\t{len(big)}\t{sum(big)}\t{max(ind, default=0)}\t{km}\n')
        o.flush(); print(f'{nm} {r}: haplotypes {len(H)}, k-mers once in GRCh38 region {len(KS)}', flush=True)
    o.close(); sys.exit(0)
if sys.argv[1] == 'projsum':   # projsum <out dir>: one row per region from the proj_*.tsv files
    import csv; D = sys.argv[2]; R = [x for f in sorted(os.listdir(D)) if f.startswith('proj_') and f.endswith('.tsv') for x in csv.DictReader(open(f'{D}/{f}'), delimiter='\t')]; print('region\thaplotypes\tcover>=0.99\tcover_min\tlength_ratio_median\tmin\tmax\tidentity_median\tp10\tmin\tsnv_per_kb_median\tindels>=50_median\tmax\thaplotypes_with_indel>=1kb\tkmers_once_median\tp10\tmin', flush=True)
    for nm in dict.fromkeys(x['region'] for x in R):
        X = [x for x in R if x['region'] == nm]; f = lambda c: np.array([float(x[c]) for x in X if x[c] != '.'])
        cv, lr, idn, sn, sv, kk = f('target_cover'), f('length_ratio'), f('identity'), f('snv_per_kb'), f('indels_ge50'), f('kmers_once'); big = sum(int(x['largest_indel']) >= 1000 for x in X)
        print(f'{nm}\t{len(X)}\t{(cv >= 0.99).sum()}\t{cv.min():.4f}\t{np.median(lr):.4f}\t{lr.min():.4f}\t{lr.max():.4f}\t{np.median(idn):.5f}\t{np.percentile(idn, 10):.5f}\t{idn.min():.5f}\t{np.median(sn):.3f}\t{np.median(sv):.0f}\t{sv.max():.0f}\t{big}' + (f'\t{np.median(kk):.4f}\t{np.percentile(kk, 10):.4f}\t{kk.min():.4f}' if len(kk) else '\t.\t.\t.'), flush=True)
    sys.exit(0)
if sys.argv[1] == 'tile':   # tile <GRCh38 chrom> <paf list> <impg index> <out tsv> [block]: per 1 Mb block and haplotype, IMPG cover, identity, indels >= 50 bp
    import subprocess, re, collections; g, PL, IX, OUT = sys.argv[2:6]; BK = int(sys.argv[6]) if len(sys.argv) > 6 else 1000000; IM = 'impg'
    L = int(next(l.split()[1] for l in open(os.environ['P'] + '/ext/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai') if l.split()[0] == g))
    with open(OUT + '.tmp', 'w') as o:
        o.write('chrom\tstart\thaplotype\tcover\tidentity\tindel_bp_ge50\tlargest_indel\tpieces\n')
        for s0 in range(0, L - BK + 1, BK):
            r = subprocess.run([IM, 'query', '--alignment-list', PL, '-i', IX, '-r', f'GRCh38#0#{g}:{s0 + 1}-{s0 + BK}', '--no-merge', '-o', 'paf', '-t', os.environ.get('SLURM_CPUS_PER_TASK', '4')], check=True, capture_output=True).stdout.decode().splitlines()
            H = collections.defaultdict(list)
            for x in r:
                x = x.split('\t')
                if not x[0].startswith('GRCh38'): H['#'.join(x[0].split('#')[:2])].append(x)
            for h, R in H.items():
                cov, p = 0, s0
                for a, b in sorted((int(x[7]), int(x[8])) for x in R): cov += max(0, min(b, s0 + BK) - max(a, p)); p = max(p, b)
                ops = [(int(n), c) for x in R for n, c in re.findall(r'(\d+)([=XID])', next((t[5:] for t in x[12:] if t.startswith('cg:Z:')), ''))]
                M = sum(n for n, c in ops if c == '='); X = sum(n for n, c in ops if c == 'X'); ind = [n for n, c in ops if c in 'ID']
                o.write(f'{g}\t{s0}\t{h}\t{cov / BK:.4f}\t{M / max(M + X + len(ind), 1):.5f}\t{sum(n for n in ind if n >= 50)}\t{max(ind, default=0)}\t{len(R)}\n')
            o.flush()
    os.replace(OUT + '.tmp', OUT); print(f'{g}: {L // BK} blocks written {OUT}', flush=True); sys.exit(0)
if sys.argv[1] == 'tilesum':   # tilesum <tile dir> <out tsv> [window blocks]: 10 Mb windows (step 1 block) ranked by the median per-haplotype bp in SVs (indels >= 50 bp + uncovered bp)
    import glob, collections; D, OUT = sys.argv[2:4]; W = int(sys.argv[4]) if len(sys.argv) > 4 else 10; B = collections.defaultdict(dict); HS = set()
    for f in glob.glob(f'{D}/*.tsv'):
        for l in list(open(f))[1:]:
            g, s0, h, cv, idn, sv, big, pc = l.split('\t'); B[(g, int(s0))][h] = (float(cv), float(idn), int(sv), int(big)); HS.add(h)
    HS = sorted(HS); rows = []
    for g in sorted({k[0] for k in B}, key=lambda x: int(x[3:]) if x[3:].isdigit() else 99):
        st = sorted(s0 for gg, s0 in B if gg == g); bk = int(np.diff(st).min()) if len(st) > 1 else 1000000   # block size: smallest step (a missing first block made chr14 steps 2 Mb)
        for s0 in st:
            ks = [(g, s0 + i * bk) for i in range(W)]
            if not all(k in B for k in ks): continue
            sv = np.array([sum(B[k].get(h, (0, 0, 0, 0))[2] + (1 - B[k].get(h, (0, 0, 0, 0))[0]) * bk for k in ks) for h in HS])   # bp in indels >= 50 + bp not covered
            idn = np.array([np.mean([B[k][h][1] for k in ks if h in B[k]] or [0]) for h in HS]); full = np.mean([all(B[k].get(h, (0,))[0] >= 0.99 for k in ks) for h in HS])
            rows.append((np.median(sv), g, s0, s0 + W * bk, np.percentile(sv, 90), sv.max(), full, np.median(idn)))
    rows.sort()
    with open(OUT, 'w') as o:
        o.write('rank\tchrom\tstart\tend\tsv_bp_median\tsv_bp_p90\tsv_bp_max\thaplotypes_full_cover\tidentity_median\n')
        for i, (m, g, a, b, p9, mx, fu, idn) in enumerate(rows, 1): o.write(f'{i}\t{g}\t{a}\t{b}\t{m:.0f}\t{p9:.0f}\t{mx:.0f}\t{fu:.4f}\t{idn:.5f}\n')
    print(f'{len(rows)} windows of {W} blocks, {len(HS)} haplotypes; written {OUT}', flush=True); sys.exit(0)
if sys.argv[1] == 'pgprep':   # pgprep <dir> [base]: sorted keys and set ids of <base>.fa (default denom)
    OD = sys.argv[2]; BS = sys.argv[3] if len(sys.argv) > 3 else 'denom'; L = open(f'{OD}/{BS}.fa').read().split('\n'); nm = [h[1:].rsplit('_', 1)[0] for h in L[0::2] if h]; q = ''.join(L[1::2]).encode(); del L
    b = LUT[np.frombuffer(q, np.uint8).reshape(-1, 31)].astype(np.uint64); k = np.zeros(len(b), np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | b[:, i]
    ST = sorted(set(nm)); sid = np.array([ST.index(x) for x in nm], np.uint8); o = np.argsort(k)
    np.save(f'{OD}/{BS}_keys.npy', k[o]); np.save(f'{OD}/{BS}_set.npy', sid[o]); open(f'{OD}/{BS}_sets.txt', 'w').write('\n'.join(ST) + '\n')
    print(f'keys {len(k)}, unique {len(np.unique(k))}, sets ' + ', '.join(f'{x} {(sid == i).sum()}' for i, x in enumerate(ST)), flush=True); sys.exit(0)
if sys.argv[1] == 'pgc':   # pgc <KMC dump> <keys npy> <out npy>: copies per key in one haplotype (uint8, 255 = 255 or more)
    K = np.load(sys.argv[3]); x = open(sys.argv[2]).read().split(); c = np.zeros(len(K), np.uint8)
    if len(x):
        b = LUT[np.frombuffer(''.join(x[0::2]).encode(), np.uint8).reshape(-1, 31)].astype(np.uint64); k = np.zeros(len(b), np.uint64)
        for i in range(31): k = (k << np.uint64(2)) | b[:, i]
        i = np.searchsorted(K, k)
        if (i >= len(K)).any() or (K[np.minimum(i, len(K) - 1)] != k).any(): sys.exit(f'ERROR {sys.argv[2]}: k-mer not in the keys')
        c[i] = np.minimum(np.array(x[1::2], np.int64), 255)
    np.save(sys.argv[4] + '.tmp.npy', c); os.replace(sys.argv[4] + '.tmp.npy', sys.argv[4]); sys.exit(0)
if sys.argv[1] == 'pgsum':   # pgsum <pg dir> <dir> [base]: copies of each set in the pangenome haplotypes; 100 kb windows of each region
    D, OD = sys.argv[2:4]; BS = sys.argv[4] if len(sys.argv) > 4 else 'denom'; P = os.environ['P']; K = np.load(f'{OD}/{BS}_keys.npy'); sid = np.load(f'{OD}/{BS}_set.npy'); ST = open(f'{OD}/{BS}_sets.txt').read().split()
    B = {l.split()[0]: l.split()[1:4] for l in list(open(f'{OD}/blocks.tsv'))[1:]}; SQ = None; KI, WK = {}, {}
    for i, t in enumerate(ST):   # keys of each set: all keys inside the region (shared by overlapping regions), else the set ids
        if t not in B and t + '_1' not in B: KI[t] = np.nonzero(sid == i)[0]; continue
        g, s0, e0 = B.get(t, B.get(t + '_1')); s0, e0 = int(s0), int(e0); SQ = SQ or dict(contigs(f'{P}/ext/ref/chm13v2.0.fa'))
        k, ok = canon(SQ[g][s0:e0 + 30]); j = np.minimum(np.searchsorted(K, k), len(K) - 1); hit = ok & (K[j] == k); KI[t] = np.unique(j[hit]); WK[t] = (j[hit], np.nonzero(hit)[0] // 100000, f'{g}:{s0}-{e0}')
    hs = sorted(f.rsplit('.', 1)[0] for f in os.listdir(D) if f.endswith(('.npy', '.npz'))); one = np.zeros(len(K), np.uint16); R = {}; W = {t: [] for t in WK}
    for h in hs:
        f = f'{D}/{h}.npz' if os.path.exists(f'{D}/{h}.npz') else f'{D}/{h}.npy'; c = np.load(f)['c'] if f.endswith('.npz') else np.load(f); one += c == 1
        R[h] = {t: ((c[x] == 1).mean(), (c[x] == 0).mean(), (c[x] >= 2).mean()) for t, x in KI.items()}
        for t, (ki, wi, _) in WK.items(): W[t].append(np.bincount(wi, c[ki].astype(np.float64)) / np.maximum(np.bincount(wi), 1))
    Y = {h: 'ysinglecopy' in KI and 1 - R[h]['ysinglecopy'][1] > 0.5 for h in hs}; print(f'haplotypes {len(hs)} (with Y by ruler presence {sum(Y.values())})', flush=True)
    print('set\tkeys\thaplotypes\tone_copy_share_median\tp10\tmin\tabsent_median\tmax\tmore_median\tmax\tkeys_one_copy_in_95pct_haplotypes\tHG002_pat_one\tHG002_mat_one\tCHM13_one', flush=True)
    for t, x in KI.items():
        hh = [h for h in hs if t != 'ysinglecopy' or Y[h]]; a, z, u = (np.array([R[h][t][j] for h in hh]) for j in (0, 1, 2)); kk = (one[x] >= 0.95 * len(hh)).mean() if t != 'ysinglecopy' else float('nan')
        sp = lambda p: next((f'{R[h][t][0]:.4f}' for h in hs if p in h.lower()), '.')
        print(f'{t}\t{len(x)}\t{len(hh)}\t{np.median(a):.4f}\t{np.percentile(a, 10):.4f}\t{a.min():.4f}\t{np.median(z):.4f}\t{z.max():.4f}\t{np.median(u):.4f}\t{u.max():.4f}\t{kk:.4f}\t{sp("hg002v1.2#1")}\t{sp("hg002v1.2#2")}\t{sp("chm13")}', flush=True)
    for t, M in W.items():
        M = np.array(M); lo, hi = (M < 0.75).sum(0), (M > 1.25).sum(0)
        print(f'windows {t} {WK[t][2]}: {M.shape[1]}; haplotypes per window with mean copies < 0.75 median {np.median(lo):.0f} max {lo.max()}; > 1.25 median {np.median(hi):.0f} max {hi.max()}; haplotypes with any window < 0.75 {(M < 0.75).any(1).sum()}, > 1.25 {(M > 1.25).any(1).sum()}', flush=True)
    sys.exit(0)
if sys.argv[1] == 'gwsel':   # gwsel <dir> <pg dir> <counts dir> <rank tsv> <out tsv> [workers]: every 10 Mb window of the IMPG ranking from 1 Mb blocks; best = lowest sum of 3 ranks: conservation, uniqueness (one copy in each haplotype), copy-number CV in the samples; Y dose CV reported only
    import subprocess, collections, re; from multiprocessing import get_context; D, PG, CD, RT, OUT = sys.argv[2:7]; NW = int(sys.argv[7]) if len(sys.argv) > 7 else 8
    K = np.load(f'{D}/cand20_keys.npy'); sid = np.load(f'{D}/cand20_set.npy'); ST = open(f'{D}/cand20_sets.txt').read().split(); B = {l.split()[0]: (l.split()[1], int(l.split()[2])) for l in list(open(f'{D}/blocks.tsv'))[1:]}
    hs = sorted(x for x in os.listdir(PG) if x.endswith('.npz')); one = np.zeros(len(K), np.uint16)
    POP = os.environ.get('POP'); SP = {l.split()[1]: l.split()[6] for l in open(os.environ['P'] + '/yhaplo/ped3202.txt') if not l.startswith('FamilyID')}   # POP=1: uniqueness and conservation in the worst superpopulation, superpopulation spread of the truth males
    grp = lambda x: SP.get((re.match(r'(HG|NA)\d+', os.path.basename(x)) or re.match('', '')).group(0), '')
    ng = collections.Counter(grp(x) for x in hs); GH = sorted(g for g in ng if g); oneg = {g: np.zeros(len(K), np.uint16) for g in GH} if POP else {}
    for x in hs:   # g: one copy inside the IMPG projection of its block (pgloc)
        z = np.load(f'{PG}/{x}'); b = np.unpackbits(z['g'])[:len(K)].astype(bool) if 'g' in z.files else z['c'] == 1; one += b
        if grp(x) in oneg: oneg[grp(x)] += b
    if POP: print('haplotypes per superpopulation ' + ', '.join(f'{g} {ng[g]}' for g in GH) + f'; without one {ng[""]}', flush=True)
    TH = [0.95, 0.5, 0.6, 0.7, 0.8, 0.9]; OK = [np.logical_and.reduce([oneg[g] >= t * ng[g] for g in GH]) if POP else one >= t * len(hs) for t in TH]; ok = OK[0]; GC = np.zeros(len(K), np.uint8)   # first threshold = the rule, others = sensitivity
    for q in range(31): b = (K >> np.uint64(2 * q)) & np.uint64(3); GC += ((b == 1) | (b == 2)).astype(np.uint8)   # G+C bases per key
    QK = orderkeys(f'{D}/tiling.fa.order'); o = len(QK); LY = []
    with open(f'{D}/tiling.fa') as fa:
        for h in fa:
            if not h.startswith(('>ysinglecopy_', '>chr2depth_')): break
            q = next(fa)
            if h[1] == 'y': LY.append(q[:31])
    LN = np.flatnonzero(np.isin(QK, rowkeys(np.frombuffer(''.join(LY).encode(), np.uint8).reshape(-1, 31)))); GCR = float(np.mean([(q.count('G') + q.count('C')) / 31 for q in LY])); del LY   # GCR: mean G+C of the Y ruler keys
    so = np.argsort(QK); j = np.minimum(np.searchsorted(QK[so], K), o - 1); assert (QK[so][j] == K).all(), 'keys without a count row'; line = so[j]; del QK, so, j
    bo = np.argsort(sid, kind='stable'); bnd = np.searchsorted(sid[bo], np.arange(len(ST) + 1)); NK = np.diff(bnd); XT = [[line[bo[bnd[b]:bnd[b + 1]]][m[bo[bnd[b]:bnd[b + 1]]]] for b in range(len(ST))] for m in OK]; XF = XT[0]; NT = [np.array([len(x) for x in X]) for X in XT]; NF = NT[0]; GBT = [np.array([GC[bo[bnd[b]:bnd[b + 1]]][m[bo[bnd[b]:bnd[b + 1]]]].sum() / 31 for b in range(len(ST))]) for m in OK]; GB = GBT[0]   # G+C bases of kept keys per block, per threshold
    print(f'haplotypes {len(hs)}, keys {len(K)}, kept at thresholds {TH}: {[round(float(m.mean()), 4) for m in OK]}; count rows {o}, Y ruler rows {len(LN)}; blocks {len(ST)}', flush=True); del K, sid, one, ok, OK, line, bo, GC, oneg
    def dep(a):
        d0, s0 = a; v = np.array(subprocess.run(['zstd', '-dcq', f'{d0}/{s0}/counts.col.zst'], capture_output=True, check=True).stdout.split(), np.float64)
        if len(v) != o: return s0, None
        tm = lambda x: x[x <= 10 * (x.mean() + 1)].mean() if len(x) else np.nan; return s0, (tm(v[LN]), np.array([[tm(v[x]) for x in X] for X in XT]))
    SS = sorted(x for x in os.listdir(CD) if os.path.exists(f'{CD}/{x}/counts.col.zst')); RC = os.environ.get('REPC')   # REPC=<replicate counts dir> with REPS: replicate CV as one more rank
    with get_context('fork').Pool(NW) as pool:
        Z = dict(pool.map(dep, [(CD, x) for x in SS])); ZR = dict(pool.map(dep, [(RC, x) for x in sorted(os.listdir(RC)) if os.path.exists(f'{RC}/{x}/counts.col.zst')])) if RC else {}
    ZR = {x: v for x, v in ZR.items() if v is not None}
    for x in [x for x in Z if Z[x] is None]: print(f'ERROR {x}: counts length differs from {o} rows', flush=True); del Z[x]
    TY = os.environ.get('TRUTHY')   # TRUTHY=cn|ydose: samples = truth males of <counts dir> (window at 2 copies and one Y copy each), Element runs left out; third criterion = window copy-number CV (cn) or Y dose CV (ydose)
    sm = sorted(x for x in Z if not (TY and x.endswith('_element'))); Y = np.array([Z[x][0] for x in sm]); BT = np.array([Z[x][1] for x in sm]); BD = BT[:, 0]; BI = {B[t]: i for i, t in enumerate(ST)}; W = []
    for l in list(open(RT))[1:]:
        rk, g, s0, e0 = l.split('\t')[:4]; ib = [BI.get((g, a)) for a in range(int(s0), int(e0), 1000000)]
        if None not in ib: W.append((int(rk), g, int(s0), int(e0), float(l.split('\t')[4]), ib, l.rstrip('\n').split('\t')[5:]))
    wd = np.array([np.nansum(BD[:, w[5]] * NF[w[5]], 1) / max(NF[w[5]].sum(), 1) for w in W]).T   # samples x windows
    U = np.array([NF[w[5]].sum() / max(NK[w[5]].sum(), 1) for w in W]); SV = np.array([w[4] for w in W])
    def crit(ix):   # copy-number CV and Y dose CV per window on samples ix
        cn = 2 * wd[ix] / np.median(wd[ix], 1, keepdims=True); yd = 2 * Y[ix, None] / wd[ix]; yd = np.where((yd > 0.5) & (yd < 1.5), yd, np.nan)
        return cn.std(0) / cn.mean(0), np.nanstd(yd, 0) / np.nanmean(yd, 0)
    rank = lambda a: np.argsort(np.argsort(a, kind='stable'), kind='stable') + 1
    PL = collections.defaultdict(lambda: collections.defaultdict(int))   # projected length per block and haplotype (IMPG, blockproj)
    if os.path.exists(f'{D}/proj.tsv'):
        for l in open(f'{D}/proj.tsv'): t, ctg, a0, b0 = l.split('\t'); PL[t][ctg.rsplit('#', 1)[0]] += int(b0) - int(a0)
    def lencv(w):   # window length per haplotype: sum over its blocks; CV over haplotypes with all blocks projected
        bl = [ST[k] for k in w[5]]; hh = set.intersection(*[set(PL[t]) for t in bl]) if all(t in PL for t in bl) else set(); v = np.array([sum(PL[t][h] for t in bl) for h in hh])
        return (len(hh), float(np.median(v)) if len(v) else np.nan, float(v.std() / v.mean()) if len(v) else np.nan)
    LC = [lencv(w) for w in W]
    GW = np.array([GB[w[5]].sum() / max(NF[w[5]].sum(), 1) for w in W]); cn0 = 2 * wd / np.median(wd, 1, keepdims=True); yd0 = 2 * Y[:, None] / wd; yd0 = np.where((yd0 > 0.5) & (yd0 < 1.5), yd0, np.nan)
    REP = os.environ.get('REPS'); GI = []   # REPS=<run, male tsv>: third criterion = median over males of the within-male Y dose CV across replicate runs (same male, same Y dose), instead of the copy-number CV
    if REP:
        RM = dict(l.split()[:2] for l in open(REP)); g = collections.defaultdict(list)
        for i, s in enumerate(sm):
            if s in RM: g[RM[s]].append(i)
        GI = [ix for ix in g.values() if len(ix) >= 2]; print(f'replicates: {len(GI)} males with >= 2 runs, {sum(map(len, GI))} runs: ' + ', '.join(f'{m} {len(ix)}' for m, ix in g.items()), flush=True)
    repcv = lambda w: np.median([(2 * Y[ix, None] / w[ix]).std(0) / (2 * Y[ix, None] / w[ix]).mean(0) for ix in GI], 0)
    CS, CY = crit(np.arange(len(sm)))
    if REP and not RC: CS = repcv(wd)
    EX, XC = 0, {}   # extra ranks and their columns
    if RC:   # replicate CV on the replicate runs, same windows
        smr = sorted(ZR); YR = np.array([ZR[x][0] for x in smr]); BRT = np.array([ZR[x][1] for x in smr]); wdrt = lambda t: np.array([np.nansum(BRT[:, t, w[5]] * NT[t][w[5]], 1) / max(NT[t][w[5]].sum(), 1) for w in W]).T
        RM = dict(l.split()[:2] for l in open(REP)); gr = collections.defaultdict(list)
        for i, s0 in enumerate(smr):
            if s0 in RM: gr[RM[s0]].append(i)
        GR = [ix for ix in gr.values() if len(ix) >= 2]; rcv = lambda wdr: np.median([(2 * YR[ix, None] / wdr[ix]).std(0) / (2 * YR[ix, None] / wdr[ix]).mean(0) for ix in GR], 0); XC['rep_cv'] = rcv(wdrt(0))
        print(f'replicate rank: {len(GR)} males with >= 2 runs, {sum(map(len, GR))} runs', flush=True)
    if POP:   # superpopulation spread of the truth males: max - min of the median window copy number
        gs = np.array([grp(x) for x in sm]); G2 = [g for g in GH if (gs == g).sum() >= 5]; psp = lambda cn: np.max([np.median(cn[gs == g], 0) for g in G2], 0) - np.min([np.median(cn[gs == g], 0) for g in G2], 0); XC['pop_spread'] = psp(cn0)
        print('truth males per superpopulation ' + ', '.join(f'{g} {int((gs == g).sum())}' for g in G2), flush=True)
        TB = collections.defaultdict(dict); HT = set()   # IMPG tiles: bp in SVs per block and haplotype
        for fn in os.listdir(os.environ['TILES']):
            for l in list(open(f'{os.environ["TILES"]}/{fn}'))[1:]:
                x = l.rstrip('\n').split('\t')
                if len(x) >= 6: TB[(x[0], int(x[1]))][x[2]] = (float(x[3]), int(x[5])); HT.add(x[2])
        HT = sorted(h for h in HT if grp(h)); hg = np.array([grp(h) for h in HT])
        def svw(w):   # worst superpopulation median of the window bp in SVs (indels >= 50 + uncovered)
            v = np.array([sum(TB[B[ST[k]]].get(h, (0, 0))[1] + (1 - TB[B[ST[k]]].get(h, (0, 0))[0]) * 1000000 for k in w[5]) for h in HT]); return max(np.median(v[hg == g]) for g in GH if (hg == g).any())
        XC['sv_worst_pop'] = np.array([svw(w) for w in W]); SV = XC['sv_worst_pop'].copy()
    for k in ('rep_cv', 'pop_spread'):
        if k in XC: EX = EX + rank(XC[k])
    if TY == 'ydose': CS = CY.copy()
    GM = os.environ.get('GCMATCH')   # GCMATCH=1: fourth rank = |G+C of the kept window keys - G+C of the Y ruler| (shared G+C bias)
    RS = rank(SV) + rank(-U) + rank(CS) + (rank(np.abs(GW - GCR)) if GM else 0) + EX; best = int(np.argmin(RS)); print(f'ruler G+C {GCR:.4f}; criteria: conservation, uniqueness, third CV' + (', |window G+C - ruler G+C|' if GM else ''), flush=True)
    M = np.stack([SV, -U, CS], 1); PF = np.array([not ((M <= M[i]).all(1) & (M < M[i]).any(1)).any() for i in range(len(W))])
    for t in range(1, len(TH)):   # sensitivity: same rule with another one-copy threshold
        wdt = np.array([np.nansum(BT[:, t, w[5]] * NT[t][w[5]], 1) / max(NT[t][w[5]].sum(), 1) for w in W]).T; cn = 2 * wdt / np.median(wdt, 1, keepdims=True); yd = 2 * Y[:, None] / wdt; yd = np.where((yd > 0.5) & (yd < 1.5), yd, np.nan)
        Ut = np.array([NT[t][w[5]].sum() / max(NK[w[5]].sum(), 1) for w in W]); GWt = np.array([GBT[t][w[5]].sum() / max(NT[t][w[5]].sum(), 1) for w in W])   # every key-dependent rank at threshold t
        EXt = (rank(rcv(wdrt(t))) if RC else 0) + (rank(psp(cn)) if POP else 0)
        rs = rank(SV) + rank(-Ut) + rank(repcv(wdt) if REP and not RC else np.nanstd(yd, 0) / np.nanmean(yd, 0) if TY == 'ydose' else cn.std(0) / cn.mean(0)) + (rank(np.abs(GWt - GCR)) if GM else 0) + EXt; bt = int(np.argmin(rs))
        print(f'threshold {TH[t]}: selected {W[bt][1]}:{W[bt][2]}-{W[bt][3]} (kept {Ut[bt]:.4f}); the 0.95 choice {W[best][1]}:{W[best][2]}-{W[best][3]} is at position {int((rs < rs[best]).sum()) + 1}; top 5 {[f"{W[q][1]}:{W[q][2] // 1000000}" for q in np.argsort(rs, kind="stable")[:5]]}', flush=True)
    with open(OUT, 'w') as w:
        w.write(('impg_rank\tchrom\tstart\tend\tsv_bp_median\tsv_bp_p90\tsv_bp_max\thaplotypes_full_cover\tidentity_median\tproj_haplotypes\tproj_length_median\tproj_length_cv\tkeys\tkeys_kept\tkept_share\tgc_kept\tcn_median\tcn_cv\tydose_median\tydose_cv\trank_conservation\trank_uniqueness\trank_cn_cv\trank_sum\tpareto\tselected' + ''.join(f'\t{k}' for k in XC) + '\n').replace('rank_cn_cv', 'rank_ydose_cv' if TY == 'ydose' else 'rank_cn_cv').replace('cn_cv', 'rep_cv' if REP and not RC else 'cn_cv'))
        r1, r2, rcs = rank(SV), rank(-U), rank(CS)
        for i, x in enumerate(W): w.write(f'{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]:.0f}\t' + '\t'.join(x[6]) + f'\t{LC[i][0]}\t{LC[i][1]:.0f}\t{LC[i][2]:.4f}\t{NK[x[5]].sum()}\t{NF[x[5]].sum()}\t{U[i]:.4f}\t{GW[i]:.4f}\t{np.median(cn0[:, i]):.4f}\t{CS[i]:.4f}\t{np.nanmedian(yd0[:, i]):.4f}\t{CY[i]:.4f}\t{r1[i]}\t{r2[i]}\t{rcs[i]}\t{RS[i]}\t{int(PF[i])}\t{int(i == best)}' + ''.join(f'\t{XC[k][i]:.5g}' for k in XC) + '\n')
    print(f'samples {len(sm)}; windows {len(W)}; Pareto front {PF.sum()}', flush=True)
    for i in list(np.argsort(RS, kind='stable')[:10]) + [i for i, x in enumerate(W) if (x[1], x[2]) in (('chr2', 50000000), ('chr7', 34000000))]:
        x = W[i]; print(f'{x[1]}:{x[2]}-{x[3]} impg rank {x[0]} sv {x[4]:.0f}; projected length median {LC[i][1]:.0f} CV {LC[i][2]:.4f} ({LC[i][0]} haplotypes); kept {U[i]:.4f}; G+C {GW[i]:.4f}; cn CV {CS[i]:.4f}; Y dose CV {CY[i]:.4f}; rank sum {RS[i]} (position {int((RS < RS[i]).sum()) + 1}); Pareto {int(PF[i])}', flush=True)
    sys.exit(0)
if sys.argv[1] == 'winkeys':   # winkeys <gw dir> <chrom> <start Mb> <end Mb> <out npy>: keys of the window blocks at one copy inside the projection in >= 95% of haplotypes (Y dose denominator)
    D, g, a, b, OF = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]; K = np.load(f'{D}/cand20_keys.npy'); sid = np.load(f'{D}/cand20_set.npy'); ST = open(f'{D}/cand20_sets.txt').read().split()
    j = np.flatnonzero(np.isin(sid, [ST.index(f'{g}_{i}') for i in range(a, b)])); n = np.zeros(len(j), np.uint16); hs = sorted(os.listdir(f'{D}/pgloc'))
    for h in hs: n += np.unpackbits(np.load(f'{D}/pgloc/{h}')['g'])[j].astype(bool)
    k = K[j[n >= 0.95 * len(hs)]]; np.save(OF, k); print(f'{g}:{a}-{b} Mb: keys {len(j)}, kept {len(k)} (one copy inside the projection in >= 95% of {len(hs)} haplotypes); written {OF}', flush=True); sys.exit(0)
if sys.argv[1] == 'winfa':   # winfa <tiling dir> <pg dir> <set> <out fa>: Y ruler + the set k-mers at one copy in >= 95% of haplotypes as class chr2depth (depth reference); <out fa>.rows.npy = their rows in tiling.fa.order
    D, PG, T, OF = sys.argv[2:6]; P = os.environ['P']; K = np.load(f'{D}/cand20_keys.npy'); B = {l.split()[0]: l.split()[1:4] for l in list(open(f'{D}/blocks.tsv'))[1:]}
    g, s0, e0 = B[T]; k, v = canon(dict(contigs(f'{P}/ext/ref/chm13v2.0.fa'))[g][int(s0):int(e0) + 30]); j = np.minimum(np.searchsorted(K, k), len(K) - 1); j = np.unique(j[v & (K[j] == k)])
    hs = sorted(x for x in os.listdir(PG) if x.endswith('.npz')); one = np.zeros(len(j), np.uint16)
    for x in hs: one += np.load(f'{PG}/{x}')['c'][j] == 1
    kk = K[j[one >= 0.95 * len(hs)]]; QK = orderkeys(f'{D}/tiling.fa.order'); so = np.argsort(QK); q = np.minimum(np.searchsorted(QK[so], kk), len(QK) - 1); assert (QK[so][q] == kk).all(), 'kept keys without a count row'; np.save(OF + '.rows.npy', so[q])
    bb = np.array([(kk >> np.uint64(2 * (30 - i))) & np.uint64(3) for i in range(31)]).T.astype(np.uint8); S = np.frombuffer(b'ACGT', np.uint8)[bb]; ny = 0
    with open(OF, 'w') as w, open(f'{D}/tiling.fa') as fa:
        for h in fa:
            if not h.startswith(('>ysinglecopy_', '>chr2depth_')): break
            q = next(fa)
            if h.startswith('>ysinglecopy_'): w.write(h + q); ny += 1
        w.writelines(f'>chr2depth_{i}\n{S[i].tobytes().decode()}\n' for i in range(len(S)))
    print(f'{T} {g}:{s0}-{e0}: keys {len(j)}, kept {len(kk)} (one copy in >= 95% of {len(hs)} haplotypes); written {OF}: ysinglecopy {ny}, chr2depth = {T} kept {len(S)}', flush=True); sys.exit(0)
if sys.argv[1] == 'c2filt':   # c2filt <tiling dir> <counts dir> <pg dir> <keys npy> <out tsv> [workers]: Y dose with all chr2depth k-mers vs those at one copy in >= 95% of haplotypes
    import subprocess; from multiprocessing import get_context; D, CD, PG, KF, OUT = sys.argv[2:7]; NW = int(sys.argv[7]) if len(sys.argv) > 7 else 8
    K = np.load(KF); hs = sorted(x for x in os.listdir(PG) if x.endswith('.npz')); one = np.zeros(len(K), np.uint16)
    for x in hs: one += np.load(f'{PG}/{x}')['c'] == 1
    QK = orderkeys(f'{D}/tiling.fa.order'); o = len(QK); S = {'y': [], 'c': []}
    with open(f'{D}/tiling.fa') as fa:   # first lines: ysinglecopy and chr2depth
        for h in fa:
            if not h.startswith(('>ysinglecopy_', '>chr2depth_')): break
            S[h[1]].append(next(fa)[:31])
    ky, kc = (rowkeys(np.frombuffer(''.join(S[t]).encode(), np.uint8).reshape(-1, 31)) for t in 'yc'); del S
    so = np.argsort(QK); qs = QK[so]; R = {}
    for t, k in (('y', ky), ('c', kc)): j = np.minimum(np.searchsorted(qs, k), o - 1); assert (qs[j] == k).all(), f'{t} keys without a count row'; R[t] = so[j]
    j = np.minimum(np.searchsorted(K, kc), len(K) - 1); fk = K[j] == kc; keep = fk & (one[j] >= 0.95 * len(hs)); del QK, so, qs
    print(f'haplotypes {len(hs)}; chr2depth k-mers {len(kc)}, in the scanned keys {fk.mean():.4f}, at one copy in >= 95% of haplotypes {keep.mean():.4f}; ruler k-mers {len(ky)}', flush=True)
    RC, RK = R['c'], R['c'][keep]
    def dep(a):
        d0, s0 = a; v = np.array(subprocess.run(['zstd', '-dcq', f'{d0}/{s0}/counts.col.zst'], capture_output=True, check=True).stdout.split(), np.float64)
        if len(v) != o: return s0, None
        tm = lambda x: x[x <= 10 * (x.mean() + 1)].mean(); return s0, (tm(v[R['y']]), tm(v[RC]), tm(v[RK]))
    SS = sorted(x for x in os.listdir(CD) if os.path.exists(f'{CD}/{x}/counts.col.zst'))
    with get_context('fork').Pool(NW) as pool: Z = dict(pool.map(dep, [(CD, x) for x in SS]))
    for x in [x for x in Z if Z[x] is None]: print(f'ERROR {x}: counts length differs from {o} rows', flush=True); del Z[x]
    with open(OUT, 'w') as w:
        w.write('sample\truler_depth\tchr2_depth_all\tchr2_depth_kept\tydose_all\tydose_kept\n')
        for x, (y, a, b) in sorted(Z.items()): w.write(f'{x}\t{y:.4f}\t{a:.4f}\t{b:.4f}\t{2 * y / a:.4f}\t{2 * y / b:.4f}\n')
    for i, nm in ((1, 'all'), (2, 'kept')):
        yd = np.array([2 * z[0] / z[i] for z in Z.values()]); yd = yd[(yd > 0.5) & (yd < 1.5)]; print(f'Y dose chr2 50-60 {nm} k-mers: n {len(yd)}, median {np.median(yd):.4f}, CV {yd.std() / yd.mean():.4f}, p10-p90 {np.percentile(yd, 10):.4f}-{np.percentile(yd, 90):.4f}', flush=True)
    sys.exit(0)
if sys.argv[1] == 'tilesel':   # tilesel <dir> <pg dir> <counts dir> <out tsv> [workers]: denominator rule (user 2026-10-30): IMPG top quarter, then lowest CV of the copy number in the CRAMs, k-mers at one copy in >= 95% of haplotypes
    import subprocess; from multiprocessing import get_context; D, PG, CD, OUT = sys.argv[2:6]; NW = int(sys.argv[6]) if len(sys.argv) > 6 else 8; P = os.environ['P']
    K = np.load(f'{D}/cand20_keys.npy'); sid = np.load(f'{D}/cand20_set.npy'); ST = open(f'{D}/cand20_sets.txt').read().split(); B = {l.split()[0]: l.split()[1:4] for l in list(open(f'{D}/blocks.tsv'))[1:]}
    RK = {l.split()[0]: int(l.split()[4]) for l in list(open(f'{D}/tiling_grch38.hdr.tsv'))[1:]}; hs = sorted(x for x in os.listdir(PG) if x.endswith('.npz')); one = np.zeros(len(K), np.uint16)
    for x in hs: one += np.load(f'{PG}/{x}')['c'] == 1
    ok = one >= 0.95 * len(hs); print(f'haplotypes {len(hs)}, keys {len(K)}, at one copy in >= 95% of them {ok.mean():.4f}', flush=True)
    QK = orderkeys(f'{D}/tiling.fa.order'); o = len(QK); ck = rowkeys
    LY = []
    with open(f'{D}/tiling.fa') as fa:   # first lines: ysinglecopy and chr2depth (denfa order)
        for h in fa:
            if not h.startswith(('>ysinglecopy_', '>chr2depth_')): break
            q = next(fa)
            if h[1] == 'y': LY.append(q[:31])
    LN = np.isin(QK, ck(np.frombuffer(''.join(LY).encode(), np.uint8).reshape(-1, 31))); del LY
    so = np.argsort(QK); j = np.minimum(np.searchsorted(QK[so], K), o - 1); hit = QK[so][j] == K; line = np.where(hit, so[j], -1); print(f'count rows {o}, Y ruler rows {LN.sum()}, keys without a count row {(~hit).sum()}', flush=True)
    SQ = dict(contigs(f'{P}/ext/ref/chm13v2.0.fa')); XA, XF = {}, {}
    for t in ST:   # keys of each region from its sequence (overlapping regions share keys)
        g, s0, e0 = B[t]; k, v = canon(SQ[g][int(s0):int(e0) + 30]); j = np.minimum(np.searchsorted(K, k), len(K) - 1); j = np.unique(j[v & (K[j] == k)]); j = j[line[j] >= 0]; XA[t], XF[t] = line[j], line[j[ok[j]]]
    del SQ, K, sid, one, line, QK, so
    def dep(s0):   # trimmed mean depth per region (all keys, kept keys) and of the Y ruler
        v = np.array(subprocess.run(['zstd', '-dcq', f'{CD}/{s0}/counts.col.zst'], capture_output=True, check=True).stdout.split(), np.float64)
        if len(v) != o: return s0, None
        tm = lambda x: x[x <= 10 * (x.mean() + 1)].mean(); return s0, (tm(v[LN]), {t: (tm(v[XA[t]]), tm(v[XF[t]])) for t in ST})
    SS = sorted(x for x in os.listdir(CD) if os.path.exists(f'{CD}/{x}/counts.col.zst'))
    with get_context('fork').Pool(NW) as pool: R = dict(pool.map(dep, SS))
    for x in [x for x in R if R[x] is None]: print(f'ERROR {x}: counts length differs from {o} lines', flush=True); del R[x]
    W = [t for t in ST if t in RK]; top = set(sorted(W, key=RK.get)[:len(W) // 4]); rows = []
    for t in ST:
        r = {}
        for a, nm in ((0, 'all'), (1, 'kept')):
            cn = np.array([2 * R[x][1][t][a] / np.median([R[x][1][w][a] for w in W]) for x in R]); yd = np.array([2 * R[x][0] / R[x][1][t][a] for x in R]); yd = yd[(yd > 0.5) & (yd < 1.5)]
            r[nm] = (np.median(cn), cn.std() / cn.mean(), np.median(yd), yd.std() / yd.mean())
        rows.append((t, B[t], RK.get(t), len(XA[t]), len(XF[t]), r))
    ck = lambda q: sorted(q, key=lambda z: z[5]['kept'][1]); pick = ck([z for z in rows if z[0] in top])[0][0]; cr = {z[0]: i + 1 for i, z in enumerate(ck([z for z in rows if z[0] in RK or z[0].startswith('ref')]))}
    with open(OUT, 'w') as w:
        w.write('set\tchrom\tstart\tend\tconservation_rank\tconservation_top_quarter\tkeys\tkeys_kept\tcn_median_all\tcn_cv_all\tydose_median_all\tydose_cv_all\tcn_median_kept\tcn_cv_kept\tydose_median_kept\tydose_cv_kept\tcn_cv_kept_rank\tselected\n')
        for t, b, rk, na, nf, r in rows: w.write(f'{t}\t' + '\t'.join(b) + f'\t{rk if rk else "."}\t{int(t in top)}\t{na}\t{nf}\t' + '\t'.join(f'{x:.4f}' for x in r['all'] + r['kept']) + f'\t{cr[t]}\t{int(t == pick)}\n')
    print(f'samples {len(R)}; windows {len(W)}, conservation top quarter {len(top)}; selected {pick} {B[pick]} (copy-number CV {dict((z[0], z[5]) for z in rows)[pick]["kept"][1]:.4f})', flush=True)
    for t in [pick] + [x for x in ST if x.startswith('ref')]: z = next(q for q in rows if q[0] == t); print(f'{t} {z[1]} conservation rank {z[2]} top quarter {t in top}; keys kept {z[4]} of {z[3]}; copy-number CV kept {z[5]["kept"][1]:.4f} (rank {cr[t]} of {len(cr)}), all {z[5]["all"][1]:.4f}; Y dose CV kept {z[5]["kept"][3]:.4f}, all {z[5]["all"][3]:.4f}', flush=True)
    sys.exit(0)
if sys.argv[1] == 'denfa':   # denfa <dir> <base> <marker fa with ysinglecopy and chr2depth> <out fa>: one class per set of <base>_keys.npy (yruler needs the two classes)
    D, BS, MF, OF = sys.argv[2:6]; K = np.load(f'{D}/{BS}_keys.npy'); sid = np.load(f'{D}/{BS}_set.npy'); ST = open(f'{D}/{BS}_sets.txt').read().split(); n = 0
    with open(OF, 'w') as w:
        L = open(MF).read().split('\n'); w.writelines(f'{h}\n{q}\n' for h, q in zip(L[0::2], L[1::2]) if h.startswith(('>ysinglecopy_', '>chr2depth_'))); del L
        for i, t in enumerate(ST):
            k = K[sid == i]; bb = np.array([(k >> np.uint64(2 * (30 - j))) & np.uint64(3) for j in range(31)]).T.astype(np.uint8); S = np.frombuffer(b'ACGT', np.uint8)[bb]
            c = os.environ.get('ONECLASS', t.replace('_', '-')); w.writelines(f'>{c}_{j}\n{S[j].tobytes().decode()}\n' for j in range(len(S))); n += len(S)   # ONECLASS: one class name (yruler takes <= 16)
    print(f'written {OF}: {n} window k-mers in {len(ST)} sets + ysinglecopy and chr2depth from {MF}', flush=True); sys.exit(0)
if sys.argv[1] == 'den':   # den <counts dir> <marker fa> <out tsv>: depth per set and Y dose per denominator
    import subprocess; CD, MF, OUT = sys.argv[2:5]; L = open(MF).read().split('\n'); SET = {q: h[1:].rsplit('_', 1)[0] for h, q in zip(L[0::2], L[1::2]) if h}; del L   # k-mer -> set
    nm = np.array([SET.get(l.rstrip('\n'), 'missing') for l in open(MF + '.order')]); print(f'markers {len(nm)}, not in {MF}: {(nm == "missing").sum()}', flush=True)
    ST = sorted(set(nm.tolist())); R = []
    for s0 in sorted(os.listdir(CD)):
        f = f'{CD}/{s0}/counts.col.zst'
        if not os.path.exists(f): continue
        v = np.array(subprocess.run(['zstd', '-dcq', f], capture_output=True, check=True).stdout.split(), np.float64)
        if len(v) != len(nm): print(f'ERROR {s0}: {len(v)} counts, {len(nm)} markers', flush=True); continue
        d = {}
        for t in ST: x = v[nm == t]; d[t] = x[x <= 10 * (x.mean() + 1)].mean()   # trimmed mean depth
        R.append((s0, d))
    D = [t for t in ST if t != 'ysinglecopy']
    with open(OUT, 'w') as o:
        o.write('sample\t' + '\t'.join(f'{t}_depth' for t in ST) + '\t' + '\t'.join(f'ydose_{t}' for t in D) + '\n')
        for s0, d in R: o.write(s0 + '\t' + '\t'.join(f'{d[t]:.3f}' for t in ST) + '\t' + '\t'.join(f'{2 * d["ysinglecopy"] / d[t]:.4f}' for t in D) + '\n')
    print(f'samples {len(R)}; denominator: Y dose median, CV, p10-p90; depth / chr2 median, CV', flush=True)
    for t in D:
        y = np.array([2 * d['ysinglecopy'] / d[t] for s0, d in R]); y = y[(y > 0.5) & (y < 1.5)]; q = np.array([d[t] / d['chr2depth'] for s0, d in R])
        print(f'{t}\tn {len(y)}\tydose {np.median(y):.4f}\tCV {y.std() / y.mean():.4f}\tp10-p90 {np.percentile(y, 10):.4f}-{np.percentile(y, 90):.4f}\tvs chr2 {np.median(q):.4f}\tCV {q.std() / q.mean():.4f}', flush=True)
    sys.exit(0)
if sys.argv[1] == 'ypq':   # [UM=markers_v6] ypq <model dir> <out npy> [workers]: X-degenerate block (XDR1-XDR8 of the assembly annotation) of each Y ruler key; label = the block of >= 90% of the assemblies with the key (as the ruler presence), else 0; Yp = XDR1-2, Yq = XDR3-8; ruler order of ysc_dyz3
    P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_v6")}'; MD, OUT = sys.argv[2:4]; NW = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    Y0 = np.load(f'{MD}/model.npz')['ysc_dyz3']; K = np.load(f'{O}/u_keys.npy')[Y0]; so = np.argsort(K); K = K[so]
    def block(x):   # votes: row 0 = on an XDR contig outside the blocks, rows 1-8 = XDR1-XDR8
        a, f = x; b = f'{P}/ext/bed/{a}.t2tv2.chrY-regions.hg38-patched.bed'
        if not os.path.exists(b): return a, None, 'no annotation'
        X = [(int(y[3][3:]), y[0], int(y[1]), int(y[2])) for y in (l.rstrip('\n').split('\t') for l in open(b) if l[0] != '#') if len(y) > 6 and y[6] == 'LABEL' and y[3][:3] == 'XDR']
        if sorted(x[0] for x in X) != list(range(1, 9)): return a, None, f'XDR blocks {sorted(x[0] for x in X)}'
        V = np.zeros((9, len(K)), np.uint16)
        for ctg, q in contigs(f):
            B = [x for x in X if x[1] == ctg]
            if not B: continue
            k, ok = canon(q); pos = np.nonzero(ok)[0]; k = k[ok]; i = np.minimum(np.searchsorted(K, k), len(K) - 1); h = K[i] == k; z = np.zeros(int(h.sum()), np.int64)
            for n, g, s0, e0 in B: z[(pos[h] >= s0) & (pos[h] < e0)] = n
            np.add.at(V, (z, i[h]), 1)
        return a, V, ' '.join(f'XDR{n} {int(V[n].sum())}' for n in range(1, 9)) + f' outside {int(V[0].sum())}'
    A = [l.split()[:2] for l in open(f'{O}/assemblies.tsv')]; T = np.zeros((9, len(K)), np.int64); nv = 0
    print(f'ruler {len(K)} keys, assemblies {len(A)}', flush=True)
    with Pool(NW) as p:
        for a, V, msg in p.imap_unordered(block, A):
            print(a, msg, flush=True)
            if V is not None: T += np.minimum(V, 1); nv += 1   # one vote per assembly and block
    n = T.sum(0); b = T[1:].argmax(0) + 1; lab = np.where((n > 0) & (T[b, np.arange(len(K))] >= 0.9 * n), b, 0).astype(np.uint8); L = np.zeros(len(K), np.uint8); L[so] = lab
    np.save(OUT, L); print(f'voting assemblies {nv}; keys per block ' + ' '.join(f'XDR{x} {int((L == x).sum())}' for x in range(1, 9)) + f'; unlabelled {int((L == 0).sum())} (no vote {int((n == 0).sum())}); votes per key median {int(np.median(n))}, min {int(n.min())}; written {OUT}', flush=True)
    sys.exit(0)
if sys.argv[1] == 'xpg':   # xpg prep <ykmers dir> <model dir> | xpg sum <ykmers dir> <scan dir> <out npz>: pangenome filter of the X keys (prep: X keys + 200,000 Y ruler keys -> <ykmers>/xpg_keys.npy for pgscan); a haplotype carries Y if it holds >= 50% of the ruler keys; X-specific = one copy in >= 90% of the haplotypes without Y and absent from >= 90% of those with Y; PAR (keys shared by X and Y) = one copy in >= 90% of the haplotypes without Y and at most one in >= 90% of those with Y
    Y = sys.argv[3]; X = np.load(f'{Y}/xkeys.npy')
    if sys.argv[2] == 'prep':
        R = np.load(f'{os.environ.get("UMD", os.environ["P"] + "/markers_v6")}/u_keys.npy')[np.load(f'{sys.argv[4]}/model.npz')['ysc_dyz3']]; R = np.random.default_rng(1).choice(R, 200000, replace=False)
        K = np.unique(np.concatenate([X, R])); np.save(f'{Y}/xpg_keys.npy', K); np.save(f'{Y}/xpg_ruler.npy', np.sort(R)); print(f'X keys {len(X)}, ruler keys {len(R)}, scan keys {len(K)}', flush=True); sys.exit(0)
    SD, OF = sys.argv[4:6]; K = np.load(f'{Y}/xpg_keys.npy'); ix = np.searchsorted(K, X); ir = np.searchsorted(K, np.load(f'{Y}/xpg_ruler.npy')); hs = sorted(x for x in os.listdir(SD) if x.endswith('.npz'))
    n1x, n0y, n1y, ny = np.zeros(len(X), np.uint16), np.zeros(len(X), np.uint16), np.zeros(len(X), np.uint16), 0
    for h in hs:
        c = np.load(f'{SD}/{h}')['c']; y = (c[ir] > 0).mean() >= 0.5; cx = c[ix]; ny += y
        if y: n0y += cx == 0; n1y += cx <= 1
        else: n1x += cx == 1
    nx = len(hs) - ny; sh = np.load(f'{Y}/xshared.npy'); xs = (n1x >= 0.9 * nx) & (n0y >= 0.9 * ny); par = (n1x >= 0.9 * nx) & (n1y >= 0.9 * ny)
    np.savez_compressed(OF, x=xs, par=par); print(f'haplotypes {len(hs)}: with Y {ny}, without Y {nx}; X-specific keys kept {int((xs & ~sh).sum())} of {int((~sh).sum())}; shared X-Y keys kept {int((par & sh).sum())} of {int(sh.sum())}; written {OF}', flush=True); sys.exit(0)
if sys.argv[1] == 'xpk':   # xpk cand <haplotype list> <out dir> [workers] | xpk sum <out dir> <scan dir>: X-specific and PAR keys from high-quality HPRC r2 haplotypes (IMPG projection of GRCh38 chrX, PAR1, PAR2); X keys 1 in 32 by fixed hash; candidates in >= 90% of the haplotypes; sum: copies in all scanned haplotypes as xpg (X: one copy in >= 90% without Y, absent in >= 90% with Y; PAR: one copy in >= 90% without Y, <= 1 in >= 90% with Y)
    import subprocess; P = os.environ['P']; IM = 'impg'; PJ = f'{P}/markers_v5/denom_proj'; AGC = os.environ['SCRATCH'] + '/sylvia/hprc/HPRC_r2_assemblies_0.6.1.agc'
    OD = sys.argv[4] if sys.argv[2] == 'cand' else sys.argv[3]; os.makedirs(OD, exist_ok=True); RG = dict(x.split('=') for x in os.environ['XPK_RG'].split(',')) if os.environ.get('XPK_RG') else {'x': 'chrX:1-156030895', 'par1': 'chrX:10001-2781479', 'par2': 'chrX:155701383-156030895'}   # GRCh38 chrX and PARs; XPK_RG=win=chr11:26000001-36000000 for a denominator window (keys 1 in 20)
    pan = lambda h: h.split('_')[0] + '#' + {'pat': '1', 'mat': '2', 'hap1': '1', 'hap2': '2'}[h.split('_')[1]]   # HPRC r2 set name -> PanSN prefix
    if sys.argv[2] == 'cand':
        HQ = [l.split()[0] for l in open(sys.argv[3]) if l.strip()]; NW = int(sys.argv[5]) if len(sys.argv) > 5 else 8; PR = {}
        for r, g in RG.items():
            out = subprocess.run([IM, 'query', '--alignment-list', f'{PJ}/hprc_r2_paf.list', '-i', f'{PJ}/hprc.impg', '-r', f'GRCh38#0#{g}', '-d', '10000', '-o', 'bed', '-t', '8'], check=True, capture_output=True).stdout.decode().splitlines()
            for x in out:
                c, a0, b0 = x.split('\t')[:3]
                if not c.startswith('GRCh38#'): PR.setdefault((r, c.rsplit('#', 1)[0]), []).append((c.rsplit('#', 1)[1], int(a0), int(b0)))
            print(f'{r} {g}: haplotypes with pieces {len({h for rr, h in PR if rr == r})}', flush=True)
        def keys(h):   # canonical keys of the projected pieces of one haplotype
            fa = subprocess.run(['agc', 'getset', AGC, h], check=True, capture_output=True).stdout; p = pan(h); Q = {}
            for n, q in contigs_bytes(fa): Q[n.rsplit('#', 1)[-1]] = q
            out = {}
            for r in RG:
                K = []
                for c, a0, b0 in PR.get((r, p), []):
                    if c not in Q: continue
                    for c0 in range(a0, b0, 50000000):   # 50 Mb chunks overlapping by 30 bp: memory bound
                        k, ok = canon(Q[c][c0:min(b0, c0 + 50000030)]); pos = np.arange(len(k)) + c0; k = k[ok]; pos = pos[ok]
                        if r == 'x':   # X keys outside the projected PARs, 1 in 32
                            inpar = np.zeros(len(k), bool)
                            for rr in ('par1', 'par2'):
                                for c2, a2, b2 in PR.get((rr, p), []):
                                    if c2 == c: inpar |= (pos >= a2) & (pos < b2)
                            k = k[~inpar]; k = k[((k * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(32)) % np.uint64(32) == 0]
                        if r == 'win': k = k[((k * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(32)) % np.uint64(20) == 0]
                        K.append(np.unique(k))
                out[r] = np.unique(np.concatenate(K)) if K else np.zeros(0, np.uint64)
            return h, out
        C = {r: [] for r in RG}
        with Pool(NW) as pool:
            for h, out in pool.imap_unordered(keys, HQ):
                for r in RG: C[r].append(out[r])
                print(h, ' '.join(f'{r} {len(out[r])}' for r in RG), flush=True)
        F = {}
        for r in RG:
            u, n = np.unique(np.concatenate(C[r]), return_counts=True); F[r] = u[n >= 0.9 * len(HQ)]; print(f'{r}: distinct keys {len(u)}, in >= 90% of {len(HQ)} haplotypes {len(F[r])}', flush=True)
        if 'par2' in F: F['par2'] = np.setdiff1d(F['par2'], F['par1']); F['x'] = np.setdiff1d(F['x'], np.union1d(F['par1'], F['par2']))
        np.savez(f'{OD}/xpk_cand.npz', **F); R = np.load(f'{P}/ykmers/xpg_ruler.npy'); np.save(f'{OD}/xpk_scan_keys.npy', np.unique(np.concatenate([F[r] for r in RG] + [R]))); print(f'written {OD}/xpk_cand.npz and xpk_scan_keys.npy', flush=True); sys.exit(0)
    SD = sys.argv[4]; F = dict(np.load(f'{OD}/xpk_cand.npz')); K = np.load(f'{OD}/xpk_scan_keys.npy'); ir = np.searchsorted(K, np.load(f'{P}/ykmers/xpg_ruler.npy')); hs = sorted(x for x in os.listdir(SD) if x.endswith('.npz')); ny = 0
    A = {r: (np.searchsorted(K, F[r]), np.zeros(len(F[r]), np.uint16), np.zeros(len(F[r]), np.uint16), np.zeros(len(F[r]), np.uint16)) for r in RG}   # index, one copy without Y, absent with Y, <= 1 with Y
    SPM = {l.split()[1]: l.split()[6] for l in open(f'{P}/yhaplo/ped3202.txt') if not l.startswith('FamilyID')}; import re; gp = lambda x: SPM.get((re.match(r'(HG|NA)\d+', x) or re.match('', '')).group(0), '')
    WG = {g: [np.zeros(len(F['win']), np.uint16), 0] for g in sorted({gp(h) for h in hs} - {''})} if 'win' in F else {}   # win: one copy per superpopulation
    for h in hs:
        c = np.load(f'{SD}/{h}')['c']; y = (c[ir] > 0).mean() >= 0.5; ny += y
        for r, (ix, n1, n0, nl) in A.items():
            cx = c[ix]
            if r == 'win':
                if gp(h) in WG: WG[gp(h)][0] += cx == 1; WG[gp(h)][1] += 1
                continue
            if y: n0 += cx == 0; nl += cx <= 1
            else: n1 += cx == 1
    nx = len(hs) - ny; O = {}
    for r, (ix, n1, n0, nl) in A.items():
        O[r] = F[r][np.logical_and.reduce([v >= 0.95 * n for v, n in WG.values()])] if r == 'win' else F[r][(n1 >= 0.9 * nx) & ((n0 if r == 'x' else nl) >= 0.9 * ny)]; print(f'{r}: candidates {len(F[r])}, kept {len(O[r])}', flush=True)   # win: one copy in >= 95% of the haplotypes of every superpopulation
    np.savez(f'{OD}/xpk_final.npz', **O); print(f'haplotypes {len(hs)} (with Y {ny}); written {OD}/xpk_final.npz', flush=True); sys.exit(0)
if sys.argv[1] == 'd17':   # d17 <out tsv> [workers]: DYZ17 arrays located with seeds from the assemblies, not CHM13: v6 DYZ17 candidates (Y keys in >= 15% of the DYZ17 regions) with >= 95% of their Y copies inside the region and found in >= 50% of the assemblies; same rule as um_newtruth.py (hits merged within 100 bp, clusters within 5 kb, cluster with most assigned bases); compared with regions.tsv
    P = os.environ['P']; O = f'{P}/markers_v6'; OUT = sys.argv[2]; NW = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    U = np.load(f'{O}/u_keys.npy'); ci = np.load(f'{O}/cand_v6.npz')['dyz17']; A = [l.split()[:2] for l in open(f'{O}/assemblies.tsv')]; ni = np.zeros(len(ci)); no = np.zeros(len(ci)); pres = np.zeros(len(ci))
    for a0, _ in A:
        z = np.load(f'{O}/copies_v6/{a0}.npz'); ni += z['dyz17_in']; no += z['dyz17_out']; pres += z['dyz17_in'] > 0
    keep = (ni / np.maximum(ni + no, 1) >= 0.95) & (pres >= 0.5 * len(A)); K = np.sort(U[ci[keep]]); print(f'DYZ17 candidates {len(ci)}, seeds {len(K)} (>= 95% of copies inside, in >= 50% of {len(A)} assemblies)', flush=True)
    OLD = {l.split()[0]: l.split()[2:5] for l in open(f'{O}/regions.tsv') if l.split()[1] == 'dyz17'}
    def clusters(p, merge=100, gap=5000):
        if len(p) == 0: return []
        p = np.sort(p); seg = []; s0, e0 = p[0], p[0] + 31
        for x in p[1:]:
            if x <= e0 + merge: e0 = max(e0, x + 31)
            else: seg.append((s0, e0)); s0, e0 = x, x + 31
        seg.append((s0, e0)); cl = []; cs, ce, ab = seg[0][0], seg[0][1], seg[0][1] - seg[0][0]
        for s0, e0 in seg[1:]:
            if s0 <= ce + gap: ce = max(ce, e0); ab += e0 - s0
            else: cl.append((cs, ce, ab)); cs, ce, ab = s0, e0, e0 - s0
        cl.append((cs, ce, ab)); return cl
    def loc(x):
        a0, f = x; best = None
        for ctg, q in contigs(f):
            k, ok = canon(q); pos = np.nonzero(ok)[0]; k = k[ok]; i = np.minimum(np.searchsorted(K, k), len(K) - 1)
            for c in clusters(pos[K[i] == k]):
                if best is None or c[2] > best[3]: best = (ctg, int(c[0]), int(c[1]), int(c[2]))
        o = OLD.get(a0)
        if best is None: return f'{a0}\tnone\t0\t0\t0\t0\t' + ('\t'.join(o) if o else 'none\t0\t0') + '\tNA\tNA'
        d = (abs(best[1] - int(o[1])), abs(best[2] - int(o[2]))) if o and o[0] == best[0] else ('NA', 'NA')
        return f'{a0}\t{best[0]}\t{best[1]}\t{best[2]}\t{best[2] - best[1]}\t{best[3]}\t' + ('\t'.join(o) if o else 'none\t0\t0') + f'\t{d[0]}\t{d[1]}'
    with Pool(NW) as pool, open(OUT, 'w') as w:
        w.write('assembly\tcontig\tstart\tend\tspan\tassigned\told_contig\told_start\told_end\tstart_diff\tend_diff\n')
        for r in pool.imap_unordered(loc, A): w.write(r + '\n'); print(r, flush=True)
    sys.exit(0)
if sys.argv[1] == 'sup':   # sup <out npz>: class candidates before the outside filter, as u_keys indices
    P = os.environ['P']; U = np.load(f'{P}/markers_v5/u_keys.npy'); LT = np.full(256, 255, np.uint8); LT[[65, 67, 71, 84]] = [0, 1, 2, 3]; Z = {}
    def keys(c):
        f = f'{P}/markers_u25/sup.{c}.txt' if os.path.exists(f'{P}/markers_u25/sup.{c}.txt') else f'{P}/markers_v5/sup.{c}.txt'
        b = LT[np.frombuffer(''.join(l[:31] for l in open(f)).encode(), np.uint8).reshape(-1, 31)].astype(np.uint64); k = np.zeros(len(b), np.uint64)
        for i in range(31): k = (k << np.uint64(2)) | b[:, i]
        i = np.searchsorted(U, k); ok = (i < len(U)) & (U[np.minimum(i, len(U) - 1)] == k)
        print(f'{c}: {f}, {len(k)} k-mers, in u_keys {ok.sum()}', flush=True)
        if not ok.all(): sys.exit(f'ERROR {c}: {(~ok).sum()} k-mers not in u_keys')
        return i
    for c, p in {'dyz3': ('dyz3',), 'yq12': ('yq12',), 'dyz1': ('dyz1',), 'dyz2': ('dyz2',), 'dyz18s': ('dyz18', 'yq3k1'), 'yq2k7': ('yq2k7',), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7'), 'dyz19': ('dyz19',), 'dyz17': ('dyz17',)}.items():
        Z[f'm_{c}'] = np.unique(np.concatenate([keys(x) for x in p]))
    np.savez(sys.argv[2], **Z); print('written', sys.argv[2], {k: len(v) for k, v in Z.items()}, flush=True); sys.exit(0)
if sys.argv[1] == 'sum':   # sum <scan dir> <out.tsv>: truth males only
    P = os.environ['P']; O = f'{P}/markers_v5'; D, OUT = sys.argv[2:4]; sys.path.insert(0, f'{P}/scripts'); import truth_cl
    T = truth_cl.load(O); TR = {(l.split('\t')[0], l.split('\t')[1]): l.split('\t')[3] for l in list(open(f'{O}/truth_regions.tsv'))[1:]}
    CS = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7')}; REG = {}
    for l in open(f'{O}/regions.tsv'):
        a, c, g, s, e = l.rstrip('\n').split('\t')[:5]; REG.setdefault((a, c), []).append((g, int(s), int(e)))
    S = {}; BN = {}
    for fn in os.listdir(D):
        if fn.endswith('.bins.tsv'):
            for l in list(open(f'{D}/{fn}'))[1:]: c, g, pl, b, n = l.split('\t'); BN.setdefault((fn[:-9], c), []).append((pl, g, int(b), int(n)))
        elif fn.endswith('.tsv'):
            for l in list(open(f'{D}/{fn}'))[1:]: x = l.split('\t'); S[(x[0], x[1])] = [int(v) for v in x[2:6]]
    print(f'scan files {len(S)} assembly-class rows', flush=True)
    rows = []
    for c in ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'yq2k7', 'dyz18sub', 'dyz19', 'dyz17']:
        for m in sorted(T.get(c, {})):
            k = 'dyz3' if c == 'dyz3' else 'yq12'; a = TR.get((m, k)) if c not in ('dyz19', 'dyz17') else None
            if a is None: a = next((x for x in (m, m + '_R2') if any((x, cc) in REG for cc in CS.get(c, (c,)))), None)
            if a is None or (a, c) not in S: print(f'ERROR {c} {m}: no scanned assembly with the class region', flush=True); continue
            mk, ci, cy, co = S[(a, c)]; iv = sorted(v for cc in CS.get(c, (c,)) for v in REG.get((a, cc), [])); tot = ci + cy + co
            top = sorted(BN.get((a, c), []), key=lambda x: -x[3])[:3]
            rows.append((c, m, a, mk, ci, cy, co, cy / max(tot, 1), co / max(tot, 1), len(iv), iv[-1][2] - iv[0][1] if iv else 0, sum(e - s for g, s, e in iv), ';'.join(f'{pl}:{g}:{b / 1e6:.1f}Mb:{n}' for pl, g, b, n in top)))
    with open(OUT[:-4] + '.truthset.tsv', 'w') as o:   # truth = all marker copies on Y: annotated truth x all copies / region copies
        for r in rows:
            if r[4]: o.write(f'{r[0]}\t{r[1]}\t{T[r[0]][r[1]] * (r[4] + r[5] + r[6]) / r[4]:.0f}\n')
            else: print(f'ERROR {r[0]} {r[1]}: no marker copies in the region', flush=True)
    with open(OUT, 'w') as o:
        o.write('class\tmale\tassembly\tmarkers\tcopies_region\tcopies_yq12_other\tcopies_y_other\tshare_yq12_other\tshare_y_other\tregion_blocks\tregion_span_bp\tregion_bp\ttop_outside_bins\n')
        for r in rows: o.write('\t'.join(f'{v:.6f}' if isinstance(v, float) else str(v) for v in r) + '\n')
    print('class\tmales\tout_share_median\tp90\tmax\tmales_out>1%\t>5%\t>10%\tyq12_other_share_of_out\tblocks_median\tblocks_max\tregion_bp/span_median', flush=True)
    for c in dict.fromkeys(r[0] for r in rows):
        R = [r for r in rows if r[0] == c]; sh = np.array([r[7] + r[8] for r in R]); oy = sum(r[5] for r in R); oo = sum(r[6] for r in R)
        bl = np.array([r[9] for r in R]); fr = np.array([r[11] / r[10] for r in R if r[10]])
        print(f'{c}\t{len(R)}\t{np.median(sh):.4f}\t{np.percentile(sh, 90):.4f}\t{sh.max():.4f}\t{(sh > 0.01).sum()}\t{(sh > 0.05).sum()}\t{(sh > 0.1).sum()}\t{oy / max(oy + oo, 1):.4f}\t{np.median(bl):.0f}\t{bl.max()}\t{np.median(fr):.4f}', flush=True)
    sys.exit(0)
P = os.environ['P']; O = f'{P}/markers_v5'; MD, UK, OUT = sys.argv[1:4]; NW = int(sys.argv[4]) if len(sys.argv) > 4 else 8
CL = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'yq2k7', 'dyz18sub', 'dyz19', 'dyz17']; CS = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7')}
U = np.load(UK); fm = dict(np.load(MD)) if os.path.isfile(MD) else {k: v for c in CL if os.path.exists(f'{MD}_{c}.npz') for k, v in np.load(f'{MD}_{c}.npz').items()}   # model npz or per-class prefix
CL = [c for c in CL if f'm_{c}' in fm]; KEYS = {c: np.sort(U[fm[f'm_{c}']]) for c in CL}
REG = {}
for l in open(f'{O}/regions.tsv'):
    a, c, ctg, s, e = l.rstrip('\n').split('\t')[:5]; REG.setdefault(a, []).append((c, ctg, int(s), int(e)))
def one(x):
    a, f = x
    if os.path.exists(f'{OUT}/{a}.tsv'): return f'{a} done'
    R = REG.get(a, []); T = {c: [0, 0, 0] for c in CL}; B = []; L = 0
    for ctg, seq in contigs(f):
        k, ok = canon(seq); L += len(seq); pos = np.nonzero(ok)[0]; k = k[ok]
        y12 = [(s, e) for c, g, s, e in R if g == ctg and c == 'yq12']
        for c in CL:
            K = KEYS[c]; i = np.minimum(np.searchsorted(K, k), len(K) - 1); h = pos[K[i] == k]
            if not len(h): continue
            ir = inside(h, [(s, e) for cc, g, s, e in R if g == ctg and cc in CS.get(c, (c,))]); iy = inside(h, y12) & ~ir
            T[c][0] += int(ir.sum()); T[c][1] += int(iy.sum()); T[c][2] += int((~ir & ~iy).sum());
            for lab, m in (('yq12_other', iy), ('y_other', ~ir & ~iy)):
                if m.any():
                    u, n = np.unique(h[m] // 100000, return_counts=True); B += [(c, ctg, lab, int(v) * 100000, int(w)) for v, w in zip(u, n)]
    with open(f'{OUT}/{a}.bins.tsv.tmp', 'w') as o:
        o.write('class\tcontig\tplace\tbin_start\tcopies\n'); o.writelines(f'{c}\t{g}\t{lab}\t{s}\t{n}\n' for c, g, lab, s, n in B)
    with open(f'{OUT}/{a}.tsv.tmp', 'w') as o:
        o.write('assembly\tclass\tmarkers\tcopies_region\tcopies_yq12_other\tcopies_y_other\tassembly_bp\tregion_rows\n')
        for c in CL: o.write(f'{a}\t{c}\t{len(KEYS[c])}\t{T[c][0]}\t{T[c][1]}\t{T[c][2]}\t{L}\t{sum(cc in CS.get(c, (c,)) for cc, g, s, e in R)}\n')
    os.replace(f'{OUT}/{a}.bins.tsv.tmp', f'{OUT}/{a}.bins.tsv'); os.replace(f'{OUT}/{a}.tsv.tmp', f'{OUT}/{a}.tsv')
    return f'{a} bp {L} ' + ' '.join(f'{c}:{T[c][0]}/{T[c][1]}/{T[c][2]}' for c in CL)
if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True); A = [l.split()[:2] for l in open(f'{O}/assemblies.tsv')]
    if os.environ.get('ASM'): A = [x for x in A if x[0] in os.environ['ASM'].split(',')]   # test subset
    print(f'model {MD}, assemblies {len(A)}, markers ' + ' '.join(f'{c} {len(KEYS[c])}' for c in CL), flush=True)
    miss = [a for a, f in A if not os.path.exists(f)]
    if miss: print(f'ERROR missing assembly files: {miss}', flush=True)
    with Pool(NW) as p:
        for r in p.imap_unordered(one, [x for x in A if x[0] not in miss]): print(r, flush=True)
