# DYZ17 and pericentromeric HSat3 truth per assembly
# usage: python3 um_newtruth.py <out.tsv> [sample ...]; CHM13 Y-specific 31-mers found only in hsat3_Y_5 (A3/DYZ17) or only in the A4 blocks hsat3_Y_1-4; hits merged within 100 bp, clusters within 5 kb; DYZ17 = cluster with most assigned bases; pericentromeric = all A4 clusters between the DYZ3 array and the DYZ17 cluster on the same contig
import os, sys, gzip, numpy as np, multiprocessing as mp, subprocess as sp
P = os.environ['P']; K = f'{P}/ykmers'; F = f'{P}/ks_truth/fa'; OUT = sys.argv[1]
LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84, 97, 99, 103, 116]] = [0, 1, 2, 3, 0, 1, 2, 3]
def log(*a): print(*a, flush=True)
def kmers(seq):   # canonical 31-mer keys and start positions, windows with N left out
    v = LUT[np.frombuffer(seq, np.uint8)]; n = len(v) - 30
    if n <= 0: return np.zeros(0, np.uint64), np.zeros(0, np.int64)
    bad = np.convolve((v == 255).astype(np.int8), np.ones(31, np.int8), 'valid') > 0; v = np.where(v == 255, 0, v).astype(np.uint64)
    k = np.zeros(n, np.uint64); rc = np.zeros(n, np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | v[i:i + n]; rc |= (np.uint64(3) - v[i:i + n]) << np.uint64(2 * i)
    return np.minimum(k, rc)[~bad], np.nonzero(~bad)[0]
IV = {'DYZ17': [(11711672, 12656463)], 'periHSat3': [(10888338, 11070712), (11076842, 11124106), (11145128, 11261834), (11512198, 11525444)]}
S = np.load(f'{K}/ykeys.npy'); z = np.load(f'{K}/ypos.npz'); pos, rank = z['pos'], z['rank']
inside = {c: np.zeros(len(pos), bool) for c in IV}
for c, v in IV.items():
    for a, b in v: inside[c] |= (pos >= a) & (pos < b)
KEY = {}
for c in IV:   # keys whose every CHM13 position is inside the class
    r_in = np.unique(rank[inside[c]]); r_out = np.unique(rank[~inside[c]]); KEY[c] = np.sort(S[np.setdiff1d(r_in, r_out)])
log('class keys', {c: len(v) for c, v in KEY.items()})
def clusters(p, merge=100, gap=5000):   # hit starts -> segments (31 bp, merged within 100 bp) -> clusters within 5 kb: (start, end, assigned bases)
    if len(p) == 0: return []
    p = np.sort(p); seg = []; s, e = p[0], p[0] + 31
    for x in p[1:]:
        if x <= e + merge: e = max(e, x + 31)
        else: seg.append((s, e)); s, e = x, x + 31
    seg.append((s, e)); cl = []; cs, ce, ab = seg[0][0], seg[0][1], seg[0][1] - seg[0][0]
    for s, e in seg[1:]:
        if s <= ce + gap: ce = max(ce, e); ab += e - s
        else: cl.append((cs, ce, ab)); cs, ce, ab = s, e, e - s
    cl.append((cs, ce, ab)); return cl
def fasta(path):
    name, seq = None, []
    with (gzip.open(path, 'rb') if path.endswith('.gz') else open(path, 'rb')) as f:
        for l in f:
            if l.startswith(b'>'):
                if name: yield name, b''.join(seq)
                name, seq = l[1:].split()[0].decode(), []
            else: seq.append(l.strip())
    if name: yield name, b''.join(seq)
DZ3 = {}
for l in open(f'{P}/hmm_ours/indep.tsv.arrays.bed'):
    x = l.rstrip('\n').split('\t')
    if x[4].startswith('DYZ3'): DZ3.setdefault(x[3], []).append((x[0], int(x[1]), int(x[2])))
def one(s):
    path = f'{F}/{s}.fa.gz' if s != 'CHM13' else None; best = None; A4 = {}
    seqs = fasta(path) if path else [('chrY', ''.join(sp.run(['samtools', 'faidx', f'{P}/ext/ref/chm13v2.0.fa', 'chrY'], capture_output=True, text=True, check=True).stdout.split('\n')[1:]).encode())]
    for name, seq in seqs:
        k, p = kmers(seq)
        for c in IV:
            i = np.minimum(np.searchsorted(KEY[c], k), len(KEY[c]) - 1); h = p[KEY[c][i] == k]
            cl = clusters(h)
            if c == 'DYZ17':
                for x in cl:
                    if best is None or x[2] > best[3]: best = (name, x[0], x[1], x[2], len(h))
            else: A4[name] = cl
    if best is None: return f'{s}\tnone\t0\t0\t0\t0\tnone\t0\t0\t0'
    nm, b0, b1 = best[0], best[1], best[2]; d3 = [x for x in DZ3.get(s, []) if x[0] == nm]
    lo = max([x[2] for x in d3 if x[2] <= b0], default=0) if d3 else max(0, b0 - 3_000_000)   # pericentromeric window: DYZ3 end (or 3 Mb) to DYZ17 start
    sel = [x for x in A4.get(nm, []) if x[0] >= lo and x[1] <= b0 + 1000]
    return f'{s}\t{nm}\t{b0}\t{b1}\t{b1 - b0}\t{best[3]}\t{"dyz3" if d3 else "3Mb"}\t{len(sel)}\t{sum(x[1] - x[0] for x in sel)}\t{sum(x[2] for x in sel)}'
L = sys.argv[2:] or ['CHM13'] + sorted(f[:-6] for f in os.listdir(F) if f.endswith('.fa.gz'))
with mp.get_context('fork').Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', 8))) as pool, open(OUT, 'w') as o:
    o.write('sample\tcontig\tdyz17_start\tdyz17_end\tdyz17_span\tdyz17_assigned\tperi_window\tperi_clusters\tperi_span_sum\tperi_assigned\n')
    for r in pool.imap(one, L): o.write(r + '\n'); log(r)
log(f'written {OUT}')
