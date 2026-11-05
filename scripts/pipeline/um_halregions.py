# Truth and class regions from the Hallast 2026 HMMER filters
# usage: python3 um_halregions.py <old marker dir> <new marker dir> [filter dir, default hallast_filt]
#   <filter dir>/<assembly>.yq12.csv, .dyz19.csv, .dyz17.csv: filterYq12HMMER.py / filterDYZ19HMMER.py outputs (DYZ17: Yq12 filter on the DYZ17 table alone)
#   length = sum of (end - start) per block, as Hallast Suppl. Table 26; Yq12 classes clipped to the Yq12 truth region; DYZ19 = filterDYZ19HMMER.py array span
#   variants for the CV choice: yq2k7a (every YqHet2.7 block, <a>.yq12all.csv); DYZ17 from the DYZ17_unit nhmmer hits: dyz17 = longest run with gaps <= 10 kb, dyz17a = hit bp outside Yq12
#   writes <new>/indep_setB.tsv, <new>/dyz19_clean.tsv (QC flags of <old>/dyz19_clean.tsv), <new>/dyz17_truth.tsv (dyz17_span; arrays passing the DYZ19 QC, D17QC), dyz17_all.tsv, <new>/regions.tsv (blocks as [start - 1, end))
import sys, os, collections, gzip
P = os.environ['P']; A, B = f'{P}/{sys.argv[1]}', f'{P}/{sys.argv[2]}'; F = f'{P}/{sys.argv[3] if len(sys.argv) > 3 else "hallast_filt"}'; os.makedirs(B, exist_ok=True)
def log(*a): print(*a, flush=True)
Q = {'DYZ1_Yq': 'dyz1', 'DYZ2_Con': 'dyz2', 'DYZ18_Yq': 'dyz18', 'Yqhet_3k1bp': 'yq3k1', 'Yqhet_2k7bp': 'yq2k7'}
yq = {x[3]: (x[4], int(x[5]), int(x[6])) for x in (l.rstrip('\n').split('\t') for l in open(f'{A}/truth_regions.tsv')) if x[1] == 'yq12'}   # assembly -> Yq12 region
asm = [l.split('\t')[0] for l in open(f'{A}/assemblies.tsv')]
def blocks(a, k):   # (contig, element, start, end, hits)
    f = f'{F}/{a}.{k}.csv'
    if not os.path.exists(f): return None
    L = [l.rstrip('\n').split(',') for l in open(f)][1:]
    if k in ('dyz19', 'dyz17'): return [(x[0], x[1], int(x[3]), int(x[4]), int(x[5])) for x in L]
    return [(x[1], x[2], int(x[4]), int(x[5]), int(x[6])) for x in L]
R = [l.rstrip('\n').split('\t') for l in open(f'{A}/regions.tsv') if l.split('\t')[1] in ('dyz3', 'yq12', 'ysinglecopy')]
n = collections.Counter(); miss = collections.defaultdict(list)
with open(f'{B}/indep_setB.tsv', 'w') as o:   # Yq12 classes (bp; hits = repeat units)
    o.write('assembly\tcontig\tyq12_start\tyq12_end\tyq12_span\tdyz1_units\tdyz2_units\tdyz18_units\tq3k1_units\tq2k7_units\t' + '\t'.join(f'{c}_hits' for c in Q.values()) + '\tq2k7a_units\tq2k7a_hits\n')   # q2k7a: every YqHet2.7 block (amount)
    for a, (ctg, s0, e0) in sorted(yq.items()):
        b = blocks(a, 'yq12')
        if b is None: miss['yq12'].append(a); continue
        bp, h = collections.Counter(), collections.Counter()
        for c, q, s, e, k in b:
            s, e = max(s, s0), min(e, e0)
            if c != ctg or e <= s: continue
            bp[Q[q]] += e - s; h[Q[q]] += k; R.append([a, Q[q], c, str(s - 1), str(e)])
        b2 = blocks(a, 'yq12all'); a2 = h2 = 0
        if b2 is None: miss['yq12all'].append(a)
        for c, q, s, e, k in b2 or []:
            s, e = max(s, s0), min(e, e0)
            if c != ctg or e <= s or q != 'Yqhet_2k7bp': continue
            a2 += e - s; h2 += k; R.append([a, 'yq2k7a', c, str(s - 1), str(e)])
        o.write(f'{a}\t{ctg}\t{s0}\t{e0}\t{e0 - s0}\t' + '\t'.join(str(bp[c]) for c in Q.values()) + '\t' + '\t'.join(str(h[c]) for c in Q.values()) + f'\t{a2}\t{h2}\n'); n['yq12'] += 1
L = [l.rstrip('\n').split('\t') for l in open(f'{A}/dyz19_clean.tsv')]; qc = {x[0]: dict(zip(L[0], x)) for x in L[1:]}
with open(f'{B}/dyz19_clean.tsv', 'w') as o:   # DYZ19 array: Hallast filter span, QC flags
    o.write('assembly\tcontig\tstart\tend\tunit_bp_array\tcopies\tflank_ok\tgap\tflag4_windows\tclean\n')
    for a in asm:
        b = blocks(a, 'dyz19')
        if b is None or a not in qc: miss['dyz19'].append(a); continue
        b = [x for x in b if x[0] == qc[a]['contig']]
        if len(b) != 1: miss['dyz19_contig'].append(a); continue
        c, q, s, e, k = b[0]; r = qc[a]; o.write(f'{a}\t{c}\t{s}\t{e}\t{e - s}\t{k}\t{r["flank_ok"]}\t{r["gap"]}\t{r["flag4_windows"]}\t{r["clean"]}\n'); n['dyz19'] += 1
        if r['clean'] == 'True': R.append([a, 'dyz19', c, str(s - 1), str(e)])
d17qc = {x[0] for x in (l.rstrip('\n').split('\t') for l in list(open(os.environ.get('D17QC', f'{P}/markers_v9/dyz17_clean.tsv')))[1:]) if x[-1] == 'True'}   # DYZ17 arrays passing the array QC
def d17hits(a):   # DYZ17_unit nhmmer hits, merged per contig
    f = f'{P}/hallast_ours/{a}/{a}.ours.{os.environ.get("D17M", "DYZ17_unit")}.hmmer-tblout.txt.gz'
    if not os.path.exists(f): return None
    H = sorted((x[0], min(int(x[8]), int(x[9])), max(int(x[8]), int(x[9]))) for x in (l.split() for l in gzip.open(f, 'rt') if not l.startswith('#')))
    M = []
    for c, s, e in H:
        if M and M[-1][0] == c and s <= M[-1][2]: M[-1][2] = max(M[-1][2], e); M[-1][3] += 1
        else: M.append([c, s, e, 1])
    return M
with open(f'{B}/dyz17_truth.tsv', 'w') as o, open(f'{B}/dyz17_all.tsv', 'w') as o2:   # array: longest run, gaps <= 10 kb
    hd = 'sample\tcontig\tstart\tend\tdyz17_span\tarray_hits\tdyz17_amount\tamount_hits'; o.write(hd + '\n'); o2.write(hd + '\tdyz17_qc\n')
    for a in asm:
        M = d17hits(a)
        if not M: miss['dyz17'].append(a); continue
        ctg, s0, e0 = yq.get(a, ('', 0, 0)); C = []; out = [x for x in M if not (x[0] == ctg and x[1] < e0 and x[2] > s0)]   # hits inside Yq12: Hallast classes
        if not out: miss['dyz17'].append(a); continue
        for c, s, e, k in out:
            if C and C[-1][0] == c and s - C[-1][2] <= 10000: C[-1][2] = max(C[-1][2], e); C[-1][3] += k
            else: C.append([c, s, e, k])
        c, s, e, k = max(C, key=lambda x: x[2] - x[1])
        if c == ctg and s < e0 and e > s0: sys.exit(f'ERROR {a}: DYZ17 array overlaps Yq12')
        am = sum(x[2] - x[1] for x in out); ka = sum(x[3] for x in out)
        R.append([a, 'dyz17', c, str(s - 1), str(e)]); R += [[a, 'dyz17a', x[0], str(x[1] - 1), str(x[2])] for x in out]
        line = f'{a}\t{c}\t{s}\t{e}\t{e - s}\t{k}\t{am}\t{ka}'; o2.write(line + f'\t{a in d17qc}\n')
        if a in d17qc: o.write(line + '\n'); n['dyz17'] += 1
bp = collections.Counter(); ns = collections.Counter(); ma = collections.defaultdict(set)
with open(f'{B}/regions.tsv', 'w') as o:
    for x in sorted(R, key=lambda x: (x[0], x[1], x[2], int(x[3]))):
        if int(x[4]) <= int(x[3]): sys.exit(f'ERROR bad region {x}')
        o.write('\t'.join(x) + '\n'); bp[x[1]] += int(x[4]) - int(x[3]); ns[x[1]] += 1; ma[x[1]].add(x[0])
log('written:', dict(n)); log('missing:', {k: (len(v), v[:8]) for k, v in miss.items()})
for c in sorted(ns): log(f'region {c}: rows {ns[c]}, assemblies {len(ma[c])}, bp {bp[c]}')
