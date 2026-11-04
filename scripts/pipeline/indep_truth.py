# Independent class truth from our model hits
import sys, glob, os, collections as co
D, OUT, HMM = sys.argv[1], sys.argv[2], sys.argv[3]; EMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-10
LN = {}; nm = None  # model lengths
for l in open(HMM):
    if l.startswith('NAME'): nm = l.split()[1]
    elif l.startswith('LENG'): LN[nm] = int(l.split()[1])
YQ = ('DYZ1', 'DYZ2', 'DYZ18', 'Yqhet_3k1bp', 'Yqhet_2k7bp'); GA, GY = 5000, 100000  # array gap, Yq12 gap
EM = {'DYZ19': 1e-5, 'DYZ3-prim': 1e-5}
def clusters(L, g):  # merge sorted segments within g
    out = []
    for a, b, u in L:
        if out and a - out[-1][1] <= g: out[-1][1] = max(out[-1][1], b); out[-1][2] += u
        else: out.append([a, b, u])
    return out
with open(OUT, 'w') as o, open(OUT + '.arrays.bed', 'w') as ob:
    cols = ['dyz1', 'dyz2', 'dyz18', 'q3k1', 'q2k7']
    o.write('assembly\tcontig\tyq12_start\tyq12_end\tyq12_span\t' + '\t'.join(c + '_units' for c in cols) + '\tnon_dyz\tdyz19_array_units\tdyz19_contig\tdyz3hor_array_units\tdyz3_contig\n')
    for f in sorted(glob.glob(f'{D}/*.tbl')):
        s = os.path.basename(f)[:-4]; H = co.defaultdict(list)
        for l in open(f):
            if l[0] == '#': continue
            x = l.split(); q = x[2]
            if float(x[12]) > EM.get(q, EMAX): continue
            a, b = int(x[6]), int(x[7]); H[x[0]].append((min(a, b) - 1, max(a, b), float(x[13]), q, int(x[5]) - int(x[4]) + 1 >= 0.5 * LN[q]))
        seg = co.defaultdict(lambda: co.defaultdict(list)); full = co.defaultdict(lambda: co.defaultdict(list))
        for c, hs in H.items():  # best score per segment
            ev = sorted({p for a, b, _, _, _ in hs for p in (a, b)}); hs.sort(); act = []; j = 0
            for k in range(len(ev) - 1):
                p0, p1 = ev[k], ev[k + 1]
                while j < len(hs) and hs[j][0] <= p0: act.append(hs[j]); j += 1
                act = [h for h in act if h[1] > p0]
                if act: w = max(act, key=lambda h: h[2]); seg[c][w[3]].append((p0, p1, p1 - p0)); full[c][w[3]].append((p0, p1, (p1 - p0) * w[4]))
        best = None  # largest Yq12 block over contigs
        for c in seg:
            L = sorted(t for q in YQ for t in full[c].get(q, []) if t[2] > 0)  # boundaries from near-full units
            for a, b, u in clusters(L, GY):
                if best is None or u > best[3]: best = (c, a, b, u)
        def arr(q):  # largest array of q over contigs
            r = (0, '-', 0, 0)
            for c in seg:
                for a, b, u in clusters(sorted(seg[c].get(q, [])), GA):
                    if u > r[0]: r = (u, c, a, b)
            if r[0]: ob.write(f'{r[1]}\t{r[2]}\t{r[3]}\t{s}\t{q}\n')
            return r
        d19, d3 = arr('DYZ19'), arr('DYZ3-sec')
        for q in LN:   # models outside Yq12 without a dedicated column (DYZ17): largest array, same rule as DYZ19
            if q not in YQ + ('DYZ19', 'DYZ3-sec', 'DYZ3-prim'): arr(q)
        if best is None: o.write(f'{s}\t-\t0\t0\t0\t' + '\t'.join('0' for _ in cols) + f'\t0\t{d19[0]}\t{d19[1]}\t{d3[0]}\t{d3[1]}\n'); continue
        c, a, b, _ = best; U = []; ob.write(f'{c}\t{a}\t{b}\t{s}\tYQ12\n')
        for q in YQ:  # DYZ18-subregion classes: main array only
            L = [(max(st, a), min(e, b), min(e, b) - max(st, a)) for st, e, _ in seg[c].get(q, []) if e > a and st < b]
            C = clusters(sorted(L), GA)
            if q in ('DYZ1', 'DYZ2'): U.append(sum(t[2] for t in L)); [ob.write(f'{c}\t{x[0]}\t{x[1]}\t{s}\t{q}\n') for x in C]
            else:
                m = max(C, key=lambda t: t[2], default=None); U.append(m[2] if m else 0)
                if m: ob.write(f'{c}\t{m[0]}\t{m[1]}\t{s}\t{q}\n')
        o.write(f'{s}\t{c}\t{a}\t{b}\t{b-a}\t' + '\t'.join(map(str, U)) + f'\t{b-a-sum(U)}\t{d19[0]}\t{d19[1]}\t{d3[0]}\t{d3[1]}\n')
        print(s, c, b - a, flush=True)
