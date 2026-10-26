# Class regions per truth assembly, one rule
import glob, os, collections
P = '/path/to/project'; W = f'{P}/ext'; O = f'{P}/markers_u'
U = {'DYZ1': 'dyz1', 'DYZ2': 'dyz2', 'DYZ18': 'dyz18', 'Yqhet_3k1bp': 'yq3k1', 'Yqhet_2k7bp': 'yq2k7', 'DYZ19': 'dyz19'}
fa = lambda a: f'{P}/ks_truth/fa/{a}.fa.gz' if os.path.exists(f'{P}/ks_truth/fa/{a}.fa.gz') else f'{P}/r2_truth/fa_bgz/{a}.fa.gz'
R = []; tr = [l.rstrip('\n').split('\t') for l in open(f'{O}/truth_regions.tsv')][1:]
for m, c, src, a, ctg, s, e in tr: R.append((a, c, ctg, int(s), int(e)))
yq = {a for m, c, src, a, *_ in tr if c == 'yq12'}
d19 = {l.split('\t')[0] for l in open(f'{O}/dyz19_clean.tsv') if l.rstrip('\n').endswith('True')}
for l in open(f'{P}/hmm_ours/indep.tsv.arrays.bed'):
    ctg, s, e, a, q = l.rstrip('\n').split('\t')
    if q in U and ((q != 'DYZ19' and a in yq) or (q == 'DYZ19' and a in d19)): R.append((a, U[q], ctg, int(s), int(e)))
for f in glob.glob(f'{W}/bed/*.t2tv2*.bed'):
    a = os.path.basename(f).split('.')[0]
    if not os.path.exists(f'{P}/ks_truth/fa/{a}.fa.gz'): continue
    for l in open(f):
        x = l.rstrip('\n').split('\t')
        if l[0] != '#' and len(x) > 6 and x[6] == 'LABEL' and x[3].startswith('XDR'): R.append((a, 'ysinglecopy', x[0], int(x[1]), int(x[2])))
n = collections.defaultdict(set); bp = collections.Counter()
with open(f'{O}/regions.tsv', 'w') as o:
    for a, c, ctg, s, e in sorted(R):
        if e <= s: raise SystemExit(f'ERROR bad region {a} {c} {ctg}:{s}-{e}')
        o.write(f'{a}\t{c}\t{ctg}\t{s}\t{e}\n'); n[c].add(a); bp[c] += e - s
with open(f'{O}/assemblies.tsv', 'w') as o:
    for a in sorted({r[0] for r in R}):
        if not os.path.exists(fa(a)): raise SystemExit(f'ERROR no fasta {a}')
        o.write(f'{a}\t{fa(a)}\n')
for c in sorted(n): print(f'{c}: assemblies {len(n[c])}, regions bp {bp[c]}')
print('assemblies', len({r[0] for r in R}))
