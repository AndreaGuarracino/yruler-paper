# Depth and lineage-out tables of the model
# usage: [MSRC=class:cv suffix,...; yq12dyz2:<suffix> for the DYZ2 part of Yq12] python3 um_mexport.py (class keys as in the CV files) <out dir> <fractions,comma>; per class CV source below; Yq12 = Yq12 - its DYZ2 part + DYZ2 with the GC floor
import os, sys, glob, re, numpy as np
P = '/path/to/project'; O = f'{P}/markers_u25'; X = f'{P}/truth_qc/x'; sys.path.insert(0, f'{P}/scripts'); import truth_cl
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
def log(*a): print(*a, flush=True)
T = truth_cl.load(O)
for l in open(f'{X}/rule_w4.txt'): c, s = l.split()[:2]; T[c].pop(s, None)
for l in open(f'{X}/truth_all.tsv'):
    c, s, v = l.split()[:3]
    if s in T[c]: T[c][s] = float(v)
SRC = {'dyz3': 'v4', 'yq12': 'w4', 'dyz1': 'v3a', 'dyz2': 'v4', 'dyz18s': 'v3a', 'yq2k7': 'v3a', 'dyz18sub': 'v3a', 'dyz19': 'v4'}   # base DYZ2 of Yq12: w4
SRC.update(x.split(':') for x in os.environ.get('MSRC', '').split(',') if x)
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'yq2k7': 'YqHet2.7', 'dyz18sub': 'DYZ18sub', 'dyz19': 'DYZ19'}
def rd(pat, c, n):   # rows of class c: id (first n-6 fields + sample) -> row
    R = {}
    for fn in glob.glob(pat):
        for l in open(fn):
            x = l.rstrip('\n').split('\t')
            if x[n - 5] == c: R[tuple(x[:n - 5]) + (x[n - 4],)] = x
    return R
def classes(pat, n, key):   # class -> {id: row}; Yq12 rows rebuilt from its parts, matched on key(id)
    R = {c: rd(pat(SRC[c]), c, n) for c in NM}; Y, D0 = R['yq12'], {key(i): x for i, x in rd(pat(SRC.get('yq12dyz2', 'w4')), 'dyz2', n).items()}; D = {key(i): x for i, x in R['dyz2'].items()}
    for i, x in list(Y.items()):
        s = i[-1]; k = key(i)
        if k not in D0 or k not in D or s not in T['dyz2']: del Y[i]; continue
        e = lambda r, t: t * (1 + float(r[-1]) / 100); x[-1] = f'{100 * ((e(x, T["yq12"][s]) - e(D0[k], T["dyz2"][s]) + e(D[k], T["dyz2"][s])) / T["yq12"][s] - 1):.4f}'
    return R
st = lambda a: f'{np.mean(np.abs(a)):.2f}\t{np.median(np.abs(a)):.2f}\t{np.percentile(np.abs(a), 90):.2f}\t{np.max(np.abs(a)):.2f}'
def dep(f, s):   # chr2 depth-marker mean count of the downsampled test run
    g = f'{P}/yruler_u25_depth/f{f}/{s}/run.log'; m = re.search(r'chr2depth_mean ([0-9.eE+-]+)', open(g).read()) if os.path.exists(g) else None
    return m.group(1) if m else 'nan'
SM = ['fraction\tclass\tn_predictions\trepeats\tmean_abs_pct\tmedian_abs_pct_mean_of_repeats\tp90_abs_pct\tmax_abs_pct\tmin_repeat_median\tmax_repeat_median']
for f in sys.argv[2].split(','):
    R = classes(lambda m: f'{O}/cv_{m}/rep*.tsv' if f == '1' else f'{O}/cv_{m}_depth_f{f}/rep*.tsv', 7, lambda i: (i[0], i[-1]))
    with open(f'{OUT}/depth_f{f}.tsv', 'w') as o:
        o.write('repeat\tfold\tclass\tsample\tchoice_presence_outside_ruler_filter\tmarkers\terror_pct\treference_bp\testimate_bp\tchr2_depth_test\n')
        for c in NM:
            for i, x in sorted(R[c].items(), key=lambda z: (int(z[0][0]), int(z[0][1]), z[0][2])): o.write('\t'.join([x[0], x[1], c, x[3], x[4], x[5], x[6], f'{T[c][x[3]]:.1f}', f'{T[c][x[3]] * (1 + float(x[6]) / 100):.0f}', dep(f, x[3])]) + '\n')
    for c in NM:
        a = np.array([float(x[6]) for x in R[c].values()]); reps = sorted({int(x[0]) for x in R[c].values()})
        if not len(a): log(f'WARNING f{f} {c}: no rows'); continue
        m = [np.median(np.abs([float(x[6]) for x in R[c].values() if int(x[0]) == q])) for q in reps]
        SM.append(f'{f}\t{NM[c]}\t{len(a)}\t{len(reps)}\t{np.mean(np.abs(a)):.2f}\t{np.mean(m):.2f}\t{np.percentile(np.abs(a), 90):.2f}\t{np.max(np.abs(a)):.2f}\t{min(m):.2f}\t{max(m):.2f}')
    log(f'f{f}: ' + ', '.join(f'{NM[c]} {len(R[c])}' for c in NM))
open(f'{OUT}/depth_summary.tsv', 'w').write('\n'.join(SM) + '\n'); log('\n'.join(SM))
L = classes(lambda m: f'{O}/lolo_{m}/*.tsv', 6, lambda i: i)
with open(f'{OUT}/lineage_out.tsv', 'w') as o, open(f'{OUT}/lineage_out_summary.tsv', 'w') as q:
    o.write('lineage\tclass\tsample\tchoice_presence_outside_ruler_filter\tmarkers\terror_pct\treference_bp\testimate_bp\n'); q.write('class\tlineage\tn\tmean_abs_pct\tmedian_abs_pct\tp90_abs_pct\tmax_abs_pct\n')
    for c in NM:
        for i, x in sorted(L[c].items()): o.write('\t'.join([x[0], c, x[2], x[3], x[4], x[5], f'{T[c][x[2]]:.1f}', f'{T[c][x[2]] * (1 + float(x[5]) / 100):.0f}']) + '\n')
        for g in ['all'] + sorted({x[0] for x in L[c].values()}):
            a = np.array([float(x[5]) for x in L[c].values() if g in ('all', x[0])]); q.write(f'{NM[c]}\t{g}\t{len(a)}\t{st(a)}\n')
log(''.join(l for l in open(f'{OUT}/lineage_out_summary.tsv') if '\tall\t' in l or l.startswith('class')))
