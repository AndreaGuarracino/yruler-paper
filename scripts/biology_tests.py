#!/usr/bin/env python3
"""Tests on the final estimates: lineages, arrays, father-son ratio, composition, HGDP/SGDP DYZ19, quantiles of one CV score per man, sister pairs, annotation floor.

usage: python3 scripts/biology_tests.py <results dir> <ours indep.tsv> <hallast indep.tsv> <out.tsv> [Dunn pairs out.tsv]
inputs: <results dir>/{cohort_1kgp,cohort_hgdp_sgdp,cv_units,holdout_21_units}.tsv, results/estimates_1kgp.relatedness_samples.tsv, manuscript/figures/figureS8_pedigree_pairs.tsv, data/method/yq12_truth_hallast.tsv, data/figure_metadata/hallast2021_supplement.xlsx (downloaded if missing), both Y trees
"""
import sys, csv, re, math, collections as co, itertools
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr, kruskal, norm, rankdata, wilcoxon
sys.path.insert(0, str(Path(__file__).resolve().parent)); import figures_1kgp as fs
D, OURS, HAL, OUT = sys.argv[1:5]; DP = sys.argv[5] if len(sys.argv) > 5 else None; D = Path(D); N = []; rng = np.random.default_rng(20261102); rng2 = np.random.default_rng(20261103)
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
def put(k, v): N.append((k, v)); print(f'{k}\t{v}', flush=True)
X = rd(D / 'cohort_1kgp.tsv'); C = {x['sample']: x for x in X}; b = lambda x, c: float(x[f'{c}_bp']); U1, U2 = 3478.9, 2404.4
keep = ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain'); inc = [x for x in X if x['sex'] == 'male' and x['estimate_basis'] == 'ruler' and x['y_state'] in keep]
kin = {r['sample']: r['kinship_component_id'] for r in rd('results/estimates_1kgp.relatedness_samples.tsv')}; rep = {}
for x in sorted(inc, key=lambda r: r['sample']): rep.setdefault(kin.get(x['sample'], x['sample']), x)
R1 = list(rep.values()); put('one_per_component_men', f'{len(R1)} of {len(inc)} (PAR2-gain men left out)')
CL = ['DYZ3', 'Yq12', 'DYZ1', 'DYZ2', 'DYZ18S', 'DYZ19', 'DYZ17']; val = lambda x, c: b(x, 'DYZ1') / b(x, 'DYZ2') * U2 / U1 if c == 'ratio' else b(x, c)
G = co.defaultdict(list)
for x in R1: G[x['lineage']].append(x)
LN = sorted(k for k, v in G.items() if len(v) >= 20 and k not in ('.', '?', '')); put('lineages_ge20', f"{len(LN)}: " + ', '.join(f'{k} {len(G[k])}' for k in LN)); P = []
for c in CL + ['ratio']:
    g = [np.array([val(x, c) for x in G[k]]) for k in LN]; H, p = kruskal(*g); put(f'kruskal_{c}', f'H {H:.1f}, p {p:.3g}')
    A = np.concatenate(g); r = rankdata(A); n = len(A); o = 0; m = {}
    for k, v in zip(LN, g): m[k] = (r[o:o + len(v)].mean(), len(v)); o += len(v)
    T = [(a, q, 2 * norm.sf(abs(m[a][0] - m[q][0]) / math.sqrt(n * (n + 1) / 12 * (1 / m[a][1] + 1 / m[q][1])))) for a, q in itertools.combinations(LN, 2)]; T.sort(key=lambda t: t[2]); adj = []; run = 0
    for i, t in enumerate(T): run = max(run, min(1, t[2] * (len(T) - i))); adj.append((c, t[0], t[1], f'{t[2]:.3g}', f'{run:.3g}'))
    P += adj; put(f'dunn_holm_{c}', f'{sum(float(a[4]) < .05 for a in adj)} of {len(adj)} pairs below 0.05')
if DP:
    with open(DP, 'w') as o: o.write('target\tlineage_a\tlineage_b\tp\tholm_p\n'); o.writelines('\t'.join(map(str, r)) + '\n' for r in P)
def hl(a, c):   # Hodges-Lehmann ratio of two groups
    a, c = np.log(a), np.log(c); e = lambda p, q: np.median(np.subtract.outer(p, q)); B = [e(rng.choice(a, len(a)), rng.choice(c, len(c))) for _ in range(2000)]; return np.exp(e(a, c)), np.exp(np.percentile(B, 2.5)), np.exp(np.percentile(B, 97.5))
v = hl(np.array([b(x, 'DYZ3') for x in R1 if x['lineage'] == 'R1b']), np.array([b(x, 'DYZ3') for x in R1 if x['lineage'] != 'R1b'])); put('hl_dyz3_r1b_over_other', '%.3f (%.3f-%.3f)' % v)
J = [x for x in R1 if x['lineage'] == 'J']; ja = np.array([b(x, 'DYZ19') for x in J if x['haplogroup'].startswith('J2a1a1')]); jo = np.array([b(x, 'DYZ19') for x in J if not x['haplogroup'].startswith('J2a1a1')]); put('hl_dyz19_j2a1a1_over_other_j', '%.3f (%.3f-%.3f)' % hl(ja, jo) + f'; {len(ja)} vs {len(jo)}')
M = np.array([[math.log(b(x, c)) for c in CL] for x in inc]); L = np.array([x['lineage'] for x in inc]); W = M.copy()
for k in set(L): W[L == k] -= W[L == k].mean(0)
for i, j in itertools.combinations(range(len(CL)), 2): put(f'corr_{CL[i]}_{CL[j]}', f'all {spearmanr(M[:, i], M[:, j])[0]:+.3f}, within lineage {spearmanr(W[:, i], W[:, j])[0]:+.3f}')
pairs = {(r['father'], r['son']): r['batch_group'] for r in rd('manuscript/figures/figureS8_pedigree_pairs.tsv')}; I = [p for p in pairs if all(C[s]['y_state'] == 'normal' for s in p)]; u = lambda s: abs(math.log(val(C[s], 'ratio')))
for g, S in (('all', I), ('cross', [p for p in I if pairs[p] == 'Cross batch']), ('same', [p for p in I if pairs[p] == 'Same batch'])):
    f, n = np.array([u(p[0]) for p in S]), np.array([u(p[1]) for p in S]); put(f'father_son_ratio_{g}', f'{len(S)} pairs: son closer to 1 in {int((n < f).sum())}, median |log ratio| father {np.median(f):.4f} son {np.median(n):.4f}, Wilcoxon p {wilcoxon(f, n).pvalue:.3g}')
def within(v, l):   # subtract the lineage means
    v = np.array(v, float)
    for k in set(l): v[l == k] -= v[l == k].mean()
    return v
def sdratio(d1, d2, l, label):   # SD log DYZ1 / SD log DYZ2 within lineage; Pitman-Morgan test
    a, c = within(np.log(d1), l), within(np.log(d2), l); n = len(a); B = [a[i].std() / c[i].std() for i in (rng2.integers(0, n, n) for _ in range(2000))]; r = pearsonr(a - c, a + c)
    put(f'composition_sd_{label}', f'{n}: SD {a.std():.3f} and {c.std():.3f}, ratio {a.std() / c.std():.3f} (95% {np.percentile(B, 2.5):.3f} to {np.percentile(B, 97.5):.3f}), Pitman-Morgan r {r[0]:.3f}, P {r[1]:.2g}'); return a, c
sc = [x for x in R1 if x['reference_male'] != 'yes']; lr = np.array([math.log(b(x, 'DYZ1') / b(x, 'DYZ2')) for x in sc]); ls = np.array([math.log(b(x, 'Yq12')) for x in sc]); l2 = np.array([x['lineage'] for x in sc])
d1, d2, dy = (np.array([b(x, k) for x in sc]) for k in ('DYZ1', 'DYZ2', 'Yq12')); sdratio(d1, d2, l2, '1kgp'); rest = dy - d1 - d2; NL = []
for _ in range(1000):   # DYZ2 moved between men of a lineage
    p = np.arange(len(sc))
    for k in sorted(set(l2)): w = np.where(l2 == k)[0]; p[w] = rng2.permutation(w)
    z = d2[p]; x_ = within(np.log(d1 + z + rest), l2); y_ = within(np.log(d1 / z), l2); NL.append((x_ * y_).sum() / (x_ * x_).sum())
put('composition_arithmetic_null', f'slope on Yq12 with DYZ2 permuted within lineage: median {np.median(NL):.3f} (central 95% {np.percentile(NL, 2.5):.3f} to {np.percentile(NL, 97.5):.3f}), 1,000 permutations')
for k in set(l2): lr[l2 == k] -= lr[l2 == k].mean(); ls[l2 == k] -= ls[l2 == k].mean()
sl = np.polyfit(ls, lr, 1)[0]; B = [np.polyfit(ls[i], lr[i], 1)[0] for i in (rng.integers(0, len(ls), len(ls)) for _ in range(2000))]; put('composition_slope_cohort', f'{sl:.3f} (95% {np.percentile(B, 2.5):.3f} to {np.percentile(B, 97.5):.3f}); {len(ls)} men, within lineage, outside the panel')
T = rd('data/method/yq12_truth_hallast.tsv'); lin = [(C[r['assembly']]['lineage'], math.log(float(r['dyz1_units']) / float(r['dyz2_units'])), math.log(float(r['yq12_span']))) for r in T if r['assembly'] in C and r['assembly'] != 'HG01261']; A = np.array([[a, q] for _, a, q in lin]); l3 = np.array([k for k, _, _ in lin])
Tc = [r for r in T if r['assembly'] in C and r['assembly'] != 'HG01261']; sdratio(np.array([float(r['dyz1_units']) for r in Tc]), np.array([float(r['dyz2_units']) for r in Tc]), l3, 'assemblies')
for k in set(l3): A[l3 == k] -= A[l3 == k].mean(0)
B = [np.polyfit(A[i, 1], A[i, 0], 1)[0] for i in (rng.integers(0, len(A), len(A)) for _ in range(2000))]; put('composition_slope_assemblies', f'{np.polyfit(A[:, 1], A[:, 0], 1)[0]:.3f} (95% {np.percentile(B, 2.5):.3f} to {np.percentile(B, 97.5):.3f}); {len(A)} assemblies, within lineage')
E = co.defaultdict(lambda: co.defaultdict(list))
for r in rd(D / 'cv_units.tsv'): E[r['class']][r['male']].append(float(r['error_pct']))
s = {c: np.std([np.mean(np.log1p(np.array(v) / 100)) for v in E[c].values()], ddof=1) for c in ('dyz1', 'dyz2')}; put('composition_error_sd', f"CV SD of log error DYZ1 {s['dyz1']:.4f}, DYZ2 {s['dyz2']:.4f}")
sys.path.insert(0, 'scripts'); import prepare_cohort_pies as pc, openpyxl
tb = list(pc.hallast_book()['Suppl2'].values); H = {}
for f in tb[2:]:
    r = dict(zip(tb[1], f))
    for k in ('Sample_ID', 'vcf_ID'):
        if pc.valid_key(r[k]): H[pc.norm(r[k])] = (str(r['hg_ISOGG_v03.10.19']), str(r['terminal marker']))
AL = co.defaultdict(list)
for r in pc.read_tsv(pc.META / 'sgdp_metadata.txt', encoding='latin-1'):
    for k in pc.ALIAS_FIELDS:
        if pc.valid_key(r[k]): AL[pc.norm(r[k])].append(r)
Q = rd(D / 'cohort_hgdp_sgdp.tsv'); seen = set(); got = []; nh = 0
def call(s):
    id0 = s.replace('SGDP_', ''); ks = [s, id0] + [r[k] for r in AL.get(pc.norm(id0), []) for k in pc.ALIAS_FIELDS if pc.valid_key(r[k])]
    for k in ks:
        if pc.norm(k) in H: return H[pc.norm(k)]
for x in sorted(Q, key=lambda r: r['sample']):
    c = call(x['sample']); nh += c is not None
    if c is None or x['library'] != 'PCR-free' or x['estimate_basis'] != 'ruler' or x['y_state'] not in keep or x['sex'] == 'female': continue
    i = x['sample'].replace('SGDP_', '')
    if i in seen or i in C: continue
    seen.add(i); got.append((x, c))
put('hgsg_calls', f'{nh} of {len(Q)} with a published call; {len(got)} PCR-free individuals after removing 1kGP aliases and duplicates'); JJ = [(x, c) for x, c in got if c[0].startswith('J')]
sh = [(x, c) for x, c in JJ if b(x, 'DYZ19') < 1e5]; put('hgsg_j', f"{len(JJ)} J men; DYZ19 below 100 kb in {len(sh)}: {min(b(x, 'DYZ19') for x, _ in sh) / 1e3:.1f}-{max(b(x, 'DYZ19') for x, _ in sh) / 1e3:.1f} kb; others min {min(b(x, 'DYZ19') for x, _ in JJ if b(x, 'DYZ19') >= 1e5) / 1e3:.1f} kb")
for co_ in ('HGDP', 'SGDP'):
    G = [(x, c) for x, c in got if x['cohort'] == co_]; sdratio(np.array([b(x, 'DYZ1') for x, _ in G]), np.array([b(x, 'DYZ2') for x, _ in G]), np.array([c[0][:1] for _, c in G]), co_.lower())
with open(D / 'hgsg_calls.tsv', 'w') as o:   # input of Fig. 2
    o.write('sample\tcohort\tcall\tDYZ3_bp\tDYZ19_bp\n'); o.writelines(f"{x['sample']}\t{x['cohort']}\t{c[0]}\t{x['DYZ3_bp']}\t{x['DYZ19_bp']}\n" for x, c in got)
M3 = lambda f: [b(x, 'DYZ3') for x, c in got if f(c[0])]; put('hgsg_dyz3', f"R1b {len(M3(lambda h: h.startswith('R1b')))} men, median {np.median(M3(lambda h: h.startswith('R1b'))):.0f}; others {len(M3(lambda h: not h.startswith('R1b')))}, median {np.median(M3(lambda h: not h.startswith('R1b'))):.0f}; O2a2b1a1 {len(M3(lambda h: h.startswith('O2a2b1a1')))}, median {np.median(M3(lambda h: h.startswith('O2a2b1a1'))):.0f}; O2a2b1a2 {len(M3(lambda h: h.startswith('O2a2b1a2')))}, median {np.median(M3(lambda h: h.startswith('O2a2b1a2'))):.0f}")
put('hgsg_j_labels_short', str(co.Counter(c[0] for _, c in sh).most_common())); put('hgsg_j_labels_other', str(co.Counter(c[0] for x, c in JJ if b(x, 'DYZ19') >= 1e5).most_common(14)))
put('hgsg_short_outside_j', str([(x['sample'], c[0], round(b(x, 'DYZ19'))) for x, c in got if b(x, 'DYZ19') < 1e5 and not c[0].startswith('J')]))
H21 = co.defaultdict(list)
for r in rd(D / 'holdout_21_units.tsv'): H21[r['class']].append(float(r['error_pct']))
for c in ('dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17'):
    sc2 = np.sort([np.median(np.abs(np.log1p(np.array(v) / 100))) for v in E[c].values()]); n = len(sc2); o = []
    for al in (.2, .1, .05):
        k = math.ceil((n + 1) * (1 - al)); q = sc2[min(k, n) - 1]; o.append(f'{100 * (math.exp(q) - 1):.2f}%')
    q90 = sc2[min(math.ceil((n + 1) * .9), n) - 1]; h = np.abs(np.log1p(np.array(H21[c]) / 100)); put(f'cv_score_quantile_{c}', f'80/90/95% {" / ".join(o)}; holdout within 90%: {int((h <= q90).sum())} of {len(h)}; {n} men')
fam = {x['sample']: x['lineage'] for x in X}
IN = {r['sample']: float(r['overlap_share']) for r in rd('results/checks/insert_1kgp.tsv')}; PX = [p for p, g in pairs.items() if g == 'Cross batch' and p[0] in IN and p[1] in IN]
for c in ('DYZ3', 'Yq12', 'DYZ1', 'DYZ2', 'DYZ18S', 'DYZ19', 'DYZ17'):
    r = spearmanr([math.log(b(C[q], c) / b(C[f_], c)) for f_, q in PX], [IN[q] - IN[f_] for f_, q in PX]); put(f'overlap_father_son_{c}', f'{len(PX)} cross-batch pairs, Spearman {r[0]:+.3f}, P {r[1]:.2g}')
v = sorted(IN.values()); put('overlap_share', f'median {np.median(v):.4f}, 1st to 99th percentile {np.percentile(v, 1):.4f} to {np.percentile(v, 99):.4f}, {len(v)} samples')
def ch(n, out):
    if not n['children']: return
    k = n['children']
    if len(k) == 2 and not k[0]['children'] and not k[1]['children'] and k[0]['name'] in C and k[1]['name'] in C and all(C[y['name']] in inc for y in k): out.append((k[0]['name'], k[1]['name']))
    for z in k: ch(z, out)
for t, f in (('v4', 'results/tree/tree.rooted.nwk'), ('v3', 'results/tree_biallelic/tree.rooted.nwk')):
    S = []; sys.setrecursionlimit(20000); ch(fs.parse_newick(fs.ROOT / f), S); put(f'sister_pairs_{t}', f'{len(S)} pairs: ' + ', '.join(f"{c} {sum(abs(math.log(b(C[a], c) / b(C[q], c))) > math.log(1.25) for a, q in S)}" for c in ('DYZ19', 'DYZ3', 'Yq12')))
O = {r['assembly']: r for r in rd(OURS)}; Hh = {r['assembly']: r for r in rd(HAL)}; MM = sorted(set(O) & set(Hh))
for c, n in (('dyz1_units', 'DYZ1'), ('dyz2_units', 'DYZ2'), ('dyz18_units', 'DYZ18'), ('yq12_span', 'Yq12 span')):
    d = np.array([100 * (float(Hh[m][c]) / float(O[m][c]) - 1) for m in MM]); r = d - np.median(d); put(f'annotation_floor_{n}', f'{len(d)} assemblies: offset {np.median(d):+.2f}%, median |difference around it| {np.median(np.abs(r)):.2f}%, p90 {np.percentile(np.abs(r), 90):.2f}%, max {np.abs(r).max():.2f}%')
with open(OUT, 'w') as o: o.write('item\tvalue\n'); o.writelines(f'{k}\t{v}\n' for k, v in N)
