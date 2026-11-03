#!/usr/bin/env python3
"""Relatedness checks: HGDP/SGDP relatives, Y-tree distance in CV, third-degree 1kGP links, male siblings.

usage: python3 scripts/relatedness_checks.py results/model results/checks/relatedness/gnomad_related_pairs.tsv results/checks/relatedness.tsv results/checks/relatedness/sgdp_king_pairs.tsv
(the ID lists gnomad_hgdp_1kg_ids.txt and sgdp_plink_ids.txt sit next to the pair files)
inputs: <results dir>/{cohort_1kgp,cohort_hgdp_sgdp,hgsg_calls}.tsv, <results dir>/depth/depth_f1.tsv, results/estimates_1kgp.relatedness_samples.tsv,
results/tree/tree.rooted.nwk, data/catalog_d376/ped.txt, data/figure_metadata/sgdp_metadata.txt, results/figures/pedigree_pairs.tsv
"""
import sys, csv, math, collections as co
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, kruskal
sys.path.insert(0, str(Path(__file__).resolve().parent)); import figures_1kgp as fs, prepare_cohort_pies as pc
D, GN, OUT = Path(sys.argv[1]), sys.argv[2], sys.argv[3]; OWN = sys.argv[4] if len(sys.argv) > 4 else None
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t')); N = []; rng = np.random.default_rng(20261103)
def put(k, v): N.append((k, v)); print(f'{k}\t{v}', flush=True)
b = lambda x, c: float(x[f'{c}_bp']); CL = ['DYZ3', 'Yq12', 'DYZ1', 'DYZ2', 'DYZ18S', 'DYZ19', 'DYZ17']; DEG = (('3rd', 0.0442), ('2nd', 0.0884), ('1st', 0.177))
X = rd(D / 'cohort_1kgp.tsv'); C = {x['sample']: x for x in X}; keep = ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain')
inc = [x for x in X if x['sex'] == 'male' and x['estimate_basis'] == 'ruler' and x['y_state'] in keep]; IN = {x['sample'] for x in inc}
GP = [(r['s1'], r['s2'], float(r['kin']), float(r['ibd0'])) for r in rd(GN)]
if OWN: GP += [(r['s1'], r['s2'], float(r['kin']), float('nan')) for r in rd(OWN)]
def comps(nodes, links):   # union-find, one id per component
    p = {n: n for n in nodes}
    def f(a):
        while p[a] != a: p[a] = p[p[a]]; a = p[a]
        return a
    for a, q in links:
        if a in p and q in p: p[f(a)] = f(q)
    return {n: f(n) for n in nodes}
def one_per(rows, comp):   # first sample of each component
    r = {}
    for x in sorted(rows, key=lambda y: y['sample']): r.setdefault(comp[x['sample']], x)
    return list(r.values())
def within(v, l):
    v = np.array(v, float)
    for k in set(l): v[l == k] -= v[l == k].mean()
    return v
def sdr(rows, lab):   # SD ratio log DYZ1 / log DYZ2 within lineage
    l = np.array(lab); a, c = within(np.log([b(x, 'DYZ1') for x in rows]), l), within(np.log([b(x, 'DYZ2') for x in rows]), l); n = len(a)
    B = [a[i].std() / c[i].std() for i in (rng.integers(0, n, n) for _ in range(2000))]; return f'{n}: ratio {a.std() / c.std():.3f} (95% {np.percentile(B, 2.5):.3f} to {np.percentile(B, 97.5):.3f})'
def hl(a, c):   # Hodges-Lehmann ratio of two groups
    a, c = np.log(a), np.log(c); return math.exp(np.median(np.subtract.outer(a, c)))

# 3: third-degree links in 1kGP
K = {r['sample']: r['kinship_component_id'] for r in rd('results/estimates_1kgp.relatedness_samples.tsv')}
byk = co.defaultdict(list)
for s in IN: byk[K.get(s, s)].append(s)
base = [(v[0], w) for v in byk.values() for w in v[1:]]; c0 = comps(IN, base)
for name, thr in DEG:
    new = [(a, q) for a, q, k, _ in GP if k >= thr and a in IN and q in IN and c0[a] != c0[q]]; put(f'kgp_new_links_{name}', f'{len(new)} pairs above kinship {thr} join different current components: {new[:6]}')
c3 = comps(IN, base + [(a, q) for a, q, k, _ in GP if k >= 0.0442]); R0, R3 = one_per(inc, c0), one_per(inc, c3)
for lab, R in (('current', R0), ('with_3rd_degree', R3)):
    nr = [x for x in R if x['reference_male'] != 'yes' and x['y_state'] != 'par2_gain']
    G = co.defaultdict(list)
    for x in R: G[x['lineage']].append(x)
    LN = sorted(k for k, v in G.items() if len(v) >= 20 and k not in ('.', '?', ''))
    J = [x for x in R if x['lineage'] == 'J']
    put(f'kgp_{lab}', f"{len(R)} representatives; SD ratio outside panel {sdr(nr, [x['lineage'] for x in nr])}; HL DYZ3 R1b/other {hl([b(x, 'DYZ3') for x in R if x['lineage'] == 'R1b'], [b(x, 'DYZ3') for x in R if x['lineage'] != 'R1b']):.3f}; "
        f"HL DYZ19 J2a1a1/other J {hl([b(x, 'DYZ19') for x in J if x['haplogroup'].startswith('J2a1a1')], [b(x, 'DYZ19') for x in J if not x['haplogroup'].startswith('J2a1a1')]):.3f}; Kruskal H DYZ3 {kruskal(*[[b(x, 'DYZ3') for x in G[k]] for k in LN])[0]:.1f} over {len(LN)} lineages")

# 1: relatives among HGDP and SGDP individuals
Q = {x['sample']: x for x in rd(D / 'cohort_hgdp_sgdp.tsv')}; HS = rd(D / 'hgsg_calls.tsv'); AL = co.defaultdict(list)
for r in pc.read_tsv(pc.META / 'sgdp_metadata.txt', encoding='latin-1'):
    for k in pc.ALIAS_FIELDS:
        if pc.valid_key(r[k]): AL[pc.norm(r[k])].append(r)
def keys(s):
    i = s.replace('SGDP_', ''); return {pc.norm(s), pc.norm(i)} | {pc.norm(r[k]) for r in AL.get(pc.norm(i), []) for k in pc.ALIAS_FIELDS if pc.valid_key(r[k])}
own = {x['sample']: x['sample'] for x in HS}; kmap = {k: x['sample'] for x in HS for k in keys(x['sample'])}
gids = {pc.norm(l.strip()) for l in open(Path(GN).parent / 'gnomad_hgdp_1kg_ids.txt')}
oids = {pc.norm(l.strip()) for l in open(Path(OWN).parent / 'sgdp_plink_ids.txt') if l.strip()} if OWN else set()
cov = co.Counter((x['cohort'], 'gnomAD' if keys(x['sample']) & gids else 'own KING' if keys(x['sample']) & oids else 'none') for x in HS); put('hgsg_relatedness_coverage', str(dict(cov)))
L1 = [(kmap[pc.norm(a)], kmap[pc.norm(q)], k) for a, q, k, _ in GP if pc.norm(a) in kmap and pc.norm(q) in kmap and kmap[pc.norm(a)] != kmap[pc.norm(q)]]
for name, thr in DEG: put(f'hgsg_pairs_{name}', f"{len([1 for *_, k in L1 if k >= thr])} pairs: {sorted({(a, q, round(k, 3)) for a, q, k in L1 if k >= thr})[:12]}")
ch = comps(list(own), [(a, q) for a, q, k in L1 if k >= 0.0442])
for lab, R in (('all', HS), ('one_per_component', one_per(HS, ch))):
    rows = [dict(Q[x['sample']], call=x['call']) for x in R]
    for co_ in ('HGDP', 'SGDP'):
        g = [x for x in rows if x['cohort'] == co_]; put(f'hgsg_sd_{co_}_{lab}', sdr(g, [x['call'][:1] for x in g]))
    JJ = [x for x in rows if x['call'].startswith('J')]; sh = [x for x in JJ if b(x, 'DYZ19') < 1e5]
    put(f'hgsg_dyz19_{lab}', f"{len(rows)} individuals; {len(sh)} of {len(JJ)} J below 100 kb, all J2a1a1: {all(x['call'].startswith('J2a1a') for x in sh)}; other lineages below 100 kb: {sum(b(x, 'DYZ19') < 1e5 for x in rows if not x['call'].startswith('J'))}")
    m = lambda f: (len(v := [b(x, 'DYZ3') for x in rows if f(x['call'])]), np.median(v))
    put(f'hgsg_dyz3_{lab}', 'R1b %d %.0f; other %d %.0f; O2a2b1a1 %d %.0f; O2a2b1a2 %d %.0f' % (*m(lambda h: h.startswith('R1b')), *m(lambda h: not h.startswith('R1b')), *m(lambda h: h.startswith('O2a2b1a1')), *m(lambda h: h.startswith('O2a2b1a2'))))

# 2: CV error against Y-tree distance to the nearest training male
T = fs.parse_newick(fs.ROOT / 'results/tree/tree.rooted.nwk'); dep, anc = {}, {}
sys.setrecursionlimit(20000)
def walk(n, d, path):
    d += n['length']; path = path + [id(n)]; dd[id(n)] = d
    if not n['children']: dep[n['name']] = d; anc[n['name']] = path
    for c in n['children']: walk(c, d, path)
dd = {}; walk(T, 0.0, [])
def dist(a, q):
    s = set(anc[a]); l = [x for x in anc[q] if x in s][-1]; return dep[a] + dep[q] - 2 * dd[l]
F = rd(D / 'depth/depth_f1.tsv'); miss = sorted({r['sample'] for r in F} - set(dep)); put('cv_males_not_in_tree', f'{len(miss)}: {miss}')
for c in ('dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17'):
    R = [r for r in F if r['class'] == c and r['sample'] in dep]; fo = co.defaultdict(set)
    for r in R: fo[(r['repeat'], r['fold'])].add(r['sample'])
    e, dm = [], []
    for r in R:
        tr = set().union(*[v for (rp, f_), v in fo.items() if rp == r['repeat'] and f_ != r['fold']]); dm.append(min(dist(r['sample'], t) for t in tr)); e.append(abs(float(r['error_pct'])))
    pm = co.defaultdict(list)   # one point per male: mean over repeats
    for r, x, d in zip(R, e, dm): pm[r['sample']].append((x, d))
    e = np.array([np.mean([x for x, _ in v]) for v in pm.values()]); dm = np.array([np.mean([d for _, d in v]) for v in pm.values()]); q = np.percentile(dm, 10); rho = spearmanr(dm, e)
    put(f'cv_tree_distance_{c}', f'{len(R)} predictions, {len(e)} males (mean over repeats); nearest training male at median {np.median(dm):.2e} substitutions per site (10th percentile {q:.2e}); Spearman |error| vs distance {rho[0]:+.3f} (P {rho[1]:.2g}); '
        f'median of male mean |error| {np.median(e[dm <= q]):.2f}% in the closest 10% and {np.median(e[dm > q]):.2f}% in the rest')

# 4: male siblings in 1kGP
ped = [l.split() for l in open('data/catalog_d376/ped.txt')][1:]; fa = co.defaultdict(list)
for f_, s, father, mother, sex, *_ in ped:
    if father != '0' and s in IN: fa[(father, mother)].append(s)
sib = [(v[i], v[j]) for v in fa.values() for i in range(len(v)) for j in range(i + 1, len(v))]
po = {tuple(sorted((r['father'], r['son']))) for r in rd('results/figures/pedigree_pairs.tsv')}
inf = [(a, q) for a, q, k, z in GP if a in IN and q in IN and 0.177 <= k < 0.354 and z > 0.1 and tuple(sorted((a, q))) not in po]
allsib = sorted({tuple(sorted(p)) for p in sib + inf}); put('sibling_pairs', f'{len(sib)} recorded full-sibling male pairs, {len(inf)} inferred first-degree non-parent pairs, {len(allsib)} unique')
if allsib:
    put('sibling_abs_diff', '; '.join(f"{c} {np.median([abs(100 * (b(C[q], c) / b(C[a], c) - 1)) for a, q in allsib]):.2f}%" for c in CL) + f' ({len(allsib)} pairs)')
with open(OUT, 'w') as o: o.write('item\tvalue\n'); o.writelines(f'{k}\t{v}\n' for k, v in N)
