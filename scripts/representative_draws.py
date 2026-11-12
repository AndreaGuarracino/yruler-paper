#!/usr/bin/env python3
"""Main tests with a random individual per kinship group, repeated: range of P values (1kGP, and HGDP and SGDP separately).

usage: python3 scripts/representative_draws.py <cohort_1kgp.tsv> <out.tsv> [draws, default 1000]
Same males and tests as biology_tests.py (all y_found males; lineages with >= 20 males),
but each draw picks one male per component of results/estimates_1kgp.relatedness_samples.tsv at random (seed 2026).
Tests: two-sided Mann-Whitney (DYZ3 R1b vs other, DYZ3 O2a2b1a1 vs O2a2b1a2, DYZ19 J2a1a1 vs other J), Kruskal-Wallis per class and unit ratio,
Pitman-Morgan (SD of log DYZ1 vs log DYZ2 within lineage, non-reference males); ED5b: Kruskal-Wallis of Y dosage by lineage
(males with status ok or low_y_dose, any Y state, NO and R2 split as in Fig. 2a).
HGDP, SGDP: results/model/hgsg_calls.tsv individuals, groups of results/model/hgsg_kinship_representatives.tsv, HaploGrouper calls; the three contrasts and Pitman-Morgan within first-letter lineages, as biology_tests.py.
"""
import csv, re, sys, collections as co
import numpy as np
from scipy.stats import kruskal, mannwhitneyu, pearsonr

COH, OUT = sys.argv[1:3]; ND = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
LIN = lambda h: h[:3] if h.startswith('R1') else (lambda t: t if len(t) > 1 else t[0])(re.match('[A-Z]+', h).group(0))   # lineage of a haplogroup call, as um_apply.py (R1a/R1b split, stem clades kept)
X = rd(COH); keep = ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain'); U1, U2 = 3478.9, 2404.6
inc = [x for x in X if x['sex'] == 'male' and x['y_status'] == 'y_found']
kin = {r['sample']: r['kinship_component_id'] for r in rd('results/estimates_1kgp.relatedness_samples.tsv')}
comp, c5 = co.defaultdict(list), co.defaultdict(list)
for x in inc:
    if x['reference_male'] != 'yes' and x['y_state'] not in ('partial', 'yq_loss', 'yp_deletion', 'par2_gain'): comp[kin.get(x['sample'], x['sample'])].append(x)   # non-reference males only
for x in X:
    if x['sex'] == 'male' and x['estimate_status'] in ('ok', 'low_y_dose') and x['reference_male'] != 'yes' and x['y_state'] not in ('partial', 'yq_loss', 'yp_deletion', 'par2_gain'): c5[kin.get(x['sample'], x['sample'])].append(x)
lin5 = lambda x: 'NO' if x['haplogroup'].split()[0] == 'NO' else 'R2' if x['lineage'] == 'R' and x['haplogroup'].startswith('R2') else x['lineage']
print(f'ED5b set: {sum(map(len, c5.values()))} males, {len(c5)} components', flush=True)
print(f'{len(inc)} males, {len(comp)} components, {sum(len(v) > 1 for v in comp.values())} with more than one male, {ND} draws', flush=True)
b = lambda x, c: float(x[f'{c}_bp']); CL = ['DYZ3', 'Yq12', 'DYZ1', 'DYZ2', 'DYZ18S', 'DYZ19', 'DYZ17']
val = lambda x, c: b(x, 'DYZ1') / b(x, 'DYZ2') * U2 / U1 if c == 'ratio' else b(x, c)
def tests(R):
    P = {}; mw = lambda a, c: mannwhitneyu(a, c, alternative='two-sided').pvalue; hg = lambda x, p: x['haplogroup'].startswith(p)
    P['mw_dyz3_r1b'] = mw([b(x, 'DYZ3') for x in R if x['lineage'] == 'R1b'], [b(x, 'DYZ3') for x in R if x['lineage'] != 'R1b'])
    P['mw_dyz3_o2a2b1a1'] = mw([b(x, 'DYZ3') for x in R if hg(x, 'O2a2b1a1')], [b(x, 'DYZ3') for x in R if hg(x, 'O2a2b1a2')])
    J = [x for x in R if x['lineage'] == 'J']; P['mw_dyz19_j2a1a1'] = mw([b(x, 'DYZ19') for x in J if hg(x, 'J2a1a1')], [b(x, 'DYZ19') for x in J if not hg(x, 'J2a1a1')])
    G = co.defaultdict(list)
    for x in R: G[x['lineage']].append(x)
    LN = sorted(k for k, v in G.items() if len(v) >= 20 and k not in ('.', '?', ''))
    for c in CL + ['ratio']: P[f'kruskal_{c}'] = kruskal(*[[val(x, c) for x in G[k]] for k in LN]).pvalue
    sc = [x for x in R if x['reference_male'] != 'yes']; l = np.array([x['lineage'] for x in sc])
    a, d = (np.log([b(x, k) for x in sc]) for k in ('DYZ1', 'DYZ2'))
    for k in set(l): a[l == k] -= a[l == k].mean(); d[l == k] -= d[l == k].mean()
    P['pitman_morgan_dyz1_dyz2'] = pearsonr(a - d, a + d)[1]; return P, len(R), len(LN)
def ed5b(R):
    G = co.defaultdict(list)
    for x in R: G[lin5(x)].append(float(x['y_dose']))
    return kruskal(*[v for k, v in G.items() if len(v) >= 20 and k not in ('.', '?', '')]).pvalue
Q = {r['sample']: r for r in rd('results/model/cohort_hgdp_sgdp.tsv')}; HC = {r['sample']: r['call'] for r in rd('results/model/hgsg_calls.tsv')}; HG = co.defaultdict(list)
for r in rd('results/model/hgsg_kinship_representatives.tsv'):
    if Q[r['sample']]['y_state'] not in ('partial', 'yq_loss', 'yp_deletion', 'par2_gain'): HG[r['kinship_group']].append(dict(Q[r['sample']], call=HC[r['sample']]))   # no structural Y profiles
print(f"HGDP/SGDP: {sum(map(len, HG.values()))} individuals, {len(HG)} groups, {sum(len(v) > 1 for v in HG.values())} with more than one", flush=True)
def hgsg(R):   # per cohort, never pooled
    P = {}; mw = lambda a, c: mannwhitneyu(a, c, alternative='two-sided').pvalue
    for co_ in ('HGDP', 'SGDP'):
        G = [x for x in R if x['cohort'] == co_]; h = lambda x, p: x['call'].startswith(p); c = co_.lower()
        P[f'{c}_mw_dyz3_r1b'] = mw([b(x, 'DYZ3') for x in G if h(x, 'R1b')], [b(x, 'DYZ3') for x in G if not h(x, 'R1b')])
        P[f'{c}_mw_dyz3_o2a2b1a1'] = mw([b(x, 'DYZ3') for x in G if h(x, 'O2a2b1a1')], [b(x, 'DYZ3') for x in G if h(x, 'O2a2b1a2')])
        P[f'{c}_mw_dyz19_j2a1a1'] = mw([b(x, 'DYZ19') for x in G if h(x, 'J2a1a1')], [b(x, 'DYZ19') for x in G if h(x, 'J') and not h(x, 'J2a1a1')])
        l = np.array([LIN(x['call']) for x in G]); a, d = (np.log([b(x, k) for x in G]) for k in ('DYZ1', 'DYZ2'))
        for k in set(l): a[l == k] -= a[l == k].mean(); d[l == k] -= d[l == k].mean()
        P[f'{c}_pitman_morgan_dyz1_dyz2'] = pearsonr(a - d, a + d)[1]
    return P
rng, rng2 = np.random.default_rng(2026), np.random.default_rng(2027); S = co.defaultdict(list); ns = set()   # rng2: HGDP/SGDP draws, 1kGP draws unchanged
for i in range(ND):
    P, n, nl = tests([v[rng.integers(len(v))] for v in comp.values()]); ns.add((n, nl)); P['ed5b_kruskal_y_dose'] = ed5b([v[rng.integers(len(v))] for v in c5.values()]); P.update(hgsg([v[rng2.integers(len(v))] for v in HG.values()]))
    for k, p in P.items(): S[k].append(p)
    if i % 200 == 0: print(f'draw {i}: ' + ', '.join(f'{k} {p:.3g}' for k, p in P.items()), flush=True)
print('males and lineages per draw', sorted(ns), flush=True)
with open(OUT, 'w') as o:
    o.write('test\tdraws\tmin_p\tmedian_p\tmax_p\tdraws_p_below_0.05\n')
    for k, v in S.items():
        v = np.array(v); r = f'{k}\t{len(v)}\t{v.min():.3g}\t{np.median(v):.3g}\t{v.max():.3g}\t{int((v < .05).sum())}'; o.write(r + '\n'); print(r, flush=True)
