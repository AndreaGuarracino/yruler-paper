#!/usr/bin/env python3
"""DYZ1:DYZ2 ratio on Yq12 length with a phylogenetic model.

usage: python3 scripts/pgls_composition.py <cohort_1kgp.tsv> <out.tsv>
GLS of log2 ratio on log2 Yq12 on tree v4 with Pagel lambda (ML), one man per kinship component, including flagged Y profiles; also with lineage terms; t test on the slope
"""
import sys, csv, math
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import t as T
sys.path.insert(0, 'scripts'); import figures_1kgp as f
from biology_phylogeny import covariance
COH, OUT = sys.argv[1:3]; rd = lambda p: list(csv.DictReader(open(p), delimiter='\t')); N = []
def put(k, v): N.append((k, v)); print(f'{k}\t{v}', flush=True)
X = rd(COH); keep = ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain'); inc = [x for x in X if x['sex'] == 'male' and x['y_status'] == 'y_found']
kin = {r['sample']: r['kinship_component_id'] for r in rd('results/estimates_1kgp.relatedness_samples.tsv')}; rep = {}
for x in sorted((x for x in inc if x['reference_male'] != 'yes' and x['y_state'] not in ('partial', 'yq_loss', 'yp_deletion', 'par2_gain')), key=lambda r: r['sample']): rep.setdefault(kin.get(x['sample'], x['sample']), x)  # one non-reference male per kinship group (user 2026-11-11: males with reference assemblies out of every population test)
tree = f.parse_newick(f.TREE_SOURCE)
def tips(n, o):
    if not n['children']: o.add(n['name'])
    for c in n['children']: tips(c, o)
    return o
TP = tips(tree, set())
gm = lambda r: (math.log2(float(r['DYZ1_bp'])) + math.log2(float(r['DYZ2_bp']))) / 2
l2 = lambda k: lambda r: math.log2(float(r[k]))
def gls(rows, label, lineage=False, size=l2('Yq12_bp')):
    rows = sorted([r for r in rows if r['sample'] in TP], key=lambda r: r['sample']); n = len(rows); diag, eig, q, _ = covariance(tree, tuple(r['sample'] for r in rows))
    y = np.array([math.log2(float(r['DYZ1_bp']) / float(r['DYZ2_bp'])) for r in rows]); x = np.array([size(r) for r in rows]); M = [np.ones(n), x]
    if lineage:
        L = sorted({r['lineage'] for r in rows}); c = {k: sum(r['lineage'] == k for r in rows) for k in L}
        M += [np.array([r['lineage'] == k for r in rows], float) for k in L if c[k] >= 10][1:]
    A = np.column_stack(M); s = 1 / np.sqrt(diag); Zy = q.T @ (y * s); ZA = q.T @ (A * s[:, None]); ld0 = np.log(diag).sum()
    def fit(lam, full=False):
        d = 1 + lam * (eig - 1); w = 1 / d; XtX = (ZA * w[:, None]).T @ ZA; b = np.linalg.solve(XtX, (ZA * w[:, None]).T @ Zy); r = Zy - ZA @ b; v = (r * r * w).sum() / n; ll = -.5 * (n * (math.log(2 * math.pi) + 1 + math.log(v)) + np.log(d).sum() + ld0)
        return (ll, b, v, XtX) if full else -ll
    g = np.linspace(0, .9999, 41); i = int(np.argmin([fit(a) for a in g])); o = minimize_scalar(fit, bounds=(g[max(0, i - 1)], g[min(40, i + 1)]), method='bounded'); lam = min([(fit(0), 0.), (fit(.9999), .9999), (o.fun, o.x)])[1]
    ll, b, v, XtX = fit(lam, True); se = math.sqrt(v * n / (n - A.shape[1]) * np.linalg.inv(XtX)[1, 1]); df = n - A.shape[1]; q975 = T.ppf(.975, df); p = 2 * T.sf(abs(b[1] / se), df)
    put(label, f'{n} men, lambda {lam:.3f}, slope {b[1]:.3f} (95% {b[1] - q975 * se:.3f} to {b[1] + q975 * se:.3f}), P {p:.3g}, terms {A.shape[1]}')
R = list(rep.values()); gls(R, 'pgls_one_per_component'); gls([r for r in R if r['reference_male'] != 'yes'], 'pgls_outside_panel'); gls(R, 'pgls_one_per_component_lineage_terms', True); gls([r for r in R if r['reference_male'] != 'yes'], 'pgls_outside_panel_lineage_terms', True)
T2 = {r['assembly'].replace('_R2', ''): r for r in rd('data/method/yq12_truth_hallast.tsv')}; A = [dict(sample=s, DYZ1_bp=T2[s]['dyz1_units'], DYZ2_bp=T2[s]['dyz2_units'], Yq12_bp=T2[s]['yq12_span'], lineage=next(x['lineage'] for x in X if x['sample'] == s)) for s in T2 if s in TP and s != 'HG01261']; gls(A, 'pgls_assemblies'); gls(A, 'pgls_assemblies_geomean', size=gm)   # HG01261: Yq12 and DYZ2 targets excluded by the four-track rule; release 2 replacements keyed by donor
O = [r for r in R if r['reference_male'] != 'yes']; gls(O, 'pgls_outside_panel_geomean', size=gm); gls(O, 'pgls_outside_panel_dyz3', size=l2('DYZ3_bp')); gls(O, 'pgls_outside_panel_rest', size=lambda r: math.log2(float(r['Yq12_bp']) - float(r['DYZ1_bp']) - float(r['DYZ2_bp'])))
tree = f.parse_newick(f.ROOT / 'results/tree_biallelic/tree.rooted.nwk'); TP = tips(tree, set()); gls(O, 'pgls_outside_panel_geomean_tree_v3', size=gm)
with open(OUT, 'w') as o: o.write('item\tvalue\n'); o.writelines(f'{k}\t{v}\n' for k, v in N)
