#!/usr/bin/env python3
"""Numbers of the Online Methods from the final tables.

usage: python3 scripts/methods_numbers.py <results dir> <gc_curves.tsv> <out.tsv>
inputs: <results dir>/{cv,cv_units,lolo,lolo_units,cv_females,cohort_1kgp,cohort_hgdp_sgdp,same_male_pairs,father_son,x_par_doses}.tsv, results/estimates_1kgp.sensitivity.tsv,
results/checks/xpar_kmer.summary.tsv, manuscript/figures/figure1_depth_summary.tsv, figureS8_pedigree_pairs.tsv, data/method/dyz19_clean.tsv, data/method/yq12_truth_hallast.tsv
"""
import sys, csv, collections as co
import numpy as np
from scipy.stats import spearmanr
D, GCF, OUT = sys.argv[1:4]; N = []
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
def put(k, v): N.append((k, v)); print(f'{k}\t{v}', flush=True)
CL = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17']; NM = dict(zip(CL, ['DYZ3', 'Yq12', 'DYZ1', 'DYZ2', 'DYZ18S', 'DYZ19', 'DYZ17']))
st = {r['class']: r for r in rd(f'{D}/cv.tsv')}; lo = {r['class']: r for r in rd(f'{D}/lolo.tsv')}; E = co.defaultdict(dict)
for r in rd(f'{D}/cv_units.tsv'): E[r['class']][(r['repeat'], r['male'])] = float(r['error_pct'])
def rmed(d): g = co.defaultdict(list); [g[k[0]].append(abs(v)) for k, v in d.items()]; return np.mean([np.median(v) for v in g.values()])
for c in CL: t = st[c]; put(f'table_{c}', f"{t['n_males']} & {t['mean_abs']} & {t['median_abs']} & {t['p90_abs']} & {t['max_abs']} & {lo[c]['median_abs']}")
R = {k: 100 * ((1 + E['dyz1'][k] / 100) / (1 + E['dyz2'][k] / 100) - 1) for k in E['dyz1'] if k in E['dyz2']}; a = np.abs(list(R.values()))
L = co.defaultdict(dict)
for r in rd(f'{D}/lolo_units.tsv'): L[r['class']][r['male']] = float(r['error_pct'])
lr = [abs(100 * ((1 + L['dyz1'][m] / 100) / (1 + L['dyz2'][m] / 100) - 1)) for m in L['dyz1'] if m in L['dyz2']]
put('table_ratio', f"{len({k[1] for k in R})} & {a.mean():.2f} & {np.median(a):.2f} & {np.percentile(a, 90):.2f} & {a.max():.2f} & {np.median(lr):.2f}; predictions {len(R)}")
P = rd('manuscript/figures/figure1_depth_summary.tsv'); put('depth_levels', ', '.join(sorted({f"{float(r['median_autosomal_depth']):.1f}" for r in P}, key=float)))
Dp = {(r['depth'], r['class']): r for r in rd(f'{D}/depth.tsv')}; put('depth_low', ', '.join(f"{NM[c]} {Dp[('d1x', c)]['median_abs']}" for c in CL))
F = {r['class']: r for r in rd(f'{D}/cv_females.tsv')}; put('female_bg', '; '.join(f"{c} median {F[c]['median_bp']} max {F[c]['max_bp']}" for c in CL) + f"; predictions {F['dyz3']['n_predictions']}")
X = rd(f'{D}/cohort_1kgp.tsv'); C = {x['sample']: x for x in X}; b = lambda x, c: float(x[f'{c}_bp']); men = [x for x in X if x['sex'] == 'male' and x['estimate_basis'] == 'ruler']
fem = [x for x in X if x['sex'] == 'female']; put('female_max_y', f"{max(float(x['y_dose']) for x in fem):.4f} ({len(fem)})")
fmax = {c: float(F[c]['max_bp']) for c in CL}; ex = [(min(b(x, NM[c]) / fmax[c] for c in CL), x['sample']) for x in men]; ex.sort(); put('male_vs_female_max', f'lowest ratio {ex[0]}, next {ex[1]}; HG02966 classes at or below: ' + str(sum(b(C['HG02966'], NM[c]) <= fmax[c] for c in CL)))
Q = {r['sample']: r for r in rd('results/checks/xpar_kmer.summary.tsv')}; ref = [s for s, r in Q.items() if r['role'] == 'reference']
g = lambda s, c: float(C[s][c]); rg = lambda v: f'{min(v):.3f}-{max(v):.3f}'
put('x_ref', rg([g(s, 'x_dose') for s in ref])); put('x_female_median', f"{np.median([float(x['x_dose']) for x in fem]):.3f}")
for s in ('HG01683', 'HG02869', 'HG03456', 'NA19226'): put(f'sex_{s}', f"X {g(s, 'x_dose'):.3f} Y {g(s, 'y_dose'):.3f} PAR1 {g(s, 'par1_dose'):.3f} PAR2 {g(s, 'par2_dose'):.3f}")
nx = [x for x in men if x['y_state'] == 'normal' and x['sex_chromosomes'] == 'XY']; put('x_normal_xy', f"{len(nx)}: {rg([float(x['x_dose']) for x in nx])}"); put('par_ref', f"PAR1 {rg([g(s, 'par1_dose') for s in ref])}, PAR2 {rg([g(s, 'par2_dose') for s in ref])}")
pc = [(x['sample'], x['par1_dose'], x['par2_dose']) for x in men if x['sex_chromosomes'] == 'XY' and (not .75 <= float(x['par2_dose']) <= 1.25 or not .75 <= float(x['par1_dose']) <= 1.19)]; put('par_changes', str(pc))
inc = [x for x in men if x['y_state'] in ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain', 'par2_gain')]; q = [x for x in inc if x['y_state'] != 'par2_gain']
put('medians', f"n {len(inc)}/{len(q)}: DYZ3 {np.median([b(x, 'DYZ3') for x in inc]):.0f}, DYZ19 {np.median([b(x, 'DYZ19') for x in inc]):.0f}, DYZ17 {np.median([b(x, 'DYZ17') for x in inc]):.0f}, DYZ18S {np.median([b(x, 'DYZ18S') for x in inc]):.0f}, Yq12 {np.median([b(x, 'Yq12') for x in q]):.0f}, DYZ1 {np.median([b(x, 'DYZ1') for x in q]):.0f}, DYZ2 {np.median([b(x, 'DYZ2') for x in q]):.0f}")
med = lambda k, c, S: np.median([b(x, c) for x in S if x['lineage'] == k]); nl = lambda k, S: sum(x['lineage'] == k for x in S)
LN = [k for k in {x['lineage'] for x in inc} if nl(k, inc) >= 20]; h = max(LN, key=lambda k: med(k, 'DYZ3', inc)); put('dyz3_largest_lineage', f'{h} {med(h, "DYZ3", inc):.0f} ({nl(h, inc)}); R1a {med("R1a", "DYZ3", inc):.0f} ({nl("R1a", inc)}); O Yq12 {med("O", "Yq12", q):.0f}; E DYZ3 {med("E", "DYZ3", inc):.0f}, DYZ19 {med("E", "DYZ19", inc):.0f} ({nl("E", inc)})')
ra = np.array([b(x, 'DYZ1') / b(x, 'DYZ2') for x in q]); put('ratio_bp', f"{np.median(ra):.4f} ({np.percentile(ra, 1):.4f}-{np.percentile(ra, 99):.4f}); rho {spearmanr([b(x, 'DYZ1') for x in q], [b(x, 'DYZ2') for x in q])[0]:.4f}")
for s in ('HG04006', 'HG02684', 'HG02686', 'HG01527', 'HG01529', 'HG02965', 'HG02966', 'NA18910', 'NA18911', 'HG01967'): x = C[s]; put(f'man_{s}', f"DYZ3 {b(x, 'DYZ3'):.0f} Yq12 {b(x, 'Yq12'):.0f} DYZ1 {b(x, 'DYZ1'):.0f} DYZ2 {b(x, 'DYZ2'):.0f} DYZ18S {b(x, 'DYZ18S'):.0f} DYZ19 {b(x, 'DYZ19'):.0f} y {g(s, 'y_dose'):.4f} ratio {b(x, 'DYZ1') / b(x, 'DYZ2'):.2f}")
J = [x for x in inc if x['lineage'] == 'J']; ja = [x for x in J if x['haplogroup'].startswith('J2a1a1')]; put('j_split', f"a {sum(x['haplogroup'].startswith('J2a1a1a') for x in ja)} b {sum(x['haplogroup'].startswith('J2a1a1b') for x in ja)}")
Dz = rd('data/method/dyz19_clean.tsv'); cl = [r for r in Dz if r.get('clean', 'True') == 'True']; sz = {r['assembly']: float(r['unit_bp_array']) for r in cl}; put('dyz19_refs', f"{len(cl)}, median {np.median(list(sz.values())):.0f}; " + ', '.join(f'{s} {sz[s]:.0f}' for s in ('HG02492', 'HG02735', 'NA20805', 'NA20809', 'NA20905') if s in sz))
T = rd('data/method/yq12_truth_hallast.tsv'); fr = [float(r['dyz1_units']) / float(r['yq12_span']) for r in T]; put('dyz1_fraction', f'{len(T)}: {100 * min(fr):.1f}-{100 * max(fr):.1f}%')
G = rd(f'{D}/cohort_hgdp_sgdp.tsv'); LIB = co.defaultdict(lambda: co.defaultdict(list))
for x in G:
    if x['cohort'] == 'HGDP' and x['sex'] == 'male' and x['estimate_basis'] == 'ruler' and x['y_state'] == 'normal': LIB[x['population']][x['library']].append(x)
pp = [p for p, d in LIB.items() if len([k for k in d if 'PCR-free' in k]) and sum(len(v) for k, v in d.items() if 'PCR-free' in k) >= 3 and sum(len(v) for k, v in d.items() if 'PCR-free' not in k) >= 3]; put('hgdp_libraries', str({k: sum(len(d[k]) for d in LIB.values()) for k in sorted({k for d in LIB.values() for k in d})}) + f'; populations {len(pp)}')
for c in ('DYZ2', 'DYZ18S', 'Yq12'):
    v = []
    for p in pp:
        d = LIB[p]; pf = [b(x, c) for k, S in d.items() if 'PCR-free' in k for x in S]; pc = [b(x, c) for k, S in d.items() if 'PCR-free' not in k for x in S]; v.append(100 * (np.median(pc) / np.median(pf) - 1))
    put(f'pcr_vs_free_{c}', f'median {np.median(v):+.1f}%, negative in {sum(a < 0 for a in v)} of {len(v)}')
S = rd(f'{D}/same_male_pairs.tsv')
for g2 in ('SGDP vs 1kGP NYGC', 'SGDP vs HGDP'):
    put(f'pairs_{g2}'.replace(' ', '_'), ', '.join(f"{NM[c]} {np.median([abs(float(r['diff_pct'])) for r in S if r['class'] == NM[c] and r['projects'] == g2]):.2f}" for c in CL) + f" ({sum(r['class'] == 'DYZ3' and r['projects'] == g2 for r in S)} pairs)")
Y = co.defaultdict(dict); B = {}
for r in rd('manuscript/figures/figureS8_pedigree_pairs.tsv'): Y[(r['father'], r['son'])][r['target']] = float(r['signed_difference_pct']); B[(r['father'], r['son'])] = r['batch_group']
for g2 in ('Cross batch', 'Same batch'):
    I = [p for p in Y if B[p] == g2]; put(f'signed_{g2}'.replace(' ', '_'), f"{len(I)} pairs: " + ', '.join(f"{c} {np.median([Y[p][c] for p in I]):+.2f}" for c in ('dyz3', 'yq12')))
AL = list(Y); put('fs_median', ', '.join(f"{c} {np.median([abs(Y[p][c]) for p in AL]):.2f}" for c in ('dyz3', 'yq12', 'dyz19', 'dyz1', 'dyz2', 'dyz18s', 'dyz17')) + f"; ratio {np.median([abs(100 * ((1 + Y[p]['dyz1'] / 100) / (1 + Y[p]['dyz2'] / 100) - 1)) for p in AL]):.2f}")
W = {r['sample']: r for r in rd(GCF)}; H = open(GCF).readline().rstrip('\n').split('\t')
col = 'curve_gc8'
cv = {s: float(r[col]) for s, r in W.items() if col and r[col] not in ('', 'nan', '.')} if col else {}
if cv:
    I = [p for p in AL if p[0] in cv and p[1] in cv and all(C[s]['y_state'] == 'normal' for s in p)]
    put('fs_lowgc', f'{len(I)} pairs: ' + ', '.join(f"{c} {spearmanr([cv[p[1]] - cv[p[0]] for p in I], [Y[p][c] for p in I])[0]:+.3f}" for c in ('dyz18s', 'dyz1', 'dyz2', 'dyz17', 'dyz3', 'dyz19')))
with open(OUT, 'w') as o: o.write('item\tvalue\n'); o.writelines(f'{k}\t{v}\n' for k, v in N)
