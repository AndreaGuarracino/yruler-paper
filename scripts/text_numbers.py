#!/usr/bin/env python3
"""Numbers of the brief text from the final-model tables.

usage: python3 scripts/text_numbers.py <results dir> <out.tsv>
inputs: <results dir>/{cv,cv_units,lolo,holdout_21,depth,father_son,cv_females,cohort_1kgp,yprofile_1kgp,dyz1_dyz2_lineages,dyz1_dyz2_truth_vs_yruler}.tsv,
manuscript/figures/figureS8_pedigree_pairs.tsv, results/ytree_v4_2026-11-01/ytree_v4.rooted.nwk, data/method/yq12_truth_hallast.tsv
"""
import sys, csv, collections
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent)); import figures_1kgp as fs

D, OUT = Path(sys.argv[1]), sys.argv[2]; ROOT = fs.ROOT; U1, U2 = 3478.9, 2404.4
rows = lambda f: list(csv.DictReader(open(f), delimiter='\t')); N = []
def put(k, v): N.append((k, v)); print(f'{k}\t{v}', flush=True)
st = {(r['class']): r for r in rows(D / 'cv.tsv')}
for c in ('dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17'): put(f'cv_median_{c}', st[c]['median_abs'])
E = collections.defaultdict(dict)
for r in rows(D / 'cv_units.tsv'): E[r['class']][(r['repeat'], r['male'])] = float(r['error_pct'])
R = {k: 100 * ((1 + E['dyz1'][k] / 100) / (1 + E['dyz2'][k] / 100) - 1) for k in E['dyz1'] if k in E['dyz2']}
put('cv_median_ratio', f'{np.median(np.abs(list(R.values()))):.2f}')
pm = lambda d: {m: np.mean([v for (r, mm), v in d.items() if mm == m]) for m in {k[1] for k in d}}   # per-man mean signed error
pma = lambda d: {m: np.mean([abs(v) for (r, mm), v in d.items() if mm == m]) for m in {k[1] for k in d}}
rr = pma(R); put('ratio_per_man_abs_median_p90', f'{np.median(list(rr.values())):.2f} / {np.percentile(list(rr.values()), 90):.2f}')
d19 = pma(E['dyz19']); put('dyz19_short_truth_men', ', '.join(f'{m} {d19[m]:.2f}' for m in ('HG02492', 'HG02735', 'NA20809') if m in d19))
for f, k in (('lolo.tsv', 'lolo'), ('holdout_21.tsv', 'h21')):
    t = {r['class']: r for r in rows(D / f)}; put(f'{k}_median', ', '.join(f"{c} {t[c]['median_abs']} ({t[c]['n_males']})" for c in ('dyz3', 'yq12', 'dyz19')))
t = {(r['depth'], r['class']): r for r in rows(D / 'depth.tsv')}; put('depth_1x_median', ', '.join(f"{c} {t[('d1x', c)]['median_abs']}" for c in ('dyz3', 'yq12', 'dyz19')))
t = {r['class']: r for r in rows(D / 'father_son.tsv') if r['pairs'] == 'all'}; put('father_son_median_signed', ', '.join(f"{c} {t[c]['median_abs']} / {t[c]['signed_median']} ({t[c]['n_pairs']})" for c in ('DYZ3', 'Yq12', 'DYZ19')))
t = {r['class']: r for r in rows(D / 'cv_females.tsv')}; put('female_bg_pct_median', ', '.join(f"{c} {t[c]['median_pct_of_median_truth']}" for c in ('dyz3', 'yq12', 'dyz19')))
put('female_bg_pct_max_any', f"{max(float(r['max_pct_of_median_truth']) for c, r in t.items() if c in ('dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17')):.2f}")
X = rows(D / 'cohort_1kgp.tsv'); ref = {x['sample']: x['reference_male'] == 'yes' for x in X}
inc = [x for x in X if x['estimate_basis'] == 'ruler' and x['sex'] != 'female' and x['y_state'] in ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain', 'par2_gain')]; I = {x['sample']: x for x in inc}; incq = [x for x in inc if x['y_state'] != 'par2_gain']   # incq: Yq12, DYZ1, DYZ2 summaries
val = [x for x in X if x['estimate_status'] in ('ok', 'low_y_dose') and x['sex'] == 'male']; put('men_estimated_included', f'{len(val)} / {len(inc)}')
b = lambda x, c: float(x[f'{c}_bp']); r1 = [b(x, 'DYZ3') for x in inc if x['lineage'] == 'R1b']; ot = [b(x, 'DYZ3') for x in inc if x['lineage'] != 'R1b']
put('dyz3_r1b_vs_other', f'{np.median(r1) / 1e3:.1f} kb ({len(r1)}) vs {np.median(ot) / 1e3:.1f} kb ({len(ot)})')
big = [x for x in val if b(x, 'DYZ3') > 1673237]; put('dyz3_above_longest_assembly_1673237', f'{len(big)} men, {sum(not ref[x["sample"]] for x in big)} outside the panel, max {max(b(x, "DYZ3") for x in big) / 1e6:.2f} Mb ({max(big, key=lambda x: b(x, "DYZ3"))["sample"]})')
J = [x for x in inc if x['lineage'] == 'J']; ja = [x for x in J if x['haplogroup'].startswith('J2a1a1')]; jo = [x for x in J if not x['haplogroup'].startswith('J2a1a1')]
put('dyz19_j2a1a1_vs_other_j', f'{len(ja)}: {min(b(x, "DYZ19") for x in ja) / 1e3:.1f}-{max(b(x, "DYZ19") for x in ja) / 1e3:.1f} kb ({sum(not ref[x["sample"]] for x in ja)} outside the panel) vs {len(jo)}: {min(b(x, "DYZ19") for x in jo) / 1e3:.1f}-{max(b(x, "DYZ19") for x in jo) / 1e3:.1f} kb')
put('dyz19_na18910_na18911', ', '.join(f'{s} {b(I[s], "DYZ19") / 1e3:.1f} kb ({I[s]["haplogroup"]})' for s in ('NA18910', 'NA18911')))
for h in ('O2a2b1a1', 'O2a2b1a2'):
    v = [b(x, 'DYZ3') for x in inc if x['haplogroup'].startswith(h)]; put(f'dyz3_{h}', f'median {np.median(v) / 1e3:.1f} kb, {len(v)} men, range {min(v) / 1e3:.1f}-{max(v) / 1e3:.1f} kb')
def cherries(n, out):   # sister tip pairs among included men
    if not n['children']: return
    k = n['children']
    if len(k) == 2 and not k[0]['children'] and not k[1]['children'] and k[0]['name'] in I and k[1]['name'] in I: out.append((k[0]['name'], k[1]['name']))
    for c in k: cherries(c, out)
P = []; cherries(fs.parse_newick(ROOT / 'results/ytree_v4_2026-11-01/ytree_v4.rooted.nwk'), P)
put('sister_pairs_over_25pct', f'{len(P)} pairs: ' + ', '.join(f"{c} {sum(abs(np.log(b(I[a], c) / b(I[s], c))) > np.log(1.25) for a, s in P if c != 'Yq12' or (I[a]['y_state'] != 'par2_gain' and I[s]['y_state'] != 'par2_gain'))}" for c in ('DYZ19', 'DYZ3', 'Yq12')))
fsr = collections.defaultdict(dict)
for r in rows(ROOT / 'manuscript/figures/figureS8_pedigree_pairs.tsv'): fsr[(r['father'], r['son'])][r['target']] = abs(float(r['signed_difference_pct']))
yq = sorted(((v['yq12'], p, v) for p, v in fsr.items() if 'HG02966' not in p and v['yq12'] > 25), reverse=True); st = {x['sample']: x['y_state'] for x in X}
for nm, g in (('with_par2_gain', [y for y in yq if 'par2_gain' in (st[y[1][0]], st[y[1][1]])]), ('without_par2_gain', [y for y in yq if 'par2_gain' not in (st[y[1][0]], st[y[1][1]])])):
    put(f'father_son_yq12_over_25pct_{nm}', f'{len(g)} pairs: ' + '; '.join(f'{p[0]}-{p[1]} Yq12 {y:.1f}%, DYZ3 {v["dyz3"]:.1f}%, DYZ19 {v["dyz19"]:.1f}%' for y, p, v in g))
put('par2_gain_men', ', '.join(f'{x["sample"]} PAR2 {x["par2_dose"]}, Yq12 {b(x, "Yq12") / 1e6:.2f} Mb' for x in X if x['y_state'] == 'par2_gain'))
Y = {r['sample']: r for r in rows(D / 'yprofile_1kgp.tsv')}; h = Y['HG02966']
put('hg02966_xdr_doses', f"blocks 1-5 {min(float(h[f'xdr{i}_dose']) for i in range(1, 6)):.2f}-{max(float(h[f'xdr{i}_dose']) for i in range(1, 6)):.2f}, blocks 6-8 {min(float(h[f'xdr{i}_dose']) for i in range(6, 9)):.4f}-{max(float(h[f'xdr{i}_dose']) for i in range(6, 9)):.4f}")
put('dyz1_dyz2_spearman_included', f"{spearmanr([b(x, 'DYZ1') for x in incq], [b(x, 'DYZ2') for x in incq])[0]:.2f}")
u = np.array([(b(x, 'DYZ1') / U1) / (b(x, 'DYZ2') / U2) for x in incq]); put('unit_ratio_median_p1_p99', f'{np.median(u):.3f} ({np.percentile(u, 1):.3f}-{np.percentile(u, 99):.3f})')
T = rows(ROOT / 'data/method/yq12_truth_hallast.tsv'); hr = lambda S: np.mean([float(r['dyz1_hits']) / float(r['dyz2_hits']) for r in S])
put('truth_hit_ratio_mean', f"{hr([r for r in T if not r['assembly'].endswith('_R2')]):.4f} ({sum(not r['assembly'].endswith('_R2') for r in T)} Hallast assemblies), {hr(T):.4f} (all {len(T)})")
L = {r['lineage']: r for r in rows(D / 'dyz1_dyz2_lineages.tsv')}; e = L['E']; low = [L[k] for k in ('D', 'I', 'C', 'O', 'R1b')]
lo, hi = min(float(r['unit_ratio_median']) for r in low), max(float(r['unit_ratio_median']) for r in low)
put('lineage_unit_ratio', f"E {float(e['unit_ratio_median']):.3f} ({float(e['ci95_low']):.3f}-{float(e['ci95_high']):.3f}; {e['men']} men); D, I, C, O, R1b {lo:.3f}-{hi:.3f} ({min(int(r['men']) for r in low)}-{max(int(r['men']) for r in low)} men), {100 * (1 - hi / float(e['unit_ratio_median'])):.1f}-{100 * (1 - lo / float(e['unit_ratio_median'])):.1f}% lower")
sg = pm(R); lin = {x['sample']: x['lineage'] for x in X}
put('ratio_signed_error_by_lineage', ', '.join(f'{k} {np.median([v for m, v in sg.items() if lin.get(m) == k]):+.2f}% ({sum(lin.get(m) == k for m in sg)})' for k in ('E', 'O', 'R1b')))
t = {r['measure']: r for r in rows(D / 'dyz1_dyz2_truth_vs_yruler.tsv')}; put('truth_lineage_unit_ratio', ', '.join(f"{k} {t[f'lineage_{k}_truth_unit_ratio']['median']} ({t[f'lineage_{k}_truth_unit_ratio']['men']})" for k in ('E', 'R1b', 'O')))
F = rows(ROOT / 'results/model/depth/depth_f1.tsv'); NMC = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
for c, n in NMC.items():   # included men outside the assembly range
    t = {r['sample']: float(r['reference_bp']) for r in F if r['class'] == c}; lo, hi = min(t.values()), max(t.values()); G = incq if n in ('Yq12', 'DYZ1', 'DYZ2') else inc; v = [b(x, n) for x in G]
    put(f'outside_range_{c}', f'assemblies {lo / 1e3:.1f}-{hi / 1e3:.1f} kb ({len(t)}); {sum(a < lo for a in v)} men below, {sum(a > hi for a in v)} above, of {len(v)}')
d19 = sorted(float(l.split('\t')[4]) for l in list(open(ROOT / 'data/method/dyz19_clean.tsv'))[1:]); put('dyz19_assemblies_below_100kb', f"{sum(a < 1e5 for a in d19)} of {len(d19)} ({min(d19) / 1e3:.1f}-{max(a for a in d19 if a < 1e5) / 1e3:.1f} kb)")
for c in ('dyz18s',):
    g = collections.defaultdict(list)
    for (r, m), e in E[c].items(): g[m].append(e)
    A = [np.mean(v) for m, v in g.items() if m in ref and {x['sample']: x['lineage'] for x in X}[m] == 'A']; put('dyz18s_lineage_A_signed', f'{np.median(A):+.2f}% ({len(A)} men)')
SP = rows(D / 'same_male_pairs.tsv')
for c in ('DYZ19', 'Yq12', 'DYZ2', 'DYZ3'):
    for g in ('SGDP vs 1kGP NYGC', 'SGDP vs HGDP'):
        v = [float(r['diff_pct']) for r in SP if r['class'] == c and r['projects'] == g]; put(f'same_male_{c}_{g}'.replace(' ', '_'), f'{len(v)} pairs, median signed {np.median(v):+.2f}%, median |diff| {np.median(np.abs(v)):.2f}%, max |diff| {np.max(np.abs(v)):.2f}%, {sum(a > 0 for a in v)} positive')
with open(OUT, 'w') as o: o.write('item\tvalue\n'); o.writelines(f'{k}\t{v}\n' for k, v in N)
print('written', OUT, flush=True)
