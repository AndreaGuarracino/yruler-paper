# Other technologies and same-male pairs
# usage: python3 um_xtech.py [model dir]; needs cohort_1kgp.tsv, cohort_hgdp_sgdp.tsv, other_runs.tsv from um_apply.py; Element runs left out; father-son pairs of the 1kGP pedigree
import os, sys, numpy as np
P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_v9")}'; X = f'{P}/truth_qc/x'; sys.path.insert(0, f'{P}/scripts'); import truth_cl
os.environ.setdefault('D17T', f'{O}/dyz17_truth.tsv'); os.environ.setdefault('D17C', 'dyz17_span')   # same truth as this model's CV
D = sys.argv[1] if len(sys.argv) > 1 else f'{P}/model'
T = truth_cl.load(O); print(f'truth from {O}, D17T {os.environ["D17T"]}, D17C {os.environ["D17C"]}: ' + ', '.join(f'{c} {len(v)}' for c, v in T.items()), flush=True)
for l in open(f'{X}/rule_w4.txt'): c, s = l.split()[:2]; T[c].pop(s, None)
if os.environ.get('TRUTHSET'):   # TRUTHSET: class, male, bp (um_prodtruth.py)
    for l in open(os.environ['TRUTHSET']):
        c, s, v = l.split()[:3]
        if c in T and s in T[c]: T[c][s] = float(v)
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'yq2k7': 'YqHet2.7', 'dyz18sub': 'DYZ18sub', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
H0 = open(f'{D}/cohort_1kgp.tsv').readline().rstrip('\n').split('\t'); NM = {c: n for c, n in NM.items() if f'{n}_bp' in H0}   # classes present in the tables
def rd(f):   # sample -> (y_dose, {class: bp}, y_status)
    R = [l.rstrip('\n').split('\t') for l in open(f'{D}/{f}')]; H = R[0]
    return {r[0]: (float(r[H.index('y_dose')]), {c: float(r[H.index(f'{n}_bp')]) for c, n in NM.items()}, r[H.index('y_status')] if 'y_status' in H else 'y_found') for r in R[1:] if r[H.index(f'{next(iter(NM.values()))}_bp')] != 'nan'}
st = lambda a: f'{np.mean(np.abs(a)):.2f}\t{np.median(np.abs(a)):.2f}\t{np.percentile(np.abs(a), 90):.2f}\t{np.max(np.abs(a)):.2f}\t{np.median(a):+.2f}'
K, H, L = rd('cohort_1kgp.tsv'), rd('cohort_hgdp_sgdp.tsv'), rd('other_runs.tsv')
CEPH = {'2281': 'NA12889', '2214': 'NA12891', '2209': 'NA12877', '2212': 'NA12882', '2189': 'NA12886'}; GR = lambda r: 'phase 3 low coverage (PCR)' if 'phase3' in r else 'other technologies (no Element, no phase 3)'; E = {g: {c: [] for c in NM} for g in ('other technologies (no Element, no phase 3)', 'phase 3 low coverage (PCR)')}   # PCR runs reported apart
with open(f'{D}/other_technologies.tsv', 'w') as o:
    o.write('run\tmale\tclass\ttruth_bp\testimate_bp\terror_pct\n'); runs = sorted(r for r in L if 'element' not in r.lower())
    print(f'other-technology runs {len(L)}, without Element {len(runs)}: {" ".join(runs)}', flush=True)
    for run in runs:
        m = CEPH.get(run.rsplit('_', 1)[0], run.rsplit('_', 1)[0])
        for c in NM:
            if m in T[c] and T[c][m] > 0: e = 100 * (L[run][1][c] / T[c][m] - 1); E[GR(run)][c].append(e); o.write(f'{run}\t{m}\t{NM[c]}\t{T[c][m]:.0f}\t{L[run][1][c]:.0f}\t{e:+.2f}\n')
PR = [(s, s[5:], 'SGDP vs 1kGP NYGC', K) for s in H if s.startswith('SGDP_') and s[5:] in K] + [(s, s[5:], 'SGDP vs HGDP', H) for s in H if s.startswith('SGDP_') and s[5:] in H]
with open(f'{D}/same_male_pairs.tsv', 'w') as o:
    o.write('sample_a\tsample_b\tprojects\tclass\tbp_a\tbp_b\tdiff_pct\n')
    for a, b, g, B in PR:
        if min(H[a][0], B[b][0]) < 0.9: print(f'pair {a} {b}: y_dose < 0.9, left out', flush=True); continue
        for c in NM: o.write(f'{a}\t{b}\t{g}\t{NM[c]}\t{H[a][1][c]:.0f}\t{B[b][1][c]:.0f}\t{100 * (H[a][1][c] / B[b][1][c] - 1):+.2f}\n')
with open(f'{D}/xtech_summary.tsv', 'w') as o:
    o.write('test\tclass\tn\tmean_abs_pct\tmedian_abs_pct\tp90_abs_pct\tmax_abs_pct\tsigned_median_pct\n')
    for g in E:
        for c in NM:
            if E[g][c]: o.write(f'{g}\t{NM[c]}\t{len(E[g][c])}\t{st(np.array(E[g][c]))}\n')
    R = [l.rstrip('\n').split('\t') for l in list(open(f'{D}/same_male_pairs.tsv'))[1:]]
    for g in ('SGDP vs 1kGP NYGC', 'SGDP vs HGDP'):
        for c in NM.values():
            a = np.array([float(r[6]) for r in R if r[2] == g and r[3] == c])
            if len(a): o.write(f'{g}\t{c}\t{len(a)}\t{st(a)}\n')
print(open(f'{D}/xtech_summary.tsv').read(), flush=True)
FS = [l.split() for l in list(open(f'{P}/lists/1kGP.3202_samples.pedigree_info.txt'))[1:]]; FS = [(x[0], x[1]) for x in FS if x[1] != '0' and x[0] in K and x[1] in K and K[x[0]][2] == K[x[1]][2] == 'y_found']; FD = [(s, f) for s, f in FS if min(K[s][0], K[f][0]) >= 0.9]   # father-son pairs: all, both y_dose >= 0.9
with open(f'{D}/father_son.tsv', 'w') as o:
    o.write('pairs\tclass\tn_pairs\tmean_abs\tmedian_abs\tp90_abs\tmax_abs\tsigned_median\n')
    for (pn, PS), c in [(q, c) for q in (('all', FS), ('y_dose_ge_0.9', FD)) for c in list(NM) + ['dyz1:dyz2']]:
        v = np.array([100 * ((K[s][1]['dyz1'] / K[s][1]['dyz2']) / (K[f][1]['dyz1'] / K[f][1]['dyz2']) - 1) if c == 'dyz1:dyz2' else 100 * (K[s][1][c] / K[f][1][c] - 1) for s, f in PS])
        o.write(f'{pn}\t{NM.get(c, "DYZ1:DYZ2")}\t{len(v)}\t{st(v)}\n')
print(open(f'{D}/father_son.tsv').read(), flush=True)
