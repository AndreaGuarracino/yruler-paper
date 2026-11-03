#!/usr/bin/env python3
"""ONT with v9b markers: same CV folds as the short reads, K per class on training males; comparison with short-read v9b estimates.

usage: python3 scripts/ont_v9b_eval.py <sums.tsv> <depth_f1.tsv> <cohort_1kgp.tsv> <lrsc list> <out.tsv>
ONT estimate = K x class sum / ruler sum, K = sum truth / sum relative value on training males; Yq12 = DYZ1 + DYZ2 + rest; groups: hprc R10, lrsc R9, lrsc R10, vienna
"""
import sys, csv, collections as co
import numpy as np
from scipy.stats import spearmanr
SUM, F1, COH, LST, OUT = sys.argv[1:6]; rd = lambda p: list(csv.DictReader(open(p), delimiter='\t')); CH = {l.split('\t')[0]: l.split('\t')[2] for l in open(LST)}
REP = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17']; Rl = co.defaultdict(lambda: co.defaultdict(dict))
def group(r): return 'hprc_R10' if r['cohort'] == 'hprc' else f"lrsc_{CH.get(r['sample'], '?')}" if r['cohort'] == 'lrsc' else 'vienna'
Gc = {}
for r in rd(SUM):
    g = group(r); Rl[g][r['sample']][r['class']] = (float(r['rel_raw']), float(r['rel_gc'])); Gc[g, r['sample']] = int(r['gc_curve'])
for g, d in Rl.items(): print(g, len(d), 'samples; with a G+C curve', sum(Gc[g, s] for s in d), flush=True)
T = co.defaultdict(dict); Fd = {}; SE = {}
for r in rd(F1): T[r['class']][r['sample']] = float(r['reference_bp']); Fd[r['repeat'], r['fold'], r['class'], r['sample']] = 1; SE[r['repeat'], r['class'], r['sample']] = (r['fold'], abs(float(r['error_pct'])))
T['yq12rest'] = {s: T['yq12'][s] - T['dyz1'][s] - T['dyz2'][s] for s in T['yq12'] if s in T['dyz1'] and s in T['dyz2']}; reps = sorted({k[0] for k in SE})
st = lambda a: (len(a), np.mean(a), np.median(a), np.percentile(a, 90), np.max(a)) if len(a) else (0, 0, 0, 0, 0)
rows = []
def cv(g, idx, c):   # one class of one group, all repeats and folds
    d = Rl[g]; base = c if c != 'yq12' else 'yq12'; men = sorted(s for s in d if s in T[base] and (c != 'yq12' or all(k in T and s in T[k] for k in ('dyz1', 'dyz2', 'yq12rest'))))
    if len(men) < 10: return [], [], men
    E, S_ = [], []
    for rp in reps:
        fold = {s: SE[rp, c, s][0] for s in men if (rp, c, s) in SE}; ms = [s for s in men if s in fold]
        for f in sorted(set(fold.values())):
            tr = [s for s in ms if fold[s] != f]; te = [s for s in ms if fold[s] == f]
            if len(tr) < 5: continue
            K = {k: sum(T[k][s] for s in tr) / sum(d[s][k][idx] for s in tr) for k in ('dyz1', 'dyz2', 'yq12rest')} if c == 'yq12' else sum(T[c][s] for s in tr) / sum(d[s][c][idx] for s in tr)
            for s in te:
                est = sum(K[k] * d[s][k][idx] for k in K) if c == 'yq12' else K * d[s][c][idx]
                E.append(100 * abs(est / T[c][s] - 1)); S_.append(SE[rp, c, s][1])
    return E, S_, ms
for g in Rl:
    for idx, nm in ((0, 'raw'), (1, 'gc')):
        for c in REP:
            E, S_, ms = cv(g, idx, c); a = st(E); b = st(S_); rows.append(('cv_vs_truth', g, nm, c, len(ms), f'{a[1]:.2f}', f'{a[2]:.2f}', f'{a[3]:.2f}', f'{a[4]:.2f}', f'{b[1]:.2f}', f'{b[2]:.2f}', f'{b[3]:.2f}', f'{b[4]:.2f}'))
X = {r['sample']: r for r in rd(COH)}; N = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}; keep = ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain')
for g, d in Rl.items():
    men = [s for s in d if s in X and X[s]['sex'] == 'male' and X[s]['y_state'] in keep and X[s]['estimate_basis'] == 'ruler']
    for idx, nm in ((0, 'raw'), (1, 'gc')):
        for c in REP:
            sr = np.array([float(X[s][N[c] + '_bp']) for s in men])
            if c == 'yq12':
                K = {k: np.median([float(X[s]['DYZ1_bp' if k == 'dyz1' else 'DYZ2_bp']) / d[s][k][idx] for s in men]) for k in ('dyz1', 'dyz2')}; rest = np.median([(float(X[s]['Yq12_bp']) - float(X[s]['DYZ1_bp']) - float(X[s]['DYZ2_bp'])) / d[s]['yq12rest'][idx] for s in men]); o = np.array([K['dyz1'] * d[s]['dyz1'][idx] + K['dyz2'] * d[s]['dyz2'][idx] + rest * d[s]['yq12rest'][idx] for s in men])
            else: o = np.array([d[s][c][idx] for s in men]) * np.median([float(X[s][N[c] + '_bp']) / d[s][c][idx] for s in men])
            e = 100 * np.abs(o / sr - 1); rows.append(('vs_short_read', g, nm, c, len(men), f'{e.mean():.2f}', f'{np.median(e):.2f}', f'{np.percentile(e, 90):.2f}', f'{e.max():.2f}', f'{spearmanr(o, sr)[0]:.3f}', '.', '.', '.'))
with open(OUT, 'w') as o:
    o.write('test\tgroup\tsums\tclass\tmen\tmean_abs\tmedian_abs\tp90_abs\tmax_abs\tshort_mean_or_spearman\tshort_median\tshort_p90\tshort_max\n'); o.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
for r in rows: print(*r, flush=True)
