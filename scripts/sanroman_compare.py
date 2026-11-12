#!/usr/bin/env python3
"""San Roman 2024 DYZ1 depth (Table S6A, https://doi.org/10.1016/j.xgen.2023.100462) vs assembly truth and yruler, same men: CV, and final model vs one K on all truth men.

usage: python3 scripts/sanroman_compare.py <S6A tsv> <truth regions tsv> <indep tsv> <truth drops> <ped> <yruler cv_units.tsv> <cohort_1kgp.tsv> <out.tsv> [their measure on 30x tsv] [yruler phase 3 table]
DYZ1 depth signal = DYZ1_norm_coverage (DYZ1 depth / depth of a 1 Mb single-copy Y region). Calibration as yruler: K = sum truth / sum signal
on the training men of each fold; 10 repeats of 5-fold CV grouped by family (rng seed 7 + repeat). Errors in % of the assembly length.
"""
import sys, collections
import numpy as np
from scipy.stats import spearmanr

def log(*a): print(*a, flush=True)
S6, REG, IND, DROP, PED, CVU, COH, OUT = sys.argv[1:9]; D30 = sys.argv[9] if len(sys.argv) > 9 and sys.argv[9] != '-' else None   # D30: their measure on our 30x CRAMs (um_srdepth.sh)
P3 = sys.argv[10] if len(sys.argv) > 10 else None   # P3: yruler final model on the 1kGP phase 3 low-coverage runs (um_apply.py table, samples <id>_phase3)
SIG = {'San Roman DYZ1 depth': {l.split('\t')[0]: float(l.split('\t')[3]) for l in list(open(S6))[1:]}}
if D30: SIG['San Roman measure on 30x CRAMs'] = {l.split('\t')[0]: float(l.split('\t')[5]) for l in open(D30)}
drop = {tuple(l.split()) for l in open(DROP) if l.strip()}
T = collections.defaultdict(dict); A12 = {}
for x in (l.rstrip('\n').split('\t') for l in open(REG)):
    if x[1] == 'yq12' and ('yq12', x[0]) not in drop: T['yq12'][x[0]] = float(int(x[6]) - int(x[5])); A12[x[0]] = x[3]
L = [l.rstrip('\n').split('\t') for l in open(IND)]; H = L[0]
for r in (dict(zip(H, x)) for x in L[1:]):
    m = r['assembly'].replace('_R2', '')
    if A12.get(m) == r['assembly'] and ('dyz1', m) not in drop: T['dyz1'][m] = float(r['dyz1_units'])
fam = {x.split()[1]: x.split()[0] for x in open(PED) if x.strip() and not x.startswith('FamilyID')} | {s: 'GIAB_AJ' for s in ('HG002', 'HG003', 'HG004')} | {s: 'GIAB_HAN' for s in ('HG005', 'HG006', 'HG007')}
E = collections.defaultdict(dict)
for l in list(open(CVU))[1:]:
    c, r, m, e = l.rstrip('\n').split('\t'); E[c][(r, m)] = float(e)
X = [l.rstrip('\n').split('\t') for l in open(COH)]; HC = X[0]; X = [dict(zip(HC, x)) for x in X[1:]]
inc = {x['sample']: x for x in X if x['y_status'] == 'y_found' and x['sex'] != 'female' and x['y_state'] in ('normal', 'mosaic_loss', 'mosaic_gain', 'Y_gain')}
def st(v): v = np.abs(np.array(v)); return np.mean(v), np.median(v), np.percentile(v, 90), np.max(v)
rows = []
for (c, col), (sname, sr) in [(x, y) for x in (('yq12', 'Yq12_bp'), ('dyz1', 'DYZ1_bp')) for y in SIG.items()]:
    M = sorted(m for m in T[c] if m in sr); fams = sorted({fam.get(m, m) for m in M}); P = collections.defaultdict(list)
    log(f'{c} / {sname}: truth men {len(T[c])}, with the signal {len(M)}, families {len(fams)}')
    for rep in range(10):   # grouped 5-fold CV
        g = list(fams); np.random.default_rng(7 + rep).shuffle(g); fold = {f: i % 5 for i, f in enumerate(g)}
        for k in range(5):
            tr = [m for m in M if fold[fam.get(m, m)] != k]; K = sum(T[c][m] for m in tr) / sum(sr[m] for m in tr)
            for m in M:
                if fold[fam.get(m, m)] == k: P[m].append(100 * (K * sr[m] / T[c][m] - 1))
    sp = [np.mean(np.abs(P[m])) for m in M]; allp = [e for m in M for e in P[m]]
    K = sum(T[c][m] for m in M) / sum(sr[m] for m in M); full = [100 * (K * sr[m] / T[c][m] - 1) for m in M]   # San Roman K on all truth men (100%)
    FM = [m for m in M if m in inc]; fy = [100 * (float(inc[m][col]) / T[c][m] - 1) for m in FM]   # yruler final model (trained on all truth men)
    fs = [100 * (K * sr[m] / T[c][m] - 1) for m in FM]
    for name, v, n in ((f'{sname}, one K on all truth men', full, len(M)), (f'{sname}, one K on all truth men, men with a yruler final estimate', fs, len(FM)), (f'yruler final model, trained on all truth men (men with {sname})', fy, len(FM))):
        a = st(v); rows.append((c, name, n, *a)); log(f'  {name}: men {n}, mean {a[0]:.2f}% median {a[1]:.2f}% p90 {a[2]:.2f}% max {a[3]:.2f}%; signed median {np.median(v):+.2f}%')
    rf = spearmanr([float(inc[m][col]) for m in FM], [T[c][m] for m in FM])[0]; rs = spearmanr([sr[m] for m in FM], [T[c][m] for m in FM])[0]
    rows.append((c, f'Spearman vs truth: yruler final model (men with {sname})', len(FM), rf, np.nan, np.nan, np.nan))
    for nm, v in ((sname, [T[c][m] / sr[m] for m in FM]), ('yruler final model', [T[c][m] / float(inc[m][col]) for m in FM])):   # scale-free: spread of truth / relative value
        v = np.array(v); rows.append((c, f'scale-free CV of truth / relative value: {nm} (men with {sname})', len(FM), 100 * v.std(ddof=1) / v.mean(), np.nan, np.nan, np.nan)); log(f'  scale-free CV of truth / relative value, {nm}: {100 * v.std(ddof=1) / v.mean():.2f}% (n {len(FM)})')
    rows.append((c, f'Spearman vs truth: {sname}, men with a yruler final estimate', len(FM), spearmanr([sr[m] for m in FM], [T[c][m] for m in FM])[0], np.nan, np.nan, np.nan))
    log(f'  Spearman vs truth on the same {len(FM)} men: yruler final {rf:.4f}, San Roman {rs:.4f}')
    yp = {m: np.mean([e for (r, mm), e in E[c].items() if mm == m]) for m in M if any(mm == m for (r, mm) in E[c])}
    yall = [e for (r, mm), e in E[c].items() if mm in yp]
    rho_sr = spearmanr([sr[m] for m in M], [T[c][m] for m in M])[0]; ym = sorted(yp)
    rho_y = spearmanr([T[c][m] * (1 + yp[m] / 100) for m in ym], [T[c][m] for m in ym])[0]
    co = sorted(m for m in inc if m in sr); rho_c = spearmanr([sr[m] for m in co], [float(inc[m][col]) for m in co])[0]
    for name, v, n in ((f'{sname}, K per fold, per-man mean |error|', sp, len(sp)), (f'{sname}, K per fold, all predictions', allp, len(M)),
                       (f'yruler CV, per-man mean |error| (men with {sname})', [np.mean([abs(e) for (r, mm), e in E[c].items() if mm == m]) for m in yp], len(yp)), (f'yruler CV, all predictions (men with {sname})', yall, len(yp))):
        a = st(v); rows.append((c, name, n, *a)); log(f'  {name}: men {n}, mean {a[0]:.2f}% median {a[1]:.2f}% p90 {a[2]:.2f}% max {a[3]:.2f}%')
    rows.append((c, f'Spearman vs truth: {sname}', len(M), rho_sr, np.nan, np.nan, np.nan)); rows.append((c, f'Spearman vs truth: yruler CV prediction (men with {sname})', len(ym), rho_y, np.nan, np.nan, np.nan))
    rows.append((c, f'Spearman {sname} vs yruler {col} (1kGP men, final model)', len(co), rho_c, np.nan, np.nan, np.nan))
    log(f'  Spearman vs truth: San Roman {rho_sr:.4f} (n {len(M)}), yruler CV {rho_y:.4f} (n {len(ym)}); San Roman vs yruler cohort {rho_c:.4f} (n {len(co)})')
PM, SUM = [], []
if P3:   # same data as the published San Roman values: phase 3 low coverage
    X3 = [l.rstrip('\n').split('\t') for l in open(P3)]; H3 = X3[0]; X3 = {r[0].replace('_phase3', ''): dict(zip(H3, r)) for r in X3[1:] if dict(zip(H3, r))['estimate_status'] == 'ok'}
    sr0 = SIG['San Roman DYZ1 depth']
    for c, col in (('yq12', 'Yq12_bp'), ('dyz1', 'DYZ1_bp')):
        M = sorted(m for m in T[c] if m in X3 and m in sr0); K = sum(T[c][m] for m in M) / sum(sr0[m] for m in M)
        fy = [100 * (float(X3[m][col]) / T[c][m] - 1) for m in M]; fs = [100 * (K * sr0[m] / T[c][m] - 1) for m in M]
        cv = lambda r: 100 * np.std(r, ddof=1) / np.mean(r)   # spread of truth / estimate: no K
        SUM.append(f'# {c}: {len(M)} males; K-free spread of truth/estimate (CV %): San Roman {cv([T[c][m] / sr0[m] for m in M]):.2f}, yruler final {cv([T[c][m] / float(X3[m][col]) for m in M]):.2f}; Spearman vs truth: San Roman {spearmanr([sr0[m] for m in M], [T[c][m] for m in M])[0]:.3f}, yruler final {spearmanr([float(X3[m][col]) for m in M], [T[c][m] for m in M])[0]:.3f}; San Roman bp use one K fitted on these males, yruler bp use the released model\n')
        Ky = sum(T[c][m] for m in M) / sum(float(X3[m][col]) for m in M); fyk = [100 * (Ky * float(X3[m][col]) / T[c][m] - 1) for m in M]   # same treatment as San Roman
        PM += [(m, c, T[c][m], sr0[m], K * sr0[m], fs[i], float(X3[m][col]), fy[i], Ky * float(X3[m][col]), fyk[i]) for i, m in enumerate(M)]   # per male, published 30
        SUM.append(f'# {c}: one K on these {len(M)} males: San Roman {K:.6g} bp per unit, yruler {Ky:.6f} times the final-model bp; all are yruler calibration males (in-sample)\n')
        for name, v in (('yruler final model on phase 3 runs (same men)', fy), ('yruler on phase 3 runs, one K on these men', fyk), ('San Roman DYZ1 depth (phase 3), one K on these men', fs)):
            a = st(v); rows.append((c, name, len(M), *a)); log(f'{c} phase 3: {name}: men {len(M)}, mean {a[0]:.2f}% median {a[1]:.2f}% p90 {a[2]:.2f}% max {a[3]:.2f}%')
        rows.append((c, 'Spearman vs truth: yruler final model on phase 3 runs', len(M), spearmanr([float(X3[m][col]) for m in M], [T[c][m] for m in M])[0], np.nan, np.nan, np.nan))
        rows.append((c, 'Spearman vs truth: San Roman DYZ1 depth (phase 3), same men', len(M), spearmanr([sr0[m] for m in M], [T[c][m] for m in M])[0], np.nan, np.nan, np.nan))
        log(f'  Spearman vs truth: yruler phase 3 {rows[-2][3]:.4f}, San Roman {rows[-1][3]:.4f}')
with open(OUT, 'w') as o:   # phase 3 only: published San Roman vs yruler final, every male
    o.writelines(SUM); o.write('male\tclass\ttruth_bp\tsan_roman_relative\tsan_roman_bp\tsan_roman_error_pct\tyruler_final_bp\tyruler_final_error_pct\tyruler_onek_bp\tyruler_onek_error_pct\n')
    o.writelines(f'{m}\t{c}\t{t:.0f}\t{r:.4f}\t{a:.0f}\t{e:+.2f}\t{y:.0f}\t{f:+.2f}\t{yk:.0f}\t{fk:+.2f}\n' for m, c, t, r, a, e, y, f, yk, fk in sorted(PM, key=lambda x: (x[1], x[0])))
log('written', OUT, len(PM), 'rows (male x class)')
