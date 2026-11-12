#!/usr/bin/env python3
"""Length-extrapolation test (um_eval.py range): per test male, absolute error and the same male's mean CV absolute error.

usage: python3 scripts/range_summary.py <range dir> <range Yq12 dir> <cv_units.tsv> <out units.tsv> <out summary.tsv>
range dir: <class>.tsv rows design, class, male, truth of the sorting class, truth of the class, error %, sorting class (estimate = class truth x (1 + error/100)).
Yq12: the four components fitted on the Yq12-length split (RANGEBY=yq12), summed per male and design, against the Yq12 truth.
Designs: training two thirds of the males sorted by true length, test outside that range.
"""
import csv, sys, glob, os, collections as co, numpy as np
RD, RY, CVU, OUT, SUM = sys.argv[1:6]
NM = {'dyz3': 'DYZ3', 'dyz17': 'DYZ17', 'dyz19': 'DYZ19', 'dyz18s': 'DYZ18S', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2'}
DES = ['short_to_long', 'long_to_short', 'middle_to_extremes']
CV = co.defaultdict(list)
for r in csv.DictReader(open(CVU), delimiter='\t'): CV[(r['class'], r['male'])].append(abs(float(r['error_pct'])))
rows = []
for f in sorted(glob.glob(f'{RD}/*.tsv')):   # single classes
    for l in open(f):
        d, c, s, t, e, err = l.rstrip('\n').split('\t')[:6]
        if c in NM and (c, s) in CV: rows.append((d, NM[c], s, float(e), float(e) * (1 + float(err) / 100), abs(float(err)), float(np.mean(CV[(c, s)]))))
        elif c in NM: print(f'WARNING {c} {s}: no CV units', flush=True)
comp = co.defaultdict(dict)   # Yq12 from its components
for f in sorted(glob.glob(f'{RY}/*.tsv')):
    for l in open(f):
        d, c, s, t, e, err = l.rstrip('\n').split('\t')[:6]; comp[(d, s)][c] = (float(e) * (1 + float(err) / 100), float(t))   # (estimate, Yq12 truth)
need = ('dyz1', 'dyz2', 'dyz18s', 'yq12rest2')
for (d, s), v in sorted(comp.items()):
    if not all(k in v for k in need): print(f'WARNING Yq12 {d} {s}: components {sorted(v)}', flush=True); continue
    est = sum(v[k][0] for k in need); t = v['dyz1'][1]
    if ('yq12', s) not in CV: print(f'WARNING Yq12 {s}: no CV units', flush=True); continue
    rows.append((d, 'Yq12', s, t, est, abs(100 * (est / t - 1)), float(np.mean(CV[('yq12', s)]))))
with open(OUT, 'w') as o:
    o.write('design\tclass\tmale\ttruth_bp\testimate_bp\tabs_error_pct\tcv_mean_abs_error_pct\n')
    o.writelines(f'{d}\t{c}\t{s}\t{t:.0f}\t{e:.0f}\t{a:.4f}\t{m:.4f}\n' for d, c, s, t, e, a, m in rows)
with open(SUM, 'w') as o:
    o.write('class\tdesign\ttest_males\tmedian_abs\tp90_abs\tmax_abs\tcv_median_abs\tcv_p90_abs\n')
    for c in ['DYZ3', 'DYZ17', 'DYZ19', 'Yq12', 'DYZ18S', 'DYZ1', 'DYZ2']:
        for d in DES:
            x = [r for r in rows if r[0] == d and r[1] == c]
            if not x: continue
            a = np.array([r[5] for r in x]); m = np.array([r[6] for r in x])
            o.write(f'{c}\t{d}\t{len(x)}\t{np.median(a):.2f}\t{np.percentile(a, 90):.2f}\t{a.max():.2f}\t{np.median(m):.2f}\t{np.percentile(m, 90):.2f}\n')
            print(f'{c} {d}: {len(x)} test males, |error| median {np.median(a):.2f} p90 {np.percentile(a, 90):.2f} max {a.max():.2f}; same males CV {np.median(m):.2f} / {np.percentile(m, 90):.2f}', flush=True)
print('written', OUT, SUM, flush=True)
