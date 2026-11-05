# Duplicate reads: correlation and removal effect
# usage: python3 scripts/duplicate_check.py results/estimates_1kgp.tsv results/checks/duprate_1kgp.tsv results/checks/dupq_bx.tsv results/checks/duplicate_check.tsv
import sys, csv, collections, numpy as np
from scipy.stats import spearmanr
est, dup, nodup, out = sys.argv[1:5]
E = {r['sample']: r for r in csv.DictReader(open(est), delimiter='\t')}
D = {r['sample']: float(r['dup_fraction']) for r in csv.DictReader(open(dup), delimiter='\t')}
C = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
YQ = {'Yq12', 'DYZ1', 'DYZ2'}   # PAR2 gains out of these
rows = []
# population-summary males, lineage-centred Spearman
for c in C.values():
    S = [r for r in E.values() if r['comparable'] == 'yes' or (r['comparable'] == 'no_Yq12' and c not in YQ)]
    miss = [r['sample'] for r in S if r['sample'] not in D]
    if miss: raise SystemExit(f'ERROR {len(miss)} males without a duplicate fraction, first {miss[:3]}')
    d, v, g = np.array([D[r['sample']] for r in S]), np.log([float(r[f'{c}_bp']) for r in S]), np.array([r['lineage'] for r in S])
    for k in set(g): m = g == k; d[m] -= d[m].mean(); v[m] -= v[m].mean()
    rho = spearmanr(d, v)[0]; rows.append(('within_lineage_spearman', c, len(S), f'{rho:+.3f}', '.', '.'))
    print(f'{c}: {len(S)} males, {len(set(g))} lineages, within-lineage rho {rho:+.3f}', flush=True)
S = [r['sample'] for r in E.values() if r['comparable'] in ('yes', 'no_Yq12')]; f = 100 * np.array([D[s] for s in S])
rows.append(('duplicate_fraction_pct', 'all', len(S), f'{np.median(f):.2f}', f'{np.percentile(f, 10):.2f}', f'{np.percentile(f, 90):.2f}'))
print(f'duplicate fraction: {len(S)} males, median {np.median(f):.2f}%, p10-p90 {np.percentile(f, 10):.2f}-{np.percentile(f, 90):.2f}%', flush=True)
# same reads without duplicates, same model
N = collections.defaultdict(list)
for r in csv.DictReader(open(nodup), delimiter='\t'):
    if r['class'] in C: N[C[r['class']]].append(100 * (float(r['estimate_bp']) / float(E[r['sample']][f"{C[r['class']]}_bp"]) - 1))
for c, x in N.items():
    rows.append(('duplicates_removed_change_pct', c, len(x), f'{np.median(x):+.2f}', f'{min(x):+.2f}', f'{max(x):+.2f}'))
    print(f'{c} without duplicates: {len(x)} males, median {np.median(x):+.2f}% ({min(x):+.2f} to {max(x):+.2f})', flush=True)
with open(out, 'w') as w:
    w.write('# value with low/high: duplicate fraction median p10 p90; removal change median min max; rho has no interval\n')
    w.write('statistic\tclass\tmales\tvalue\tlow\thigh\n'); w.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
print(f'wrote {out}: {len(rows)} rows', flush=True)
