# Parent-offspring differences across the 1kGP sequencing batches and their association with read base quality (Methods, parent-offspring pairs)
# usage: python3 scripts/batch_pairs.py results/figures/pedigree_pairs.tsv results/checks/base_quality_1kgp.tsv results/checks/batch_pairs.tsv
import sys, csv, numpy as np
from scipy import stats
R = list(csv.DictReader(open(sys.argv[1]), delimiter='\t')); Q = {r['sample']: float(r['mean_base_quality']) for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}
rows = []
for g in ('Cross batch', 'Same batch'):
    for lab in ('DYZ3', 'DYZ17', 'DYZ19', 'Yq12', 'DYZ18S', 'DYZ1', 'DYZ2'):
        X = [r for r in R if r['batch_group'] == g and r['label'] == lab]
        d = np.array([float(r['signed_difference_pct']) for r in X]); Y = [r for r in X if r['son'] in Q and r['father'] in Q]
        if len(d) < 20 or len(Y) < 20: sys.exit(f'ERROR {g} {lab}: {len(d)} pairs, {len(Y)} with base quality')
        y = np.array([float(r['signed_difference_pct']) for r in Y]); x = np.array([Q[r['son']] - Q[r['father']] for r in Y]); rho = stats.spearmanr(x, y)
        rows.append((g, lab, len(d), f'{np.median(d):+.2f}', f'{stats.wilcoxon(d).pvalue:.2g}', len(Y), f'{rho.statistic:+.3f}', f'{rho.pvalue:.2g}'))
        print('\t'.join(map(str, rows[-1])), flush=True)
with open(sys.argv[3], 'w') as o:
    o.write('pairs\tclass\tn_pairs\tmedian_offspring_minus_parent_pct\twilcoxon_p_two_sided\tn_with_base_quality\tspearman_rho_vs_base_quality_difference\tspearman_p\n'); o.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
print('written', sys.argv[3], len(rows), 'rows')
