# Bootstrap intervals, San Roman vs Yruler
# usage: python3 scripts/sanroman_bootstrap.py results/model/sanroman_comparison.tsv results/model/sanroman_bootstrap.tsv
import sys, csv, numpy as np
from scipy.stats import spearmanr
src, out = sys.argv[1:3]; N = 4000
R = list(csv.DictReader((l for l in open(src) if not l.startswith('#')), delimiter='\t'))
cv = lambda a: 100 * np.std(a, ddof=1) / np.mean(a)
rows = []
for c in ('yq12', 'dyz1'):
    x = [r for r in R if r['class'] == c]; t, s, y = (np.array([float(r[k]) for r in x]) for k in ('truth_bp', 'san_roman_relative', 'yruler_final_bp'))
    if len(t) < 10: raise SystemExit(f'ERROR {c}: only {len(t)} males in {src}')
    f = lambda i: (spearmanr(t[i], y[i])[0] - spearmanr(t[i], s[i])[0], cv(t[i] / s[i]) / cv(t[i] / y[i]))
    o = f(np.arange(len(t))); rng = np.random.default_rng(2026); B = np.array([f(rng.integers(0, len(t), len(t))) for _ in range(N)])
    for j, k in enumerate(('spearman_difference', 'cv_ratio')):
        lo, hi = np.percentile(B[:, j], [2.5, 97.5]); rows.append((c, k, len(t), f'{o[j]:.3f}', f'{lo:.3f}', f'{hi:.3f}'))
        print(f'{c} {k}: {o[j]:.3f} (95% {lo:.3f} to {hi:.3f}), {len(t)} males, {N} resamples', flush=True)
with open(out, 'w') as w:
    w.write(f'# Yruler minus San Roman Spearman with assembled length; San Roman over Yruler CV of length/estimate; {N} resamples of males\n')
    w.write('class\tstatistic\tmales\tvalue\tci95_low\tci95_high\n'); w.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
print(f'wrote {out}: {len(rows)} rows', flush=True)
