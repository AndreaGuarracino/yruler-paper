# Marker order of counts.tsv -> index into superset keys (once)
import numpy as np, sys
sys.argv = ['x']; exec(open('/path/to/project/scripts/um_eval.py').read().split('fa = open(')[0])
order = [l.strip() for l in open(f'{O}/counts_order.txt')]; k = enc_lines(order); oi = np.searchsorted(U, k)
if not (U[oi] == k).all(): raise SystemExit('ERROR counts order not in superset')
np.save(f'{O}/counts_order_idx.npy', oi); print('order', len(oi))
