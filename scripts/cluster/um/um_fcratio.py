# Per-flowcell estimates over the whole-library estimate, current model
# usage: python3 um_fcratio.py <sample> ...; counts from fc_split/<sample>/<flowcell> and pass_1kgp/<sample>/whole
import sys, glob, subprocess as sp, numpy as np
P = '/path/to/project'; O = f'{P}/markers_u25'; sys.path.insert(0, f'{P}/scripts'); import um_model
Z = '/path/to/conda/zstd/1.5.7/bin/zstd'
U = np.load(f'{O}/u_keys.npy'); oi = np.load(f'{O}/counts_order_idx.npy'); gc = um_model.gcount(U); fm, M = um_model.load(); y0 = fm[f'ysc_{M[0]["class"]}']
def est(f):
    x = np.array(sp.run([Z, '-dcq', f], check=True, capture_output=True).stdout.split(), np.float64)
    if len(x) != len(oi): sys.exit(f'ERROR {f}: {len(x)} counts')
    c = np.zeros(len(U)); c[oi] = x; e = um_model.lengths(c, fm, M, gc, y0)
    if e is None: sys.exit(f'ERROR {f}: no G+C curve')
    return e[0], e[2]
for s in sys.argv[1:]:
    W, YW = est(f'{P}/pass_1kgp/{s}/whole/counts.col.zst')
    for d in sorted(glob.glob(f'{P}/fc_split/{s}/*/counts.col.zst')):
        L, Y = est(d); fc = d.split('/')[-2]
        print(f'{s}\t{fc}\truler_sum {Y:.0f}\t' + '\t'.join(f'{c} {L[c] / W[c]:.3f}' for c in ('dyz3', 'yq12', 'dyz19')), flush=True)
