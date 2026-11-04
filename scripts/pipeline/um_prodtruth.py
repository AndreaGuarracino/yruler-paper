# All-copies truth for the production markers
# usage: UM=markers_v6 python3 um_prodtruth.py <model dir> <cand npz> <copies dir> <out.tsv>; truth = annotated truth x (in + out copies) / in copies of the model markers in the male's assembly (as um_eval.py tcopies)
import os, sys, numpy as np
P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_v6")}'; sys.path.insert(0, f'{P}/scripts'); import truth_cl
MD, CA, CO, OUT = sys.argv[1:5]; REG = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7')}   # merged classes, as um_eval.py
A = {}
for l in open(f'{O}/regions.tsv'): a, c = l.split('\t')[:2]; A.setdefault(a.replace('_R2', ''), {})[c] = a
T = truth_cl.load(O); M = np.load(f'{MD}/model.npz'); CAND = np.load(CA); n = 0
with open(OUT, 'w') as o:
    for c in T:
        if f'm_{c}' not in M.files: print(f'{c}: no markers in {MD}, skipped', flush=True); continue
        m = M[f'm_{c}']; p = np.searchsorted(CAND[c], m)
        if not (CAND[c][np.minimum(p, len(CAND[c]) - 1)] == m).all(): sys.exit(f'ERROR {c}: model markers not in the candidates of {CA}')
        for s, v in sorted(T[c].items()):
            a = A.get(s, {}).get(REG.get(c, (c,))[0])
            if not a or not os.path.exists(f'{CO}/{a}.npz'): print(f'WARNING {c} {s}: no copies file, skipped', flush=True); continue
            z = np.load(f'{CO}/{a}.npz'); i, u = z[f'{c}_in'][p].astype(np.int64).sum(), z[f'{c}_out'][p].astype(np.int64).sum()
            if i == 0: print(f'WARNING {c} {s}: none of {len(m)} markers in the region, annotated truth kept', flush=True)
            o.write(f'{c}\t{s}\t{v * (i + u) / i if i else v:.0f}\n'); n += 1
print(f'written {OUT}: {n} truths', flush=True)
