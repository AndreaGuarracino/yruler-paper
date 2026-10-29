# Apply the model to saved counts
# usage: python3 um_apply.py <pass dir> <out.tsv> [model dir]; lengths, Y dose, flags and Y state per sample
import os, sys, glob, subprocess as sp, numpy as np, multiprocessing as mp
P = '/path/to/project'; O = f'{P}/markers_u25'; sys.path.insert(0, f'{P}/scripts'); import um_model, truth_cl
D, OUT = f'{P}/{sys.argv[1]}', sys.argv[2]; MD = sys.argv[3] if len(sys.argv) > 3 else f'{P}/model'; Z = '/path/to/conda/zstd/1.5.7/bin/zstd'
U = np.load(f'{O}/u_keys.npy'); oi = np.load(f'{O}/counts_order_idx.npy'); gc = um_model.gcount(U); fm, M = um_model.load(MD); y0 = fm[f'ysc_{M[0]["class"]}']
LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]
fa = open(f'{O}/markers_u.fa').read().split('\n'); c2 = [fa[i + 1] for i in range(0, len(fa) - 1, 2) if fa[i].startswith('>chr2depth_')]
v = LUT[np.frombuffer(''.join(c2).encode(), np.uint8).reshape(-1, 31)]; k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
kc = np.minimum(k, rc); C2 = np.minimum(np.searchsorted(U, kc), len(U) - 1)
if not np.array_equal(U[C2], kc): sys.exit(f'ERROR {(U[C2] != kc).sum()} chr2 markers not in u_keys')
del fa, c2, v, k, rc, kc
def one(s):   # Y dose and reported lengths of one sample
    try:
        x = np.array(sp.run([Z, '-dcq', f'{D}/{s}/whole/counts.col.zst'], check=True, capture_output=True).stdout.split(), np.float64)
        if len(x) != len(oi): raise ValueError(f'{len(x)} counts, expected {len(oi)}')
        c = np.zeros(len(U)); c[oi] = x; e = um_model.lengths(c, fm, M, gc, y0)
        return s, 2 * c[y0].mean() / c[C2].mean(), e[0] if e else None
    except Exception as e: print(f'WARNING {s}: {e}', flush=True); return s, None, None
S = sorted(f.split('/')[-3] for f in glob.glob(f'{D}/*/whole/counts.col.zst')); print(f'{D}: {len(S)} samples, model {MD}', flush=True)
with mp.Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', 8))) as p: R = p.map(one, S)
ref = set().union(*[set(x) for x in truth_cl.load(O).values()])
ped = {l.split()[1]: l.split() for l in open(f'{P}/yhaplo/ped3202.txt') if not l.startswith('FamilyID')}
hg = dict(l.rstrip('\n').split('\t')[:2] for l in open(f'{O}/lineage.tsv') if l.strip()); bat = {l.split('\t')[0]: l.split('\t')[3] for l in open(f'{P}/lists/pass1kgp.tsv')}
YS = {x[0]: x[2:6] for x in (l.rstrip('\n').split('\t') for l in list(open(f'{P}/results/ystate_1kgp.tsv'))[1:])}   # Yp, Yq, state, ruler bias
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'yq2k7': 'YqHet2.7', 'dyz18sub': 'DYZ18sub', 'dyz19': 'DYZ19'}
cls = [m['class'] for m in M if m['length'] != 'none']; lin = lambda h: '.' if h == '.' else h[:3] if h.startswith('R1') else h[0]; n = {}
with open(OUT, 'w') as o:
    o.write('\t'.join(['sample', 'sex_ped', 'population', 'superpopulation', 'batch', 'haplogroup', 'lineage', 'y_dose'] + [f'{NM[c]}_bp' for c in cls] + ['reference_male', 'flag', 'yp_ruler_dose', 'yq_ruler_dose', 'y_state', 'yp_deletion_bias_pct']) + '\n')
    for s, dose, e in R:
        pd = ped.get(s, ['.'] * 7); h = hg.get(s, '.')
        f = 'female' if pd[4] == '2' else 'no_counts' if dose is None else 'no_estimate' if e is None else 'y_dose_below_0.9' if dose < 0.9 else 'ok'; n[f] = n.get(f, 0) + 1
        o.write('\t'.join([s, pd[4], pd[5], pd[6], bat.get(s, '.'), h, lin(h), 'nan' if dose is None else f'{dose:.4f}'] + [f'{e[c]:.0f}' if e else 'nan' for c in cls] + ['yes' if s in ref else 'no', f] + (YS[s][:3] + [YS[s][3] if YS[s][2] == 'yp_deletion' else '.'] if s in YS else ['.'] * 4)) + '\n')
print(f'written {OUT}: {len(R)} samples, flags {n}', flush=True)
ok = [e for s, d, e in R if e and d >= 0.9 and ped.get(s, ['.'] * 7)[4] != '2']
for c in cls: q = np.percentile([e[c] for e in ok], [25, 50, 75]); print(f'{NM[c]}: {len(ok)} males y_dose >= 0.9, median {q[1]:,.0f} bp (IQR {q[0]:,.0f}-{q[2]:,.0f})', flush=True)
