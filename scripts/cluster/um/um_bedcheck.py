# BED vs whole-file estimates of the model
# usage: [MDIR=<model dir>] python3 um_bedcheck.py <pass dir>; writes <model dir>/bedcheck_<pass>.tsv
import os, sys, glob, subprocess as sp, numpy as np, multiprocessing as mp
P = '/path/to/project'; O = f'{P}/markers_u25'; D = f'{P}/{sys.argv[1] if len(sys.argv) > 1 else "pass_1kgp"}'; Z = '/path/to/conda/zstd/1.5.7/bin/zstd'
U = np.load(f'{O}/u_keys.npy'); oi = np.load(f'{O}/counts_order_idx.npy'); sys.path.insert(0, f'{P}/scripts'); import um_model
MD = os.environ.get('MDIR', f'{P}/model'); fm, FM = um_model.load(MD); gc = um_model.gcount(U); y0 = fm[f'ysc_{FM[0]["class"]}']; M = [m for m in FM if m['length'] != 'none']   # load before fork: a shared npz handle corrupts reads
def est(f):   # reported lengths from one counts file
    v = np.array(sp.run([Z, '-dcq', f], check=True, capture_output=True).stdout.split(), np.float64)
    if len(v) != len(oi): raise SystemExit(f'ERROR {f}: {len(v)} counts')
    c = np.zeros(len(U)); c[oi] = v; e = um_model.lengths(c, fm, FM, gc, y0)
    if e is None: raise ValueError(f'{f}: no G+C curve (no Y)')
    return e[0]
def one(s):
    try: w, b = est(f'{D}/{s}/whole/counts.col.zst'), est(f'{D}/{s}/bed/counts.col.zst'); return s, {k: 100 * (b[k] / w[k] - 1) for k in w}
    except Exception as e: print(f'WARNING {s}: {e}', flush=True); return s, None
S = sorted(os.path.basename(os.path.dirname(os.path.dirname(f))) for f in glob.glob(f'{D}/*/bed/counts.col.zst') if os.path.exists(f.replace('/bed/', '/whole/')))
print(f'samples with whole and BED counts: {len(S)}', flush=True)
with mp.Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', 8))) as p: R = [x for x in p.map(one, S) if x[1]]
with open(f'{MD}/bedcheck_{os.path.basename(D)}.tsv', 'w') as o:
    o.write('sample\t' + '\t'.join(m['class'] for m in M) + '\n')
    for s, e in R: o.write(s + '\t' + '\t'.join(f'{e[m["class"]]:.4f}' for m in M) + '\n')
for m in M:
    x = np.array([e[m['class']] for _, e in R]); ok = np.isfinite(x); a = np.abs(x[ok]); nm = [R[i][0] for i in np.flatnonzero(ok)]   # no-Y samples give 0 / 0
    print(f'{m["class"]}: n {ok.sum()} (not finite {(~ok).sum()}), BED/whole - 1 median {np.median(x[ok]):+.3f}%, median |d| {np.median(a):.3f}%, p99 |d| {np.percentile(a, 99):.3f}%, max |d| {a.max():.3f}% ({nm[int(a.argmax())]})', flush=True)
