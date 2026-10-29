# Per assembly: superset k-mers present in whole Y and in each class region
import sys, os, subprocess as sp, numpy as np, tempfile
P = '/path/to/project'; O = f'{P}/{os.environ.get("UM", "markers_u")}'; D = '/path/to/scratch/umk/db'; K = '/path/to/conda/kmc/3.2.4/bin'
LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]
def enc(f):   # KMC dump (31-mer, tab, count 1) -> sorted uint64
    b = np.fromfile(f, np.uint8)
    if b.size == 0: return np.zeros(0, np.uint64)
    if b.size % 34: raise SystemExit(f'ERROR {f}: lines not 34 bytes')
    v = LUT[b.reshape(-1, 34)[:, :31]]
    if (v == 255).any(): raise SystemExit(f'ERROR {f}: non-ACGT')
    k = np.zeros(len(v), np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
    return k
U = np.load(f'{O}/u_keys.npy'); a = sys.argv[1] if len(sys.argv) > 1 else open(f'{O}/assemblies.tsv').read().split('\n')[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1].split('\t')[0]
T = tempfile.mkdtemp(dir='/path/to/scratch/umk'); out = {}
for db in sorted(f[:-8] for f in os.listdir(D) if f.startswith(a + '.') and f.endswith('.kmc_pre')):
    name = db[len(a) + 1:]
    sp.run([f'{K}/kmc_tools', '-t2', 'simple', f'{D}/{db}', f'{O}/u_db', 'intersect', f'{T}/x', '-ocleft'], check=True, stdout=sp.DEVNULL)
    sp.run([f'{K}/kmc_tools', '-t2', 'transform', f'{T}/x', 'set_counts', '1', f'{T}/y'], check=True, stdout=sp.DEVNULL)
    sp.run([f'{K}/kmc_tools', '-t2', 'transform', f'{T}/y', 'dump', f'{T}/d.txt'], check=True, stdout=sp.DEVNULL)
    k = enc(f'{T}/d.txt'); i = np.searchsorted(U, k)
    if len(k) and (i.max() >= len(U) or (U[i] != k).any()): raise SystemExit(f'ERROR {a}.{name}: k-mer not in superset')
    m = np.zeros(len(U), bool); m[i] = True; out[name] = np.packbits(m); print(f'{a} {name}: superset k-mers present {m.sum()}', flush=True)
os.makedirs(f'{O}/bits', exist_ok=True); np.savez(f'{O}/bits/{a}.npz', **out); sp.run(['rm', '-rf', T])
