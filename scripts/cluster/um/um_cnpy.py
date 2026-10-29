# Per male: per-marker counts (marker order) -> uint32 array in superset key order
import os, sys, subprocess as sp, numpy as np
P = '/path/to/project'; O = f'{P}/{os.environ.get("UM", "markers_u")}'; YU = os.environ.get('OUTD', 'yruler_u'); Z = '/path/to/conda/zstd/1.5.7/bin/zstd'
s = sys.argv[1] if len(sys.argv) > 1 else sorted(os.listdir(f'{P}/{YU}'))[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]
oi = np.load(f'{O}/counts_order_idx.npy'); U = len(np.load(f'{O}/u_keys.npy', mmap_mode='r'))
v = np.array(sp.run([Z, '-dcq', f'{P}/{YU}/{s}/counts.col.zst'], check=True, capture_output=True).stdout.split(), np.uint32)
if len(v) != len(oi): raise SystemExit(f'ERROR {s}: {len(v)} counts vs {len(oi)} markers')
c = np.zeros(U, np.uint32); c[oi] = v; CO = os.environ.get('CNPYOUT', 'cnpy'); os.makedirs(f'{O}/{CO}', exist_ok=True); np.save(f'{O}/{CO}/{s}.npy', c); print(s, 'markers', len(v), 'sum', int(v.sum(dtype=np.uint64)))
