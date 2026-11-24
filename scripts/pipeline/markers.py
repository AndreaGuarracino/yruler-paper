# Marker candidates per assembly
# usage: python3 markers.py outside <regions.tsv> <assembly> <class> <fasta .fai>   (Y sequence outside the class regions, as samtools regions)
#        python3 markers.py bits <universe keys .npy> <universe KMC db> <KMC db dir> <assembly> <out .npz>   (presence of each universe key in whole Y, class and outside dbs)
import os, sys, subprocess as sp, tempfile, numpy as np
MODE = sys.argv[1] if len(sys.argv) > 1 else ''
if MODE == 'outside':   # composite classes: dyz18s = dyz18 + yq3k1, dyz18sub adds yq2k7, dyz18suba adds yq2k7a
    reg, a, c, fai = sys.argv[2:6]; cs = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7'), 'dyz18suba': ('dyz18', 'yq3k1', 'yq2k7a')}.get(c, (c,)); iv = {}
    for l in open(reg):
        x = l.rstrip('\n').split('\t')
        if x[0] == a and x[1] in cs: iv.setdefault(x[2], []).append((int(x[3]), int(x[4])))
    if not iv: sys.exit(0)
    for l in open(fai):
        ctg, n = l.split('\t')[:2]; n = int(n); p = 0
        for s, e in sorted(iv.get(ctg, [])):
            if s > p: print(f'{ctg}:{p + 1}-{s}')
            p = max(p, e)
        if p < n: print(f'{ctg}:{p + 1}-{n}')
elif MODE == 'bits':
    UK, UDB, D, a, OUT = sys.argv[2:7]; K = os.environ['KMC_BIN']
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
    U = np.load(UK); T = tempfile.mkdtemp(dir=os.environ['TMPD']); out = {}
    for db in sorted(f[:-8] for f in os.listdir(D) if f.startswith(a + '.') and f.endswith('.kmc_pre')):
        name = db[len(a) + 1:]
        sp.run([f'{K}/kmc_tools', '-t2', 'simple', f'{D}/{db}', UDB, 'intersect', f'{T}/x', '-ocleft'], check=True, stdout=sp.DEVNULL)
        sp.run([f'{K}/kmc_tools', '-t2', 'transform', f'{T}/x', 'set_counts', '1', f'{T}/y'], check=True, stdout=sp.DEVNULL)
        sp.run([f'{K}/kmc_tools', '-t2', 'transform', f'{T}/y', 'dump', f'{T}/d.txt'], check=True, stdout=sp.DEVNULL)
        k = enc(f'{T}/d.txt'); i = np.searchsorted(U, k)
        if len(k) and (i.max() >= len(U) or (U[i] != k).any()): raise SystemExit(f'ERROR {a}.{name}: k-mer not in the universe')
        m = np.zeros(len(U), bool); m[i] = True; out[name] = np.packbits(m); print(f'{a} {name}: universe keys present {m.sum()}', flush=True)
    if not out: raise SystemExit(f'ERROR {a}: no KMC db in {D}')
    np.savez(OUT, **out); print(f'written {OUT}: {len(out)} sets', flush=True)
else: sys.exit(f'ERROR mode {MODE}: outside or bits')
