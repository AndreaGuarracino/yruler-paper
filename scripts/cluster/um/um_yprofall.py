# Yp/Yq ruler dose for all males
# usage: python3 um_yprofall.py <pass dir> <out.tsv>
import sys, os, gzip, numpy as np, multiprocessing as mp
D, OUT = sys.argv[1], sys.argv[2]
def prof(s):   # ruler hits per segment / chr2 hits
    f = f'{D}/{s}/whole/bins.tsv.gz'
    if not os.path.exists(f): return s, None
    yp = yq = c2 = 0
    with gzip.open(f, 'rt') as g:
        h = g.readline().rstrip('\n').split('\t'); iy, ic = h.index('ysinglecopy_hits'), h.index('chr2depth_hits')
        for l in g:
            x = l.split('\t'); c2 += int(x[ic])
            if x[0] == 'chrY':
                p = int(x[1])
                if p < 10_300_000: yp += int(x[iy])
                elif 10_600_000 <= p < 26_600_000: yq += int(x[iy])
    return s, (yp / c2, yq / c2) if c2 else None
S = sorted(os.listdir(D))
with mp.Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', 8))) as p: R = {s: v for s, v in p.map(prof, S) if v}
# Reference: males with both ratios > 0.5 x median
mp_, mq = np.median([v[0] for v in R.values()]), np.median([v[1] for v in R.values()])
ref = [s for s, v in R.items() if v[0] > 0.5 * mp_ and v[1] > 0.5 * mq]
mp_, mq = np.median([R[s][0] for s in ref]), np.median([R[s][1] for s in ref])
with open(OUT, 'w') as o:
    o.write('sample\typ_ruler_dose\tyq_ruler_dose\n')
    for s in sorted(R): o.write(f'{s}\t{R[s][0] / mp_:.4f}\t{R[s][1] / mq:.4f}\n')
print(f'samples {len(R)} of {len(S)}; reference {len(ref)}', file=sys.stderr)
