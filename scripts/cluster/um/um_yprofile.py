# Ruler (single-copy Y) coverage on Yp vs Yq, relative to chr2
import sys, gzip, random, numpy as np, multiprocessing as mp
P = '/path/to/project'; D = f'{P}/pass_1kgp'
def prof(s):   # ruler hits per segment / chr2 depth-marker hits
    yp = yq = c2 = 0
    with gzip.open(f'{D}/{s}/whole/bins.tsv.gz', 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); iy, ic = h.index('ysinglecopy_hits'), h.index('chr2depth_hits')
        for l in f:
            x = l.split('\t'); c2 += int(x[ic])
            if x[0] == 'chrY':
                p = int(x[1])
                if p < 10_300_000: yp += int(x[iy])
                elif 10_600_000 <= p < 26_600_000: yq += int(x[iy])
    return s, yp / c2, yq / c2
hi = sys.argv[1].split(','); allm = [l.split('\t')[0].strip() for l in open(sys.argv[2]) if l.strip()]
random.seed(5); ref = random.sample([s for s in allm if s not in hi], 100)
with mp.Pool(8) as p: R = dict((s, (a, b)) for s, a, b in p.map(prof, hi + ref))
mp_, mq = np.median([R[s][0] for s in ref]), np.median([R[s][1] for s in ref])
sp_ = np.percentile([R[s][0] / mp_ for s in ref], [1, 99]); sq = np.percentile([R[s][1] / mq for s in ref], [1, 99])
print(f'reference males 100: Yp ruler dose 1-99% {sp_[0]:.3f}-{sp_[1]:.3f}; Yq ruler dose 1-99% {sq[0]:.3f}-{sq[1]:.3f}')
for s in hi: print(f'{s}: Yp ruler dose {R[s][0] / mp_:.3f}, Yq ruler dose {R[s][1] / mq:.3f}')
