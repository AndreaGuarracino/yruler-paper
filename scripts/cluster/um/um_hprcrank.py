# Rank HPRC X haplotypes; quality subsets
import csv, os
P = '/path/to/project'; Q = f'{P}/openhgl/hprcqc'; O = f'{P}/markers_u25'
q = {r['haplotype']: r for r in csv.DictReader(open(f'{Q}/hprc_x_quality.tsv'), delimiter='\t')}
t = {r['haplotype']: r for r in csv.DictReader(open(f'{Q}/hprc_x_t2t_check.tsv'), delimiter='\t')}
F = [r for r in q.values() if r['sex'] == 'F']; fmax = {k: max(float(r[k]) for r in F) for k in ('dyz3_present', 'yq12_present', 'dyz1_present', 'dyz19_present')}
R = []
for h, r in q.items():
    if h not in t: continue
    yok = all(float(r[k]) <= v for k, v in fmax.items())   # no more Y satellite than any female
    R.append((h, r['sex'], int(t[h]['t2t_ks']), int(r['n_gaps']), int(r['flagger_hifi_bp']) + int(r['flagger_ont_bp']), int(r['nucflag_bp']), yok))
R.sort(key=lambda x: (-x[2], x[4], x[5]))   # verified T2T, then Flagger, then NucFlag bp
with open(f'{Q}/hprc_x_rank.tsv', 'w') as o:
    o.write('rank\thaplotype\tsex\tt2t_verified\tn_gaps\tflagger_bp\tnucflag_bp\ty_satellite_ok\n')
    for i, x in enumerate(R, 1): o.write('\t'.join(map(str, (i,) + x)) + '\n')
print('male X haplotypes failing the Y-satellite check:', sum(1 for x in R if x[1] == 'M' and not x[6]), 'of', sum(1 for x in R if x[1] == 'M'), '| female max present', {k: round(100 * v, 3) for k, v in fmax.items()})
for s in ('F', 'M'):
    v = sorted(x[2] for x in R if x[1] == s); print(s, 'verified T2T chromosomes: median', v[len(v) // 2], 'min', v[0], 'max', v[-1], '| >= 20:', sum(1 for x in v if x >= 20), '| >= 18:', sum(1 for x in v if x >= 18))
T = int(os.environ.get('T2TMIN', 12))   # top quality: verified T2T chromosomes
S = {'pg_topF': [x[0] for x in R if x[1] == 'F' and x[2] >= T], 'pg_topFM': [x[0] for x in R if x[2] >= T and x[6]]}
for d, hs in S.items():
    os.makedirs(f'{O}/{d}', exist_ok=True)
    for h in hs:
        if not os.path.exists(f'{O}/{d}/{h}.txt'): os.symlink(f'{O}/pg_nony/{h}.txt', f'{O}/{d}/{h}.txt')
    print(d, 'haplotypes', len(hs))
