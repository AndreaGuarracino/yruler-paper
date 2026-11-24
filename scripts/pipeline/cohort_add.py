# Cohort tables, second stage: rows from whole-file recounts, pangenome X/PAR doses, merge
# usage: python3 cohort_add.py check <cohort tsv> <n> <out.tsv> | apply <list: sample, cohort, URL> <out.tsv> (MD=<model dir>, YPQ=<ruler Yp/Yq blocks>, from $P)
#        python3 cohort_add.py xpar <xs72_doses.tsv> <1kGP table> <in.tsv> <out.tsv> [<in.tsv> <out.tsv> ...] | merge <table> <rows> <out.tsv>
#   check: shared samples, recount values against the table; apply: rows in the cohort_hgdp_sgdp.tsv format for CRAMs counted whole (final1 recount)
#   xpar: x_dose, par1_dose, par2_dose from the pangenome X and PAR rulers over the chr8:91-101 Mb window (xs72_doses.tsv), divided by the median of the 1kGP XY males with y_state normal
#   merge: table rows plus new rows, sorted by sample; a sample in both stops the run
import os, sys, csv, subprocess as sp, numpy as np, multiprocessing as mp
P = os.environ['P']; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import model as um_model
MODE, A, OUT = sys.argv[1], sys.argv[2], sys.argv[-1]; IN = os.environ['IN']; RC = f'{IN}/counts/whole_recount'; Z = 'zstd'; MD = os.environ['MD']
YDN, NOY = 0.9617, 0.025; XN = {'x': 0.96768, 'par1': 0.98332, 'par2': 0.95907}   # 1kGP normal-male medians of the raw doses (xpar mode checks them)
def log(*a): print(*a, flush=True)
if MODE == 'xpar':
    X = {r['sample']: r for r in csv.DictReader(open(A), delimiter='\t')}
    def raw(s): w = float(X[s]['depth1_mean']); return np.array([2 * float(X[s]['xruler_mean']) / w, float(X[s]['parone_mean']) / w, float(X[s]['partwo_mean']) / w])
    ref = [r['sample'] for r in csv.DictReader(open(sys.argv[3]), delimiter='\t') if r['sex'] == 'male' and r['y_status'] == 'y_found' and r['y_state'] == 'normal']
    med = np.median([raw(s) for s in ref], axis=0); log(f'reference 1kGP XY males with y_state normal {len(ref)}; raw medians X {med[0]:.5f} PAR1 {med[1]:.5f} PAR2 {med[2]:.5f}')
    if [round(float(v), 5) for v in med] != list(XN.values()): sys.exit(f'ERROR medians {med} differ from XN {XN} (apply rows use XN)')
    for i in range(4, len(sys.argv), 2):
        R = list(csv.reader(open(sys.argv[i]), delimiter='\t')); h = R[0]; ix = [h.index(c) for c in ('x_dose', 'par1_dose', 'par2_dose')]; miss = [r[0] for r in R[1:] if r[0] not in X]
        if miss: sys.exit(f'ERROR {sys.argv[i]}: {len(miss)} samples without pangenome doses: {miss[:5]}')
        for r in R[1:]:
            for j, v in zip(ix, raw(r[0]) / med): r[j] = f'{v:.4f}'
        with open(sys.argv[i + 1], 'w') as o: o.write(''.join('\t'.join(r) + '\n' for r in R))
        log(f'{sys.argv[i]} -> {sys.argv[i + 1]}: {len(R) - 1} rows, X and PAR doses replaced')
    sys.exit(0)
if MODE == 'merge':
    T, N = [list(csv.reader(open(f), delimiter='\t')) for f in (A, sys.argv[3])]
    if T[0] != N[0]: sys.exit('ERROR headers differ')
    d = {r[0] for r in T[1:]} & {r[0] for r in N[1:]}
    if d: sys.exit(f'ERROR samples in both: {sorted(d)[:5]}')
    with open(OUT, 'w') as o: o.write(''.join('\t'.join(r) + '\n' for r in [T[0]] + sorted(T[1:] + N[1:], key=lambda r: r[0])))
    log(f'merge: {len(T) - 1} + {len(N) - 1} rows -> {OUT}'); sys.exit(0)
K = np.load(f'{IN}/counts/whole_recount_keys.npy', mmap_mode='r'); U = np.load(f'{IN}/markers/u_keys.npy')
def idx(X, name):   # positions of keys in final1, all required
    X = np.asarray(X, np.uint64); j = np.searchsorted(K, X); j[j >= len(K)] = 0
    if not (K[j] == X).all(): sys.exit(f'ERROR {name}: {int((K[j] != X).sum())} keys missing from final1')
    return j
iu = np.searchsorted(K, U); iu[iu >= len(K)] = 0; hit = K[iu] == U
fm, M = um_model.load(MD); y0 = fm[f'ysc_{M[0]["class"]}']; gc = um_model.gcount(U); YPQ = np.load(os.environ['YPQ'])
if not hit[y0].all(): sys.exit('ERROR ruler keys missing from final1')
iw = idx(np.load(f'{IN}/windows/chr8_keys.npy'), 'g72 window'); RZ = np.load(f'{IN}/windows/xpar_rulers.npz'); IX = {k: idx(RZ[k], k) for k in ('x', 'par1', 'par2')}
log(f'final1 {len(K)} keys; universe in final1 {int(hit.sum())}; ruler {len(y0)}; window {len(iw)}; X {len(IX["x"])}, PAR1 {len(IX["par1"])}, PAR2 {len(IX["par2"])}')
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz18s': 'DYZ18S', 'dyz18sub': 'DYZ18sub', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17', 'dyz2': 'DYZ2', 'yq2k7': 'YqHet2.7'}
def one(s):   # doses and lengths of one sample
    x = np.frombuffer(sp.run([Z, '-dcq', f'{RC}/{s}/counts.u32.zst'], check=True, capture_output=True).stdout, np.uint32).astype(np.float64)
    if len(x) != len(K): raise SystemExit(f'ERROR {s}: {len(x)} counts vs {len(K)} keys')
    c = np.zeros(len(U)); c[hit] = x[iu[hit]]; d = x[iw].mean(); yd = 2 * c[y0].mean() / d / YDN
    v = {'y_dose': yd, 'yp_ruler_dose': 2 * c[y0[(YPQ == 1) | (YPQ == 2)]].mean() / d / YDN, 'yq_ruler_dose': 2 * c[y0[YPQ >= 3]].mean() / d / YDN,
         'x_dose': 2 * x[IX['x']].mean() / d / XN['x'], 'par1_dose': x[IX['par1']].mean() / d / XN['par1'], 'par2_dose': x[IX['par2']].mean() / d / XN['par2']}
    e = um_model.lengths(c, fm, M, gc, y0, flat=True); e = e[0] if e else {}
    for k, n in NM.items(): v[f'{n}_bp'] = e.get(k, float('nan')) if round(yd, 4) > NOY else float('nan')
    return s, v
def state(v):   # um_ystate.py rules
    y, p, q, a, b = v['y_dose'], v['yp_ruler_dose'], v['yq_ruler_dose'], v['par1_dose'], v['par2_dose']
    if round(y, 4) <= NOY: return 'no_Y'
    if p < 0.5 and abs(q - 1) <= 0.1: return 'yp_deletion'
    if abs(p - 1) <= 0.1 and q < p - 0.1: return 'yq_loss'
    if abs(p - q) > 0.1: return 'partial'
    if b >= 1.25 and abs(a - 1) <= 0.1: return 'par2_gain'
    return 'Y_gain' if y >= 1.75 else 'mosaic_gain' if y >= 1.25 else 'mosaic_loss' if y < 0.9 else 'normal'
F = ['y_dose', 'yp_ruler_dose', 'yq_ruler_dose', 'x_dose', 'par1_dose', 'par2_dose'] + [f'{n}_bp' for n in NM.values()]
NP = int(os.environ.get('SLURM_CPUS_PER_TASK', 8))
if MODE == 'check':
    T = [r for r in csv.DictReader(open(A), delimiter='\t') if r['y_dose'] not in ('nan', 'NA', '') and os.path.exists(f'{RC}/{r["sample"]}/done')]
    T = T[:: max(1, len(T) // int(sys.argv[3]))][:int(sys.argv[3])]; log(f'check {len(T)} shared samples: {T[0]["sample"]} ... {T[-1]["sample"]}')
    with mp.Pool(NP) as p: Q = dict(p.map(one, [r['sample'] for r in T]))
    with open(OUT, 'w') as o:
        o.write('sample\tfield\ttable\trecount\tdiff_pct\n')
        for r in T:
            for f in F:
                a, b = r.get(f, 'nan'), Q[r['sample']][f]
                if a not in ('nan', 'NA', '.', '') and b == b and float(a) != 0: o.write(f"{r['sample']}\t{f}\t{a}\t{b:.4f}\t{100 * (b / float(a) - 1):.3f}\n")
    D = {}
    for l in list(open(OUT))[1:]: x = l.rstrip('\n').split('\t'); D.setdefault(x[1], []).append(float(x[4]))
    for f in F:
        if f in D: z = np.array(D[f]); log(f'{f}: {len(z)} samples, recount vs table median {np.median(z):+.3f}%, max |.| {np.abs(z).max():.3f}%')
    st = sum(state(Q[r['sample']]) == r['y_state'] for r in T if r.get('y_state') not in (None, '.', '')); log(f'y_state equal in {st} of {sum(1 for r in T if r.get("y_state") not in (None, ".", ""))}')
else:
    L = [l.rstrip('\n').split('\t')[:3] for l in open(A) if l.strip()]; L = [x for x in L if os.path.exists(f'{RC}/{x[0]}/done')]; log(f'apply {len(L)} samples with a final1 recount')
    SX = {x[2]: x[3] for f in ('PRJEB6463_sex.tsv', 'PRJEB9586_sex.tsv') for x in (l.rstrip('\n').split('\t') for l in list(open(f'{IN}/lists/{f}'))[1:])}
    with mp.Pool(NP) as p: Q = dict(p.map(one, [x[0] for x in L]))
    H = ['sample', 'cohort', 'sex', 'population', 'superpopulation', 'batch', 'haplogroup', 'lineage', 'y_dose', 'x_dose', 'DYZ3_bp', 'Yq12_bp', 'DYZ1_bp', 'DYZ18S_bp', 'DYZ18sub_bp', 'DYZ19_bp', 'DYZ17_bp', 'DYZ2_bp', 'YqHet2.7_bp', 'estimate_status', 'reference_male', 'yp_ruler_dose', 'yq_ruler_dose', 'y_state', 'yp_deletion_bias_pct', 'sex_chromosomes', 'par1_dose', 'par2_dose', 'library', 'y_status']
    with open(OUT, 'w') as o:
        o.write('\t'.join(H) + '\n')
        for s, co, u in L:
            v = Q[s]; st = state(v); sex = SX.get(s.replace('SGDP_', ''), '.'); pop = u.split('/data/')[1].split('/')[0] if '/data/' in u else '.'
            if co == 'HGDP' and not os.path.exists(f'{RC}/{s}/header.sam.gz'): sys.exit(f'ERROR {s}: no {RC}/{s}/header.sam.gz for the read-group rule')
            rg = sp.run(f'zcat {RC}/{s}/header.sam.gz | grep -c "^@RG"', shell=True, capture_output=True, text=True).stdout.strip()   # HGDP: PCR-free with >= 2 read groups, also when the published type differs (Methods)
            lib = ('PCR-free' if int(rg or 0) >= 2 else 'PCR') if co == 'HGDP' else 'PCR-free' if 'simons' in u else '.'
            nx, ny = round(v['x_dose']), round(v['y_dose']); sc = 'X' * nx + 'Y' * ny if abs(v['x_dose'] - nx) <= 0.15 and abs(v['y_dose'] - ny) <= 0.15 and nx > 0 else 'mosaic'
            ys = 'no_y' if round(v['y_dose'], 4) <= NOY else 'y_found'; es = 'no_y' if ys == 'no_y' else 'ok'
            row = {'sample': s, 'cohort': co, 'sex': sex, 'population': pop, 'superpopulation': '.', 'batch': '.', 'haplogroup': '.', 'lineage': '.', 'estimate_status': es, 'reference_male': 'no',
                   'y_state': st, 'yp_deletion_bias_pct': f"{100 * (v['yq_ruler_dose'] / v['y_dose'] - 1):+.2f}" if ys == 'y_found' else '.', 'sex_chromosomes': sc, 'library': lib, 'y_status': ys}
            for f in F: row[f] = 'nan' if v[f] != v[f] else f'{v[f]:.0f}' if f.endswith('_bp') else f'{v[f]:.4f}'
            o.write('\t'.join(row[h] for h in H) + '\n'); log(s, co, sex, pop, lib, f"y_dose {v['y_dose']:.4f}", st, sc, f"Yq12 {row['Yq12_bp']}")
log('written', OUT)
