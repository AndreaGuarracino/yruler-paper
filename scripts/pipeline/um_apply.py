# Apply the model to saved counts
# usage: python3 um_apply.py <pass dir> <out.tsv> [model dir]; per sample: cohort, sex, population, Y and X dose, lengths, estimate_status, Y state
import os, sys, glob, subprocess as sp, numpy as np, multiprocessing as mp
P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_u25")}'; A25 = f'{P}/markers_u25'; sys.path.insert(0, f'{P}/scripts'); import um_model, truth_cl
D, OUT = f'{P}/{sys.argv[1]}', sys.argv[2]; MD = sys.argv[3] if len(sys.argv) > 3 else f'{P}/model'; Z = 'zstd'
U = np.load(f'{O}/u_keys.npy'); oi = np.load(f'{A25}/counts_order_idx.npy'); gc = um_model.gcount(U); fm, M = um_model.load(MD); y0 = fm[f'ysc_{M[0]["class"]}']
LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]
fa = open(f'{O}/markers_u.fa').read().split('\n'); c2 = [fa[i + 1] for i in range(0, len(fa) - 1, 2) if fa[i].startswith('>chr2depth_')]
v = LUT[np.frombuffer(''.join(c2).encode(), np.uint8).reshape(-1, 31)]; k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
kc = np.minimum(k, rc); C2 = np.minimum(np.searchsorted(U, kc), len(U) - 1)
if not np.array_equal(U[C2], kc): sys.exit(f'ERROR {(U[C2] != kc).sum()} chr2 markers not in u_keys')
del fa, c2, v, k, rc, kc
V5 = O != A25   # v5, v6: markers_u25 whole counts plus saved-read counts of the other keys
if V5:
    U0 = np.load(f'{A25}/u_keys.npy'); o2n = np.searchsorted(U, U0); N0 = len(o2n); DD = f'{O}/{os.environ.get("DYZ17DIR", "dyz17c")}'   # markers_u25 whole counts -> UM order; DYZ17DIR/DFA: saved-read counts of the keys not in markers_u25 (v5: dyz17c/dyz17.fa, v6: deltac/delta.fa)
    if not (U[np.minimum(o2n, len(U) - 1)] == U0).all(): sys.exit(f'ERROR {O} keys do not contain markers_u25 keys')
    sq = open(f'{O}/{os.environ.get("DFA", "dyz17.fa")}.order').read().split('\n')[:-1]; v = LUT[np.frombuffer(''.join(sq).encode(), np.uint8).reshape(-1, 31)]; k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
    for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
    kd = np.minimum(k, rc); DJ = np.searchsorted(U, kd); assert (U[DJ] == kd).all(); NEW = ~np.isin(DJ, o2n); del sq, v, k, rc, kd
    print(f'{O}: {len(DJ)} markers in the saved-read count file, {int(NEW.sum())} keys not in markers_u25, counts from {DD}', flush=True)
DN = os.environ.get('DENOMC')   # DENOMC=<slice counts dir> DENOMFA=<its marker fa> DENOMK=<window keys npy>[,<npy>...]: Y dose denominator = median over windows of the mean window key count in the region slice
if DN:
    sq = open(os.environ['DENOMFA'] + '.order').read().split('\n')[:-1]; v = LUT[np.frombuffer(''.join(sq).encode(), np.uint8).reshape(-1, 31)]; k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
    for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
    WN = os.environ['DENOMK'].split(','); kc = np.minimum(k, rc); NDR = len(sq); DR = []; del sq, v, k, rc
    WL = os.environ['DENOMN'].split(',') if os.environ.get('DENOMN') else [os.path.basename(f)[:-4] for f in WN]   # DENOMN: window labels
    SEL = os.environ.get('DENOMSEL', 'median'); print(f'denominator rule: {"first window, the others as checks" if SEL == "first" else "median of the windows"}', flush=True)   # DENOMSEL=first: denominator = first window of DENOMK
    if len(set(WL)) != len(WN): sys.exit(f'ERROR denominator window labels not unique: {WL} (set DENOMN)')
    for f in WN:
        WK = np.load(f); DR.append(np.isin(kc, WK)); print(f'denominator window {os.path.basename(f)}: {int(DR[-1].sum())} of {len(WK)} keys in {os.environ["DENOMFA"]}, counts from {DN}', flush=True)
        if DR[-1].sum() < 0.99 * len(WK): sys.exit(f'ERROR denominator: keys of {f} missing from the slice marker file')
    del kc
def den(s, c):   # denominator depth (median over windows) and per-window means: window keys in the slice, else the chr2 50-60 markers
    if not DN: return c[C2].mean(), [c[C2].mean()]
    d = np.array(sp.run([Z, '-dcq', f'{DN}/{s}/counts.col.zst'], check=True, capture_output=True).stdout.split(), np.float64)
    if len(d) != NDR: raise ValueError(f'{len(d)} slice counts, expected {NDR}')
    w = [d[m].mean() for m in DR]; return (w[0] if SEL == 'first' else float(np.median(w))), w
YPQ = np.load(os.environ['YPQ']) if os.environ.get('YPQ') else None   # YPQ=<um_v5pos.py ypq XDR block labels>: raw Yp (XDR1-2), Yq (XDR3-8) and per-block ruler doses -> <out>.ypq.tsv
if YPQ is not None and len(YPQ) != len(y0): sys.exit(f'ERROR {len(YPQ)} Yp/Yq labels, ruler {len(y0)}')
NOY = float(os.environ['NOY']) if os.environ.get('NOY') else None   # NOY=<no-Y cut on the scaled y_dose>: basis rule of the release
YDN = float(os.environ.get('YDN', 1))   # YDN=<reference median of the raw 1kGP y_dose>: y_dose relative to the median 1kGP male (scale B, user 2026-10-31)
GB = [] if YPQ is None else [y0[(YPQ == 1) | (YPQ == 2)], y0[YPQ >= 3]] + [y0[YPQ == b] for b in range(1, 9)]
def one(s):   # Y dose and reported lengths of one sample
    try:
        x = np.array(sp.run([Z, '-dcq', f'{D}/{s}/whole/counts.col.zst'], check=True, capture_output=True).stdout.split(), np.float64)
        if len(x) != len(oi): raise ValueError(f'{len(x)} counts, expected {len(oi)}')
        if V5:
            d = np.array(sp.run([Z, '-dcq', f'{DD}/{s}/counts.col.zst'], check=True, capture_output=True).stdout.split(), np.float64)
            if len(d) != len(DJ): raise ValueError(f'{len(d)} saved-read counts, expected {len(DJ)}')
            c0 = np.zeros(N0); c0[oi] = x; c = np.zeros(len(U)); c[o2n] = c0; c[DJ[NEW]] = d[NEW]
        else: c = np.zeros(len(U)); c[oi] = x
        d, wm = den(s, c); pq = [2 * c[g].mean() / d for g in GB]; yd = 2 * c[y0].mean() / d / YDN
        e = um_model.lengths(c, fm, M, gc, y0, flat=NOY is not None and round(yd, 4) > NOY)   # NOY: ruler basis above the no-Y cut, no G+C curve needed
        if e and NOY is not None and round(yd, 4) <= NOY: e = None   # no Y signal: haploid basis
        if e and not um_model.gcbins(c[y0].astype(np.float64), gc[y0]): print(f'{s}: y_dose {yd:.4f} above the no-Y cut without a G+C curve: ruler basis, no G+C correction', flush=True)
        return s, yd, (e[0], 'ruler') if e else (um_model.haploid(c, fm, M, d), 'haploid_depth'), pq, wm
    except Exception as e: print(f'WARNING {s}: {e}', flush=True); return s, None, None, [], []
S = sorted(f.split('/')[-3] for f in glob.glob(f'{D}/*/whole/counts.col.zst')); print(f'{D}: {len(S)} samples, model {MD}', flush=True)
with mp.Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', 8))) as p: R = p.map(one, S)
ref = set().union(*[set(x) for x in truth_cl.load(O).values()])
ped = {l.split()[1]: l.split() for l in open(f'{P}/yhaplo/ped3202.txt') if not l.startswith('FamilyID')}
for l in list(open(f'{P}/lists/1kGP.3202_samples.pedigree_info.txt'))[1:]: x = l.split(); ped[x[0]][4] = x[3]   # sex from the corrected 1kGP pedigree (Byrska-Bishop 2022): HG02300 female
hg = dict(l.rstrip('\n').split('\t')[:2] for l in open(f'{O}/lineage.tsv') if l.strip()); bat = {l.split('\t')[0]: l.split('\t')[3] for l in open(f'{P}/lists/pass1kgp.tsv')}
YS = {x[0]: x[2:10] for x in (l.rstrip('\n').split('\t') for l in list(open(os.environ.get('YSF', f'{P}/results/ystate_1kgp.tsv')))[1:])}   # Yp, Yq, state, ruler bias; YSF=<um_ystate.py output> (v6: model_v5_v6/ystate_1kgp.tsv)
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'yq2k7': 'YqHet2.7', 'dyz18sub': 'DYZ18sub', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
SX = {x[2]: x[3] for f in ('PRJEB6463_sex.tsv', 'PRJEB9586_sex.tsv') for x in (l.rstrip('\n').split('\t') for l in list(open(f'{P}/lists/{f}'))[1:])}   # HGDP, SGDP sex by sample alias
PO = {x[0]: x[1].split('/data/')[1].split('/')[0] for x in (l.split('\t') for l in open(f'{P}/lists/hgdp_sgdp_urls.tsv')) if '/data/' in x[1]}   # HGDP, SGDP population from the URL
LB = {x[0]: ('PCR-free' if int(x[1]) >= 2 else 'PCR') if x[0].startswith('HGDP') else 'PCR-free' if x[2].startswith('LP60') else 'unknown' for x in (l.split('\t') for l in list(open(f'{P}/model/cram_meta_hgdp_sgdp.tsv'))[1:])}   # library: HGDP >= 2 read groups, SGDP LP60* plates; 1kGP all PCR-free
def meta(s):   # cohort, sex, population, superpopulation
    if s in ped: pd = ped[s]; return '1kGP', {'1': 'male', '2': 'female'}.get(pd[4], '.'), pd[5], pd[6]
    if s.startswith('SGDP_'): return 'SGDP', SX.get(s[5:], '.'), PO.get(s, '.'), '.'
    if s.startswith('HGDP'): return 'HGDP', SX.get(s, '.'), PO.get(s, '.'), '.'
    return 'other_run:' + s.split('_', 1)[-1], 'male', '.', '.'
cls = [m['class'] for m in M if m['length'] != 'none']; lin = lambda h: '.' if h == '.' else h[:3] if h.startswith('R1') else h[0]; n = {}
with open(OUT, 'w') as o:
    o.write('\t'.join(['sample', 'cohort', 'sex', 'population', 'superpopulation', 'batch', 'haplogroup', 'lineage', 'y_dose', 'x_dose'] + [f'{NM[c]}_bp' for c in cls] + ['estimate_status', 'reference_male', 'yp_ruler_dose', 'yq_ruler_dose', 'y_state', 'yp_deletion_bias_pct', 'sex_chromosomes', 'par1_dose', 'par2_dose', 'library', 'estimate_basis']) + '\n')
    for s, dose, e, pq, wm in R:
        co, sx, po, su = meta(s); h = hg.get(s, '.'); y = YS.get(s, []); y = y + ['.'] * (8 - len(y))
        f = 'no_counts' if dose is None else 'female' if sx == 'female' else 'haploid_depth' if e[1] != 'ruler' else 'low_y_dose' if dose < 0.9 else 'ok'; n[f] = n.get(f, 0) + 1
        o.write('\t'.join([s, co, sx, po, su, bat.get(s, '.'), h, lin(h), 'nan' if dose is None else f'{dose:.4f}', y[4]] + [f'{e[0][c]:.0f}' if e else 'nan' for c in cls] + [f, 'yes' if s in ref else 'no'] + y[:3] + [y[3] if y[2] == 'yp_deletion' else '.'] + y[5:8] + ['PCR-free' if co == '1kGP' else LB.get(s, '.'), e[1] if e else '.']) + '\n')
print(f'written {OUT}: {len(R)} samples, flags {n}', flush=True)
if YPQ is not None:
    with open(f'{OUT}.ypq.tsv', 'w') as o: o.write('\t'.join(['sample', 'yp_raw', 'yq_raw'] + [f'xdr{b}_raw' for b in range(1, 9)]) + '\n'); o.writelines(s + ''.join(f'\t{v:.5f}' for v in pq) + '\n' for s, d, e, pq, wm in R if pq)   # every sample, sex-blind
    print(f'written {OUT}.ypq.tsv: keys Yp {len(GB[0])}, Yq {len(GB[1])}, per block ' + ' '.join(str(len(g)) for g in GB[2:]), flush=True)
if DN and len(DR) > 1:   # per-window depth / median: a window off by > 10% flags a local copy change in that sample
    with open(f'{OUT}.den.tsv', 'w') as o: o.write('sample\t' + '\t'.join(WL) + '\tmedian\tdenominator\t' + '\t'.join(f'{w}_ratio' for w in WL) + '\n'); o.writelines(f'{s}\t' + '\t'.join(f'{x:.4f}' for x in wm) + f'\t{np.median(wm):.4f}\t{(wm[0] if SEL == "first" else np.median(wm)):.4f}\t' + '\t'.join(f'{x / np.median(wm):.4f}' for x in wm) + '\n' for s, d, e, pq, wm in R if wm)
    fl = [(s, [w for w, x in zip(WL, wm) if abs(x / np.median(wm) - 1) > 0.1]) for s, d, e, pq, wm in R if wm]; fl = [x for x in fl if x[1]]
    print(f'written {OUT}.den.tsv; samples with a window off the median by > 10%: {len(fl)}' + ''.join(f'; {s} {",".join(w)}' for s, w in fl[:40]), flush=True)
ok = [e[0] for s, d, e, pq, wm in R if e and e[1] == 'ruler' and d >= 0.9 and meta(s)[1] != 'female']
for c in cls if ok else []: q = np.percentile([e[c] for e in ok], [25, 50, 75]); print(f'{NM[c]}: {len(ok)} males y_dose >= 0.9, median {q[1]:,.0f} bp (IQR {q[0]:,.0f}-{q[2]:,.0f})', flush=True)
