# Model: per-class GC floors, Yq12 composite
# build: python3 um_model.py <base model> <floor model> <class:floor,...> <out dir>; library: load(), lengths()
import os, re, sys, hashlib, numpy as np
P = '/path/to/project'; O = f'{P}/markers_u25'
def log(*a): print(*a, flush=True)
def rtsv(f): R = [l.rstrip('\n').split('\t') for l in open(f)]; return R[0], [dict(zip(R[0], r)) for r in R[1:]]
def gcount(U):   # G+C bases of each 31-mer key
    g = np.zeros(len(U), np.uint8)
    for i in range(31): b = (U >> np.uint64(2 * i)) & np.uint64(3); g += ((b == 1) | (b == 2)).astype(np.uint8)
    return g
def load(D=f'{P}/model'): return dict(np.load(f'{D}/model.npz').items()), rtsv(f'{D}/model.tsv')[1]   # load before fork
def gcbins(y, g):   # G+C curve: ruler counts pooled per G+C bin, counts above 10 x (mean + 1) dropped; needs 1e6 ruler counts, bins 1,000 markers and 1e4 counts
    m = y.mean(); k = y <= 10 * (m + 1)
    if y.sum() < 1e6 or not m > 0: return None
    al = y[k].mean(); d = [(b, y[k & (g == b)].mean() / al) for b in range(32) if (g == b).sum() >= 1000 and y[k & (g == b)].sum() >= 1e4]
    return d or None
def lengths(c, fm, M, gc, y0):   # counts in u_keys order -> (reported bp, sum_gc per class, ruler sum); None without a G+C curve
    y = c[y0].astype(np.float64); g = gc[y0]; Y = y.sum(); d = gcbins(y, g)
    if not d: return None
    r = np.array([min(d, key=lambda e: abs(e[0] - b))[1] for b in range(32)]); fr = (np.bincount(g, minlength=32) * r).sum() / len(g)
    S = {m['class']: (np.bincount(gc[fm[f'm_{m["class"]}']], weights=c[fm[f'm_{m["class"]}']], minlength=32) / r).sum() * fr for m in M}
    L = {m['class']: float(m['K']) * S[m['class']] / Y for m in M}
    return {m['class']: sum(-L[t[1:]] if t[0] == '-' else L[t.lstrip('+')] for t in re.findall(r'[+-]?[^+-]+', m['length'])) for m in M if m['length'] != 'none'}, S, Y
if __name__ == '__main__':
    B, F, FL, D = sys.argv[1:5]; FL = {c: int(g) for c, g in (x.split(':') for x in FL.split(','))}; os.makedirs(D, exist_ok=True)
    gc = gcount(np.load(f'{O}/u_keys.npy')); fb, ff, fa = (dict(np.load(f'{O}/{m}.npz').items()) for m in (B, F, 'final_model'))
    (H, mb), (H2, mf) = rtsv(f'{O}/{B}.tsv'), rtsv(f'{O}/{F}.tsv'); mb = {m['class']: m for m in mb}; c0 = mf[0]['class']; y0 = ff[f'ysc_{c0}']
    if H != H2: sys.exit(f'ERROR {B} and {F}: different columns')
    for m in mf:
        c = m['class']
        if not all(np.array_equal(z[f'ysc_{c}'], y0) for z in (fb, ff)): sys.exit(f'ERROR ruler of {c} differs')
        if c in FL:
            g = gc[ff[f'm_{c}']]; log(f'{c}: GC >= {100 * FL[c] / 31:.1f}% ({FL[c]} of 31): markers {len(g)} (base {len(fb[f"m_{c}"])}), lowest GC {g.min()} of 31, K {float(m["K"]):.1f} (base {float(mb[c]["K"]):.1f})')
            if g.min() < FL[c]: sys.exit(f'ERROR {c}: marker below the floor')
        elif m != mb[c] or not np.array_equal(ff[f'm_{c}'], fb[f'm_{c}']): sys.exit(f'ERROR {c}: no floor but differs between {B} and {F}')
        else: log(f'{c}: no floor, markers {m["markers"]}, equal in both models')
    log(f'ruler: {len(y0)} markers, equal to the active final_model ruler: {np.array_equal(fa[f"ysc_{c0}"], y0)}')
    np.savez(f'{D}/model.npz', **ff, m_yq12dyz2=fb['m_dyz2'], ysc_yq12dyz2=fb['ysc_dyz2'])   # DYZ2 part of Yq12: base DYZ2
    with open(f'{D}/model.tsv', 'w') as o:
        o.write('\t'.join(H + ['min_gc', 'length']) + '\n')
        for m in mf + [dict(mb['dyz2'], **{'class': 'yq12dyz2'})]:
            c = m['class']; o.write('\t'.join([m[h] for h in H] + [str(FL.get(c, 0)), 'yq12-yq12dyz2+dyz2' if c == 'yq12' else 'none' if c == 'yq12dyz2' else c]) + '\n')
    with open(f'{D}/source.txt', 'w') as o:
        for m in (B, F): o.write(f'{m}.npz md5 {hashlib.md5(open(f"{O}/{m}.npz", "rb").read()).hexdigest()}; {m}.tsv md5 {hashlib.md5(open(f"{O}/{m}.tsv", "rb").read()).hexdigest()}\n')
        o.write(f'GC floors (of 31 bases): {",".join(f"{c}:{g}" for c, g in FL.items())}; yq12dyz2 = dyz2 of {B}\n')
    log(open(f'{D}/model.tsv').read() + open(f'{D}/source.txt').read() + f'written {D}')
