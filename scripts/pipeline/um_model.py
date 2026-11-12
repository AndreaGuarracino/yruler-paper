# Model: per-class GC floors, Yq12 = DYZ1 + DYZ2 + rest
# build: (TRUTHX, TRUTHSET, EXCL as the model) python3 um_model.py <base model or per-class prefix> <floor model or prefix> <class:floor,...> <out dir>; library: load(), lengths()
import os, re, sys, hashlib, numpy as np
P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_u25")}'
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
def lengths(c, fm, M, gc, y0, flat=False):   # counts in u_keys order -> (reported bp, sum_gc per class, ruler sum); None without a G+C curve unless flat (no G+C correction, as the CV)
    y = c[y0].astype(np.float64); g = gc[y0]; Y = y.sum(); d = gcbins(y, g)
    if not d and not (flat and Y > 0): return None
    r = np.array([min(d, key=lambda e: abs(e[0] - b))[1] for b in range(32)]) if d else np.ones(32); fr = (np.bincount(g, minlength=32) * r).sum() / len(g)
    S = {m['class']: (np.bincount(gc[fm[f'm_{m["class"]}']], weights=c[fm[f'm_{m["class"]}']], minlength=32) / r).sum() * fr for m in M}
    return report({m['class']: float(m['K']) * S[m['class']] / Y for m in M}, M), S, Y
def report(L, M):   # reported bp per class from model.tsv 'length': a class, a constant, or their sum (Yq12 = DYZ1 + DYZ2 + rest)
    v = lambda t: float(t) if t[0].isdigit() else L[t]
    return {m['class']: sum(-v(t[1:]) if t[0] == '-' else v(t.lstrip('+')) for t in re.findall(r'[+-]?[^+-]+', m['length'])) for m in M if m['length'] != 'none'}
def haploid(c, fm, M, d2):   # no Y ruler signal (no G+C curve): bp at haploid autosomal depth d2 / 2 (as the CV female background), reported as in model.tsv
    return report({m['class']: float(m['K']) * c[fm[f'm_{m["class"]}']].sum() / (len(fm[f'ysc_{m["class"]}']) * d2 / 2) for m in M}, M)
CL = ['dyz3', 'yq12', 'dyz1', 'dyz18s', 'dyz18sub', 'dyz19', 'dyz17', 'dyz2', 'yq2k7', 'yq12rest', 'yq12rest2']
def joined(pres):   # per-class production files <prefix>_<class> -> (npz dict, header, rows)
    Z, H, R = {}, None, []
    for c, p in pres:
        z = np.load(f'{O}/{p}_{c}.npz'); h, rows = rtsv(f'{O}/{p}_{c}.tsv'); r = [m for m in rows if m['class'] == c]
        if (H and h != H) or len(r) != 1: sys.exit(f'ERROR {p}_{c}: columns or rows')
        H = h; Z[f'm_{c}'], Z[f'ysc_{c}'] = z[f'm_{c}'], z[f'ysc_{c}']; R.append(r[0])
    return Z, H, R
def model(m, base, floors):   # a model name, or a prefix of per-class files (classes in floors from it, others from base)
    if os.path.exists(f'{O}/{m}.npz'): return (dict(np.load(f'{O}/{m}.npz').items()),) + rtsv(f'{O}/{m}.tsv')
    return joined([(c, m if (base is None or c in floors) else base) for c in CL if os.path.exists(f'{O}/{(m if (base is None or c in floors) else base)}_{c}.npz')])
if __name__ == '__main__':
    B, F, FL, D = sys.argv[1:5]; FL = {} if FL == 'none' else {c: int(g) for c, g in (x.split(':') for x in FL.split(','))}; os.makedirs(D, exist_ok=True)   # FL none: no G+C floors
    gc = gcount(np.load(f'{O}/u_keys.npy')); (fb, H, mb), (ff, H2, mf) = model(B, None, FL), model(F, B, FL); fa = dict(np.load(f'{O}/final_model.npz').items()) if os.path.exists(f'{O}/final_model.npz') else None
    mb = {m['class']: m for m in mb}; c0 = mf[0]['class']; y0 = ff[f'ysc_{c0}']
    if H != H2: sys.exit(f'ERROR {B} and {F}: different columns')
    for m in mf:
        c = m['class']
        if not all(np.array_equal(z[f'ysc_{c}'], y0) for z in (fb, ff)): sys.exit(f'ERROR ruler of {c} differs')
        if c in FL:
            g = gc[ff[f'm_{c}']]; log(f'{c}: GC >= {100 * FL[c] / 31:.1f}% ({FL[c]} of 31): markers {len(g)} (base {len(fb[f"m_{c}"])}), lowest GC {g.min()} of 31, K {float(m["K"]):.1f} (base {float(mb[c]["K"]):.1f})')
            if g.min() < FL[c]: sys.exit(f'ERROR {c}: marker below the floor')
        elif m != mb[c] or not np.array_equal(ff[f'm_{c}'], fb[f'm_{c}']): sys.exit(f'ERROR {c}: no floor but differs between {B} and {F}')
        else: log(f'{c}: no floor, markers {m["markers"]}, equal in both models')
    log(f'ruler: {len(y0)} markers, equal to the active final_model ruler: {np.array_equal(fa[f"ysc_{c0}"], y0) if fa else "no final_model in " + O}')
    sys.path.insert(0, f'{P}/scripts'); import truth_cl; T = truth_cl.load(O)   # Yq12 rest: training truth males
    for l in open(os.environ['TRUTHX']) if os.environ.get('TRUTHX') else (): c, s = l.split()[:2]; T[c].pop(s, None)
    for l in open(os.environ['TRUTHSET']) if os.environ.get('TRUTHSET') else ():
        c, s, v = l.split()[:3]
        if s in T[c]: T[c][s] = float(v)
    EX = {l.split()[0] for l in open(os.environ['EXCL'])} if os.environ.get('EXCL') else set(); CM = {f[:-4] for f in os.listdir(f'{O}/cnpy') if f.endswith('.npy')}
    np.savez(f'{D}/model.npz', **ff)
    with open(f'{D}/model.tsv', 'w') as o:
        o.write('\t'.join(H + ['min_gc', 'length']) + '\n')
        for m in mf:
            c = m['class']; o.write('\t'.join([m[h] for h in H] + [str(FL.get(c, 0)), ('dyz1+dyz2+dyz18s+yq12rest2' if c == 'yq12' else 'none' if c in ('yq12rest', 'yq12rest2') else c)]) + '\n')
    with open(f'{D}/source.txt', 'w') as o:
        for m in sorted({B, F}): o.write(''.join(f'{x}.npz md5 {hashlib.md5(open(f"{O}/{x}.npz", "rb").read()).hexdigest()}; {x}.tsv md5 {hashlib.md5(open(f"{O}/{x}.tsv", "rb").read()).hexdigest()}\n' for x in ([m] if os.path.exists(f'{O}/{m}.npz') else [f'{m}_{c}' for c in CL if os.path.exists(f'{O}/{m}_{c}.npz')])))
        o.write(f'GC floors (of 31 bases): {",".join(f"{c}:{g}" for c, g in FL.items())}; Yq12 = DYZ1 + DYZ2 + DYZ18S + estimated rest (class yq12rest2)' + '\n')
    log(open(f'{D}/model.tsv').read() + open(f'{D}/source.txt').read() + f'written {D}')
