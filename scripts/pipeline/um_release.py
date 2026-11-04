# Release markers and model.tsv for yruler
# usage: [UM=markers_v6 DENOMK=<window keys npy>[,<npy>...] YDN=<reference raw y_dose> NOY=<no-Y cut, scaled y_dose>] python3 um_release.py [model dir]; writes <model dir>/release/markers.fa and model.tsv; denominator classes depth1..depthN = the window keys (default: chr2depth = chr2 50-60 markers); YDN -> #y_dose_scale in model.tsv
import os, sys, numpy as np
P = os.environ['P']; O = f'{P}/{os.environ.get("UM", "markers_u25")}'; sys.path.insert(0, f'{P}/scripts'); import um_model
MD = sys.argv[1] if len(sys.argv) > 1 else f'{P}/model'; R = f'{MD}/release'; os.makedirs(R, exist_ok=True)
def log(*a): print(*a, flush=True)
U = np.load(f'{O}/u_keys.npy'); fm, M = um_model.load(MD)
if any('_' in r['class'] for r in M): raise SystemExit('ERROR class names with _ (yruler reads the header name before _)')
LUT = np.full(256, 255, np.uint8); LUT[[65, 67, 71, 84]] = [0, 1, 2, 3]
def enc(lines):   # canonical keys
    v = LUT[np.frombuffer(''.join(lines).encode(), np.uint8).reshape(-1, 31)]; k = np.zeros(len(v), np.uint64); rc = np.zeros(len(v), np.uint64)
    for i in range(31): k = (k << np.uint64(2)) | v[:, i].astype(np.uint64)
    for i in range(30, -1, -1): rc = (rc << np.uint64(2)) | (3 - v[:, i]).astype(np.uint64)
    return np.minimum(k, rc)
def dec(keys):   # uint64 keys -> 31-mers
    b = np.zeros((len(keys), 31), np.uint8)
    for i in range(31): b[:, 30 - i] = ((keys >> np.uint64(2 * i)) & np.uint64(3)).astype(np.uint8)
    return np.frombuffer(b'ACGT', np.uint8)[b]
if os.environ.get('DENOMK'): DW = [(f'depth{i + 1}', np.unique(np.load(f).astype(np.uint64))) for i, f in enumerate(os.environ['DENOMK'].split(','))]
else: fa = open(f'{O}/markers_u.fa').read().split('\n'); DW = [('chr2depth', U[np.unique(np.searchsorted(U, enc([fa[i + 1] for i in range(0, len(fa) - 1, 2) if fa[i].startswith('>chr2depth_')])))])]
y0 = fm[f'ysc_{M[0]["class"]}']
for r in M:
    if not np.array_equal(fm[f'ysc_{r["class"]}'], y0): raise SystemExit(f'ERROR ruler of {r["class"]} differs')
    if len(fm[f'm_{r["class"]}']) != int(r['markers']) or len(np.unique(fm[f'm_{r["class"]}'])) != int(r['markers']): raise SystemExit(f'ERROR marker count {r["class"]}')
S = [(r['class'], U[np.sort(fm[f'm_{r["class"]}'])]) for r in M] + [('ysinglecopy', U[np.sort(y0)])] + DW   # sorted canonical keys
with open(f'{R}/markers.fa', 'wb') as o:
    for c, ix in S:
        s = dec(ix); o.write(f'>{c}\n'.encode()); o.write(np.hstack([s, np.full((len(s), 1), 10, np.uint8)]).tobytes())
        if not np.array_equal(np.sort(enc([s.tobytes().decode()])), ix): raise SystemExit(f'ERROR round trip {c}')   # all k-mers
        log(f'{c}: {len(ix)} markers')
for (a, x), (b, y) in [(p, q) for i, p in enumerate(S) for q in S[i + 1:]]:
    n = len(np.intersect1d(x, y))
    if n: log(f'shared {a} {b}: {n}')
with open(f'{R}/model.tsv', 'w') as o:
    o.write('# length_bp = K * sum_gc(class) / sum(ysinglecopy); K per technology (-x); length: reported length as a signed sum of class lengths, none = not reported\n' + (f'#y_dose_scale={float(os.environ["YDN"]):.6f}\n' if os.environ.get('YDN') else '') + (f'#no_y_dose={float(os.environ["NOY"]):.6f}\n' if os.environ.get('NOY') else '') + 'class\tK\ttech\tlength\n')
    for r in M: o.write(f'{r["class"]}\t{r["K"]}\tillumina\t{r["length"]}\n')
log('written', f'{R}/markers.fa', f'{R}/model.tsv')
