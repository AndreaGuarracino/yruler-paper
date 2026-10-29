# Release markers and model.tsv for yruler
# usage: python3 um_release.py [model dir]; writes <model dir>/release/markers.fa and model.tsv
import os, sys, numpy as np
P = '/path/to/project'; O = f'{P}/markers_u25'; sys.path.insert(0, f'{P}/scripts'); import um_model
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
fa = open(f'{O}/markers_u.fa').read().split('\n'); c2 = [fa[i + 1] for i in range(0, len(fa) - 1, 2) if fa[i].startswith('>chr2depth_')]
C2 = np.unique(np.searchsorted(U, enc(c2))); y0 = fm[f'ysc_{M[0]["class"]}']
for r in M:
    if not np.array_equal(fm[f'ysc_{r["class"]}'], y0): raise SystemExit(f'ERROR ruler of {r["class"]} differs')
    if len(fm[f'm_{r["class"]}']) != int(r['markers']) or len(np.unique(fm[f'm_{r["class"]}'])) != int(r['markers']): raise SystemExit(f'ERROR marker count {r["class"]}')
S = [(r['class'], fm[f'm_{r["class"]}']) for r in M] + [('ysinglecopy', y0), ('chr2depth', C2)]
with open(f'{R}/markers.fa', 'wb') as o:
    for c, ix in S:
        s = dec(U[np.sort(ix)]); o.write(f'>{c}\n'.encode()); o.write(np.hstack([s, np.full((len(s), 1), 10, np.uint8)]).tobytes())
        if not np.array_equal(np.sort(enc([s.tobytes().decode()])), U[np.sort(ix)]): raise SystemExit(f'ERROR round trip {c}')   # all k-mers
        log(f'{c}: {len(ix)} markers')
for (a, x), (b, y) in [(p, q) for i, p in enumerate(S) for q in S[i + 1:]]:
    n = len(np.intersect1d(x, y))
    if n: log(f'shared {a} {b}: {n}')
with open(f'{R}/model.tsv', 'w') as o:
    o.write('# length_bp = K * sum_gc(class) / sum(ysinglecopy); K per technology (-x); length: reported length as a signed sum of class lengths, none = not reported\nclass\tK\ttech\tlength\n')
    for r in M: o.write(f'{r["class"]}\t{r["K"]}\tillumina\t{r["length"]}\n')
log('written', f'{R}/markers.fa', f'{R}/model.tsv')
