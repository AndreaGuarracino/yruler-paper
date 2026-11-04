#!/usr/bin/env python3
# Repair a download from corrupt copies
# usage: blockfix.py OUT MD5 URL COPY1 COPY2 [COPY3 ...]
# per 1 MiB block: majority of copies; else HTTP range fetch until two fetches agree
# MD5 '-' = unknown: accept the vote result
import sys, hashlib, subprocess, collections as co
out, md5, url, copies = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
B = 1 << 20; fh = [open(c, 'rb') for c in copies]; o = open(out, 'wb'); h = hashlib.md5()
log = lambda *a: print(*a, file=sys.stderr, flush=True)
def fetch(i, n):
    r = subprocess.run(['curl', '-sS', '--fail', '--retry', '5', '-r', f'{i * B}-{i * B + n - 1}', url], capture_output=True)
    if r.returncode or len(r.stdout) != n: log(f'block {i}: range fetch exit {r.returncode} bytes {len(r.stdout)}'); return None
    return r.stdout
i = disputed = voted = fetched = 0
while True:
    bl = [f.read(B) for f in fh]
    if not bl[0]: break
    c = co.Counter(bl); best, n = c.most_common(1)[0]
    if len(c) > 1:
        disputed += 1
        if n >= 2: voted += 1
        else:
            got = [];
            for t in range(6):
                x = fetch(i, len(bl[0]))
                if x is not None: got.append(x)
                if len(got) >= 2 and co.Counter(got).most_common(1)[0][1] >= 2: break
            if not got or co.Counter(got).most_common(1)[0][1] < 2: log(f'ERROR block {i}: no two fetches agree'); sys.exit(1)
            best = co.Counter(got).most_common(1)[0][0]; fetched += 1
    o.write(best); h.update(best); i += 1
o.close(); m = h.hexdigest()
log(f'blocks {i}, disputed {disputed}, by vote {voted}, by range fetch {fetched}; md5 {m} expected {md5}')
sys.exit(0 if md5 == '-' or m == md5 else 2)
