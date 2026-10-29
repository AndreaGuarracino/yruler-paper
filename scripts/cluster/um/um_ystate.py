# Y state per male: 4-class rule
# usage: python3 um_ystate.py <yprofile.tsv> <pass dir> <out.tsv>
import sys, re, os
YP, D, OUT = sys.argv[1:4]; n = {}
with open(OUT, 'w') as o:
    o.write('sample\ty_dose\typ_ruler_dose\tyq_ruler_dose\ty_state\truler_bias_pct\n')
    for l in list(open(YP))[1:]:
        s, p, q = l.split('\t'); p, q = float(p), float(q); f = f'{D}/{s}/whole/run.log'
        m = re.search(r'y_dose (\S+)', open(f).read()) if os.path.exists(f) else None
        if not m: print('no y_dose', s, file=sys.stderr); continue
        y = float(m.group(1).rstrip(';'))
        if p < 0.5 and q < 0.5: st = 'no_Y'
        elif p < 0.5 and abs(q - 1) <= 0.1: st = 'yp_deletion'   # germline AMELY-type deletion
        elif abs(p - q) > 0.1: st = 'partial'
        elif y >= 1.75: st = 'XYY'
        elif y >= 1.25: st = 'mosaic_gain'
        elif y < 0.9: st = 'mosaic_loss'
        else: st = 'normal'
        n[st] = n.get(st, 0) + 1; o.write(f'{s}\t{y:.4f}\t{p:.4f}\t{q:.4f}\t{st}\t{100 * (q / y - 1):+.2f}\n')   # estimate bias if ruler lost
print(n, file=sys.stderr)
