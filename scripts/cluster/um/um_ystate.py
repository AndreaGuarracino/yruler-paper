# Y state per male: 4-class rule
# usage: python3 um_ystate.py <yprofile.tsv> <pass dir> <out.tsv> <xdose.tsv>; X dose = chrX (no PAR) / chr2 MQ20 reads over the median of males with Y dose 0.95-1.05
import sys, re, os
YP, D, OUT, XD = sys.argv[1:5]; n = {}
XR = {x[0]: int(x[2]) / int(x[1]) for x in (l.split('\t') for l in list(open(XD))[1:]) if int(x[1]) > 0}
YD = {}
for l in list(open(YP))[1:]:
    s = l.split('\t')[0]; f = f'{D}/{s}/whole/run.log'; m = re.search(r'y_dose (\S+)', open(f).read()) if os.path.exists(f) else None
    if m: YD[s] = float(m.group(1).rstrip(';'))
import statistics; XREF = statistics.median(XR[s] for s, y in YD.items() if 0.95 <= y <= 1.05 and s in XR)
print(f'X reference ratio {XREF:.6f} from {sum(1 for s, y in YD.items() if 0.95 <= y <= 1.05 and s in XR)} samples', file=sys.stderr)
with open(OUT, 'w') as o:
    o.write('sample\ty_dose\typ_ruler_dose\tyq_ruler_dose\ty_state\truler_bias_pct\tx_dose\tsex_chromosomes\n')
    for l in list(open(YP))[1:]:
        s, p, q = l.split('\t'); p, q = float(p), float(q); f = f'{D}/{s}/whole/run.log'
        m = re.search(r'y_dose (\S+)', open(f).read()) if os.path.exists(f) else None
        if not m: print('no y_dose', s, file=sys.stderr); continue
        y = float(m.group(1).rstrip(';'))
        if y <= 0.0222: st = 'no_Y'   # at most the highest heldout female dose
        elif p < 0.5 and abs(q - 1) <= 0.1: st = 'yp_deletion'   # germline AMELY-type deletion
        elif abs(p - 1) <= 0.1 and q < p - 0.1: st = 'yq_loss'   # Yp intact, part of Yq lost
        elif abs(p - q) > 0.1: st = 'partial'
        elif y >= 1.75: st = 'Y_gain'
        elif y >= 1.25: st = 'mosaic_gain'
        elif y < 0.9: st = 'mosaic_loss'
        else: st = 'normal'
        x = XR[s] / XREF if s in XR else float('nan'); nx, ny = round(x) if x == x else 0, round(y)   # X and Y copies when both near integers
        sc = 'X' * nx + 'Y' * ny if x == x and abs(x - nx) <= 0.15 and abs(y - ny) <= 0.15 and nx > 0 else 'mosaic' if x == x else '.'
        n[st] = n.get(st, 0) + 1; n[sc] = n.get(sc, 0) + 1; o.write(f'{s}\t{y:.4f}\t{p:.4f}\t{q:.4f}\t{st}\t{100 * (q / y - 1):+.2f}\t{x:.4f}\t{sc}\n')   # estimate bias if ruler lost
print(n, file=sys.stderr)
