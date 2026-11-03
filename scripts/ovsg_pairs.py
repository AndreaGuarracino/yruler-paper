#!/usr/bin/env python3
"""SGDP vs 1kGP same-individual pairs, with mate overlap counted once on both sides.

usage: python3 scripts/ovsg_pairs.py <dir with <sample>/estimates.tsv or lengths.tsv> <same_male_pairs.tsv> <out.tsv>
"""
import sys, csv, os
import numpy as np
D, SMP, OUT = sys.argv[1:4]; CL = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17']
L = {}
for s in os.listdir(D):
    f = os.path.join(D, s, 'estimates.tsv'); f = f if os.path.exists(f) else os.path.join(D, s, 'lengths.tsv')   # new name, else old
    e = os.path.join(D, s, 'view.err')   # a CRAM read error leaves a non-empty view.err
    if os.path.exists(f) and not (os.path.exists(e) and os.path.getsize(e)):
        L[s] = {r['class']: float(r.get('estimate_bp', r.get('length_bp'))) for r in csv.DictReader(open(f), delimiter='\t')}
P = sorted(s[5:] for s in L if s.startswith('SGDP_') and s[5:] in L); print(f'{len(P)} complete pairs: {P}', flush=True)
B = [r for r in csv.DictReader(open(SMP), delimiter='\t')]; print('same_male_pairs columns', list(B[0])[:12], flush=True)
rows = []
for c in CL:
    d = np.array([100 * (L['SGDP_' + s][c] / L[s][c] - 1) for s in P])
    rows.append((c, len(d), f'{np.median(d):+.2f}', f'{np.median(np.abs(d)):.2f}', int((d > 0).sum())))
    print(f'{c}: {len(d)} pairs, mate overlap once on both sides: median signed {np.median(d):+.2f}%, median |diff| {np.median(np.abs(d)):.2f}%, SGDP higher in {int((d > 0).sum())}', flush=True)
with open(OUT, 'w') as o: o.write('class\tpairs\tmedian_signed_pct\tmedian_abs_pct\tsgdp_higher\n'); o.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
