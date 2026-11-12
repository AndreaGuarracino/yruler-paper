#!/usr/bin/env python3
"""All females of the three cohorts, released binary and model at the no-Y cut 0.025: one table per female and the background per class.

usage: python3 scripts/female_all.py <held-out estimates> <held-out females> <cohort estimates> <cohort table> <truth tsv> <out females> <out background> <SGDP metadata>
held-out estimates: results/checks/female_heldout_haploid.tsv (sample, class, estimate_bp, relative, y_dose, y_status; um_femhap.sh: every female on the haploid basis, the basis of all females below the cut, so HGDP00247 above the cut is comparable)
held-out females: results/checks/female_heldout_ydose.tsv (cohort; *_control rows are cohort females run again, taken from the cohort table)
cohort estimates: results/checks/female_cohort_estimates.tsv (the 346 model-cohort females counted again, um_femre2.sh; also the SGDP females without ENA sex, lists/femre_sgdp_ena.tsv)
cohort table: results/model/cohort_fem.tsv (Y dosage of the 346 model-cohort females; the recount gives the same y_dose to 4 decimals, 189 of 189 checked)
truth tsv: results/checks/truth_long_v9bBx.tsv (median assembly length per class)
SGDP metadata: data/sgdp_metadata.279public.21signedLetter.samples.txt; a recounted SGDP sample outside the cohort table enters only if its Gender is F (set sgdp_metadata_sex); Gender U (NOR111) is left out
Background: estimate as % of the median assembly length, for all females and without the females in EXPLAINED.
"""
import sys, collections
import numpy as np

HE, HF, CE, CF, TRUTH, OUTF, OUTB, SGM = sys.argv[1:9]
CLASSES = ['dyz3', 'yq12', 'dyz1', 'dyz2', 'dyz18s', 'dyz19', 'dyz17']
EXPLAINED = {'HG02772': 'Yq12 (DYZ1 and DYZ2) without other Y sequence; mother of HG02773'}   # the one female off the dosage trend (32 robust SD above, next 8)
NOY = 0.025   # no-Y cut of the released model, above the highest of the 1,970 females (0.0221; user 2026-11-06, was 0.02)
def log(*a): print(*a, flush=True)
def rows(p):
    L = [l.rstrip('\n').split('\t') for l in open(p)]
    return [dict(zip(L[0], x)) for x in L[1:]]

SG = {}   # SGDP name -> Gender
for x in [l.rstrip('\n').split('\t') for l in open(SGM, encoding='latin-1')][1:]:
    for a in {x[2], *(a.strip() for a in x[3].split(','))}: SG.setdefault(a, set()).add(x[10])
F = {}; ROUND = []; SKIP = set()   # sample -> {cohort, set, y_dose, status, class: bp}; y_dose one unit apart in the 4th decimal
for r in rows(CF):   # every model-cohort female: Y dosage from the cohort table
    d = float(r['y_dose']); F[r['sample']] = {'cohort': r['cohort'], 'set': 'model_cohort', 'y_dose': d, 'status': 'y_found' if d > NOY else 'no_y'}
for r in rows(CE):   # class estimates from the recount (released binary)
    f = F.get(r['sample']); g = SG.get(r['sample'].removeprefix('SGDP_')) if r['cohort'] == 'SGDP' else None
    if f is None and g == {'F'}:   # recorded female only in the SGDP metadata
        d = float(r['y_dose']); f = F[r['sample']] = {'cohort': 'SGDP', 'set': 'sgdp_metadata_sex', 'y_dose': d, 'status': 'y_found' if d > NOY else 'no_y'}
    elif f is None and g:
        SKIP.add(f'{r["sample"]} (SGDP Gender {"".join(sorted(g))})'); continue
    if f is None: sys.exit(f'ERROR recounted {r["sample"]} not in the cohort table')
    dd = abs(f['y_dose'] - float(r['y_dose']))
    if dd > 0.00011: sys.exit(f'ERROR {r["sample"]}: y_dose {r["y_dose"]} in the recount, {f["y_dose"]} in the cohort table')   # both printed to 4 decimals: one unit apart is rounding
    if dd > 0.00005 and r['class'] == 'dyz3': ROUND.append(f'{r["sample"]} {r["y_dose"]} vs {f["y_dose"]}')
    f[r['class']] = float(r['estimate_bp'])
coh = {r['sample']: r['cohort'] for r in rows(HF)}
for r in rows(HE):
    c = coh.get(r['sample'], 'NA')
    if c.endswith('_control') or r['sample'] in F and F[r['sample']]['set'] == 'model_cohort':
        continue
    if c == 'NA':
        sys.exit(f'ERROR {r["sample"]}: no cohort')
    d = float(r['y_dose']); f = F.setdefault(r['sample'], {'cohort': c, 'set': 'held_out', 'y_dose': d, 'status': 'y_found' if d > NOY else 'no_y'})
    f[r['class']] = float(r['estimate_bp'])
part = [s for s, f in F.items() if 0 < sum(c in f for c in CLASSES) < len(CLASSES)]
if part: sys.exit(f'ERROR {len(part)} females with some classes only: {part[:5]}')
est = {s for s, f in F.items() if all(c in f for c in CLASSES)}
log(f'recounted but not female by recorded sex, left out: {sorted(SKIP)}')
log(f'recount y_dose one unit apart in the 4th decimal from the cohort table: {len(ROUND)} {ROUND}')
log(f'females with class estimates {len(est)} of {len(F)}; without (recount not done): {len(F) - len(est)}')
n = collections.Counter((f['cohort'], f['set']) for f in F.values())
log('females', len(F), dict(sorted(n.items())), 'status', dict(collections.Counter(f['status'] for f in F.values())))
with open(OUTF, 'w') as o:
    o.write('sample\tcohort\tset\ty_dose\ty_status\t' + '\t'.join(f'{c}_bp' for c in CLASSES) + '\n')
    for s in sorted(F):
        f = F[s]
        o.write(f'{s}\t{f["cohort"]}\t{f["set"]}\t{f["y_dose"]}\t{f["status"]}\t' + '\t'.join(f'{f[c]:.0f}' if c in f else 'NA' for c in CLASSES) + '\n')

T = collections.defaultdict(list)
for r in rows(TRUTH): T[r['class']].append(float(r['truth_bp']))
med = {c: float(np.median(T[c])) for c in CLASSES}
with open(OUTB, 'w') as o:
    o.write('set\tclass\tfemales\tpct_mean\tpct_median\tpct_p90\tpct_max\tmax_sample\tmax_cohort\tmedian_truth_bp\n')
    for name, drop in (('all', set()), ('without explained', set(EXPLAINED))):
        for c in CLASSES:
            S = sorted(s for s in est if s not in drop)
            p = np.array([100 * F[s][c] / med[c] for s in S]); top = S[int(np.argmax(p))]
            o.write(f'{name}\t{c}\t{len(p)}\t{p.mean():.4f}\t{np.median(p):.4f}\t{np.percentile(p, 90):.4f}\t{p.max():.4f}\t{top}\t{F[top]["cohort"]}\t{med[c]:.0f}\n')
            log(f'{name:17s} {c:7s} females {len(p)} % of median length: mean {p.mean():.3f} median {np.median(p):.3f} p90 {np.percentile(p, 90):.3f} max {p.max():.3f} ({top}, {F[top]["cohort"]})')
for s, why in EXPLAINED.items():
    if s in F: log(f'explained {s} ({F[s]["cohort"]}, y_dose {F[s]["y_dose"]}): {why}; ' + ', '.join(f'{c} {100 * F[s][c] / med[c]:.1f}%' for c in CLASSES))
log(f'wrote {OUTF} and {OUTB}')
