#!/usr/bin/env python3
"""Recorded females of the three cohorts, one row per individual, for ED2a and the no-Y specificity numbers.

usage: python3 scripts/females_merge.py [out.tsv]   (default results/model/cohort_fem_all.tsv)
sources: results/model/cohort_fem.tsv (saved whole-file counts), results/checks/female_whole_file.tsv (final release on whole files,
reads checked against the 1kGP index or md5-checked downloads), results/checks/female_recount_final1.tsv (final release applied to the
whole-file recounts of md5-checked CRAMs), recorded females of results/model/cohort_hgdp_sgdp.tsv.
An individual sequenced in two cohorts counts once (SGDP aliases), kept in the first of 1kGP, HGDP, SGDP; within a cohort the first source above wins.
"""
import csv, sys, collections
from pathlib import Path
import prepare_cohort_pies as pc

ROOT = Path(__file__).resolve().parents[1]; rd = lambda p: list(csv.DictReader(open(ROOT / p), delimiter='\t'))
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'results/model/cohort_fem_all.tsv'
AL = collections.defaultdict(list)
for r in pc.read_tsv(pc.META / 'sgdp_metadata.txt', encoding='latin-1'):
    for k in pc.ALIAS_FIELDS:
        if pc.valid_key(r[k]): AL[pc.norm(r[k])].append(r)
def ident(s):   # one key per individual: the smallest normalized alias
    i = s.replace('SGDP_', ''); return min({pc.norm(i)} | {pc.norm(r[k]) for r in AL.get(pc.norm(i), []) for k in pc.ALIAS_FIELDS if pc.valid_key(r[k])})
rows = [dict(sample=r['sample'], cohort=r['cohort'], y_dose=r['y_dose'], y_status=r.get('y_status', '.'), source='cohort_fem') for r in rd('results/model/cohort_fem.tsv')]
rows += [dict(sample=r['sample'], cohort=r['cohort'], y_dose=r['y_dose'], y_status=r['y_status'], source=f"whole_file_{r['source']}") for r in rd('results/checks/female_whole_file.tsv')]
rows += [dict(sample=r['sample'], cohort=r['cohort'], y_dose=r['y_dose'], y_status=r['y_state'], source='recount_final1') for r in rd('results/checks/female_recount_final1.tsv')]
rows += [dict(sample=r['sample'], cohort=r['cohort'], y_dose=r['y_dose'], y_status=r['y_status'], source='cohort_hgdp_sgdp') for r in rd('results/model/cohort_hgdp_sgdp.tsv') if r['sex'] == 'female' and r['y_dose'] not in ('nan', 'NA', '.')]
order = {'1kGP': 0, 'HGDP': 1, 'SGDP': 2}; seen, keep, dup = {}, [], []
for r in sorted(rows, key=lambda r: order[r['cohort']]):   # stable: source order kept within a cohort
    k = ident(r['sample'])
    if k in seen: dup.append((r['sample'], r['source'], seen[k])); continue
    seen[k] = r['sample']; keep.append(r)
print(f"{len(rows)} rows, {len(keep)} individuals; dropped {len(dup)} repeats (same sample in two sources, or one individual in two cohorts)", flush=True)
print('individuals in two cohorts:', [d for d in dup if d[0].replace('SGDP_', '') != d[2].replace('SGDP_', '')], flush=True)
print('by cohort and source:', dict(collections.Counter((r['cohort'], r['source']) for r in keep)), flush=True)
y = [(float(r['y_dose']), r['sample']) for r in keep]; print('highest Y dosage:', sorted(y)[-4:], flush=True)
with OUT.open('w') as o:
    o.write('sample\tcohort\ty_dose\ty_status\tsource\n'); o.writelines(f"{r['sample']}\t{r['cohort']}\t{r['y_dose']}\t{r['y_status']}\t{r['source']}\n" for r in keep)
print('written', OUT, flush=True)
