#!/usr/bin/env python3
"""Final-model cohort table (cohort_1kgp.tsv) to the snapshot table read by figures_1kgp.py.

usage: python3 scripts/cohort_to_snapshot.py <cohort_1kgp.tsv> <out estimates_1kgp.tsv> <CV depth_f1.tsv>
flag: ok -> ok, low_y_dose -> y_dose_below_0.9, female -> female, no_y (male at or below the no-Y cut) -> no_y; sex_ped: male 1, female 2.
comparable: no = arm-specific profile (yp_deletion, yq_loss, partial), no_Yq12 = par2_gain (Yq12, DYZ1, DYZ2 not comparable), yes otherwise, . = no estimate.
outside_assembly_range: classes whose estimate lies outside the assembly lengths of the CV men.
"""
import csv, sys, collections

COH, OUT, CV = sys.argv[1:4]
COLS = ['sample', 'sex_ped', 'population', 'superpopulation', 'batch', 'haplogroup', 'lineage', 'y_dose', 'DYZ3_bp', 'Yq12_bp', 'DYZ1_bp', 'DYZ2_bp', 'DYZ18S_bp', 'DYZ19_bp', 'DYZ17_bp', 'reference_male', 'flag', 'yp_ruler_dose', 'yq_ruler_dose', 'y_state', 'yp_deletion_bias_pct', 'x_dose', 'sex_chromosomes', 'par1_dose', 'par2_dose', 'library', 'y_status', 'comparable', 'outside_assembly_range']
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
RG = collections.defaultdict(list)
for r in csv.DictReader(open(CV), delimiter='\t'): RG[NM[r['class']]].append(float(r['reference_bp']))
RG = {c: (min(v), max(v)) for c, v in RG.items()}; print('assembly ranges', {c: (round(a), round(b)) for c, (a, b) in RG.items()}, flush=True)
FLAG = {'ok': 'ok', 'low_y_dose': 'y_dose_below_0.9', 'female': 'female', 'no_y': 'no_y'}
rows = list(csv.DictReader(open(COH), delimiter='\t')); bad = [r['sample'] for r in rows if r['estimate_status'] not in FLAG or r['sex'] not in ('male', 'female')]
if bad: sys.exit(f'ERROR unmapped status or sex: {len(bad)} samples, first {bad[:5]}')
for r in rows:
    r['sex_ped'] = {'male': '1', 'female': '2'}[r['sex']]; r['flag'] = FLAG[r['estimate_status']]; est = r['y_status'] == 'y_found' and r['sex'] == 'male'
    r['comparable'] = '.' if not est else 'no' if r['y_state'] in ('yp_deletion', 'yq_loss', 'partial') else 'no_Yq12' if r['y_state'] == 'par2_gain' else 'yes'
    r['outside_assembly_range'] = ';'.join(c for c, (a, b) in RG.items() if not a <= float(r[f'{c}_bp']) <= b) or '.' if est else '.'
with open(OUT, 'w', newline='') as o:
    w = csv.DictWriter(o, COLS, delimiter='\t', extrasaction='ignore', lineterminator='\n'); w.writeheader(); w.writerows(rows)
print(f'{len(rows)} samples, flags {dict(collections.Counter(r["flag"] for r in rows))}, comparable {dict(collections.Counter(r["comparable"] for r in rows))}, men with a class outside the range {sum(r["outside_assembly_range"] != "." for r in rows)} -> {OUT}', flush=True)
