#!/usr/bin/env python3
"""SGDP samples without ENA BioSample sex in the cohort table: sex from the public SGDP metadata (Mallick 2016).

usage: python3 scripts/sgdp_sex_fill.py <cohort_hgdp_sgdp.tsv> <SGDP metadata> <out tsv>
SGDP metadata: data/sgdp_metadata.279public.21signedLetter.samples.txt (Gender column: M, F, U)
Samples named SGDP_NA_<BioSample> get their name from the EBI BioSamples API.
"""
import sys, csv, json, urllib.request
from collections import Counter

COH, META, OUT = sys.argv[1:4]
def log(*a): print(*a, flush=True)
L = [l.rstrip('\n').split('\t') for l in open(META, encoding='latin-1')]
gi = L[0].index('Gender'); G = {}
for x in L[1:]:
    for a in {x[2], *(a.strip() for a in x[3].split(','))}:
        G.setdefault(a, set()).add(x[gi])

def biosample_name(s):   # NA_<BioSample> aliases
    if not s.startswith('NA_SAMEA'): return s
    return json.load(urllib.request.urlopen(f'https://www.ebi.ac.uk/biosamples/samples/{s[3:]}.json', timeout=60))['name']

rows = [r for r in csv.DictReader(open(COH), delimiter='\t') if r['cohort'] == 'SGDP' and r['sex'] in ('', '.')]
log(f'SGDP samples without sex in {COH}: {len(rows)}')
out = []
for r in rows:
    n = biosample_name(r['sample'][5:]); g = G.get(n, set())
    if len(g) != 1: sys.exit(f'ERROR {r["sample"]} ({n}): SGDP Gender {sorted(g)}')
    g = g.pop(); y = float(r['y_dose']); xd = float(r['x_dose'])
    dosage = 'XY' if y > 0.025 else 'XX' if xd > 1.5 else 'X0'
    agree = {'M': 'XY', 'F': 'XX'}.get(g) == dosage
    out.append((r['sample'], n, g, r['y_dose'], r['x_dose'], dosage, 'yes' if agree else 'no', r['estimate_status']))
    log(f'{r["sample"]}\t{n}\tGender {g}\ty_dose {y}\tx_dose {xd}\t{dosage}\tagree {agree}\t{r["estimate_status"]}')
with open(OUT, 'w') as o:
    o.write('sample\tsgdp_name\tsgdp_gender\ty_dose\tx_dose\tdosage_karyotype\tgender_agrees\testimate_status\n')
    o.writelines('\t'.join(x) + '\n' for x in out)
log('Gender', dict(Counter(x[2] for x in out)), 'agree', dict(Counter(x[6] for x in out)), f'wrote {OUT}')
