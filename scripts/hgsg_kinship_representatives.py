#!/usr/bin/env python3
"""Choose one HGDP/SGDP individual per connected kinship group."""
import collections
import csv
from pathlib import Path
import prepare_cohort_pies as pc

ROOT = Path(__file__).resolve().parent.parent
rd = lambda p: list(csv.DictReader(open(ROOT / p), delimiter='\t'))
rows = rd('results/model/hgsg_calls.tsv')
ids = {r['sample'] for r in rows}
aliases = collections.defaultdict(list)
for r in pc.read_tsv(pc.META / 'sgdp_metadata.txt', encoding='latin-1'):
    for k in pc.ALIAS_FIELDS:
        if pc.valid_key(r[k]): aliases[pc.norm(r[k])].append(r)
def keys(s):
    i = s.replace('SGDP_', '')
    return {pc.norm(s), pc.norm(i)} | {pc.norm(r[k]) for r in aliases.get(pc.norm(i), []) for k in pc.ALIAS_FIELDS if pc.valid_key(r[k])}
lookup = {k: s for s in sorted(ids) for k in keys(s)}
parent = {s:s for s in ids}
def root(s):
    while parent[s] != s:
        parent[s] = parent[parent[s]]; s = parent[s]
    return s
for path in ('results/checks/relatedness/gnomad_related_pairs.tsv', 'results/checks/relatedness/sgdp_king_pairs.tsv'):
    for r in rd(path):
        a, b = lookup.get(pc.norm(r['s1'])), lookup.get(pc.norm(r['s2']))
        if a and b and a != b and float(r['kin']) > .05:
            parent[root(a)] = root(b)
groups = collections.defaultdict(list)
for s in sorted(ids): groups[root(s)].append(s)
output = ROOT / 'results/model/hgsg_kinship_representatives.tsv'
with output.open('w') as f:
    f.write('sample\tkinship_group\trepresentative\n')
    for v in sorted(groups.values(), key=lambda x:x[0]):
        for s in v: f.write(f'{s}\t{v[0]}\t{int(s == v[0])}\n')
print(f'{len(ids)} individuals; {len(groups)} groups; {output}')
