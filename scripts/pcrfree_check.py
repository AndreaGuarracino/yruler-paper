# PCR-free labels against published library types
# usage: python3 scripts/pcrfree_check.py <cram_meta_hgdp_sgdp.tsv> <hgdp_wgs.20190516.metadata.txt> <ENA PRJEB9586 read_run tsv> results/checks/pcrfree_labels.tsv
#   rule as um_apply.py: HGDP >= 2 read groups, SGDP LP60 plates
import sys, csv, collections
cram, meta, ena, out = sys.argv[1:5]
M = {r['sample']: r['library_type'] for r in csv.DictReader(open(meta), delimiter='\t')}
E = collections.defaultdict(set)
for r in csv.DictReader(open(ena), delimiter='\t'): E[r['sample_alias']].add(r['library_construction_protocol'] or 'blank')
rows = []
for r in csv.DictReader(open(cram), delimiter='\t'):
    s = r['sample']
    if s.startswith('HGDP'): rule, pub = ('PCR-free' if int(r['read_groups']) >= 2 else 'PCR'), M.get(s, 'missing')
    else: rule, pub = ('PCR-free' if r['libraries'].startswith('LP60') else 'unknown'), '|'.join(sorted(E.get(s.replace('SGDP_', ''), {'missing'})))
    rows.append((s, 'HGDP' if s.startswith('HGDP') else 'SGDP', r['read_groups'], r['libraries'].split(',')[0][:4], rule, pub))
for co in ('HGDP', 'SGDP'):
    for k, v in sorted(collections.Counter((x[4], x[5]) for x in rows if x[1] == co).items()): print(f'{co}: rule {k[0]}, published {k[1]}: {v}', flush=True)
with open(out, 'w') as w:
    w.write('sample\tcohort\tread_groups\tplate_prefix\trule_label\tpublished_library\n'); w.writelines('\t'.join(x) + '\n' for x in rows)
print(f'wrote {out}: {len(rows)} CRAM files', flush=True)
