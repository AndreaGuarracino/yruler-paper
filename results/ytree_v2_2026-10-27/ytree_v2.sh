#!/bin/bash
# 1kGP male Y tree v2: callable mask, singletons kept, full search, supports
set -euo pipefail
P=/path/to/project; C=/path/to/conda; B=$C/bcftools/1.22/bin/bcftools; IQ=$(ls $C/iqtree/*/bin/iqtree* | head -1); PY=$C/ktrk-bench/1.0/bin/python
V=$P/yhaplo/vcf/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz; MASK=$P/yhaplo/hg/data/chrY_callable_regions.20130802_b38.bed
W=$P/ytree_v2; mkdir -p $W; cd $W; echo "iqtree $IQ ($($IQ --version | head -1)); bcftools $($B --version | head -1); mask $MASK ($(awk '{s+=$3-$2} END{print s}' $MASK) bp); start $(date)"
# biallelic PASS SNPs in the callable mask, males only; positions kept with genotype codes
$B view -S $P/yhaplo/males.txt --force-samples -R $MASK -v snps -m2 -M2 -f PASS -Ou $V | $B query -f '%CHROM\t%POS[\t%GT]\n' > gt.tsv
$B view -S $P/yhaplo/males.txt --force-samples -h $V | $B query -l > ids.txt
echo "sites $(wc -l < gt.tsv), males $(wc -l < ids.txt)"
$PY -I - <<'PYEOF'
import numpy as np
ids = [l.strip() for l in open('ids.txt')]; code = {'0': 0, '0/0': 0, '0|0': 0, '1': 1, '1/1': 1, '1|1': 1}; rows, pos = [], []
for l in open('gt.tsv'):
    x = l.rstrip('\n').split('\t'); g = np.array([code.get(v, 2) for v in x[2:]], np.uint8); n = g != 2; a = (g == 1).sum()
    if n.mean() >= 0.95 and 1 <= a <= n.sum() - 1: rows.append(g); pos.append((x[0], x[1]))   # call rate 95%, polymorphic (singletons kept)
M = np.array(rows).T; ok = (M == 2).mean(1) <= 0.5
print('dropped (more than 50% missing):', [s for s, k in zip(ids, ok) if not k], flush=True); ids = [s for s, k in zip(ids, ok) if k]; M = M[ok]
keep = ((M == 0).any(0) & (M == 1).any(0)); print('sites invariant after dropping males:', int((~keep).sum()), flush=True); M = M[:, keep]; pos = [p for p, k in zip(pos, keep) if k]   # +ASC needs variable sites only
mac = np.minimum((M == 1).sum(0), (M == 0).sum(0)); print('males', M.shape[0], 'sites', M.shape[1], 'singletons', int((mac == 1).sum()), 'missing', round(float((M == 2).mean()), 5), flush=True)
with open('aln.phy', 'w') as o:
    o.write(f'{M.shape[0]} {M.shape[1]}\n')
    for s, r in zip(ids, M): o.write(s + ' ' + ''.join('01?'[x] for x in r) + '\n')
with open('sites.tsv', 'w') as o:
    o.write('column\tCHROM\tPOS\tbuild\tmac\n')
    for i, ((c, p), m) in enumerate(zip(pos, mac), 1): o.write(f'{i}\t{c}\t{p}\tGRCh38\t{m}\n')
with open('missing_per_sample.tsv', 'w') as o:
    o.write('sample\tmissing_fraction\n')
    for s, r in zip(ids, M): o.write(f'{s}\t{(r == 2).mean():.5f}\n')
PYEOF
$IQ -s aln.phy -st BIN -m GTR2+FO+ASC -B 1000 -alrt 1000 -seed 20261027 -T ${SLURM_CPUS_PER_TASK:-32} --prefix ytree -redo > iqtree.log 2>&1 || { echo "ERROR iqtree"; tail -5 iqtree.log; exit 1; }
grep -E "Log-likelihood|Total tree length|Wall-clock" ytree.iqtree | head -5; echo "end $(date)"
