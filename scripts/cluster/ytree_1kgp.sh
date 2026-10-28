#!/bin/bash
# 1kGP male Y tree from genotype codes
set -euo pipefail
P=/path/to/project; C=/path/to/conda; B=$C/bcftools/1.22/bin/bcftools; IQ=$(ls $C/iqtree/*/bin/iqtree* | head -1)
V=$P/yhaplo/vcf/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz; W=$P/ytree; mkdir -p $W; cd $W
echo "iqtree $IQ; start $(date)"
# biallelic PASS SNPs, male-specific region, males only
$B view -S $P/yhaplo/males.txt --force-samples -r chrY:2781480-56887902 -v snps -m2 -M2 -f PASS -Ou $V | $B query -f '[%GT\t]\n' > gt.tsv
$B query -l -S $P/yhaplo/males.txt --force-samples $V > ids.txt 2>/dev/null || $B view -S $P/yhaplo/males.txt --force-samples -h $V | $B query -l > ids.txt
echo "sites $(wc -l < gt.tsv), males $(wc -l < ids.txt)"
/path/to/conda/ktrk-bench/1.0/bin/python -I - <<'PY'
import numpy as np
ids = [l.strip() for l in open('ids.txt')]; code = {'0': 0, '0/0': 0, '0|0': 0, '1': 1, '1/1': 1, '1|1': 1}
rows = []
for l in open('gt.tsv'):
    g = np.array([code.get(x, 2) for x in l.rstrip('\t\n').split('\t')], np.uint8)
    n = g != 2; a = (g == 1).sum()
    if n.mean() >= 0.95 and 2 <= a <= n.sum() - 2: rows.append(g)   # call rate 95%, informative
M = np.array(rows).T; ok = (M == 2).mean(1) <= 0.5   # drop males without Y calls
print('dropped (more than 50% missing):', [s for s, k in zip(ids, ok) if not k], flush=True); ids = [s for s, k in zip(ids, ok) if k]; M = M[ok]
print('males', M.shape[0], 'informative sites', M.shape[1], 'missing', round(float((M == 2).mean()), 5), flush=True)
with open('aln.phy', 'w') as o:
    o.write(f'{M.shape[0]} {M.shape[1]}\n')
    for s, r in zip(ids, M): o.write(s + ' ' + ''.join('01?'[x] for x in r) + '\n')
PY
$IQ -s aln.phy -st BIN -m GTR2+FO+ASC -fast -T ${SLURM_CPUS_PER_TASK:-32} --prefix ytree -redo > iqtree.log 2>&1 || { echo "ERROR iqtree"; tail -5 iqtree.log; exit 1; }
grep -E "Log-likelihood|Total tree length|Wall-clock" ytree.iqtree | head -5; echo "end $(date)"
