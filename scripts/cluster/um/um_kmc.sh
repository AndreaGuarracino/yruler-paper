#!/bin/bash
# per assembly: whole-Y and class k-mer presence dbs
set -euo pipefail
P=/path/to/project; C=/path/to/conda; K=$C/kmc/3.2.4/bin; SAM=$C/samtools/1.24/bin/samtools
O=$P/markers_u; D=/path/to/scratch/umk/db; mkdir -p $D; read a f < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $O/assemblies.tsv)
T=$(mktemp -d /path/to/scratch/umk/t.XXXX); trap "rm -rf $T" EXIT; echo "host $(hostname) $a start $(date)"
db(){ $K/kmc -k31 -ci1 -cs255 -m4 -t4 -fm $1 $T/raw $T > /dev/null; $K/kmc_tools -t4 transform $T/raw set_counts 1 $2 > /dev/null; }
db $f $D/$a.whole
for c in $(awk -F"\t" -v a=$a '$1==a{print $2}' $O/regions.tsv | sort -u); do
  awk -F"\t" -v a=$a -v c=$c '$1==a && $2==c{print $3":"$4+1"-"$5}' $O/regions.tsv > $T/rg; $SAM faidx -r $T/rg $f > $T/$c.fa; db $T/$c.fa $D/$a.$c
  echo "$c regions $(wc -l < $T/rg) bp $(grep -v ">" $T/$c.fa | tr -d "\n" | wc -c)"; done
echo "done $a $(date)"
