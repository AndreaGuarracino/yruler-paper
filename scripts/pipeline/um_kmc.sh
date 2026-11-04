#!/bin/bash
# per assembly: whole-Y and class k-mer presence dbs
# usage: [UM=<marker dir>] [ONLY=<class>] sbatch -a 1-N um_kmc.sh; the whole db is kept when present
set -euo pipefail
P=${P:?set P to the project directory}; K=${KMC_BIN:?set KMC_BIN}; SAM=samtools
O=$P/${UM:-markers_u}; D=${KDB:-${SCRATCH:?set SCRATCH}/umk/db}; mkdir -p $D; read a f < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $O/assemblies.tsv)
T=$(mktemp -d ${SCRATCH:?set SCRATCH}/umk/t.XXXX); trap "rm -rf $T" EXIT; echo "host $(hostname) $a start $(date)"
db(){ $K/kmc -k31 -ci1 -cs255 -m4 -t4 -fm $1 $T/raw $T > /dev/null; $K/kmc_tools -t4 transform $T/raw set_counts 1 $2 > /dev/null; }
[ -e $D/$a.whole.kmc_pre ] || db $f $D/$a.whole
for c in ${ONLY:-$(awk -F"\t" -v a=$a '$1==a{print $2}' $O/regions.tsv | sort -u)}; do
  grep -q -P "^$a\t$c\t" $O/regions.tsv || continue
  awk -F"\t" -v a=$a -v c=$c '$1==a && $2==c{print $3":"$4+1"-"$5}' $O/regions.tsv > $T/rg; $SAM faidx -r $T/rg $f > $T/$c.fa; db $T/$c.fa $D/$a.$c
  echo "$c regions $(wc -l < $T/rg) bp $(grep -v ">" $T/$c.fa | tr -d "\n" | wc -c)"; done
echo "done $a $(date)"
