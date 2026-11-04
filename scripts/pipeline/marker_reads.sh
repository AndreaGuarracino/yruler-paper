#!/bin/bash
# marker reads: separate untimed pass, permissive markers
set -euo pipefail
s=$1; f=$2; O=$3/reads; P=${P:?set P to the project directory}; 
Y=${YRULER:-yruler}; R=$P/ext/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; Z=zstd
M=${READS_MARKERS:?READS_MARKERS not set}; mkdir -p $O; echo "marker reads $s markers $M"
/usr/bin/time -f "%e s %M kB %P cpu" $Y -i $f -r $R -m $M -n $s -t ${SLURM_CPUS_PER_TASK:-8} -o $O 2> $O/run.log > /dev/null || { echo "ERROR yruler reads $s"; tail -3 $O/run.log; exit 1; }
[ -s ${M%.fa}.order ] || { cut -f1 $O/counts.tsv > ${M%.fa}.order.$$ && mv -n ${M%.fa}.order.$$ ${M%.fa}.order; [ -e ${M%.fa}.order.$$ ] && rm -f ${M%.fa}.order.$$ || echo "marker order written by $s"; }   # first run sets the order
[ "$(cut -f1 $O/counts.tsv | md5sum | cut -c1-32)" = "$(md5sum < ${M%.fa}.order | cut -c1-32)" ] || { echo "ERROR marker order $s"; exit 1; }
n=$(grep -c '^>' $O/any.fa || true); cut -f2 $O/counts.tsv | $Z -q -3 -o $O/counts.col.zst; rm -f $O/counts.tsv; gzip -f $O/bins.tsv; $Z -q -10 -T4 --rm $O/any.fa
echo "reads $(tail -1 $O/run.log | cut -c1-200); records $n, any.fa.zst $(stat -c %s $O/any.fa.zst) bytes, counts $(stat -c %s $O/counts.col.zst) bytes"
