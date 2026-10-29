#!/bin/bash
# truth males: yruler with superset markers, per-marker counts kept
set -euo pipefail
P=/path/to/project; W=$P/ext; Y=$P/tools/yruler-src/target-curl13/release/yruler; ZS=/path/to/conda/zstd
Z=$(ls $ZS/*/bin/zstd 2>/dev/null | head -1); [ -x "$Z" ] || { echo "ERROR no zstd under $ZS"; exit 1; }
read s f < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${LIST:-$P/lists/ascairn_males.txt}); O=$P/${OUTD:-yruler_u}/$s
[ -s $O/counts.col.zst ] && { echo "$s exists"; exit 0; }; mkdir -p $O; echo "host $(hostname) $s $f start $(date)"
/usr/bin/time -f "%e s %M kB %P cpu" $Y -i $f -r $W/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -m $P/${UM:-markers_u}/markers_u.fa -n $s -t ${SLURM_CPUS_PER_TASK:-16} ${BED:+-b $BED} -o $O 2> $O/run.log > /dev/null || { echo "ERROR $s"; tail -3 $O/run.log; exit 1; }
tail -1 $O/run.log; [ -s $P/${UM:-markers_u}/counts_order.txt ] || { cut -f1 $O/counts.tsv > $O/order.tmp && mv $O/order.tmp $P/${UM:-markers_u}/counts_order.txt; }
cut -f2 $O/counts.tsv | $Z -q -3 -o $O/counts.col.zst; rm -f $O/any.fa $O/counts.tsv; gzip -f $O/bins.tsv; echo end $(date)
