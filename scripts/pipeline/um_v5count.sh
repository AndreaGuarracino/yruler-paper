#!/bin/bash
# Marker counts per sample (default DYZ17: markers_v5/dyz17.fa)
# usage: [MK=<markers fa> OC=<out dir> BED=<regions, BAM/CRAM only>] sbatch --array=1-K -c 4 --mem=8G um_v5count.sh <list: sample, saved-reads dir or BAM/CRAM> [samples per task]; -> markers_v5/dyz17c/<sample>/ (counts.col.zst, summary.tsv, run.log)
set -euo pipefail
P=${P:?set P to the project directory}; Z=zstd; Y=${YRULER:-yruler}; M=${MK:-$P/markers_v5/dyz17.fa}; O=${OC:-$P/markers_v5/dyz17c}; R=$P/ext/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
N=${2:-1}; i0=$(( (SLURM_ARRAY_TASK_ID - 1) * N + 1 )); mkdir -p $O; e=0; echo "host $(hostname) task $SLURM_ARRAY_TASK_ID lines $i0-$(( i0 + N - 1 )) start $(date)"
while IFS=$'\t' read -r s f; do
  [ -s $O/$s/counts.col.zst ] && { echo "$s done"; continue; }; T=$(mktemp -d ${SCRATCH:?set SCRATCH}/v5c.XXXX); r=()
  if [[ $f == *.bam || $f == *.cram ]]; then i=$f; r=(-r $R ${BED:+-b $BED}); else
    [ -s $f/any.fa.zst ] || { echo "ERROR $s: no $f/any.fa.zst"; e=1; rm -rf $T; continue; }; grep -q "markers_perm2.fa" $f/run.log || { echo "ERROR $s: saved reads not made with markers_perm2.fa ($f/run.log)"; e=1; rm -rf $T; continue; }; $Z -dcq $f/any.fa.zst > $T/r.fa; i=$T/r.fa; fi
  [[ $i == /* || $i == *://* || $i == *,* ]] || i=$PWD/$i   # absolute input: yruler runs in $T, so remote indexes land there and not in $P
  if (cd $T && $Y -i $i "${r[@]}" -m $M -n $s -t ${SLURM_CPUS_PER_TASK:-4} -o $T/o > /dev/null 2> $T/run.log); then
    [ -s $M.order ] || { cut -f1 $T/o/counts.tsv > $M.order.$$ && mv -n $M.order.$$ $M.order; rm -f $M.order.$$; }
    [ "$(cut -f1 $T/o/counts.tsv | md5sum | cut -c1-32)" = "$(md5sum < $M.order | cut -c1-32)" ] || { echo "ERROR $s: marker order"; e=1; rm -rf $T; continue; }
    mkdir -p $O/$s.tmp; cut -f2 $T/o/counts.tsv | $Z -q -3 -o $O/$s.tmp/counts.col.zst; cp $T/o/summary.tsv $O/$s.tmp/; { echo "input $f"; cat $T/run.log; } > $O/$s.tmp/run.log; rm -rf $O/$s; mv $O/$s.tmp $O/$s
    echo "$s $(tail -1 $T/run.log | cut -c1-150)"
  else echo "ERROR $s: $(tail -2 $T/run.log)"; e=1; fi; rm -rf $T
done < <(sed -n "${i0},$(( i0 + N - 1 ))p" $1)
echo "end $(date)"; exit $e
