#!/bin/bash
# one CRAM: download once, yruler whole (+ BED), per-marker counts
# env: LIST, OUT, BED, UM (marker dir), DLX (aria2 connections), POST (extra scripts on the local file, colon-separated), DLMAX (concurrent downloads), LOCKDIR (token dir), PASSMAX_GB (disk guard)
set -euo pipefail
P=${P:?set P to the project directory}; R=$P/ext/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
Y=${YRULER:-yruler}; A=aria2c; SAM=samtools; Z=zstd
M=$P/${UM:-markers_u}/markers_u.fa; OM=$(md5sum < $P/${UM:-markers_u}/counts_order.txt | cut -c1-32); T=${SLURM_CPUS_PER_TASK:-16}
IFS=$'\t' read -r s url md5 _ bedurl _ < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST); O=$P/$OUT/$s   # column 5: file for the BED run (DRAGEN BAM URL), else the local file
[ -s $O/whole/counts.col.zst ] && { [ -z "${BED:-}" ] || [ -s $O/bed/counts.col.zst ]; } && { [ -z "${POST:-}" ] || [ -e $O/post.ok ]; } && { echo "$s done"; exit 0; }
rm -rf $O/whole $O/bed $O/ychr $O/reads $O/srk $O/post.ok   # partial outputs of a failed run
D=${SCRATCH:?set SCRATCH}/yp_pass/$s; mkdir -p $D $O; TOK=; trap 'rm -rf $D; [ -z "$TOK" ] || rm -rf $TOK' EXIT; trap 'exit 143' TERM; echo "host $(hostname) $s start $(date)"
# disk guard: wait while yp_pass >= PASSMAX_GB
wait_disk() { local u; [ -n "${PASSMAX_GB:-}" ] || return 0; while u=$(du -s --block-size=1G ${SCRATCH:?set SCRATCH}/yp_pass | cut -f1); [ "$u" -ge "$PASSMAX_GB" ]; do echo "disk wait: yp_pass $u GB >= $PASSMAX_GB GB"; sleep 60; done; }
# download URL: PROTO=ftp swaps EBI https for ftp
dlurl() { [ "${PROTO:-}" = ftp ] && echo "${url/https:\/\/ftp.1000genomes.ebi.ac.uk/ftp://ftp.1000genomes.ebi.ac.uk}" || echo "$url"; }
# download token: at most DLMAX at once
take_token() { local i o st; LK=${LOCKDIR:-$P/locks/dl}; mkdir -p $LK; while :; do for i in $(seq 1 $DLMAX); do
  if mkdir $LK/$i 2>/dev/null; then echo ${SLURM_JOB_ID:-$$} > $LK/$i/owner; TOK=$LK/$i; echo "download token $i"; return; fi
  o=$(cat $LK/$i/owner 2>/dev/null) || continue; [ -n "$o" ] || continue; st=$(sacct -n -X -j $o -o State 2>/dev/null | head -1 | tr -d " ")
  case $st in COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_ME*|NODE_FAIL|PREEMPTED) echo "stale token $i of job $o ($st), removed"; rm -rf $LK/$i;; esac; done; sleep 30; done; }
t0=$(date +%s); pf=$(ls ${SCRATCH:?set SCRATCH}/yp_prefetch/*/$s.*.cram $P/prefetch/*/$s.*.cram 2>/dev/null | head -1 || true)
F=$D/$s.cram; if [[ $url == /* ]]; then F=$url; echo "local file $url (read in place, not deleted)"
elif [ -n "$pf" ] && [ -s "$pf.crai" ]; then mv "$pf" $D/$s.cram; mv "$pf.crai" $D/$s.cram.crai; echo "prefetched $pf (md5 checked at download)"   # kept by keep_dl.sh
else wait_disk; [ -z "${DLMAX:-}" ] || take_token
  # md5 mismatch: second copy, then blockfix.py repairs from both copies
  dl() { $A -q -x ${DLX:-8} -s ${DLX:-8} -k 20M -c --max-tries=0 --retry-wait=10 --file-allocation=none -d $D -o $1 "$(dlurl)" >&2 || { echo "ERROR aria2 $s"; exit 1; }; md5sum < $D/$1 | cut -c1-32; }
  m1=$(dl c1); if [ "$md5" = "-" ] || [ "$m1" = "$md5" ]; then mv $D/c1 $D/$s.cram
  else echo "md5 mismatch $s copy 1 on $(hostname) ($(stat -c %s $D/c1) bytes)"; m2=$(dl c2)
    if [ "$m2" = "$md5" ]; then mv $D/c2 $D/$s.cram; rm -f $D/c1
    else echo "md5 mismatch $s copy 2 on $(hostname)"; python3 $P/scripts/blockfix.py $D/$s.cram $md5 "$url" $D/c1 $D/c2 || { echo "ERROR md5 $s: blockfix failed"; exit 1; }
      [ "$(md5sum < $D/$s.cram | cut -c1-32)" = "$md5" ] || { echo "ERROR md5 $s after blockfix"; exit 1; }; echo "repaired $s with blockfix"; rm -f $D/c1 $D/c2; fi; fi
  curl -sS --fail --retry 5 -o $D/$s.cram.crai "$url.crai" || $SAM index $D/$s.cram
  [ -s $D/$s.cram.crai ] || { echo "index missing or empty for $s: samtools index"; rm -f $D/$s.cram.crai; $SAM index $D/$s.cram; }
  [ -z "$TOK" ] || { rm -rf $TOK; TOK=; }; fi
$SAM quickcheck $F || { echo "ERROR quickcheck $s"; exit 1; }
$SAM idxstats $F > /dev/null 2>&1 || { echo "index unreadable for $s: samtools index"; rm -f $F.crai $F.bai; $SAM index $F; }
echo "download $(( $(date +%s) - t0 )) s, $(stat -c %s $F) bytes, md5 ok"
for m in whole ${BED:+bed}; do o=$O/$m; mkdir -p $o; b=""; in=$F; [ $m = bed ] && { b="-b $BED"; in=${bedurl:-$F}; }
  for try in 1 2 3; do   # remote reads can break: retry
    /usr/bin/time -f "%e s %M kB %P cpu" $Y -i $in -r $R -m $M -n $s -t $T $b -o $o 2> $o/run.log > /dev/null && break
    echo "yruler $s $m try $try failed: $(grep -m1 ERROR $o/run.log | cut -c1-160)"
    [ $m = whole ] && grep -q "index" $o/run.log && { echo "rebuilding the index of $s"; rm -f $F.crai $F.bai; $SAM index $F; }   # local index unusable
    [ $try = 3 ] && { echo "ERROR yruler $s $m after 3 tries"; tail -3 $o/run.log; exit 1; }; sleep 60; done
  [ "$(cut -f1 $o/counts.tsv | md5sum | cut -c1-32)" = "$OM" ] || { echo "ERROR marker order $s $m"; exit 1; }
  cut -f2 $o/counts.tsv | $Z -q -3 -o $o/counts.col.zst; rm -f $o/counts.tsv $o/any.fa; gzip -f $o/bins.tsv
  echo "$m input $in $(tail -1 $o/run.log) $(grep -o "records read (primary) [0-9]*" $o/run.log) $(head -1 $o/summary.tsv | cut -f1)"; done
# chrY genotypes at ISOGG 2019 sites (haploid) for HaploGrouper
G=$P/yhaplo/chrY_isogg2019-decode1_loci_b38.reg; BC=bcftools; mkdir -p $O/ychr
t1=$(date +%s); $BC mpileup -f $R -r chrY -T $G -a AD,DP -Ou $F 2> $O/ychr/mpileup.log | $BC call -m --ploidy 1 -Oz -o $O/ychr/$s.chrY.isogg.vcf.gz 2> $O/ychr/call.log || { echo "ERROR chrY genotypes $s"; exit 1; }
echo "chrY genotypes $(( $(date +%s) - t1 )) s, records $($BC view -H $O/ychr/$s.chrY.isogg.vcf.gz | wc -l)"
# extra untimed scripts on the local file (POST, colon-separated)
if [ -n "${POST:-}" ]; then IFS=: read -ra H <<< "$POST"; for h in "${H[@]}"; do t2=$(date +%s); bash $h $s $F $O || { echo "ERROR post $h $s"; exit 1; }; echo "post $(basename $h) $(( $(date +%s) - t2 )) s"; done; touch $O/post.ok; fi
echo "end $s $(date)"
