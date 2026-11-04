#!/bin/bash
# Download one file to keep, verified
# LIST rows: name, URL, md5 or -; OUT dir; index fetched if .cram/.bam
set -euo pipefail
P=${P:?set P to the project directory}; A=aria2c; SAM=samtools
read n u m < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST); f=$(basename $u); D=$OUT; mkdir -p $D/.work/$n; cd $D/.work/$n
echo "host $(hostname) $n $f md5 $m start $(date)"
[ -s $D/$n.$f ] && { echo "exists, skip"; exit 0; }
dl() { $A -q -x ${CONN:-16} -s ${CONN:-16} -k 20M -c --max-tries=0 --retry-wait=10 --file-allocation=none -o $1 "$u" >&2 || { echo "ERROR aria2 exit $? on $1" >&2; exit 1; }; md5sum $1 | cut -d' ' -f1; }
# check before keeping
ok() { case $f in *.cram|*.bam) $SAM quickcheck $1 || { echo "ERROR quickcheck $1"; exit 1; };; *.gz) gzip -t $1 || { echo "ERROR gzip $1"; exit 1; };; esac; mv $1 $D/$n.$f; }
a=$(dl c1) || exit 1; echo "copy 1 md5 $a"
if [ "$a" = "$m" ]; then ok c1
else b=$(dl c2) || exit 1; echo "copy 2 md5 $b"
  if [ "$b" = "$m" ] || { [ "$m" = - ] && [ "$a" = "$b" ]; }; then ok c2
  elif [ "$m" != - ]; then python3 $P/scripts/blockfix.py fixed $m "$u" c1 c2 && ok fixed
  else c=$(dl c3) || exit 1; echo "copy 3 md5 $c"; python3 $P/scripts/blockfix.py fixed - "$u" c1 c2 c3 && ok fixed; fi; fi
case $f in *.cram) curl -sS --fail --retry 5 -o $D/$n.$f.crai "$u.crai";; *.bam) curl -sS --fail --retry 5 -o $D/$n.$f.bai "$u.bai" || $SAM index $D/$n.$f;; esac
case $f in *.cram|*.bam) $SAM quickcheck $D/$n.$f || { echo "ERROR quickcheck"; exit 1; };; *.gz) gzip -t $D/$n.$f || { echo "ERROR gzip"; exit 1; };; esac
cd $D; rm -rf $D/.work/$n; echo "kept $D/$n.$f $(stat -c %s $D/$n.$f) bytes $(date)"
