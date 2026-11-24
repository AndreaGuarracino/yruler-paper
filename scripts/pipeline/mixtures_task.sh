#!/bin/bash
#SBATCH -J v1-mixtures -c 24 --mem 100G -t 03:00:00 -o logs/%x.%A_%a.log
# Male plus female read mixtures with the v1 release: one male and one female 1kGP CRAM per task, subsampled and mixed, run through yruler
# usage: OUT=<v1 dir> [BASE=5.5] [FR="1 .5 .2 .1 .05 .03 .02 .01 0"] [DP="5 1 .2"] [SEEDS="1 2 3"] [PAR=4] sbatch -a 1-6 v1/pipeline/mixtures_task.sh (from $P; $IN/lists/mixtures_pairs.tsv: male, url, cov, female, url, cov)
#   -> $OUT/mixtures/base<BASE>/<male>_<female>/f<fr>_d<dp>_s<seed>/{estimates,summary}.tsv; v1: BASE=5.5 DP="5 1 .2" and BASE=10.5 DP=10
# Pairs stay together: --subsample hashes the read name; base BAMs at BASE x per seed; depth is that of the mixture; fraction is the male share of reads; fraction 0 runs: female only
set -euo pipefail; P=${P:?}; IN=${IN:?}; . ${PIPE:?}/config.sh; export SAM; A=$ARIA2; export REF=$IN/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa; export R=${OUT:?}/release; export Y=${YBIN:-$OUT/tool_build/release/yruler}
IFS=$'\t' read -r m mu mc f fu fc < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $IN/lists/mixtures_pairs.tsv); export m f BASE=${BASE:-5.5}; export O=${MIXDIR:-$OUT/mixtures/base$BASE}/${m}_$f D=$TMPBASE/mix_base${BASE}_${SLURM_ARRAY_JOB_ID}/${m}_$f; rm -rf "${D:?}"; mkdir -p $O $D; trap 'rm -rf "${D:?}"' EXIT; echo "host $(hostname) $m $mc x $f $fc start $(date)"
get() { if ! timeout 2400 $A -q -x 8 -s 8 -k 20M -c --max-tries=20 --retry-wait=10 --file-allocation=none -d $D -o $1.cram "$2"; then echo "aria2 failed for $1, curl"; rm -f "${D:?}/${1:?}".cram*; curl -sS --retry 5 -o $D/$1.cram "$2"; fi; $SAM quickcheck $D/$1.cram; echo "$1 $(stat -c %s $D/$1.cram) bytes"; }
get $m "$mu" & a=$!; get $f "$fu" & b=$!; wait $a; wait $b; echo "downloads done $(date)"
base() { local p=$(awk -v b=$BASE -v c=$2 'BEGIN{printf "%.6f", b / c}'); $SAM view -b -1 -F 0x900 --subsample $p --subsample-seed $3 -T $REF -@ 3 -o $D/$1.$3.bam $D/$1.cram; }
ok=0; pids=(); for s in ${SEEDS:-1 2 3}; do base $m $mc $s & pids+=($!); base $f $fc $s & pids+=($!); done; for j in "${pids[@]}"; do wait $j || ok=1; done; [ $ok = 0 ] || { echo "ERROR base bam"; exit 1; }; rm -f "${D:?}/${m:?}.cram" "${D:?}/${f:?}.cram"; echo "base bams done $(date) $(du -sh $D | cut -f1)"
has() { [ -s $1/estimates.tsv ] || [ -s $1/lengths.tsv ]; }   # Bx binary: estimates.tsv
run1() { local fr=$1 dp=$2 s=$3; local o=$O/f${fr}_d${dp}_s$s; local w=$D/w_${fr}_${dp}_$s; has $o && return 0; mkdir -p $o $w
  local pm=$(awk -v f=$fr -v d=$dp -v b=$BASE 'BEGIN{printf "%.6f", f * d / b}') pf=$(awk -v f=$fr -v d=$dp -v b=$BASE 'BEGIN{printf "%.6f", (1 - f) * d / b}'); local sd=$(( s * 100003 + $(echo "$fr$dp" | cksum | cut -d" " -f1) % 99991 ))
  ( h=-h; [ "$pm" != 0.000000 ] && { $SAM view $h --subsample $pm --subsample-seed $sd $D/$m.$s.bam; h=; }; [ "$pf" != 0.000000 ] && $SAM view $h --subsample $pf --subsample-seed $(( sd + 7 )) $D/$f.$s.bam; true ) 2> $o/view.err | (cd $w && /usr/bin/time -f "%e s %M kB" $Y -i - -m $R/markers.ms.fa.zst -M $R/model.tsv -n ${m}_$f -t 5 -o $o > /dev/null 2> $o/run.log) || { echo "ERROR run $o"; tail -3 $o/run.log; return 1; }
  [ -s $o/view.err ] && { echo "ERROR samtools $o: $(head -c 200 $o/view.err)"; return 1; }; has $o || { echo "ERROR no estimates $o"; return 1; }; rm -f "${o:?}/counts.tsv" "${o:?}/any.fa"; rm -rf "${w:?}"; echo "$(basename $o) $(tail -1 $o/run.log | cut -c1-100) $(grep -oE 'total [0-9.]+ s' $o/run.log | tail -1)"; }
export -f run1 has
for fr in ${FR:-1 .5 .2 .1 .05 .03 .02 .01 0}; do for dp in ${DP:-5 1 .2}; do for s in ${SEEDS:-1 2 3}; do echo "$fr $dp $s"; done; done; done | sort -k2,2nr | xargs -P ${PAR:-4} -L 1 bash -c 'run1 "$@"' _
k=$(( $(echo ${FR:-1 .5 .2 .1 .05 .03 .02 .01 0} | wc -w) * $(echo ${DP:-5 1 .2} | wc -w) * $(echo ${SEEDS:-1 2 3} | wc -w) )); n=$(find $O -mindepth 2 -maxdepth 2 -name estimates.tsv | wc -l); echo "runs with estimates $n of $k; end $(date)"; [ $n = $k ] || { echo "ERROR $n runs of $k"; exit 1; }
