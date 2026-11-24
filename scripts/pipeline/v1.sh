#!/bin/bash
# Paper model v1: one driver, steps in run order; outputs in the layout of v1/ under OUT
# usage: OUT=<results dir, e.g. $P/v1/results> [IDX=<task indices, e.g. 0,100>] [DEP=afterok:<job>] bash v1/pipeline/v1.sh <step> (from the project dir, login node: submits only; site settings in config.sh); inputs in $P/v1/inputs, scripts in $P/v1/pipeline
#   steps in run order: model, cv, lolo, holdout, range, depth, nanopore, cohorts, release, mixtures, export
#   model, cv, lolo, holdout, range, depth, nanopore: Slurm arrays of cv_task.sh, then summary.py or model.py; cohorts, release, mixtures: cohort_task.sh, release_task.sh, mixtures_task.sh
set -euo pipefail
trap 'echo "EXIT line ${LINENO} status $?"' EXIT
log() { echo "[$(date '+%F %T %Z')] $*"; }
ST="model cv lolo holdout range depth nanopore cohorts release mixtures export"; S=${1:-}; [[ " $ST " == *" $S "* ]] || { echo "ERROR usage: bash v1/pipeline/v1.sh <$ST>"; exit 1; }
export P=$PWD; export PIPE=$P/v1/pipeline IN=$P/v1/inputs; . $PIPE/config.sh; export OUT=${OUT:?OUT (results dir) is required}; U=v1/inputs/markers
mkdir -p $OUT/{cv,lineage_out,batch_holdout,range,depth,model/markers,nanopore/cv,nanopore/model/markers} $P/logs/ids; R=$(realpath -m --relative-to=$P/$U $OUT)   # cv_eval.py paths are relative to $P/$U
# v1 settings
export UM=$U FULL=0.15 RULERS=y90 THRS=0.5,0.6,0.7,0.8,0.9,0.95 OUTFS=0,0.05,0.1,0.25 INNER=5 INNERSCORE=mean OUTLAB=1 \
  RULEX=$IN/markers/ruler_exclusion.npy CLSX=$IN/markers/class_exclusion.npy D17T=$IN/truth/dyz17_truth.tsv D17C=dyz17_span \
  TRUTHX=$IN/truth/exclusions.txt MARKX=1 CNPY=../counts/males NYSCORE=$IN/markers/nony_count.npy NYN=71 NYSTRICT=1 PED=$IN/lists/ped.txt \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset BYREP MINGC CVDIR LOLODIR HOLDDIR RANGEDIR FINAL OUTLAB_CLASSES ALTDIR
OLC=dyz1,dyz2,dyz18s,yq12rest2   # v1/results/cv, lineage_out, depth, nanopore cv: OUTLAB=1 for these classes, the other classes from the run without OUTLAB (as v1); model, holdout, range: OUTLAB=1 for all classes
G=dyz2:14; F=$IN/counts/females.npz; export EXCL=$IN/lists/external_males.txt   # females: d2 = chr8:91-101 Mb window mean
CI=(dyz3 yq12 dyz1 dyz2 dyz18s yq2k7 dyz18sub dyz19 dyz17 yq12rest yq2k7a dyz18suba dyz17a yq12rest2); NC=${#CI[@]}; export LINLIST=$IN/lists/lineages.txt; NL=$(wc -l < $LINLIST)
LV="${LEVELS:-full d30x d25x d20x d15x d10x d5x d2x d1x}"; CH8=$IN/windows/chr8_depth_levels.tsv   # depth levels; chr8 window 31-mer depth per level and sample
fi_() { local i o=""; for i in "${!CI[@]}"; do [ ${CI[$i]} = ${G%%:*} ] && o="$(( i * $1 ))-$(( i * $1 + $1 - 1 ))"; done; echo $o; }   # task range of the floored class (block size $1)
ck() { local x; x=$(sbatch --parsable "$@" 2>&1); echo "$x" | grep "^[0-9]" || { echo "SUBMIT ERROR: $(echo "$x" | grep -i error | head -1)" >&2; exit 1; }; }
ev() { local n=$1 k=${PACK:-1} f=$P/logs/ids/v1$1.$(date +%s%N).txt; echo $2 | tr ',' '\n' | awk -F- '{for (j = $1; j <= ($2 == "" ? $1 : $2); j++) print j}' > $f; shift 2   # cv_task.sh array on the indices of $2
  IDSF=$f PACK=$k ck -J v1-$n -o logs/v1-$n.%A_%a.log --array=0-$(( ($(wc -l < $f) + k - 1) / k - 1 )) -c $(( ${EVC:-3} * k )) --mem=$(( ${EVM:-40} * k ))G -t ${EVT:-00:40:00} ${DEP:+-d $DEP} "$@" --export=ALL $PIPE/cv_task.sh; }
sm() { ck -J v1-sum -o logs/v1-sum.%j.log -c 1 --mem=8G -t 01:00:00 "$@" --export=ALL --wrap "$W"; }
log "step $S OUT $OUT"; J=""; W=true
case $S in
  model) J=":$(MODE=production PRE=$R/model/markers/{c} ev prod ${IDX:-0-$(( NC - 1 ))})"; J="$J:$(MODE=production PRE=$R/model/markers/{c}_gcfloor MINGC=$G ev prodfloor $(fi_ 1))"
    W="python3 $PIPE/model.py '$R/model/markers/{c}' '$R/model/markers/{c}_gcfloor' $G $OUT/model";;
  cv) export OUTLAB_CLASSES=$OLC; J=":$(MODE=cv FEMNPZ=$F CVDIR=$R/cv/folds ev cv ${IDX:-0-$(( NC * 50 - 1 ))})"; [ -z "${IDX:-}" ] && J="$J:$(MODE=cv FEMNPZ=$F CVDIR=$R/cv/folds_dyz2_gcfloor MINGC=$G ev cvfloor $(fi_ 50))"
    W="python3 $PIPE/summary.py cv $R/cv/folds $R/cv/folds_dyz2_gcfloor $OUT/cv/cv.tsv";;
  lolo) export OUTLAB_CLASSES=$OLC; J=":$(MODE=lolo LOLODIR=$R/lineage_out/folds ev lolo ${IDX:-0-$(( NC * NL - 1 ))})"; [ -z "${IDX:-}" ] && J="$J:$(MODE=lolo LOLODIR=$R/lineage_out/folds_dyz2_gcfloor MINGC=$G ev lolofloor $(fi_ $NL))"
    W="python3 $PIPE/summary.py lolo $R/lineage_out/folds $R/lineage_out/folds_dyz2_gcfloor $OUT/lineage_out/lolo.tsv";;
  holdout) export EXCL=$IN/lists/batch_holdout_698.txt; J=":$(MODE=holdout HOLDDIR=$R/batch_holdout/folds HOLDLIST=$EXCL ev hold ${IDX:-0-$(( NC - 1 ))})"   # batch holdout: train 2,504, test 698
    [ -z "${IDX:-}" ] && J="$J:$(MODE=holdout HOLDDIR=$R/batch_holdout/folds_dyz2_gcfloor HOLDLIST=$EXCL MINGC=$G ev holdfloor $(fi_ 1))"
    W="python3 $PIPE/summary.py holdout $R/batch_holdout/folds $R/batch_holdout/folds_dyz2_gcfloor $OUT/batch_holdout/holdout.tsv";;
  range) export MINGC=$G; unset THRS OUTFS INNER NYSCORE NYN NYSTRICT   # v1 ran the length test with the cv_eval.py defaults for these: no non-Y genome filter, no inner rule choice
    J=":$(MODE=range RANGEDIR=$R/range EVT=00:50:00 ev range ${IDX:-2,3,4})"; J="$J:$(MODE=range RANGEDIR=$R/range/yq12 RANGEBY=yq12 EVT=00:50:00 ev rangeyq12 ${IDX:-2,3,4,13})";;
  depth) export OUTLAB_CLASSES=$OLC BYREP=1; W=""   # one task per class and repeat; test counts of each level in v1/inputs/counts/depth/<level>
    for l in $LV; do J="$J:$(MODE=cv FEMNPZ=$F ALTDIR=v1/inputs/counts/depth/$l CVDIR=$R/depth/folds/$l EVT=${DEVT:-02:00:00} ev d$l ${IDX:-0-$(( NC * 10 - 1 ))})"
      [ -z "${IDX:-}" ] && J="$J:$(MODE=cv FEMNPZ=$F ALTDIR=v1/inputs/counts/depth/$l CVDIR=$R/depth/folds_dyz2_gcfloor/$l MINGC=$G EVT=${DEVT:-02:00:00} ev df$l $(fi_ 10))"
      W="$W python3 $PIPE/summary.py cv $R/depth/folds/$l $R/depth/folds_dyz2_gcfloor/$l $OUT/depth/cv_depth_$l.tsv &&"; done
    W="$W python3 $PIPE/summary.py depth $OUT/depth $OUT/depth/depth.tsv $LV";;
  nanopore) export CNPY=../counts/nanopore OUTLAB_CLASSES=$OLC; J=":$(MODE=cv FEMNPZ=$F CVDIR=$R/nanopore/cv/folds EVM=28 ev ontcv ${IDX:-0-$(( NC * 50 - 1 ))})"   # nanopore counts, 82 R10 males
    [ -z "${IDX:-}" ] && J="$J:$(MODE=cv FEMNPZ=$F CVDIR=$R/nanopore/cv/folds_dyz2_gcfloor MINGC=$G EVM=28 ev ontcvfloor $(fi_ 50))"
    J="$J:$(OUTLAB_CLASSES= MODE=production PRE=$R/nanopore/model/markers/{c} EVM=28 ev ontprod 0-$(( NC - 1 ))):$(OUTLAB_CLASSES= MODE=production PRE=$R/nanopore/model/markers/{c}_gcfloor MINGC=$G EVM=28 ev ontprodfloor $(fi_ 1))"
    W="python3 $PIPE/summary.py cv $R/nanopore/cv/folds $R/nanopore/cv/folds_dyz2_gcfloor $OUT/nanopore/cv/cv.tsv && python3 $PIPE/model.py '$R/nanopore/model/markers/{c}' '$R/nanopore/model/markers/{c}_gcfloor' $G $OUT/nanopore/model";;
  cohorts) J=":$(ck ${DEP:+-d $DEP} --export=ALL $PIPE/cohort_task.sh)";;
  release) J=":$(ck ${DEP:+-d $DEP} --export=ALL $PIPE/release_task.sh):$(TECH=ont ck ${DEP:+-d $DEP} --export=ALL $PIPE/release_task.sh)";;
  mixtures) J=":$(BASE=5.5 DP='5 1 .2' ck -a 1-6 ${DEP:+-d $DEP} --export=ALL $PIPE/mixtures_task.sh):$(BASE=10.5 DP=10 ck -a 1-6 ${DEP:+-d $DEP} --export=ALL $PIPE/mixtures_task.sh)"
    W="python3 $PIPE/mixtures.py $OUT/cohorts/cohort_1kgp.tsv $OUT/mixtures/base5.5 $OUT/mixtures/base5.5.tsv $OUT/mixtures/base5.5_summary.tsv && python3 $PIPE/mixtures.py $OUT/cohorts/cohort_1kgp.tsv $OUT/mixtures/base10.5 $OUT/mixtures/base10.5.tsv $OUT/mixtures/base10.5_summary.tsv && (cat $OUT/mixtures/base5.5.tsv; tail -n +2 $OUT/mixtures/base10.5.tsv) > $OUT/mixtures/mixtures.tsv";;
  export) W="CHR8DEPTH=$CH8 python3 $PIPE/summary.py export $OUT $OUT/depth/figexport";;
esac
[ -n "${IDX:-}" ] && [ $S != model ] && W="true"   # partial runs (tests): no summary
[ "$W" = true ] && [ -n "$J" ] && { log "step $S: jobs${J//:/ }"; exit 0; }
if [ -n "$J" ]; then D2="-d afterok$J"; else D2=${DEP:+-d $DEP}; fi; S1=$(sm $D2); log "step $S: jobs${J//:/ }, then $S1"
