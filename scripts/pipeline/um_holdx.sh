#!/bin/bash
# Holdout of extra males with the final model procedure
# usage: [DEP=afterok:<job>] [UM=markers_v6] [SUMD=model_v5_v6d14_sum] bash um_holdx.sh <holdout list> <extra cnpy dir in UM> <tag>; counts dir cnpy_<tag> = cnpy + extra males; training = every CV male outside the holdout families (same env as um_v5run.sh TAG=v6, DYZ2 floor 14 as v6d14)
set -euo pipefail
P=${P:?set P to the project directory}; cd $P; HL=$1; XC=$2; T=$3; U=${UM:-markers_v6}; SD=${SUMD:-model_v5_v6d14_sum}; CD=$P/$U/cnpy_$T
mkdir -p $CD; for f in $P/$U/cnpy/*.npy $P/$U/$XC/*.npy; do ln -sf $f $CD/; done; echo "counts dir $CD: $(ls $CD | wc -l) males; holdout $(wc -l < $HL) males"
for s in $(cat $HL); do [ -e $CD/$s.npy ] || { echo "ERROR no counts for $s"; exit 1; }; done
export UM=$U FULL=0.15 RULERS=y90 THRS=0.5,0.6,0.7,0.8,0.9,0.95 OUTFS=0,0.05,0.1,0.25 TRUTHX=$P/truth_qc/x/rule_w4.txt TRUTHSET= COPIES=${COPIES-$P/$U/copies_${U#markers_}} INNER=5 NYSCORE=$P/$U/pg_topFM.noours.count.npy NYN=71 NYSTRICT=1 NYRULER=$P/$U/nyruler_none.npy OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FILT= FEMT=inf EXCL=$P/markers_v5/excl_external.txt CNPY=cnpy_$T HOLDLIST=$P/$HL MODE=holdout
unset BYREP MINGC
sb() { sbatch --parsable -c 6 --mem=40G -t 00:15:00 ${DEP:+-d $DEP} --export=ALL "$@" scripts/um_cvsplit.sh; }
B=$(HOLDDIR=hold_$T sb -J yp-h$T -o logs/yp-h$T.%A_%a.log --array=0-$(( ${NCL:-13} - 1 ))); F=$(HOLDDIR=hold_${T}ne MINGC=dyz2:14 sb -J yp-h${T}ne -o logs/yp-h${T}ne.%A_%a.log --array=3)
S=$(sbatch --parsable -J yp-h${T}sum -o logs/yp-h${T}sum.%j.log -c 1 --mem=8G -p bigmem -t 01:00:00 -d afterany:$B:$F --export=ALL --wrap "mkdir -p $SD/hold_$T; YQ12=${YQ12:-rest} python3 scripts/um_v5sum.py holdout hold_$T hold_${T}ne $SD/hold_$T/holdout.tsv")
echo "holdout $T: base classes $B, DYZ2 floor $F, summary $S"
