#!/bin/bash
# Full v5 run: model, CVs, applies, summaries
# walls measured on v9b: CV, lineage-out, holdout max 13.3 min (EVT 20 min); depth p95 54 min, max 164 min (DEVT 65 min, then rerun ONLY=depth with DEVT=04:00:00 for the timed-out tasks)
# usage: [TAG=..] [EXCLF=markers_v5/excl_external.txt] [ONLY=depth | ONLY=floors | ONLY=cv | ONLY=production | DEPTH=0] [GCFLOORS=dyz2:16,yq2k7:10] [DEP=afterok:<job>] bash um_v5run.sh (from $P); defaults: inner 5-fold CV, per-fold all-copies truth, candidates absent from the non-Y sequence of CHM13 and 71 HPRC r2 haplotypes (NYSTRICT=1), outside fraction for Y only; ONLY=depth: depth tasks without output or queue entry; truth rule B; outputs carry TAG: model_v5[_TAG], markers_v5/{cv,lolo,holdout}_v5TAG*, model_v5[_TAG]_sum
set -euo pipefail
P=${P:?set P to the project directory}; cd $P; t=${TAG:-}; M=model_v5${t:+_$t}; SM=${M}_sum; G=${GCFLOORS:-}; U=${UM:-markers_v5}; F=${FEMNPZ:-$P/$U/fem_v5.npz}; mkdir -p $SM   # UM: marker set dir (markers_v6: universe at 15%, with FULL=0.15)
CI=(dyz3 yq12 dyz1 dyz2 dyz18s yq2k7 dyz18sub dyz19 dyz17 yq12rest yq2k7a dyz18suba dyz17a yq12rest2); NC=${#CI[@]}; NL=$(wc -l < $P/lists/v5_lineages.txt)
fi_() { local c i o=""; for c in $(echo $G | tr ',' '\n' | cut -d: -f1); do for i in "${!CI[@]}"; do [ ${CI[$i]} = $c ] && o="$o,$(( i * $1 ))-$(( i * $1 + $1 - 1 ))"; done; done; echo ${o#,}; }   # task ranges of the floored classes (block size $1)
ck() { local x; x=$(sbatch --parsable "$@" 2>&1); echo "$x" | grep "^[0-9]" || { echo "SUBMIT ERROR: $(echo "$x" | grep -i error | head -1)" >&2; exit 1; }; }
ev() { local n=$1 k=${PACK:-4} f=$P/logs/ids/yp-v5$t$1.$(date +%s%N).txt; mkdir -p $P/logs/ids; echo $2 | tr ',' '\n' | awk -F- '{for (j = $1; j <= ($2 == "" ? $1 : $2); j++) print j}' > $f; shift 2; IDSF=$f PACK=$k ck -J yp-v5$t$n -o logs/yp-v5$t$n.%A_%a.log --array=0-$(( ($(wc -l < $f) + k - 1) / k - 1 )) -c $(( ${EVC:-6} * k )) --mem=$(( ${EVM:-40} * k ))G -t ${EVT:-00:20:00} --comment=$f ${DEP:+-d $DEP} "$@" --export=ALL scripts/um_cvsplit.sh; }   # um_eval.py tasks, PACK per job in parallel (the scheduler starts few jobs per minute); comment: file of task indices; EVT 09:00:00: cmedium QOS; DEP: dependency
C=("${CI[@]}")   # index order of um_cvsplit.sh
mi() { local i a="" r=" $(squeue -u $USER -h -r -n $1 -o '%K %k' | awk '{if ($2 ~ /^\//) {while ((getline x < $2) > 0) print x; close($2)} else print $1}' | tr '\n' ' ') "; for i in $(echo $2 | tr ',' '\n' | awk -F- '{for (j = $1; j <= ($2 == "" ? $1 : $2); j++) print j}'); do f=$P/$U/$3/rep$(( i % 50 / 5 ))_k$(( i % 5 ))_c${C[$(( i / 50 ))]}.tsv; [ ${4:-} = rep ] && f=$P/$U/$3/rep$(( i % 10 ))_c${C[$(( i / 10 ))]}.tsv; [ -s $f ] || [[ $r == *" $i "* ]] || a="$a,$i"; done; echo ${a#,} | tr ',' '\n' | awk 'NR == 1 {s = p = $1; next} $1 == p + 1 {p = $1; next} {printf "%s%s,", s, (p > s ? "-" p : ""); s = p = $1} END {printf "%s%s\n", s, (p > s ? "-" p : "")}'; }   # missing tasks as ranges: per fold (class*50 + repeat*5 + fold) or, with 4th argument rep, per repeat (class*10 + repeat)
sm() { ck -J yp-v5${t}sum -o logs/yp-v5${t}sum.%j.log -c 1 --mem=8G -p bigmem -t 02:00:00 "$@" --export=ALL --wrap "$W"; }   # small job on bigmem
export UM=$U FULL=${FULL:-0.25} RULERS=y90 THRS=${THRS:-0.5,0.6,0.7,0.8,0.9,0.95} OUTFS=${OUTFS:-0,0.05,0.1,0.25} TRUTHX=${TRUTHX:-$P/truth_qc/x/rule_w4.txt} TRUTHSET=${TRUTHSET:-} COPIES=${COPIES-$P/$U/copies_$(basename $U | sed 's/markers_//')} INNER=${INNER:-5} NYSCORE=${NYSCORE:-$P/$U/pg_topFM.noours.count.npy} NYN=${NYN:-71} NYSTRICT=${NYSTRICT:-1} NYRULER=${NYRULER:-$P/$U/nyruler_none.npy} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FILT=${FILT:-} FEMT=${FEMT:-inf}
if [ -n "${EXCLF:-}" ]; then export EXCL=$P/$EXCLF; else unset EXCL; fi; unset BYREP MINGC   # CV tasks per fold; floors only through GCFLOORS
if [ "${ONLY:-}" = floors ]; then   # floored classes only, after a base run with the same TAG: production, CV, lineage-out, holdout, model with floors, summary
  [ -n "$G" ] || { echo "ERROR ONLY=floors needs GCFLOORS" >&2; exit 1; }
  PF=$(MODE=production PRE=fm_v5${t}f MINGC=$G ev pf $(fi_ 1)); Q=$(squeue -u $USER -h -o :%F -n yp-v5${t}pb,yp-v5${t}cv,yp-v5${t}lolo,yp-v5${t}hold | sort -u | tr -d '\n')   # base jobs still queued
  J=":$(MODE=cv BYREP=1 FEMNPZ=$F CVDIR=cv_v5${t}ne MINGC=$G ev cvne $(fi_ 10)):$(MODE=lolo LOLODIR=lolo_v5${t}ne MINGC=$G ev lolone $(fi_ $NL))"
  W="python3 scripts/um_v5sum.py cv cv_v5$t cv_v5${t}ne $SM/cv.tsv && python3 scripts/um_v5sum.py lolo lolo_v5$t lolo_v5${t}ne $SM/lolo.tsv"
  if [ -n "${EXCL:-}" ]; then J="$J:$(MODE=holdout HOLDDIR=holdout_v5${t}ne HOLDLIST=$EXCL MINGC=$G ev holdne $(fi_ 1))"; W="$W && python3 scripts/um_v5sum.py holdout holdout_v5$t holdout_v5${t}ne $SM/holdout.tsv"; fi
  AP=$(ck -J yp-v5${t}app -o logs/yp-v5${t}app.%j.log -c 16 --mem=64G -p bigmem -t 04:00:00 -d afterany$Q:$PF --export=ALL --wrap "python3 scripts/um_model.py fm_v5${t}b fm_v5${t}f $G $P/$M && python3 scripts/um_apply.py pass_1kgp $M/cohort_1kgp.tsv $P/$M && python3 scripts/um_apply.py pass_hgsg $M/cohort_hgdp_sgdp.tsv $P/$M && DYZ17DIR=dyz17lots python3 scripts/um_apply.py lots_view $M/other_runs.tsv $P/$M && python3 scripts/um_xtech.py $P/$M")
  S1=$(sm -d afterany$Q$J); echo "tag '$t' floors $G: production $PF, model with floors + applies $AP, summary $S1"; exit 0
fi
if [ "${ONLY:-}" = production ]; then   # model only, no cohort applies (nanopore: CNPY=cnpy_ont)
  PB=$(MODE=production PRE=fm_v5${t}b ev pb 0-$(( NC - 1 ))); PF=$PB; [ -n "$G" ] && PF=$(MODE=production PRE=fm_v5${t}f MINGC=$G ev pf $(fi_ 1))
  AP=$(ck -J yp-v5${t}app -o logs/yp-v5${t}app.%j.log -c 16 --mem=64G -p bigmem -t 04:00:00 -d afterok:$PB:$PF --export=ALL --wrap "python3 scripts/um_model.py fm_v5${t}b fm_v5${t}$([ -n "$G" ] && echo f || echo b) ${G:-none} $P/$M")
  echo "tag '$t' production only: $PB $PF, model $AP"; exit 0
fi
if [ "${ONLY:-}" = cv ]; then   # CV only (e.g. PED=<other groups>), summary to $SM/cv.tsv
  J=":$(MODE=cv FEMNPZ=$F CVDIR=cv_v5$t ev cv ${CVIDX:-0-$(( NC * 50 - 1 ))})"; [ -n "$G" ] && J="$J:$(MODE=cv FEMNPZ=$F CVDIR=cv_v5${t}ne MINGC=$G ev cvne $(fi_ 50))"
  W="python3 scripts/um_v5sum.py cv cv_v5$t cv_v5${t}ne $SM/cv.tsv"; S1=$(sm -d afterok$J); echo "tag '$t' CV only, PED ${PED:-default}: summary $S1"; exit 0
fi
if [ "${ONLY:-}" != depth ]; then
PB=$(MODE=production PRE=fm_v5${t}b ev pb 0-$(( NC - 1 ))); PF=$PB; [ -n "$G" ] && PF=$(MODE=production PRE=fm_v5${t}f MINGC=$G ev pf $(fi_ 1))
AP=$(ck -J yp-v5${t}app -o logs/yp-v5${t}app.%j.log -c 16 --mem=64G -p bigmem -t 04:00:00 -d afterok:$PB:$PF --export=ALL --wrap "python3 scripts/um_model.py fm_v5${t}b fm_v5${t}$([ -n "$G" ] && echo f || echo b) ${G:-none} $P/$M && python3 scripts/um_apply.py pass_1kgp $M/cohort_1kgp.tsv $P/$M && python3 scripts/um_apply.py pass_hgsg $M/cohort_hgdp_sgdp.tsv $P/$M && DYZ17DIR=dyz17lots python3 scripts/um_apply.py lots_view $M/other_runs.tsv $P/$M && python3 scripts/um_xtech.py $P/$M")
J=":$(MODE=cv FEMNPZ=$F CVDIR=cv_v5$t ev cv 0-$(( NC * 50 - 1 ))):$(MODE=lolo LOLODIR=lolo_v5$t ev lolo 0-$(( NC * NL - 1 )))"   # one task per class, repeat and fold
[ -n "$G" ] && J="$J:$(MODE=cv FEMNPZ=$F CVDIR=cv_v5${t}ne MINGC=$G ev cvne $(fi_ 50)):$(MODE=lolo LOLODIR=lolo_v5${t}ne MINGC=$G ev lolone $(fi_ $NL))"
W="python3 scripts/um_v5sum.py cv cv_v5$t cv_v5${t}ne $SM/cv.tsv && python3 scripts/um_v5sum.py lolo lolo_v5$t lolo_v5${t}ne $SM/lolo.tsv"
if [ -n "${EXCL:-}" ]; then J="$J:$(MODE=holdout HOLDDIR=holdout_v5$t HOLDLIST=$EXCL ev hold 0-$(( NC - 1 )))"; [ -n "$G" ] && J="$J:$(MODE=holdout HOLDDIR=holdout_v5${t}ne HOLDLIST=$EXCL MINGC=$G ev holdne $(fi_ 1))"; W="$W && python3 scripts/um_v5sum.py holdout holdout_v5$t holdout_v5${t}ne $SM/holdout.tsv"; fi
S1=$(sm -d afterok$J); fi; [ "${DEPTH:-1}" = 0 ] && { echo "tag '$t' $M: production $PB $PF, model + applies + pairs $AP, CV/lolo/holdout summary $S1, no depth"; exit 0; }; L="full d30x d25x d20x d15x d10x d5x d2x d1x"; D=""
for l in $L; do   # afterany: um_v5sum.py logs missing repeats
  for k in "d $l cv_v5${t}_depth_$l 0-$(( NC * 10 - 1 ))" "dn $l cv_v5${t}ne_depth_$l $(fi_ 10)"; do set -- $k; [ $1 = dn ] && [ -z "$G" ] && continue; a=$(mi yp-v5$t$1$2 $4 $3 rep); D="$D$(squeue -u $USER -h -n yp-v5$t$1$2 -o :%F | sort -u | tr -d '\n')"   # one task per class and repeat: 5x fewer starts
    [ -z "$a" ] && continue; echo "submit yp-v5$t$1$2 tasks $a" >&2
    if [ $1 = d ]; then D="$D:$(EVT=${DEVT:-01:05:00} MODE=cv BYREP=1 FEMNPZ=$F ALTDIR=$U/depth/$2 CVDIR=$3 ev $1$2 $a)"; else D="$D:$(EVT=${DEVT:-01:05:00} MODE=cv BYREP=1 FEMNPZ=$F ALTDIR=$U/depth/$2 CVDIR=$3 MINGC=$G ev $1$2 $a)"; fi; done; done
W=""; for l in $L; do W="$W python3 scripts/um_v5sum.py cv cv_v5${t}_depth_$l cv_v5${t}ne_depth_$l $SM/cv_depth_$l.tsv &&"; done; W="$W MDIR=$SM python3 scripts/um_v5sum.py depth $SM/depth.tsv $L"; S2=$(sm -d afterany$D)
echo "tag '$t' $M: production ${PB:-} ${PF:-}, model + applies + pairs ${AP:-}, CV/lolo/holdout summary ${S1:-}, depth summary $S2"
