#!/bin/bash
# One cv_eval.py task per class (Slurm array task, submitted by v1.sh)
# usage: (model env exported) MODE=cv|production|lolo|holdout|range [IDSF=<file of indices> PACK=<indices per job>] sbatch --array=... cv_task.sh
#   OUTLAB_CLASSES=<class,...>: OUTLAB=1 for these classes, 0 for the others
#   cv: index class*50+repeat*5+fold, BYREP=1: class*10+repeat (CVDIR); production: class (PRE, '{c}' = class); lolo: class*n+lineage of LINLIST (LOLODIR); holdout: class (HOLDLIST, HOLDDIR); range: class (RANGEDIR)
set -euo pipefail
P=${P:?}; D=${PIPE:?}; . $D/config.sh; i=$SLURM_ARRAY_TASK_ID; M=${MODE:-cv}
if [ -n "${IDSF:-}" ]; then mapfile -t A < $IDSF; k=${PACK:-1}; e=0; p=()   # packed job: PACK indices of file IDSF in parallel
  for i in "${A[@]:$(( SLURM_ARRAY_TASK_ID * k )):$k}"; do l=$P/logs/$SLURM_JOB_NAME.${SLURM_ARRAY_JOB_ID}_i$i.log; SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=$(( SLURM_CPUS_PER_TASK / k )) IDSF= bash $D/cv_task.sh > $l 2>&1 & p+=($!); echo "index $i pid $! log $l"; done
  for x in "${p[@]}"; do wait $x || { echo "ERROR pid $x exit $?"; e=1; }; done; echo "end $(date)"; exit $e; fi
C=(dyz3 yq12 dyz1 dyz2 dyz18s yq2k7 dyz18sub dyz19 dyz17 yq12rest yq2k7a dyz18suba dyz17a yq12rest2); export UM=${UM:?} FULL=${FULL:?} RULERS=${RULERS:?} REP=1
case $M in
  cv) if [ -n "${BYREP:-}" ]; then c=${C[$(( i / 10 ))]}; export REPSTART=$(( i % 10 )); else f=$(( i % 50 )); c=${C[$(( i / 50 ))]}; export REPSTART=$(( f / 5 )) FOLDS=$(( f % 5 )); fi;;
  production) c=${C[$i]}; [[ ${PRE:?} == *'{c}'* ]] && export FINAL=${PRE//\{c\}/$c} || export FINAL=${PRE}_$c;;
  lolo) LL=${LINLIST:?}; n=$(wc -l < $LL); c=${C[$(( i / n ))]}; export LIN=$(sed -n "$(( i % n + 1 ))p" $LL);;
  holdout) c=${C[$i]}; mkdir -p $P/$UM/${HOLDDIR:?}; export HOLDTSV=$P/$UM/$HOLDDIR/$c.tsv;;
  range) c=${C[$i]};;
  *) echo "ERROR MODE $M"; exit 1;;
esac
if [ -n "${OUTLAB_CLASSES:-}" ]; then [[ ",$OUTLAB_CLASSES," == *",$c,"* ]] && export OUTLAB=1 || export OUTLAB=0; fi   # OUTLAB per class (v1 CV: the Yq12 classes only)
export CLASSES=ysinglecopy,$c; echo "host $(hostname) mode $M class $c OUTLAB ${OUTLAB:-} repeat ${REPSTART:-} fold ${FOLDS:-all} lineage ${LIN:-} out ${CVDIR:-}${LOLODIR:-}${HOLDDIR:-}${FINAL:-} MINGC ${MINGC:-} EXCL ${EXCL:-} start $(date)"
$PY -I $D/cv_eval.py $M; echo "end $(date)"
