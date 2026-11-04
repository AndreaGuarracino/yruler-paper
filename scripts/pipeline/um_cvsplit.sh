#!/bin/bash
# One um_eval.py task per class
# usage: (model env exported) MODE=cv|production|lolo|holdout [IDSF=<file of indices> PACK=<indices per job>] sbatch --array=... um_cvsplit.sh; cv: index class*50+repeat*5+fold, BYREP=1: class*10+repeat (CVDIR); production: class (PRE); lolo: class*n+lineage of lists/v5_lineages.txt (LOLODIR); holdout: class (HOLDLIST, HOLDDIR)
set -euo pipefail
P=${P:?set P to the project directory}; PY=python3; i=$SLURM_ARRAY_TASK_ID; M=${MODE:-cv}
if [ -n "${IDSF:-}" ]; then mapfile -t A < $IDSF; k=${PACK:-1}; e=0; p=()   # packed job: PACK indices of file IDSF in parallel
  for i in "${A[@]:$(( SLURM_ARRAY_TASK_ID * k )):$k}"; do l=$P/logs/$SLURM_JOB_NAME.${SLURM_ARRAY_JOB_ID}_i$i.log; SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=$(( SLURM_CPUS_PER_TASK / k )) IDSF= bash $P/scripts/um_cvsplit.sh > $l 2>&1 & p+=($!); echo "index $i pid $! log $l"; done
  for x in "${p[@]}"; do wait $x || { echo "ERROR pid $x exit $?"; e=1; }; done; echo "end $(date)"; exit $e; fi
C=(dyz3 yq12 dyz1 dyz2 dyz18s yq2k7 dyz18sub dyz19 dyz17 yq12rest yq2k7a dyz18suba dyz17a yq12rest2); export UM=${UM:-markers_u25} FULL=${FULL:-0.25} RULERS=y90 REP=1
case $M in
  cv) if [ -n "${BYREP:-}" ]; then c=${C[$(( i / 10 ))]}; export REPSTART=$(( i % 10 )); else f=$(( i % 50 )); c=${C[$(( i / 50 ))]}; export REPSTART=$(( f / 5 )) FOLDS=$(( f % 5 )); fi;;
  production) c=${C[$i]}; export FINAL=${PRE:?}_$c;;
  lolo) n=$(wc -l < $P/lists/v5_lineages.txt); c=${C[$(( i / n ))]}; export LIN=$(sed -n "$(( i % n + 1 ))p" $P/lists/v5_lineages.txt);;
  holdout) c=${C[$i]}; mkdir -p $P/$UM/${HOLDDIR:?}; export HOLDTSV=$P/$UM/$HOLDDIR/$c.tsv;;
  *) echo "ERROR MODE $M"; exit 1;;
esac
export CLASSES=ysinglecopy,$c; echo "host $(hostname) mode $M class $c repeat ${REPSTART:-} fold ${FOLDS:-all} lineage ${LIN:-} out ${CVDIR:-}${LOLODIR:-}${HOLDDIR:-}${FINAL:-} MINGC ${MINGC:-} FILT ${FILT:-} FEMT ${FEMT:-} EXCL ${EXCL:-} start $(date)"
$PY -I $P/scripts/um_eval.py $M; echo "end $(date)"
