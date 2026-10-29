#!/bin/bash
# CV split by fold and class: 400 tasks
set -euo pipefail
P=/path/to/project; PY=/path/to/conda/ktrk-bench/1.0/bin/python; i=$SLURM_ARRAY_TASK_ID
C=(dyz3 yq12 dyz1 dyz2 dyz18s yq2k7 dyz18sub dyz19); f=$(( i % 50 )); c=${C[$(( i / 50 ))]}
export UM=${UM:-markers_u25} FULL=${FULL:-0.25} RULERS=y90 REP=1 REPSTART=$(( f / 5 )) FOLDS=$(( f % 5 )) CLASSES=ysinglecopy,$c
echo "class $c repeat $REPSTART fold $FOLDS thresholds ${THRS:-} outside ${OUTFS:-} filters ${FILTS:-} dir $CVDIR"; $PY -I $P/scripts/um_eval.py cv
