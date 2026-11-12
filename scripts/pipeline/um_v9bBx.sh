#!/bin/bash
# v9bB with QC exclusions also in marker selection
# usage: bash um_v9bBx.sh (from $P, login node: submits only); as um_v9bBH.sh v9bB (truth Hallast + HPRC R2,) plus MARKX=1:
# the TRUTHX pairs (HG01261 yq12 and dyz2, and the yq12 rest classes) leave marker selection as well as calibration; depth CV included
set -euo pipefail
P=${P:?set P to the project directory}; cd $P; U=markers_v9; X=$P/truth_qc/x; t=v9bBx
for d in model_v5_$t model_v5_${t}_sum $U/cv_v5$t $U/lolo_v5$t $U/holdout_v5$t; do [ -e $d ] && { echo "ERROR $d exists"; exit 1; }; done
[ -s $X/rule_w4.txt ] || { echo "ERROR no $X/rule_w4.txt"; exit 1; }; echo "TRUTHX $(wc -l < $X/rule_w4.txt) lines: $(head -3 $X/rule_w4.txt | tr '\t\n' ' ')"
export D17T=$P/$U/dyz17_truth.tsv D17C=dyz17_span COPIES= TRUTHX=$X/rule_w4.txt MARKX=1
echo "== $t TRUTHX $TRUTHX MARKX $MARKX start $(date)"
TAG=$t UM=$U FULL=0.15 GCFLOORS=dyz2:14 EXCLF=markers_v5/excl_external.txt PACK=1 EVC=5 EVM=48 EVT=01:00:00 DEVT=04:00:00 bash scripts/um_v5run.sh
