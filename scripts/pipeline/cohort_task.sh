#!/bin/bash
#SBATCH -J v1-cohorts -c 16 --mem 100G -t 01:30:00 -o logs/%x.%j.log
# Cohort tables of v1: per-sample doses and lengths from the region slices (cohort.py), rows of the CRAMs counted whole (cohort_add.py apply), pangenome X/PAR doses (xpar), merge
# usage: OUT=<v1 dir> sbatch v1/pipeline/cohort_task.sh (from $P) -> $OUT/cohorts/{cohort_1kgp,cohort_hgdp_sgdp,cohort_fem}.tsv, y_dose_scale.txt, pairs/; stages in $OUT/cohorts/{slices,whole}
#   reads $OUT/model (model.npz, model.tsv); YPQ (ruler Yp/Yq blocks), YDN (raw y_dose of the reference XY males, g72 window) and Y state files as below
set -euo pipefail
trap 'echo "EXIT line ${LINENO} status $?"' EXIT
log() { echo "[$(date '+%F %T %Z')] $*"; }
P=${P:?}; cd $P; D=${PIPE:?}; . $D/config.sh; IN=${IN:?}; OUT=${OUT:?}; C=$OUT/cohorts; S=$C/slices; W=$C/whole; mkdir -p $S $W
export P UM=v1/inputs/markers DFA=extra_keys.fa DENOMC=$IN/counts/window DENOMFA=$IN/windows/chr8_window.fa DENOMSEL=first DENOMN=chr8_91_101 \
  DENOMK=$IN/windows/chr8_keys.npy YPQ=$IN/markers/ruler_ypq.npy YDN=0.961700 NOY=0.025 DYZ17DIR=../counts/extra MD=$OUT/model D17T=$IN/truth/dyz17_truth.tsv D17C=dyz17_span TRUTHX=$IN/truth/exclusions.txt
YS=$IN/lists; ADD=$IN/lists/cohort_add.tsv   # ADD: 58 HGDP/SGDP CRAMs added after the slice counts (lists/rc_add_20261109.tsv without SGDP_HG00126, whose CRAM md5 never matched)
for f in $MD/model.npz $MD/model.tsv $YPQ $DENOMK $DENOMFA.order $YS/y_states_1kgp.tsv $YS/y_states_hgdp_sgdp.tsv $ADD $IN/windows/xpar_means.tsv; do [ -s $f ] || { echo "ERROR missing $f"; exit 1; }; done
log "host $(hostname) OUT $OUT YPQ $YPQ YDN $YDN NOY $NOY window keys $DENOMK, slices counted $(ls $DENOMC | wc -l)"
echo $YDN > $C/y_dose_scale.txt
if [ -s $S/cohort_1kgp.tsv ] && [ -s $S/cohort_fem.tsv ] && [ -s $S/cohort_hgdp_sgdp.tsv ]; then log "slices: 3 tables exist, skipped"; else log "slices: 3 tables"; e=0; p=()
YSF=$YS/y_states_1kgp.tsv python3 $D/cohort.py v1/inputs/counts/slices_1kgp $S/cohort_1kgp.tsv $MD > $S/1kgp.log 2>&1 & p+=($!)
YSF=$YS/y_states_1kgp.tsv python3 $D/cohort.py v1/inputs/counts/slices_females $S/cohort_fem.tsv $MD > $S/fem.log 2>&1 & p+=($!)
YSF=$YS/y_states_hgdp_sgdp.tsv python3 $D/cohort.py v1/inputs/counts/slices_hgdp_sgdp $S/cohort_hgdp_sgdp.tsv $MD > $S/hgsg.log 2>&1 & p+=($!)
for x in "${p[@]}"; do wait $x || { echo "ERROR pid $x"; e=1; }; done; tail -n 2 $S/*.log; [ $e = 0 ] || exit 1; fi
log "whole-file recounts: $(wc -l < $ADD) CRAMs"; python3 $D/cohort_add.py apply $ADD $W/rows.tsv > $W/apply.log 2>&1 || { tail -5 $W/apply.log; exit 1; }; tail -2 $W/apply.log
n=$(( $(wc -l < $W/rows.tsv) - 1 )); [ $n = $(wc -l < $ADD) ] || { echo "ERROR $n rows for $(wc -l < $ADD) CRAMs"; exit 1; }
log "X and PAR doses"; python3 $D/cohort_add.py xpar $IN/windows/xpar_means.tsv $S/cohort_1kgp.tsv $S/cohort_1kgp.tsv $C/cohort_1kgp.tsv $S/cohort_fem.tsv $C/cohort_fem.tsv $S/cohort_hgdp_sgdp.tsv $W/slices_xpar.tsv
python3 $D/cohort_add.py merge $W/slices_xpar.tsv $W/rows.tsv $C/cohort_hgdp_sgdp.tsv
log "pairs"; python3 $D/pairs.py $C $C/pairs > $C/pairs.log 2>&1 || { tail -5 $C/pairs.log; exit 1; }; wc -l $C/pairs/*.tsv
for t in cohort_1kgp cohort_hgdp_sgdp cohort_fem; do log "$t rows $(( $(wc -l < $C/$t.tsv) - 1 )) md5 $(md5sum < $C/$t.tsv | cut -c1-32)"; done; log end
