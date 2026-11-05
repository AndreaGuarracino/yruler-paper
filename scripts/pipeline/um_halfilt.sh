#!/bin/bash
# Hallast HMMER filters on one assembly
# usage: sbatch --array=1-N um_halfilt.sh <list: assembly, table dir, kind yq12|yq12all|dyz19|dyz17|dyz17blk>; out $P/hallast_filt/<assembly>.<kind>.csv (DYZ17: DYZ19 array rule on the DYZ17 table alone)
set -euo pipefail
P=${P:?set P to the project directory}; PY=${PY:?set PY to a Python with pandas}; cd $P
read a d k < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $1); o=hallast_filt/$a.$k.csv; echo "host $(hostname) $a $k $d start $(date)"
case $k in
  yq12) $PY -I scripts/filterYq12_hal.py -i $d -o $o;;
  dyz19) $PY -I scripts/filterDYZ19_hal.py -i $d -o $o;;
  yq12all) ALL2K7=1 $PY -I scripts/filterYq12_hal.py -i $d -o $o;;   # every YqHet2.7 block (amount)
  dyz17blk) m=${D17M:-DYZ17_unit}; T=$(mktemp -d ${SCRATCH:?set SCRATCH}/hf.XXXX); trap "rm -rf $T" EXIT; ln -s $P/$d/$a.ours.$m.hmmer-tblout.txt.gz $T/
         EXTRA=$m NOMERGE=1 $PY -I scripts/filterYq12_hal.py -i $T -o $o;;   # Hallast blocks without the final merge (amount)
  dyz17) m=${D17M:-DYZ17_unit}; T=$(mktemp -d ${SCRATCH:?set SCRATCH}/hf.XXXX); trap "rm -rf $T" EXIT; ln -s $P/$d/$a.ours.$m.hmmer-tblout.txt.gz $T/   # DYZ17 table alone
         MOTIF=$m $PY -I scripts/filterDYZ19_hal.py -i $T -o $o;;   # DYZ19 array rule
  *) echo "ERROR kind $k"; exit 1;;
esac
echo "$a $k blocks $(($(wc -l < $o) - 1)), end $(date)"
