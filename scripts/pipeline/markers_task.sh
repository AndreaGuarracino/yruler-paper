#!/bin/bash
#SBATCH -J v1-asmbits -c 4 --mem 32G -t 03:00:00 -o logs/%x.%A_%a.log
# Per assembly (one array task each): KMC 31-mer sets of the whole Y, of each class region and of the Y outside each class, then their presence over the universe keys
# usage: MK=<dir with assemblies.tsv, regions.tsv, u_keys.npy, u_db.kmc_pre/suf> OUTD=<out dir> sbatch -a <line numbers of assemblies.tsv> v1/pipeline/markers_task.sh (from $P)
#   -> $OUTD/db/<assembly>.{whole,<class>,out_<class>}.kmc_*, $OUTD/bits/<assembly>.npz; outside sets for the classes of the assembly plus dyz18s, dyz18sub, dyz18suba (KMC 3.2.4, samtools 1.24)
set -euo pipefail
trap 'echo "EXIT line ${LINENO} status $?"; [ -n "${T:-}" ] && rm -rf -- "$T"' EXIT
log() { echo "[$(date '+%F %T %Z')] $*"; }
P=${P:?}; cd $P; D=${PIPE:?}; . $D/config.sh; K=$KMC_BIN
MK=${MK:?}; OUTD=${OUTD:?}; DB=$OUTD/db; mkdir -p $DB $OUTD/bits $TMPBASE; T=$(mktemp -d $TMPBASE/asm.XXXX); export TMPD=$T
read a f < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $MK/assemblies.tsv); [ -n "$a" ] && [ -s $f ] && [ -s $f.fai ] || { echo "ERROR line $SLURM_ARRAY_TASK_ID: assembly '$a' fasta '$f' or its .fai missing"; exit 1; }
log "host $(hostname) assembly $a $f, markers $MK, out $OUTD"
db() { $K/kmc -k31 -ci1 -cs255 -m4 -t4 -fm $1 $T/raw $T > /dev/null; $K/kmc_tools -t4 transform $T/raw set_counts 1 $2 > /dev/null; }
[ -e $DB/$a.whole.kmc_pre ] || db $f $DB/$a.whole; log "whole Y db"
CL=$(awk -F'\t' -v a=$a '$1 == a {print $2}' $MK/regions.tsv | sort -u); log "classes with regions: $(echo $CL | tr '\n' ' ')"
for c in $CL; do awk -F'\t' -v a=$a -v c=$c '$1 == a && $2 == c {print $3":"$4+1"-"$5}' $MK/regions.tsv > $T/rg; $SAM faidx -r $T/rg $f > $T/$c.fa; db $T/$c.fa $DB/$a.$c
  log "$c regions $(wc -l < $T/rg) bp $(grep -v '>' $T/$c.fa | tr -d '\n' | wc -c)"; done
for c in $CL dyz18s dyz18sub dyz18suba; do python3 $D/markers.py outside $MK/regions.tsv $a $c $f.fai > $T/rg; [ -s $T/rg ] || continue   # Y outside the class
  $SAM faidx -r $T/rg $f > $T/o.fa; db $T/o.fa $DB/$a.out_$c; log "$c outside: segments $(wc -l < $T/rg) bp $(grep -v '>' $T/o.fa | tr -d '\n' | wc -c)"; done
python3 $D/markers.py bits $MK/u_keys.npy $MK/u_db $DB $a $OUTD/bits/$a.npz; log end
