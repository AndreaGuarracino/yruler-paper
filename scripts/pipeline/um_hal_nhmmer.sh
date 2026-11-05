#!/bin/bash
# nhmmer as Hallast 2026 (-E $E ${ZMB:+-Z $ZMB} --dna), one task = assembly x motif
# usage: [ZMB=<database Mb for E-values>] [SUF=<tag>] sbatch --array=1-N um_hal_nhmmer.sh <tasks.tsv: assembly, fasta, motif>; out $P/hallast_ours<SUF>/<a>/<a>.ours.<motif>.hmmer-tblout.txt.gz
set -euo pipefail
P=${P:?set P to the project directory}   # nhmmer (HMMER 3.4) in PATH
read a fa m < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $1); O=$P/hallast_ours${SUF:-}/$a; mkdir -p $O; f=$O/$a.ours.$m.hmmer-tblout.txt.gz
[ -s $f ] && { echo "$a $m done"; exit 0; }
case $m in DYZ19_Yq) E=1.60E-15;; Yqhet_2k7bp|Yqhet_3k1bp) E=1.60E-200;; DYZ1_Yq|DYZ2_Con|DYZ18_Yq) E=1.60E-150;; *) E=${EVAL:?ERROR no E-value for $m (set EVAL)};; esac
T=$(mktemp -d /dev/shm/hn.XXXX); trap "rm -rf $T" EXIT; zcat -f $fa > $T/t.fa
echo "host $(hostname) $a $m E $E Z ${ZMB:-default} start $(date); target bp $(grep -v '>' $T/t.fa | tr -d '\n' | wc -c), contigs $(grep -c '>' $T/t.fa)"
nhmmer -o /dev/null --tblout $T/out.txt -E $E ${ZMB:+-Z $ZMB} --dna --cpu ${SLURM_CPUS_PER_TASK:-6} $P/hallast_motifs/$m.fasta $T/t.fa
gzip -c $T/out.txt > $f.tmp && mv $f.tmp $f; echo "$a $m hits $(grep -vc '^#' $T/out.txt), end $(date)"
