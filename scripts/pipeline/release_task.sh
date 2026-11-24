#!/bin/bash
#SBATCH -J v1-release -c 4 --mem 12G -t 01:30:00 -o logs/%x.%j.log
# Release of v1: markers and model.tsv for yruler, masked superstrings; Illumina: yruler build, regions BED, checks on HG00126
# usage: OUT=<v1 dir> [TECH=illumina|ont] [YSRC=<yruler source, relative to $P; default yruler>] sbatch v1/pipeline/release_task.sh (from $P)
#   illumina: reads $OUT/model and $OUT/cohorts (cohort_1kgp.tsv, y_dose_scale.txt) -> $OUT/release/{markers.fa,markers.ms.fa.zst,model.tsv,model_pipeline.tsv,regions.bed,pack,HG00126_{dragen_bed_plain,dragen_bed,whole}}, $OUT/tool_build
#   ont: [ONTDK=<window keys> ONTDEPTH=<window runs>] reads $OUT/nanopore/model -> $OUT/nanopore/release/{markers.fa,markers.ms.fa.zst,model.tsv,model_pipeline.tsv,pack,y_dose_scale.txt,y_dose_scale.tsv}
set -euo pipefail
trap 'echo "EXIT line ${LINENO} status $?"; [ -n "${T:-}" ] && rm -rf -- "$T"' EXIT
log() { echo "[$(date '+%F %T %Z')] $*"; }
P=${P:?}; cd $P; export P UM=v1/inputs/markers; D=${PIPE:?}; . $D/config.sh; IN=${IN:?}; OUT=${OUT:?}; TECH=${TECH:-illumina}; C=$CLANG; B=$SPSS_BIN
REF=$IN/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa; YSRC=$P/${YSRC:-yruler}; CO=$OUT/cohorts
case $TECH in illumina) MS=$OUT/model; R=$OUT/release;; ont) MS=$OUT/nanopore/model; R=$OUT/nanopore/release;; *) echo "ERROR TECH $TECH"; exit 1;; esac
for f in $MS/model.npz $MS/model.tsv; do [ -s $f ] || { echo "ERROR missing $f"; exit 1; }; done
[ ! -e $R/markers.fa ] || { echo "ERROR $R/markers.fa exists"; exit 1; }; mkdir -p $R; mkdir -p $TMPBASE; T=$(mktemp -d $TMPBASE/v1rel.XXXX); log "host $(hostname) TECH $TECH OUT $OUT"
M=$T/relmodel; mkdir -p $M; awk -F'\t' 'index(" dyz18sub yq2k7 yq12rest ", " " $1 " ") == 0' $MS/model.tsv > $M/model.tsv; ln -s $MS/model.npz $M/model.npz   # released classes: 7 lengths + parts
if [ $TECH = ont ]; then
  DK=${ONTDK:-$IN/windows/chr8_keys.npy}; python3 $D/release.py ontscale $M $IN/lists/ont_r10_82.tsv ${ONTDEPTH:-$IN/counts/nanopore_window}; YDN=$(cat $M/ydn.txt)   # 82 R10 males of the nanopore CV; ONTDEPTH: chr8 window depth of each male (ontwin72: g72 window)
else
  for f in $CO/cohort_1kgp.tsv $CO/y_dose_scale.txt; do [ -s $f ] || { echo "ERROR missing $f"; exit 1; }; done; YDN=$(cat $CO/y_dose_scale.txt); DK=$IN/windows/chr8_keys.npy
fi
log "classes $(cut -f1 $M/model.tsv | tail -n +2 | paste -sd ' '); YDN $YDN; no-Y cut 0.025; window keys $DK"
DENOMK=$DK YDN=$YDN NOY=0.025 python3 $D/release.py write $M
if [ $TECH = ont ]; then awk 'BEGIN {FS = OFS = "\t"} /^#/ || $1 == "class" {print; next} {$3 = "ont"; print}' $M/release/model.tsv > $M/release/model.tsv.ont; mv $M/release/model.tsv.ont $M/release/model.tsv; cp $M/ydn.txt $R/y_dose_scale.txt; cp $M/ydn.tsv $R/y_dose_scale.tsv; fi
cp $M/release/markers.fa $M/release/model.tsv $R/; cp $M/model.tsv $R/model_pipeline.tsv
log "masked superstrings"; mkdir -p $R/pack; (cd $R/pack && awk '/^>/{c=substr($0,2); next} {print > c".txt"}' $R/markers.fa && : > ms.fa
  for c in $(grep '^>' $R/markers.fa | cut -c2-); do awk '{print ">"NR; print}' $c.txt > $c.kfa
    $B/kmercamel -k 31 -p $c.kfa -o $c.ms.txt 2> $c.camel.log || $B/kmercamel compute -k 31 -o $c.ms.txt $c.kfa 2>> $c.camel.log
    s=$(grep -v '>' $c.ms.txt | tr -d '\n'); up=$(echo -n "$s" | tr -cd 'ACGT' | wc -c); log "$c kmers $(wc -l < $c.txt) superstring chars ${#s} uppercase $up"
    [ "$up" = "$(wc -l < $c.txt)" ] || { echo "ERROR $c uppercase starts differ from k-mers"; exit 1; }; printf '>%s\n%s\n' $c "$s" >> ms.fa; done
  $B/zstd -19 -q -f ms.fa -o $R/markers.ms.fa.zst); log "markers.ms.fa.zst md5 $(md5sum < $R/markers.ms.fa.zst | cut -c1-32), model.tsv md5 $(md5sum < $R/model.tsv | cut -c1-32)"
if [ $TECH = ont ]; then log end; exit 0; fi
log "build yruler from $YSRC"; (cd $YSRC && PATH=$BUILDDEPS:$PATH CARGO_TARGET_DIR=$OUT/tool_build LIBCLANG_PATH=$C/lib BINDGEN_EXTRA_CLANG_ARGS="-I$C/lib/clang/18/include" $CARGO build --release 2> $T/build.log) || { tail -20 $T/build.log; exit 1; }
Y=$OUT/tool_build/release/yruler; log "built $Y $($Y --version)"
( cat $IN/reference/regions_grch38_1kgp_dragen376.bed; printf 'chr8\t90990000\t101010000\n' ) | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} $1 == c && $2 <= e {if ($3 > e) e = $3; next} {if (c != "") print c, s, e; c = $1; s = $2; e = $3} END {print c, s, e}' > $R/regions.bed
log "regions (DRAGEN 3.7.6 BED plus the chr8:91-101 Mb window +-10 kb) $(wc -l < $R/regions.bed), bp $(awk '{s += $3 - $2} END {print s}' $R/regions.bed)"
chk() { python3 - $1/estimates.tsv $CO/cohort_1kgp.tsv <<'PY'
import sys
L = {x.split('\t')[1]: float(x.split('\t')[2]) for x in list(open(sys.argv[1]))[1:]}; R = [l.rstrip('\n').split('\t') for l in open(sys.argv[2])]; H = R[0]; r = next(x for x in R if x[0] == 'HG00126')
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
if set(L) != set(NM): sys.exit(f'ERROR reported classes {sorted(L)}')
for c, n in NM.items(): t = float(r[H.index(f'{n}_bp')]); print(f'{c}: yruler {L[c]:.0f} cohort table {t:.0f} diff {100 * (L[c] / t - 1):+.3f}%', flush=True)
PY
}
run() { local n=$1 f=$2; shift 2; (cd $T && $Y -i $f -r $REF "$@" -M $R/model.tsv -n HG00126 -t 4 -o $R/HG00126_$n > /dev/null 2> $R/HG00126_$n.log) || { echo "ERROR yruler $n"; tail -3 $R/HG00126_$n.log; exit 1; }
  COHORT=$CO/cohort_1kgp.tsv python3 $D/release.py check HG00126 $R/HG00126_$n $M || { echo "ERROR release check $n"; exit 1; }; chk $R/HG00126_$n; }
U=https://s3.us-west-2.amazonaws.com/1000genomes-dragen-3.7.6/data/individuals/hg38-graph-based/HG00126/HG00126.bam   # DRAGEN 3.7.6 alignment (the BED is for this aligner only), as the README demo
run dragen_bed_plain "$U" -b $R/regions.bed -m $R/markers.fa; run dragen_bed "$U" -b $R/regions.bed -m $R/markers.ms.fa.zst
cmp <(cut -f1-3 $R/HG00126_dragen_bed_plain/estimates.tsv) <(cut -f1-3 $R/HG00126_dragen_bed/estimates.tsv) && log "packed and plain markers: identical estimates"
run whole https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3240199/HG00126.final.cram -m $R/markers.ms.fa.zst   # NYGC CRAM, whole file
log end
