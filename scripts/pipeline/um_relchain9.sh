#!/bin/bash
# release, 7 classes: build, markers, checks
# usage: [T=v9bB] [TRY=yruler-try21] sbatch ... um_relchain9.sh (from $P); model model_v5_$T, final tables model_v5_$T/final_chr8_91_101, out model_v5_${T}_rel7
set -euo pipefail; P=${P:?set P to the project directory}; C=${CLANG:?set CLANG to the clang 18 prefix}; Z=zstd
T=${T:-v9b}; M=$P/model_v5_${T}_rel7; R=$M/release; T9=$P/model_v5_$T/final_chr8_91_101; REF=$P/ext/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; echo "host $(hostname) start $(date)"
mkdir -p $M; D="dyz18sub yq2k7"; [ $T = v9b ] || D="$D yq12rest"; awk -F'\t' -v d=" $D " 'index(d, " " $1 " ") == 0' $P/model_v5_$T/model.tsv > $M/model.tsv; ln -sfn $P/model_v5_$T/model.npz $M/model.npz; echo "classes: $(cut -f1 $M/model.tsv | tail -n +2 | paste -sd ' ')"
cd ${YRULER_SRC:?set YRULER_SRC to the yruler repository}; LIBCLANG_PATH=$C/lib BINDGEN_EXTRA_CLANG_ARGS="-I$C/lib/clang/18/include" cargo build --release 2> build.log || { tail -20 build.log; exit 1; }
Y=$YRULER_SRC/target/release/yruler; cd $P; echo "built $Y"
YDN=$(cat $T9/ydn.txt); NOY=$(cat $T9/noy.txt); echo "YDN $YDN, no-Y cut $NOY"; UM=markers_v9 DENOMK=$P/ykmers/wpk_chr8_91_101/win_keys.npy YDN=$YDN NOY=$NOY python3 scripts/um_release.py $M
R=$R bash scripts/release_pack.sh
( cat $YRULER_SRC/data/regions_grch38_1kgp_dragen376.bed; grep -P '^chr8\t' $P/markers_v6/win5.bed ) | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} $1 == c && $2 <= e {if ($3 > e) e = $3; next} {if (c != "") print c, s, e; c = $1; s = $2; e = $3} END {print c, s, e}' > $R/regions.bed; echo "regions $(wc -l < $R/regions.bed) bp $(awk '{s += $3 - $2} END {print s}' $R/regions.bed)"
chk() { python3 - $1/lengths.tsv $T9/cohort_1kgp.tsv <<'PY'
import sys
import os; f = sys.argv[1] if os.path.exists(sys.argv[1]) else sys.argv[1].replace('lengths.tsv', 'estimates.tsv'); L = {x.split('\t')[1]: float(x.split('\t')[2]) for x in list(open(f))[1:]}; R = [l.rstrip('\n').split('\t') for l in open(sys.argv[2])]; H = R[0]; r = next(x for x in R if x[0] == 'HG00126')
NM = {'dyz3': 'DYZ3', 'yq12': 'Yq12', 'dyz1': 'DYZ1', 'dyz2': 'DYZ2', 'dyz18s': 'DYZ18S', 'dyz19': 'DYZ19', 'dyz17': 'DYZ17'}
if set(L) != set(NM): sys.exit(f'ERROR reported classes {sorted(L)}')
for c, n in NM.items(): t = float(r[H.index(f'{n}_bp')]); print(f'{c}: yruler {L[c]:.0f} cohort table {t:.0f} diff {100 * (L[c] / t - 1):+.3f}%', flush=True)
PY
}
U=$(grep -P "^HG00126\t" lists/v6_slice_s3.tsv | cut -f2); W=$(mktemp -d ${SCRATCH:?set SCRATCH}/rel.XXXX)
for m in markers.fa markers.ms.fa.zst; do (cd $W && $Y -i $U -r $REF -b $R/regions.bed -m $R/$m -M $R/model.tsv -n HG00126 -t 16 -o $R/HG00126_$m > /dev/null 2> $R/HG00126_$m.log) || { echo "ERROR yruler $m"; tail -3 $R/HG00126_$m.log; exit 1; }; tail -2 $R/HG00126_$m.log; done
E=lengths.tsv; [ -s $R/HG00126_markers.fa/$E ] || E=estimates.tsv; cmp $R/HG00126_markers.fa/$E $R/HG00126_markers.ms.fa.zst/$E && echo "packed and plain markers: identical lengths"
echo "== slice"; UM=markers_v9 COHORT=$T9/cohort_1kgp.tsv python3 scripts/um_relcheck.py HG00126 $R/HG00126_markers.fa $M || { echo "ERROR relcheck slice"; exit 1; }; cat $R/HG00126_markers.fa/$E; chk $R/HG00126_markers.fa
CR=$(ls ${SCRATCH:?set SCRATCH}/yp_keep/truth/HG00126*.cram | head -1); B=${SCRATCH:?set SCRATCH}/yp_keep/dragen376/HG00126.HG00126.bam; echo "CRAM $CR; DRAGEN BAM $B"
for run in "whole $CR" "dragen_bed $B -b $R/regions.bed"; do set -- $run; n=$1; shift; f=$1; shift
  (cd $W && $Y -i $f -r $REF "$@" -m $R/markers.ms.fa.zst -M $R/model.tsv -n HG00126 -t 16 -o $R/HG00126_$n > /dev/null 2> $R/HG00126_$n.log) || { echo "ERROR yruler $n"; tail -3 $R/HG00126_$n.log; exit 1; }; tail -2 $R/HG00126_$n.log
  echo "== $n"; UM=markers_v9 COHORT=$T9/cohort_1kgp.tsv python3 scripts/um_relcheck.py HG00126 $R/HG00126_$n $M || { echo "ERROR relcheck $n"; exit 1; }; chk $R/HG00126_$n
done; rmdir $W 2>/dev/null || true; echo "end $(date)"
