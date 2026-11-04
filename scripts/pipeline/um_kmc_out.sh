#!/bin/bash
# per assembly and class: k-mers of the Y outside the class regions (sequence complement)
set -euo pipefail
P=${P:?set P to the project directory}; K=${KMC_BIN:?set KMC_BIN}; SAM=samtools; PY=python3
O=$P/${UM:-markers_u}; D=${KDB:-${SCRATCH:?set SCRATCH}/umk/db}; read a f < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $O/assemblies.tsv)
T=$(mktemp -d ${SCRATCH:?set SCRATCH}/umk/o.XXXX); trap "rm -rf $T" EXIT; echo "host $(hostname) $a start $(date)"
[ -s $f.fai ] || { echo "ERROR no fai for $f"; exit 1; }
for c in ${ONLY:-$(awk -F"\t" -v a=$a '$1==a{print $2}' $O/regions.tsv | sort -u) dyz18s dyz18sub}; do
  $PY -I - $O/regions.tsv $a $c $f.fai > $T/rg <<'PYEOF'
import sys
reg, a, c, fai = sys.argv[1:5]; cs = {'dyz18s': ('dyz18', 'yq3k1'), 'dyz18sub': ('dyz18', 'yq3k1', 'yq2k7'), 'dyz18suba': ('dyz18', 'yq3k1', 'yq2k7a')}.get(c, (c,)); iv = {}
for l in open(reg):
    x = l.rstrip('\n').split('\t')
    if x[0] == a and x[1] in cs: iv.setdefault(x[2], []).append((int(x[3]), int(x[4])))
if not iv: sys.exit(0)
for l in open(fai):
    ctg, n = l.split('\t')[:2]; n = int(n); p = 0
    for s, e in sorted(iv.get(ctg, [])):
        if s > p: print(f'{ctg}:{p + 1}-{s}')
        p = max(p, e)
    if p < n: print(f'{ctg}:{p + 1}-{n}')
PYEOF
  [ -s $T/rg ] || continue
  $SAM faidx -r $T/rg $f > $T/o.fa; $K/kmc -k31 -ci1 -cs255 -m4 -t4 -fm $T/o.fa $T/raw $T > /dev/null; $K/kmc_tools -t4 transform $T/raw set_counts 1 $D/$a.out_$c > /dev/null
  echo "$c outside: segments $(wc -l < $T/rg) bp $(grep -v ">" $T/o.fa | tr -d "\n" | wc -c)"; done
echo "done $a $(date)"
