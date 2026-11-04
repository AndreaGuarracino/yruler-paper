#!/bin/bash
# Release markers as masked superstrings
# usage: R=<release dir> bash release_pack.sh; markers.fa -> markers.ms.fa.zst
set -euo pipefail
P=${P:?set P to the project directory}; B=${SPSS_BIN:?set SPSS_BIN to the masked-superstring tools}; R=${R:-$P/model/release}; T=$R/pack; mkdir -p $T; cd $T; echo "host $(hostname) start $(date)"
awk '/^>/{c=substr($0,2); next} {print > c".txt"}' $R/markers.fa; rm -f ms.fa
for c in $(grep '^>' $R/markers.fa | cut -c2-); do
  awk '{print ">"NR; print}' $c.txt > $c.kfa
  $B/kmercamel -k 31 -p $c.kfa -o $c.ms.txt 2> $c.camel.log || $B/kmercamel compute -k 31 -o $c.ms.txt $c.kfa 2>> $c.camel.log
  s=$(grep -v '>' $c.ms.txt | tr -d '\n'); up=$(echo -n "$s" | tr -cd 'ACGT' | wc -c)
  echo "$c kmers $(wc -l < $c.txt) superstring chars ${#s} uppercase $up"; [ "$up" = "$(wc -l < $c.txt)" ] || { echo "ERROR $c uppercase starts differ from k-mers"; exit 1; }
  printf '>%s\n%s\n' $c "$s" >> ms.fa; done
$B/zstd -19 -q -f ms.fa -o $R/markers.ms.fa.zst; ls -l $R/markers.fa ms.fa $R/markers.ms.fa.zst; md5sum $R/markers.ms.fa.zst; echo "end $(date)"
