#!/bin/bash
# superset db and sorted keys
set -euo pipefail
P=/path/to/project; O=$P/${UM:-markers_u}; K=/path/to/conda/kmc/3.2.4/bin; PY=/path/to/conda/ktrk-bench/1.0/bin/python
T=$(mktemp -d /path/to/scratch/umk/u.XXXX); trap "rm -rf $T" EXIT
$K/kmc -k31 -ci1 -cs255 -m8 -t4 -fm $O/markers_u.fa $T/raw $T > /dev/null; $K/kmc_tools -t4 transform $T/raw set_counts 1 $O/u_db > /dev/null
$K/kmc_tools -t4 transform $O/u_db dump $T/u.txt > /dev/null
$PY -I -c "
import sys; sys.argv=['x']; exec(open('$P/scripts/um_bits.py').read().split('U = np.load')[0]); import numpy as np
k=enc('$T/u.txt'); assert (k[1:]>k[:-1]).all(); np.save('$O/u_keys.npy', k); print('superset distinct k-mers', len(k))"
