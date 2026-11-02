"""Read the delivered ten-repeat depth benchmark and its measured coverage."""
import csv
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/model/depth"
FRACTIONS = ("1", "d20x", "d15x", "d10x", "d5x", "d2x", "d1x")   # 1: full depth (the CV); dNx: test reads thinned to N x base coverage (d30x, d25x in results/model/depth, not drawn)


@lru_cache(maxsize=1)
def load():
    depth = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    files, coverage = [], {}
    for fraction in FRACTIONS:
        path = SOURCE / f"depth_f{fraction}.tsv"
        files.append(path)
        seen, sample_depth = set(), {}
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                rep, fold = int(row["repeat"]), int(row["fold"])
                target, sample = row["class"], row["sample"]
                key = (target, rep, sample)
                assert key not in seen, (path, key)
                seen.add(key)
                measured = float(row["chr2_depth_test"])
                reference, estimate = float(row["reference_bp"]), float(row["estimate_bp"])
                signed = float(row["error_pct"])
                assert np.isfinite(measured) and measured > 0
                assert reference > 0 and abs(estimate - reference * (1 + signed / 100)) <= 1
                assert sample not in sample_depth or np.isclose(sample_depth[sample], measured), (path, sample)
                sample_depth[sample] = measured
                depth[fraction][target][rep].append({
                    "sample": sample, "fold": fold, "signed": signed,
                    "reference_bp": reference, "estimate_bp": estimate,
                    "chr2_depth": measured,
                })
        assert len(depth[fraction]) == 7, path
        for target, repeats in depth[fraction].items():
            assert set(repeats) == set(range(10)), (fraction, target)
            for rows in repeats.values():
                assert len(rows) == len(depth["1"][target][0]), (fraction, target)   # same men at every depth
                assert {row["fold"] for row in rows} == set(range(5))
        # Every individual counts once, independently of its number of targets.
        coverage[fraction] = float(np.median(list(sample_depth.values()))) * 150 / 120   # chr2 31-mer depth to base coverage, 150 bp reads (31-mers per read 120)
    return depth, files, coverage
