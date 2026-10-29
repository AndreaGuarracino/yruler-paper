"""Read the delivered ten-repeat depth benchmark and its measured coverage."""
import csv
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/model/depth"
FRACTIONS = ("1", "0.667", "0.5", "0.333", "0.167", "0.067", "0.033")


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
        assert len(seen) == 7120 and len(sample_depth) == 118, path
        assert len(depth[fraction]) == 8, path
        for target, repeats in depth[fraction].items():
            assert set(repeats) == set(range(10)), (fraction, target)
            for rows in repeats.values():
                assert len(rows) == ({"dyz3": 104, "yq12": 81, "dyz2": 81, "dyz19": 118}.get(target, 82)), (fraction, target)
                assert {row["fold"] for row in rows} == set(range(5))
        # Every individual counts once, independently of its number of targets.
        coverage[fraction] = float(np.median(list(sample_depth.values())))
    return depth, files, coverage
