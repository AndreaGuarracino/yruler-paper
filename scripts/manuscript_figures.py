#!/usr/bin/env python3
"""Render the cohort map and the 1kGP snapshot figures (Fig. 1 is scripts/figure1.py).

The 1kGP table carries the grouped CV, lineage-out and depth evaluations,
with reference males marked in the population and pedigree panels.
No batch adjustment is applied.
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")


ROOT = Path(__file__).resolve().parents[1]

# target, available assembly references, measurement definition
REFERENCES = [
    ("DYZ3", 129, "Hallast 2026 array coordinates"),
    ("Yq12", 91, "Hallast 2026 region span"),
    ("DYZ1", 91, "sum of Hallast 2026 HMMER block lengths"),
    ("DYZ2", 91, "sum of Hallast 2026 HMMER block lengths"),
    ("DYZ18S", 91, "sum of Hallast 2026 HMMER block lengths (DYZ18 and YqHet3.1)"),
    ("DYZ19", 140, "Hallast 2026 array span, assembly QC"),
    ("DYZ17", 140, "longest run of DYZ17 unit hits (gaps <= 10 kb), assembly QC"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=ROOT / "data/catalog_d376/out_rule/catalog.tsv")
    parser.add_argument("--pairs", type=Path, default=ROOT / "data/catalog_d376/fs_pairs.tsv")
    parser.add_argument("--snapshot", type=Path, help="Dated frozen-model population table")
    parser.add_argument("--dyz19", type=Path, default=ROOT / "data/method/dyz19_clean.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "manuscript/figures")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with (ROOT / "results/figures/reference_targets.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["target", "reference_arrays", "length_definition"])
        writer.writerows(REFERENCES)
    from figures_1kgp import DEFAULT_SNAPSHOT, load_snapshot, validation_summary, render_snapshot
    snapshot = args.snapshot or DEFAULT_SNAPSHOT
    frozen = snapshot.exists()
    if frozen:
        raw, valid, pedigree = load_snapshot(snapshot)
        benchmarks, depth = validation_summary()
        import validation_plots as vp; vp.CV, vp.LO, vp.DEP, _ = vp.load()   # ED4, ED5, ED8 read these
    from figure2_tree import make_figure
    from prepare_cohort_pies import prepare
    prepare(snapshot if frozen else None)
    if frozen:
        render_snapshot(snapshot, raw, valid, pedigree, benchmarks, depth, args.out)
        return
    metadata = ROOT / "data/figure_metadata"
    make_figure(args.catalog, ROOT / "data/method/haplogroups_1kgp_isogg2019_decode1.tsv",
                metadata / "isogg2019_tree.tsv", metadata / "cohort_population_locations.tsv",
                metadata / "world_land.geojson", args.dyz19, args.out)


if __name__ == "__main__":
    main()
