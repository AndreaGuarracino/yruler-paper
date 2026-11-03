#!/usr/bin/env python3
"""Render Figure 1 and the cohort map from the 1kGP table.

The 1kGP table carries the grouped CV, lineage-out and depth evaluations,
with reference males marked in the population and pedigree panels.
No batch adjustment is applied.
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
DRAFT = "#914119"

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


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def save(fig, out, name):
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def panel(ax, label, title):
    ax.set_title(title, loc="left", fontsize=11, pad=12)
    ax.text(-0.10, 1.08, label, transform=ax.transAxes, fontsize=14,
            fontweight="bold", va="bottom")


def clean_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e8ebed", linewidth=0.7)


def method_figure(catalog, pairs_path, out, benchmarks=None, depth=None, snapshot_label=None):
    from figure2_tree import father_son_panels

    fig = plt.figure(figsize=(10.5, 9.6))
    diagram = fig.add_axes([0.055, 0.86, 0.90, 0.13])
    diagram.set_xlim(0, 1)
    diagram.set_ylim(0, 1)
    diagram.axis("off")
    fig.text(0.055, 1.0044, "a", fontsize=15, fontweight="bold", va="bottom")

    def box(x, y, width, height, text, color):
        diagram.add_patch(FancyBboxPatch(
            (x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.018",
            facecolor=color, edgecolor="#52646e", linewidth=0.8))
        diagram.text(x + width / 2, y + height / 2, text,
                     ha="center", va="center", fontsize=9.8, linespacing=1.4)

    def arrow(start, end):
        diagram.add_patch(FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=12,
            color="#52646e", linewidth=1.0))

    box(0.02, 0.66, 0.27, 0.24, "Y pangenome references\narray spans and unit labels", "#eef3f5")
    box(0.37, 0.66, 0.27, 0.24, "Training families\nselect markers and fit scale", "#e7f1f7")
    box(0.72, 0.66, 0.26, 0.24, "Test families\nassess prediction error", "#f9eee7")
    box(0.02, 0.15, 0.27, 0.24, "BAM/CRAM or SAM\nall eligible read sequences", "#eef3f5")
    box(0.37, 0.15, 0.27, 0.24, "Count marker k-mers\narray, single-copy Y, chr8 window", "#e7f1f7")
    box(0.72, 0.15, 0.26, 0.24, "Array lengths and composition\nY dosage", "#eaf4ef")
    arrow((0.30, 0.78), (0.355, 0.78))
    arrow((0.30, 0.27), (0.355, 0.27))
    arrow((0.65, 0.27), (0.705, 0.27))
    arrow((0.505, 0.645), (0.505, 0.405))
    arrow((0.85, 0.405), (0.85, 0.645))
    diagram.text(0.50, 0.04, "Marker selection and scale fitting use training families only",
                 ha="center", fontsize=9.4, color="#43525b")

    if benchmarks is None:
        validation = fig.add_axes([0.055, 0.40, 0.90, 0.25])
        validation.set_xlim(0, 1)
        validation.set_ylim(0, 1)
        validation.axis("off")
        panel(validation, "b", "Independent accuracy and generalization — validation plan")
        tests = [
            (0.015, 0.55, "Grouped cross-validation", "Families grouped; markers rebuilt within training folds", "Corrected errors and interval coverage pending"),
            (0.525, 0.55, "External test and Y-lineage generalization", "Untouched samples and lineages excluded from training", "External test planned; lineage evaluation pending"),
            (0.015, 0.06, "Depth requirements", "Read downsampling; performance assessed per class", "Planned"),
            (0.525, 0.06, "Across-cohort reproducibility", "Shared individuals; verify independent libraries", "Planned"),
        ]
        for x, y, title, description, status in tests:
            validation.add_patch(FancyBboxPatch((x, y), 0.45, 0.38,
                                 boxstyle="round,pad=0.012,rounding_size=0.012",
                                 facecolor="#faf7f1", edgecolor="#c7b99e", linewidth=0.7, linestyle="--"))
            validation.text(x + 0.014, y + 0.285, title, fontsize=9.5, fontweight="bold", color="#394f5b")
            validation.text(x + 0.014, y + 0.17, description, fontsize=8.3, color="#536975")
            validation.text(x + 0.014, y + 0.055, status, fontsize=8.2, color=DRAFT)
    else:
        from figures_1kgp import validation_panels
        validation_panels(fig, benchmarks, depth)

    # Keep the primary-array and composition comparisons together by target.
    from pedigree_plots import collect, scatter, differences, legend, LOOKUP
    observations=collect(catalog)
    for row,keys in enumerate((("dyz3","yq12","dyz19"),("dyz1","dyz2","dyz1_dyz2_ratio"))):
        bottom,top=(.285,.43) if row==0 else (.07,.215)
        grid=fig.add_gridspec(1,6,left=.085,right=.985,bottom=bottom,top=top,
                             width_ratios=[1.15,1.35]*3,wspace=.72)
        for col,key in enumerate(keys):
            scatter_ax=fig.add_subplot(grid[0,col*2]);delta_ax=fig.add_subplot(grid[0,col*2+1])
            scatter(scatter_ax,observations,key,compact=True)
            differences(delta_ax,observations,key,compact=True,reverse=False)
            x0=scatter_ax.get_position().x0;x1=delta_ax.get_position().x1
            fig.text((x0+x1)/2,top+.006,LOOKUP[key][2],ha="center",fontsize=9,color=LOOKUP[key][-1])
            if row==0:fig.text(x0-.03,top+.027,chr(ord('c')+col),fontsize=15,fontweight='bold')
        if row==1:fig.text(.055,top+.027,'f',fontsize=15,fontweight='bold')
    legend(fig,observations,y=.496,compact=True)
    save(fig, out, "figure1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=ROOT / "data/catalog_d376/out_rule/catalog.tsv")
    parser.add_argument("--pairs", type=Path, default=ROOT / "data/catalog_d376/fs_pairs.tsv")
    parser.add_argument("--snapshot", type=Path, help="Dated frozen-model population table")
    parser.add_argument("--dyz19", type=Path, default=ROOT / "data/method/dyz19_clean.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "manuscript/figures")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 10, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False,
    })
    with (args.out / "reference_targets.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["target", "reference_arrays", "length_definition"])
        writer.writerows(REFERENCES)
    from figures_1kgp import DEFAULT_SNAPSHOT, load_snapshot, legacy_catalog, validation_summary, render_snapshot
    snapshot = args.snapshot or DEFAULT_SNAPSHOT
    frozen = snapshot.exists()
    if frozen:
        raw, valid, pedigree = load_snapshot(snapshot)
        benchmarks, depth = validation_summary()
        method_figure(legacy_catalog(valid), args.pairs, args.out, benchmarks, depth,
                      f"preliminary {len(valid)}-male 1kGP snapshot")
    else:
        method_figure(read_tsv(args.catalog), args.pairs, args.out)
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
