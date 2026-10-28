#!/usr/bin/env python3
"""Render preliminary manuscript figures from the current local evidence.

The dated snapshot uses the frozen procedure's grouped CV, lineage-out and
depth evaluations, with calibration references marked in population/pedigree
panels. Historical development checks remain in a separate supplement.
No batch adjustment is applied. An untouched external evaluation is pending.
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BLUE = "#277da8"
GREEN = "#338565"
ORANGE = "#a86828"
RED = "#b34a45"
GRAY = "#8c969c"
DRAFT = "#914119"

# target, available assembly references, measurement definition
REFERENCES = [
    ("DYZ3", 129, "active array span"),
    ("Yq12", 91, "region span"),
    ("DYZ1", 91, "assigned unit bases"),
    ("DYZ2", 91, "assigned unit bases"),
    ("DYZ18S", 91, "assigned unit bases"),
    ("DYZ18sub", 91, "DYZ18S plus YqHet2.7 assigned bases; overlapping output"),
    ("DYZ19", 139, "main array span"),
    ("YqHet2.7", 91, "assigned unit bases"),
]


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def save(fig, out, name):
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out / f"{name}.svg", bbox_inches="tight")
    fig.savefig(out / f"{name}.png", dpi=220, bbox_inches="tight")
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
    box(0.37, 0.15, 0.27, 0.24, "Count marker k-mers\narray, single-copy Y, chr2", "#e7f1f7")
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
        from frozen_snapshot import validation_panels
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
    save(fig, out, "figure1_draft")


def boxplot(ax, values, positions, colors, rng, widths=0.58, size=6):
    plot = ax.boxplot(values, positions=positions, widths=widths,
                      showfliers=False, patch_artist=True,
                      medianprops={"color": "#17272f", "linewidth": 1.1},
                      whiskerprops={"color": "#63727b", "linewidth": 0.8},
                      capprops={"color": "#63727b", "linewidth": 0.8})
    for patch, color in zip(plot["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.23)
        patch.set_edgecolor(color)
    for x, data, color in zip(positions, values, colors):
        jitter = rng.uniform(-widths * 0.32, widths * 0.32, len(data))
        ax.scatter(x + jitter, data, s=size, alpha=0.36,
                   color=color, linewidths=0, zorder=2)


def biology_figure(catalog, pairs_path, dyz19, out):
    raw_catalog = {row["sample"]: row for row in catalog}
    catalog = [row for row in catalog if row["sample"] not in {"HG02300", "NA19226"}
               and float(row["y_dosage"]) >= 0.9]
    usable = {row["sample"]: row for row in catalog}
    rng = np.random.default_rng(20261026)
    by_lineage = defaultdict(list)
    for row in catalog:
        by_lineage[row["haplogroup"][0]].append(row)
    groups = {key: rows for key, rows in by_lineage.items() if len(rows) >= 10}
    order = sorted(groups, key=lambda key: np.median([float(r["dyz3_bp"]) for r in groups[key]]))

    fig, axes = plt.subplots(2, 2, figsize=(9.6, 8.0),
                              gridspec_kw={"hspace": 0.48, "wspace": 0.30})
    for ax, key, label, title, color in [
        (axes[0, 0], "dyz3_bp", "a", "DYZ3 by Y lineage", BLUE),
        (axes[0, 1], "yq12_bp", "b", "Yq12 by Y lineage", GREEN),
    ]:
        values = [[float(row[key]) / 1e6 for row in groups[g]] for g in order]
        positions = np.arange(len(order))
        boxplot(ax, values, positions, [color] * len(order), rng, size=3.2)
        ax.set_xticks(positions, [f"{g}\n{len(groups[g])}" for g in order], fontsize=7.6)
        ax.set_xlabel("Major haplogroup; number of individuals", fontsize=9)
        ax.set_ylabel("Estimated length (Mb)")
        panel(ax, label, title)
        clean_axes(ax)

    dyz_groups = {"J2a": [], "J2b": [], "Other references": []}
    for row in dyz19:
        if row["clean"] != "True":
            continue
        span = int(row["end"]) - int(row["start"])
        if span != int(row["unit_bp_array"]):
            raise ValueError(f"DYZ19 span mismatch for {row['assembly']}")
        haplogroup = raw_catalog.get(row["assembly"], {}).get("haplogroup", "")
        group = "J2a" if haplogroup.startswith("J2a") else "J2b" if haplogroup.startswith("J2b") else "Other references"
        dyz_groups[group].append(span / 1000)
    ax = axes[1, 0]
    labels = list(dyz_groups)
    values = list(dyz_groups.values())
    boxplot(ax, values, np.arange(len(labels)), [RED, ORANGE, GRAY], rng, size=13)
    ax.set_xticks(np.arange(len(labels)), [f"{g}\nn={len(v)}" for g, v in dyz_groups.items()])
    ax.set_ylabel("Assembly array span (kb)")
    panel(ax, "c", "DYZ19 in assembly references")
    clean_axes(ax)

    paired = defaultdict(lambda: {"dyz3_bp": [], "yq12_bp": []})
    pair_counts = Counter()
    with pairs_path.open() as stream:
        for fields in csv.reader(stream, delimiter="\t"):
            if len(fields) != 2:
                raise ValueError("Expected father and son columns")
            father, son = fields
            if father not in usable or son not in usable:
                continue
            f, s = usable[father], usable[son]
            batch_pair = (f["batch"], s["batch"])
            pair_counts[batch_pair] += 1
            if batch_pair == ("2504", "698_related"):
                group = "Cross-batch"
            elif batch_pair == ("698_related", "698_related"):
                group = "Same batch"
            else:
                continue
            for key in ("dyz3_bp", "yq12_bp"):
                if float(f[key]) <= 0:
                    raise ValueError(f"Nonpositive paternal estimate for {father}")
                paired[group][key].append(100 * (float(s[key]) / float(f[key]) - 1))
    ax = axes[1, 1]
    for index, group in enumerate(("Cross-batch", "Same batch")):
        values = [paired[group][key] for key in ("dyz3_bp", "yq12_bp")]
        boxplot(ax, values, [index - 0.18, index + 0.18],
                [BLUE, GREEN], rng, widths=0.30, size=5)
    ax.set_xticks([0, 1], [f"{g}\nn={len(paired[g]['dyz3_bp'])} pairs"
                           for g in ("Cross-batch", "Same batch")])
    ax.set_xlim(-0.55, 1.55)
    ax.axhline(0, color="#4d5960", linewidth=0.9, linestyle="--")
    ax.set_yscale("symlog", linthresh=5)
    ax.set_yticks([-50, -20, -5, 0, 5, 20, 50])
    ax.set_yticklabels(["−50", "−20", "−5", "0", "5", "20", "50"])
    ax.set_ylabel("Signed son/father difference (%)")
    panel(ax, "d", "Pedigree differences before correction")
    ax.legend(handles=[Patch(facecolor=BLUE, label="DYZ3"), Patch(facecolor=GREEN, label="Yq12")],
              frameon=False, loc="upper right", fontsize=8.5)
    clean_axes(ax)

    fig.suptitle("Development dataset and complete-array references",
                 color=DRAFT, fontsize=10, y=0.99)
    fig.text(0.5, 0.015,
             "All shown points are retained; panel d uses a symmetric log scale outside ±5%",
             ha="center", fontsize=9, color="#43525b")
    save(fig, out, "figureS1_development_checks")

    print(f"Filtered catalog: {len(catalog)} individuals; {len(order)} displayed lineages")
    print(f"Quality-passing DYZ19 references: {sum(map(len, dyz_groups.values()))}")
    print("Filtered pedigree pairs:", dict(pair_counts))
    for group, classes in paired.items():
        for key, v in classes.items():
            print(f"{group} {key}: n={len(v)}, median signed={np.median(v):+.3f}%, "
                  f"median absolute={np.median(np.abs(v)):.3f}%")

    rows = []
    for g in order:
        for key in ("dyz3_bp", "yq12_bp"):
            v = np.array([float(row[key]) for row in groups[g]])
            rows.append(("development catalog", g, key, len(v), np.median(v)))
    for g, v in dyz_groups.items():
        rows.append(("assembly span", g, "dyz19_bp", len(v), np.median(v) * 1000))
    with (out / "figure2_summary.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["source", "group", "target", "n", "median_bp"])
        writer.writerows(rows)


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
    from frozen_snapshot import DEFAULT_SNAPSHOT, load_snapshot, legacy_catalog, validation_summary, render_snapshot
    snapshot = args.snapshot or DEFAULT_SNAPSHOT
    frozen = snapshot.exists()
    if frozen:
        raw, valid, pedigree = load_snapshot(snapshot)
        benchmarks, depth = validation_summary()
        method_figure(legacy_catalog(valid), args.pairs, args.out, benchmarks, depth,
                      f"preliminary {len(valid)}-male 1kGP snapshot")
    else:
        method_figure(read_tsv(args.catalog), args.pairs, args.out)
    biology_figure(read_tsv(args.catalog), args.pairs, read_tsv(args.dyz19), args.out)
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
