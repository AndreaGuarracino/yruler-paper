#!/usr/bin/env python3
"""Figure 2: cohort geography and satellite measurements aligned to Y clades.

The tree is a pruned ISOGG 2019 haplogroup hierarchy, not an inferred sample
phylogeny. Branch lengths have no time or substitution interpretation. The
map uses reported population coordinates, before final sample selection.
"""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Wedge
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
from palette import CLASS, RATIO, TITLE, NAMES, FS, MM, LINEAGE, LAND, UNITS, LINE, MUTED
BLUE, GREEN, PURPLE = CLASS["DYZ3"], CLASS["Yq12"], CLASS["DYZ19"]   # DYZ3, Yq12, DYZ19
COHORTS = {"1kGP": ("#167dad", "D"), "HGDP": ("#d2942d", "s"), "SGDP": ("#b24c74", "o")}
LINEAGE_COLORS = LINEAGE
PIE_RADIUS_FACTOR = 0.765  # 10% smaller diameters; areas remain proportional to counts.
COMPOSITION = [
    ("dyz1", "DYZ1", CLASS["DYZ1"]),
    ("dyz2", "DYZ2", CLASS["DYZ2"]),
    ("dyz18s", "DYZ18S", CLASS["DYZ18S"]),
    ("rest", "Rest of Yq12", UNITS["Unannotated"]),   # final model parts: DYZ1 + DYZ2 + DYZ18S + rest
]
CUT = {
    "A00", "A0", "A1a", "A1b1", "B", "C", "D", "E", "F1", "F2", "F3", "F4",
    "G", "H", "I", "J1", "J2a", "J2b", "L", "T", "N", "O", "Q", "R1a", "R1b", "R2", "M", "S",
}
TIP_PRIORITY = {name: i for i, name in enumerate([
    "A00", "A0", "A1a", "A1b1", "B", "D", "E", "C", "F1", "F2", "F3", "F4",
    "G", "H", "I", "J1", "J2a", "J2b", "L", "T", "NO*", "N", "O", "M", "S", "Q", "R1a", "R1b", "R2",
])}


def read_tsv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def save(fig, out, name):
    from palette import save as put
    put(fig, out / f"{name}.pdf")
    plt.close(fig)


def world_map(ax, locations, land, label="a"):
    polygons = []
    for feature in land["features"]:
        geometry = feature["geometry"]
        rings = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
        polygons.extend(np.asarray(ring[0], dtype=float) for ring in rings)
    ax.add_collection(PolyCollection(polygons, facecolors="#e8ecee", edgecolors="#b5bfc4", linewidths=0.30))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-58, 84)
    ax.set_aspect("equal", adjustable="box")
    ax.set_facecolor("#fafcfd")
    counts = Counter(row["cohort"] for row in locations)
    # Open outer shapes keep coincident cohort locations visible.
    for cohort in ("HGDP", "1kGP", "SGDP"):
        rows = [row for row in locations if row["cohort"] == cohort]
        color, marker = COHORTS[cohort]
        ax.scatter([float(row["longitude"]) for row in rows],
                   [float(row["latitude"]) for row in rows],
                   marker=marker, s={"HGDP": 49, "1kGP": 39, "SGDP": 13}[cohort],
                   facecolors=color if cohort == "SGDP" else "none",
                   edgecolors=color, linewidths=0.85 if cohort != "SGDP" else 0.25,
                   alpha=0.82, zorder={"HGDP": 3, "1kGP": 4, "SGDP": 5}[cohort])
    handles = [Line2D([], [], linestyle="none", marker=COHORTS[c][1],
                      markerfacecolor=COHORTS[c][0] if c == "SGDP" else "none",
                      markeredgecolor=COHORTS[c][0], markersize=6,
                      label=f"{c} ({counts[c]})") for c in COHORTS]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.025), ncol=3,
              title="Populations with coordinates", title_fontproperties=TITLE(8.5),
              fontsize=9, frameon=False, columnspacing=2.0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Reported population locations of the three cohorts", loc="left", fontsize=11.5, pad=10)
    if label:
        ax.text(-0.05, 1.045, label, fontsize=15, fontweight="bold", transform=ax.transAxes)


def pie_positions(original, radii):
    """Separate overlapping circles, retaining their reported coordinate anchors."""
    positions = original.copy()
    for _ in range(260):
        delta = positions[:, None, :] - positions[None, :, :]
        distance = np.sqrt((delta ** 2).sum(axis=2))
        np.fill_diagonal(distance, np.inf)
        # Deterministic directions for coincident reported locations.
        coincident = distance < 1e-8
        if coincident.any():
            i, j = np.where(coincident)
            angle = (i * 0.731 + j * 2.17) % (2 * np.pi)
            delta[i, j, 0], delta[i, j, 1] = np.cos(angle), np.sin(angle)
            distance[i, j] = 1.0
        overlap = np.maximum(radii[:, None] + radii[None, :] + 0.65 - distance, 0)
        force = (delta / distance[:, :, None] * overlap[:, :, None] * 0.32).sum(axis=1)
        force += 0.018 * (original - positions)
        length = np.sqrt((force ** 2).sum(axis=1))
        positions += force * np.minimum(1, 1.2 / np.maximum(length, 1e-9))[:, None]
        positions[:, 0] = np.clip(positions[:, 0], -178 + radii, 178 - radii)
        positions[:, 1] = np.clip(positions[:, 1], -56 + radii, 83 - radii)
    return positions


def pie_map(ax, cohort, samples, land, size_legend=False, legend_y=-47):
    polygons = [np.asarray(feature["geometry"]["coordinates"][0], float) for feature in land["features"]]
    ax.add_collection(PolyCollection(polygons, facecolors=LAND["face"], edgecolors=LAND["edge"], linewidths=0.25))
    selected = samples if cohort is None else [row for row in samples if row["cohort"] == cohort]
    grouped = defaultdict(Counter)
    for row in selected:
        # Pool identical reported locations across cohorts after individual deduplication.
        grouped[(round(float(row["longitude"]), 6), round(float(row["latitude"]), 6))][row["lineage"]] += 1
    items = sorted(grouped.items(), key=lambda item: -sum(item[1].values()))
    original = np.array([key for key, _ in items])
    # A common radius factor makes circle area proportional to sample count.
    radii = PIE_RADIUS_FACTOR * np.sqrt([sum(counts.values()) for _, counts in items])
    centers = pie_positions(original, radii)
    for (key, counts), anchor, center, radius in zip(items, original, centers, radii):
        if np.linalg.norm(center - anchor) > 0.25:
            ax.plot([anchor[0], center[0]], [anchor[1], center[1]], color=LINE, linewidth=0.35, zorder=2)
            ax.plot(*anchor, ".", color=LINE, markersize=1.5, zorder=2)
        angle, n = 0.0, sum(counts.values())
        for lineage in LINEAGE_COLORS:
            value = counts.get(lineage, 0)
            if value:
                end = angle + 360 * value / n
                ax.add_patch(Wedge(center, radius, angle, end, facecolor=LINEAGE_COLORS[lineage],
                                   edgecolor="white", linewidth=0.25, zorder=3))
                angle = end
        ax.add_patch(Circle(center, radius, fill=False, edgecolor=MUTED, linewidth=0.3, zorder=4))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-58, 84)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    if cohort is not None:
        stage = "current catalog" if cohort == "1kGP" else "planned samples"
        ax.set_title(f"{cohort}  (n={len(selected):,})\n{stage}", fontsize=8.7, pad=4)
    if size_legend:
        longitude, previous_radius = -158.0, None
        for n in (1, 5, 10, 25, 50):
            radius = PIE_RADIUS_FACTOR * np.sqrt(n)
            if previous_radius is not None:
                longitude += previous_radius + radius + 3.0
            ax.add_patch(Circle((longitude, legend_y), radius, facecolor="white",
                                edgecolor=MUTED, linewidth=0.45, zorder=5))
            ax.text(longitude, legend_y - 10, str(n), ha="center", fontsize=FS["legend"], color=MUTED)
            previous_radius = radius
        ax.text(-159, legend_y + 5, "Individuals", fontsize=FS["legend"], ha="left", color=MUTED)


def lineage_legend(fig, samples, y, x=0.5, fontsize=7.6, title=None):
    counts = Counter(row["lineage"] for row in samples)
    handles = [Patch(facecolor=color, edgecolor="none", label=f"{name} (n={counts[name]:,})")
               for name, color in LINEAGE_COLORS.items() if counts[name]]
    # Keep the long Unassigned entry on the lowest row, below Australia.
    if len(handles)==20:
        handles.insert(-1,Patch(facecolor="none",edgecolor="none",label=""))
    return fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(x, y), ncol=7,
               fontsize=fontsize, frameon=False, labelspacing=0.6, borderaxespad=0.15, title=title, title_fontproperties=TITLE(fontsize))


def align_map_legend_titles(ax, legend):
    """Use identical black headings on one figure-coordinate baseline."""
    fig=ax.figure
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    labels=[t for t in legend.get_texts() if t.get_text()]
    font=labels[0].get_fontproperties().copy()
    top=max(t.get_window_extent(renderer).y1 for t in labels)
    bounds=legend.get_window_extent(renderer)
    heading=next(t for t in ax.texts if t.get_text()=='Individuals')
    left=fig.transFigure.inverted().transform(ax.transData.transform(heading.get_position()))[0]
    heading.remove()
    baseline=fig.transFigure.inverted().transform((0,top))[1]+4/72/fig.get_figheight()
    center=fig.transFigure.inverted().transform((bounds.x0+bounds.width/2,0))[0]+1.5*MM/fig.get_figwidth()   # over the legend columns
    for x,label,alignment in [(left,'Individuals','left'),(center,NAMES['lineage'],'center')]:
        bold=font.copy();bold.set_weight('bold')
        fig.text(x,baseline,label,fontproperties=bold,color='black',ha=alignment,va='baseline')
    for text in ax.texts:
        if text.get_text() in {'1','5','10','25','50'}:
            text.set_fontproperties(font);text.set_color('black')
    for text in labels:text.set_color('black')


def father_son_panels(axes, catalog, pairs_path, out, scatter_axes=None,
                      panel_labels=("c", "d", "e", "f"), legend_y=0.09,
                      output_prefix="figure1_pedigree", compact_legend=False, minimum_y_dosage=0.9):
    usable = {row["sample"]: row for row in catalog if row["sample"] not in {"HG02300", "NA19226"}
              and float(row["y_dosage"]) >= minimum_y_dosage}
    observations = []
    for fields in csv.reader(Path(pairs_path).open(), delimiter="\t"):
        father, son = fields
        if father not in usable or son not in usable:
            continue
        f, s = usable[father], usable[son]
        group = ("Cross batch" if (f["batch"], s["batch"]) == ("2504", "698_related") else
                 "Same batch" if f["batch"] == s["batch"] else "Reverse cross batch")
        for target in ("dyz3_bp", "yq12_bp"):
            denominator = float(f[target])
            if denominator <= 0:
                raise ValueError(f"Nonpositive paternal length: {father}, {target}")
            observations.append({"father": father, "son": son, "batch_group": group, "target": target,
                                 "father_bp": float(f[target]), "son_bp": float(s[target]),
                                 "father_reference": f.get("reference_male", "unknown"),
                                 "son_reference": s.get("reference_male", "unknown"),
                                 "signed_difference_pct": 100 * (float(s[target]) / denominator - 1)})
    with (ROOT / f"results/figures/{output_prefix}_pairs.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(observations[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(observations)
    shown = [row for row in observations if row["batch_group"] != "Reverse cross batch"]
    maximum = max(abs(row["signed_difference_pct"]) for row in shown) * 1.12
    rng = np.random.default_rng(1026)
    for target_index, (ax, target, color, title) in enumerate(zip(axes, ("dyz3_bp", "yq12_bp"), (BLUE, GREEN), ("DYZ3", "Yq12"))):
        labels = []
        for y, group in enumerate(("Cross batch", "Same batch")):
            values = np.array([row["signed_difference_pct"] for row in shown
                               if row["target"] == target and row["batch_group"] == group])
            q1, median, q3 = np.percentile(values, [25, 50, 75])
            selected = [row for row in shown if row["target"] == target and row["batch_group"] == group]
            jitter = y + rng.uniform(-0.17, 0.17, len(values))
            in_reference = np.array([r["father_reference"] == "yes" or r["son_reference"] == "yes" for r in selected])
            ax.scatter(values[~in_reference], jitter[~in_reference], s=5, color=color, alpha=0.30, linewidths=0)
            ax.scatter(values[in_reference], jitter[in_reference], s=10, facecolors="none",
                       edgecolors=color, alpha=0.6, linewidths=0.45)
            ax.plot([q1, q3], [y, y], color=color, linewidth=3, zorder=4)
            ax.plot(median, y, "|", color="#263e49", markersize=9, markeredgewidth=1.3, zorder=5)
            ax.annotate(f"{median:+.2f}%", (median, y), xytext=(4, 6), textcoords="offset points", fontsize=7.4, color="#263e49")
            labels.append(f"{group}\n(n={len(values)})")
        ax.set_yticks([0, 1], labels, fontsize=7.8)
        ax.set_ylim(1.45, -0.45)
        ax.set_xscale("symlog", linthresh=5)
        ax.set_xlim(-maximum, maximum)
        ticks = [x for x in (-100, -20, -5, 0, 5, 20, 100) if abs(x) <= maximum]
        ax.set_xticks(ticks, [str(x) for x in ticks], fontsize=7)
        ax.axvline(0, color="#697d88", linestyle="--", linewidth=0.7)
        ax.grid(axis="x")
        ax.set_axisbelow(True)
        ax.spines[["left", "right", "top"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        if scatter_axes is None:
            ax.set_title(title, color=color, fontsize=9.8, pad=12)
        else:
            ax.set_title(f"{title} difference", color=color, fontsize=9.8, pad=10)
        ax.text(-0.12, 1.065, panel_labels[2 * target_index + 1], transform=ax.transAxes,
                fontsize=15, fontweight="bold")
        ax.set_xlabel("Signed offspring/parent difference (%)", fontsize=8.2)
    if scatter_axes is not None:
        styles = {"Cross batch": ("o", None, 0.43), "Same batch": ("^", "#cf9036", 0.85),
                  "Reverse cross batch": ("D", "#435460", 0.9)}
        for target_index, (ax, target, color, title) in enumerate(zip(scatter_axes, ("dyz3_bp", "yq12_bp"), (BLUE, GREEN), ("DYZ3", "Yq12"))):
            selected = [row for row in observations if row["target"] == target]
            step = 1 if target == "dyz3_bp" else 20
            upper = step * np.ceil(1.05 * max(max(row["father_bp"], row["son_bp"]) for row in selected) / (1e6 * step))
            ax.plot([0, upper], [0, upper], color="#71838c", linestyle="--", linewidth=0.75, zorder=1)
            for group, (marker, point_color, alpha) in styles.items():
                subset = [row for row in selected if row["batch_group"] == group]
                for reference in (False, True):
                    shown_subset = [r for r in subset if (r["father_reference"] == "yes" or r["son_reference"] == "yes") == reference]
                    ax.scatter([row["father_bp"] / 1e6 for row in shown_subset],
                               [row["son_bp"] / 1e6 for row in shown_subset],
                               s=10 if group == "Cross batch" else 20, marker=marker,
                               facecolors="none" if reference else point_color or color,
                               edgecolors=point_color or color if reference else "white",
                               alpha=0.8 if reference else alpha, linewidths=0.55 if reference else 0.25, zorder=3)
            for row in sorted(selected, key=lambda row: -abs(row["signed_difference_pct"]))[:2]:
                x, y = row["father_bp"] / 1e6, row["son_bp"] / 1e6
                right = x > upper / 2
                offset = -8 if y > 0.8 * upper else 5
                ax.annotate(row["son"], (x, y), xytext=(-3 if right else 3, offset), textcoords="offset points",
                            ha="right" if right else "left", va="top" if offset < 0 else "bottom",
                            fontsize=5.8, color="#596c76")
            median = np.median([abs(row["signed_difference_pct"]) for row in selected])
            ax.text(0.04, 0.96, f"Median |Δ| {median:.2f}%", transform=ax.transAxes, va="top", fontsize=7.1, color="#324f5b")
            ax.set_xlim(0, upper)
            ax.set_ylim(0, upper)
            ax.set_aspect("equal", adjustable="box")
            ticks = np.arange(0, upper + step / 2, step)
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.tick_params(labelsize=7, length=2.5)
            ax.set_xlabel("Parent (Mb)", fontsize=8.2)
            ax.set_ylabel("Offspring (Mb)", fontsize=8.2, labelpad=2)
            ax.spines[["right", "top"]].set_visible(False)
            ax.set_title(title, color=color, fontsize=9.8, pad=10)
            ax.text(-0.25, 1.065, panel_labels[2 * target_index], transform=ax.transAxes,
                    fontsize=15, fontweight="bold")
        header = scatter_axes[0]
        handles = [Line2D([], [], linestyle="none", marker=marker, color=point_color or "#7895a1", markersize=4.5,
                          label=f"{group} (n={sum(row['batch_group'] == group and row['target'] == 'dyz3_bp' for row in observations)})")
                   for group, (marker, point_color, _) in styles.items()
                   if any(row["batch_group"] == group for row in observations)]
        if any(row["father_reference"] == "yes" or row["son_reference"] == "yes" for row in observations):
            handles.append(Line2D([], [], linestyle="none", marker="o", markerfacecolor="none",
                                  markeredgecolor="#536c79", markersize=4.5, label="Open: reference involved"))
        header.figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.53, legend_y),
                             ncol=4 if compact_legend else 3,
                             title=None if compact_legend else "Parent–offspring pairs", title_fontsize=8,
                             fontsize=6.8 if compact_legend else 7.4, frameon=False, columnspacing=1.3 if compact_legend else 1.7)
    summaries = []
    for target in ("dyz3_bp", "yq12_bp"):
        values = [abs(row["signed_difference_pct"]) for row in observations if row["target"] == target]
        print("Father-son",target,"pairs",len(values),"median absolute difference",f"{np.median(values):.2f}%")
        for group in ("All pairs", "Cross batch", "Same batch", "Reverse cross batch"):
            selected = [row for row in observations if row["target"] == target
                        and (group == "All pairs" or row["batch_group"] == group)]
            if not selected:
                continue
            signed = np.array([row["signed_difference_pct"] for row in selected])
            summaries.append({"target": target, "batch_group": group, "n_pairs": len(selected),
                              "n_fathers": len({row["father"] for row in selected}),
                              "n_sons": len({row["son"] for row in selected}),
                              "median_signed_pct": np.median(signed), "median_absolute_pct": np.median(abs(signed)),
                              "within_5pct_n": int((abs(signed) <= 5).sum()), "within_10pct_n": int((abs(signed) <= 10).sum())})
    with (ROOT / f"results/figures/{output_prefix}_summary.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    return observations


def tree_data(path):
    parents, children = {}, defaultdict(list)
    for line in Path(path).read_text().splitlines():
        node, parent, _ = line.split("\t")
        node, parent = node.split(" ")[0], parent.split(" ")[0]
        if node in parents:
            raise ValueError(f"Ambiguous normalized tree name: {node}")
        parents[node] = parent
        children[parent].append(node)
    return parents, children


def add_reference_calls(calls, parents):
    """Published calls for assemblies without a 1kGP call."""
    for row in read_tsv(ROOT / "data/method/reference_haplogroups_published.tsv"):
        h = row["haplogroup"]
        while h and h not in parents: h = h[:-1]   # deepest ISOGG 2019 node
        if row["sample"] not in calls: calls[row["sample"]] = h; print(f"reference call {row['sample']}: {row['haplogroup']} -> {h} ({row['source']})", flush=True)
    return calls


def clade_of(haplogroup, parents):
    node = haplogroup.split(" ")[0]
    if node not in parents:
        raise ValueError(f"Haplogroup absent from the reference hierarchy: {haplogroup}")
    original, seen = node, set()
    while node and node not in CUT:
        if node in seen:
            raise ValueError("Cycle in haplogroup hierarchy")
        seen.add(node)
        node = parents[node]
    # An unresolved ancestral call is retained as an explicitly marked tip.
    return node if node in CUT else original + "*"


def pruned_tree(groups, parents, original_children):
    include, group_nodes = set(), {}
    for group in groups:
        base = group.rstrip("*")
        group_nodes[group] = base
        node = base
        while node:
            include.add(node)
            node = parents[node]
    children = {node: [child for child in original_children[node] if child in include] for node in include}
    for group, base in group_nodes.items():
        if group.endswith("*"):
            children[base].append(group)
            children[group] = []

    def compact(node):
        branches = [compact(child) for child in children[node]]
        if len(branches) == 1 and node not in groups:
            return branches[0]
        return {"name": node, "children": branches}

    root = next(node for node in include if not parents[node])
    tree = compact(root)
    order = []

    def rotate(node):
        if not node["children"]:
            return TIP_PRIORITY[node["name"]]
        ranked = [(rotate(child), child) for child in node["children"]]
        ranked.sort(key=lambda item: item[0])
        node["children"] = [child for _, child in ranked]
        return ranked[0][0]

    rotate(tree)

    def walk(node):
        if not node["children"]:
            if node["name"] not in groups:
                raise ValueError(f"Unexpected displayed leaf: {node['name']}")
            order.append(node["name"])
        for child in node["children"]:
            walk(child)

    walk(tree)
    return tree, order


def place_tree(tree, order):
    ranks = {name: i for i, name in enumerate(order)}

    def visit(node, depth):
        node["depth"] = depth
        if not node["children"]:
            node["y"] = ranks[node["name"]]
            return depth
        last = max(visit(child, depth + 1) for child in node["children"])
        node["y"] = np.mean([child["y"] for child in node["children"]])
        return last

    maximum = visit(tree, 0)

    def assign(node):
        node["x"] = maximum + 0.6 if not node["children"] else node["depth"]
        for child in node["children"]:
            assign(child)

    assign(tree)
    return maximum + 0.6


def draw_tree(ax, tree, order, counts):
    tip_x = place_tree(tree, order)

    def draw(node):
        children = node["children"]
        if children:
            ax.plot([node["x"], node["x"]], [min(c["y"] for c in children), max(c["y"] for c in children)],
                    color="#52616a", linewidth=0.85)
            for child in children:
                ax.plot([node["x"], child["x"]], [child["y"], child["y"]], color="#52616a", linewidth=0.85)
                draw(child)
            if node["name"] in {"BT", "CT", "DE", "CF", "IJ", "NO", "R"}:
                ax.annotate(node["name"], (node["x"], node["y"]), xytext=(-2, 4), textcoords="offset points",
                            color="#677780", fontsize=6.8, ha="right")
        else:
            name = node["name"]
            major = "NO" if name.startswith("NO") else name[0]
            color = LINEAGE_COLORS.get(major, "#263b46")
            ax.plot(tip_x, node["y"], "o", markersize=2.5, color=color)
            ax.text(tip_x + 0.17, node["y"], f"{name} ({counts[name]})", va="center", fontsize=8.7,
                    color=color, fontweight="bold" if name in {"J2a", "R1b"} else "normal")

    draw(tree)
    ax.set_xlim(-0.15, tip_x + 3.8)
    ax.set_ylim(len(order) - 0.3, -0.8)
    ax.axis("off")
    ax.set_title("Y haplogroup (n)", loc="left", fontsize=10.5, pad=10)
    ax.text(-0.09, 1.035, "b", fontsize=15, fontweight="bold", transform=ax.transAxes)


def summaries(ax, order, groups, key, color, scale, limit, ticks, title, panel, references=False):
    rng = np.random.default_rng(426 if references else 319)
    rows = []
    for index, group in enumerate(order):
        values = np.asarray([float(row[key]) / scale for row in groups.get(group, [])])
        if len(values):
            ax.scatter(values, index + rng.uniform(-0.16, 0.16, len(values)),
                       s=8 if references else 3.5, color=color, alpha=0.65 if references else 0.21,
                       linewidths=0, zorder=3)
            q1, median, q3 = np.percentile(values, [25, 50, 75])
            ax.plot([q1, q3], [index, index], color=color, linewidth=3.0, solid_capstyle="butt", zorder=4)
            ax.plot(median, index, "|", markersize=8, markeredgewidth=1.3, color="#253c47", zorder=5)
            rows.append((group, key, len(values), median * scale, q1 * scale, q3 * scale))
        else:
            ax.text(limit * 0.48, index, "—", color="#b0bbc0", ha="center", va="center", fontsize=8)
        if references:
            ax.text(limit * 1.025, index, str(len(values)), color="#70808a", fontsize=7.5,
                    ha="left", va="center", clip_on=False)
    ax.set_xlim(0, limit)
    ax.set_ylim(len(order) - 0.3, -0.8)
    ax.set_xticks(ticks)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=8)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_title(title, loc="center", fontsize=10.5, pad=10, color=color)
    ax.text(-0.06, 1.035, panel, fontsize=15, fontweight="bold", transform=ax.transAxes, color="#182a34")
    ax.set_xlabel("Assembly span (kb)" if references else "Estimated length (Mb)", fontsize=8.7)
    if references:
        ax.text(1.025, 1.01, "n", transform=ax.transAxes, fontsize=8, color="#70808a")
    return rows


def load_composition(annotations_path, regions_path, calls, parents, displayed_groups):
    annotations = {row["assembly"]: row for row in read_tsv(annotations_path)}
    drops = {tuple(l.split()[:2]) for l in open(ROOT / "data/method/truth_qc/x/rule_w4.txt") if l.strip()}   # truth exclusions of the final model
    rows = []
    for reference in read_tsv(regions_path):
        if reference["class"] != "yq12" or ("yq12", reference["male"]) in drops:
            continue
        sample, assembly = reference["male"], reference["asm"]
        annotation = annotations[assembly]
        span = int(annotation["yq12_span"])
        if annotation["contig"] != reference["contig"] or span != int(annotation["yq12_end"]) - int(annotation["yq12_start"]):
            raise ValueError(f"Inconsistent Yq12 annotation coordinates: {assembly}")
        unannotated = int(annotation["non_dyz"]) if "non_dyz" in annotation else span - sum(int(annotation[k]) for k in ("dyz1_units", "dyz2_units", "dyz18_units", "q3k1_units", "q2k7_units"))
        components = {
            "dyz1": int(annotation["dyz1_units"]),
            "dyz2": int(annotation["dyz2_units"]),
            "dyz18s": int(annotation["dyz18_units"]) + int(annotation["q3k1_units"]),   # DYZ18 and the 3.1-kb repeat
            "rest": int(annotation["q2k7_units"]) + unannotated,   # 2.7-kb repeat and unannotated sequence   # Hallast table: span minus annotated bases
        }
        if span <= 0 or min(components.values()) < 0 or sum(components.values()) != span:
            raise ValueError(f"Yq12 components do not partition the annotated span: {assembly}")
        group = clade_of(calls[sample], parents) if sample in calls else "Unassigned"
        if group not in displayed_groups:
            group = "Unassigned"
        row = {"sample": sample, "assembly": assembly, "reference_source": reference["source"],
               "contig": annotation["contig"], "haplogroup_2019": calls.get(sample, ""),
               "display_group": group, "reference_span_bp": int(reference["end"]) - int(reference["start"]),
               "annotation_start": int(annotation["yq12_start"]), "annotation_end": int(annotation["yq12_end"]),
               "annotated_span_bp": span}
        row.update({key + "_bp": value for key, value in components.items()})
        row.update({key + "_fraction": value / span for key, value in components.items()})
        rows.append(row)
    if len({row["sample"] for row in rows}) != len(rows):
        raise ValueError("More than one retained Yq12 assembly per individual")
    return rows


def composition_panel(ax, order, rows, out):
    grouped = {group: [row for row in rows if row["display_group"] == group] for group in order}
    summaries = []
    rng = np.random.default_rng(2701)
    for y, group in enumerate(order):
        individuals = grouped[group]
        if individuals:
            fractions = np.array([[row[key + "_fraction"] for key, _, _ in COMPOSITION] for row in individuals])
            means = 100 * fractions.mean(axis=0)
            left = 0.0
            for (key, label, color), mean, column in zip(COMPOSITION, means, fractions.T):
                ax.barh(y, mean, left=left, height=0.53, color=color, edgecolor="white", linewidth=0.2, zorder=2)
                left += mean
                summaries.append((group, key, len(individuals), mean, 100 * np.median(column),
                                  100 * np.percentile(column, 25), 100 * np.percentile(column, 75)))
            ax.scatter(100 * fractions[:, 0], y + rng.uniform(-0.14, 0.14, len(individuals)),
                       s=7, facecolors="white", edgecolors="#26444b", linewidths=0.4, alpha=0.85, zorder=3)
        else:
            ax.text(50, y, "—", ha="center", va="center", fontsize=8, color="#b0bbc0")
        ax.text(103, y, str(len(individuals)), va="center", fontsize=7.5, color="#70808a", clip_on=False)
    ax.set_xlim(0, 100)
    ax.set_ylim(len(order) - 0.3, -0.8)
    ax.set_xticks([0, 50, 100])
    ax.tick_params(axis="x", labelsize=8)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.set_xlabel("Annotated bases (%)", fontsize=8.7)
    ax.set_title("Yq12 composition", fontsize=10, pad=10, color="#315b60")
    ax.text(-0.06, 1.035, "f", transform=ax.transAxes, fontsize=15, fontweight="bold", color="#182a34")
    ax.text(1.03, 1.01, "n", transform=ax.transAxes, fontsize=8, color="#70808a")
    ax.legend(handles=[Patch(facecolor=color, label=label) for _, label, color in COMPOSITION], title=NAMES["class"], title_fontproperties=TITLE(6.7),
              loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=6.7, frameon=False)
    with (ROOT / "results/figures/figure2_composition_samples.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (ROOT / "results/figures/figure2_composition_summary.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["haplogroup", "component", "n", "mean_pct", "median_pct", "q25_pct", "q75_pct"])
        writer.writerows(summaries)


def composition_supplement(rows, order, out, pop=()):
    """ED2: a, parts of Yq12 in the references; b, DYZ1 share against Yq12 length within the five largest lineages (1kGP, Yruler)."""
    from scipy.stats import spearmanr
    ranks = {group: i for i, group in enumerate(order + ["Unassigned"])}
    ordered = sorted(rows, key=lambda row: (ranks[row["display_group"]], row["dyz1_fraction"], row["sample"]))
    fig = plt.figure(figsize=(180 * MM, 120 * MM)); a = fig.add_axes([.07, .47, .92, .44])
    gi = {g: k for k, g in enumerate(dict.fromkeys(row["display_group"] for row in ordered))}   # lineage index of each bar
    x = np.arange(len(ordered)) + 0.6 * np.array([gi[row["display_group"]] for row in ordered]); bottom = np.zeros(len(ordered))   # small gap between lineages
    for key, label, color in COMPOSITION:
        values = np.array([100 * row[key + "_fraction"] for row in ordered])
        a.bar(x, values, bottom=bottom, width=0.92, color=color, edgecolor="none", label=label); bottom += values
    centers, labels = [], []
    for group in order + ["Unassigned"]:
        indices = [i for i, row in enumerate(ordered) if row["display_group"] == group]
        if indices:
            centers.append(np.mean(x[indices])); labels.append(f"{group} ({len(indices)})" if group != "Unassigned" else f"No lineage ({len(indices)})")
            if len(centers) > 1: centers[-1] = max(centers[-1], centers[-2] + 1.6)   # labels of adjacent single bars apart
    a.set_xticks(centers, labels, rotation=90); a.set_xlim(x[0] - 1, x[-1] + 1); a.set_ylim(0, 100); a.set_ylabel("Annotated Yq12 bases (%)")
    leg = a.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=4, borderaxespad=0.2, columnspacing=0.8)   # one line at the top right
    fig.canvas.draw(); bb = leg.get_window_extent().transformed(a.transAxes.inverted())
    a.text(bb.x0, (bb.y0 + bb.y1) / 2, NAMES["class"], transform=a.transAxes, ha="right", va="center", fontsize=FS["legend"], fontweight="bold")   # legend title on the same line
    a.text(.33, 1.0 + 4 / 72 / (fig.get_figheight() * .44), f"Parts of Yq12 in the {len(ordered)} reference assemblies", transform=a.transAxes, ha="center", va="baseline", fontsize=FS["title"], fontweight="bold")   # title left of the legend
    major = lambda g: g if g in ("R1a", "R1b") else ("NO" if g.startswith("NO") else g[0])
    lins = ["E", "R1b", "O", "R1a", "J"]
    for k, lin in enumerate(lins):
        ax = fig.add_axes([.07 + k * .188, .07, .155, .2]); c = LINEAGE_COLORS.get(lin[0], LINE)
        P = [r for r in pop if r["lineage"] == lin]; y = np.array([r["Yq12_bp"] for r in P], float) / 1e6; sh = 100 * np.array([r["DYZ1_bp"] for r in P], float) / 1e6 / y
        ax.scatter(y, sh, s=3, color=c, alpha=.45, lw=0, zorder=2)
        A = [r for r in ordered if major(r["display_group"]) == lin]
        ax.scatter([r["annotated_span_bp"] / 1e6 for r in A], [100 * r["dyz1_fraction"] for r in A], marker="D", s=10, facecolor="white", edgecolor="k", lw=.6, zorder=4)
        sl, ic = np.polyfit(y, sh, 1); xs = np.array([np.percentile(y, 2), np.percentile(y, 98)]); ax.plot(xs, ic + sl * xs, color="k", lw=.9, zorder=3)
        rho = spearmanr(y, sh).statistic; ax.set_title(f"{lin} ({len(P)}), ρ = {rho:.2f}", color=c, fontsize=FS["title"] - 1.5)
        ax.set_xlim(np.percentile(y, .5) - 1, np.percentile(y, 99.5) + 1); ax.set_ylim(47, 68)
        ax.set_xlabel("Yq12 length (Mb)")
        if k == 0: ax.set_ylabel("DYZ1 fraction of Yq12 (%)")
        else: ax.set_yticklabels([])
        ax.grid(True); ax.set_axisbelow(True)
        print(f"ED2b {lin}: {len(P)} males, Spearman {rho:.3f}, slope {sl:.3f} points per Mb, assemblies {len(A)}", flush=True)
    nt = sum(1 for r in pop if r["lineage"] in lins)   # males in the five panels
    fig.text(.45, .335, f"DYZ1 fraction against Yq12 length within lineages in 1000 Genomes ({nt:,})", fontsize=FS["title"], fontweight="bold", ha="center")
    fig.add_artist(plt.Line2D([.865], [.343], marker="D", ms=3.5, mfc="white", mec="k", mew=.6, ls="none", transform=fig.transFigure)); fig.text(.873, .336, "= assemblies", fontsize=FS["legend"])
    ta = .91 + 4 / 72 / fig.get_figheight()   # baseline of the panel a title
    for s_, x_, y_ in (("a", .005, ta), ("b", .005, .335)): fig.text(x_, y_, s_, fontsize=FS["letter"], fontweight="bold", va="baseline")   # letters on the title baselines
    for axx in fig.axes: axx.spines[["right", "top"]].set_visible(False)
    save(fig, out, "ed/ed2_yq12_composition")


def make_figure(catalog_path, calls_path, tree_path, locations_path, land_path, dyz19_path, out):
    catalog = read_tsv(catalog_path)
    calls = {row["ID"]: row["Haplogroup"] for row in read_tsv(calls_path)}
    parents, children = tree_data(tree_path)
    add_reference_calls(calls, parents)
    groups = defaultdict(list)
    for row in catalog:
        if row["sample"] in {"HG02300", "NA19226"} or float(row["y_dosage"]) < 0.9:
            continue
        group = clade_of(calls[row["sample"]], parents)
        groups[group].append(row)
    counts = {group: len(rows) for group, rows in groups.items()}
    tree, order = pruned_tree(groups, parents, children)
    eligible = sum(row["sample"] not in {"HG02300", "NA19226"} and float(row["y_dosage"]) >= 0.9
                   for row in catalog)
    assert len(order) == len(groups) and sum(counts.values()) == eligible

    reference_groups, unmapped = defaultdict(list), []
    for row in read_tsv(dyz19_path):
        if row["clean"] != "True":
            continue
        if int(row["end"]) - int(row["start"]) != int(row["unit_bp_array"]):
            raise ValueError("DYZ19 target is not the recorded array span")
        sample = row["assembly"]
        if sample not in calls:
            unmapped.append(sample)
            continue
        group = clade_of(calls[sample], parents)
        if group not in groups:
            unmapped.append(sample)
            continue
        reference_groups[group].append(row)

    locations = read_tsv(locations_path)
    land = json.loads(Path(land_path).read_text())
    from cohort_lineages import load_combined_samples
    pie_samples = load_combined_samples()
    lineage_counts = Counter(row["lineage"] for row in pie_samples)
    with (ROOT / "results/figures/figure2_lineage_counts.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["lineage", "n_unique_individuals", "percent"])
        for lineage in LINEAGE_COLORS:
            if lineage_counts[lineage]:
                writer.writerow([lineage, lineage_counts[lineage], 100 * lineage_counts[lineage] / len(pie_samples)])
    composition = load_composition(ROOT / "data/method/indep_setB.tsv", ROOT / "data/markers_u/truth_regions.tsv",
                                   calls, parents, groups)
    fig = plt.figure(figsize=(9.8, 9.2))
    map_ax = fig.add_axes([0.045, 0.61, 0.91, 0.37])
    pie_map(map_ax, None, pie_samples, land, size_legend=True)
    map_ax.text(-0.025, 1.035, "a", transform=map_ax.transAxes, fontsize=15, fontweight="bold")
    lineage_legend(fig, pie_samples, 0.595, title=NAMES["lineage"])
    grid = fig.add_gridspec(1, 5, left=0.055, right=0.96, bottom=0.12, top=0.49,
                            width_ratios=[1.98, 1.03, 1.03, 1.12, 1.45], wspace=0.22)
    tree_ax = fig.add_subplot(grid[0, 0])
    axes = [tree_ax] + [fig.add_subplot(grid[0, column]) for column in (1, 2, 3, 4)]
    for index, group in enumerate(order):
        if group == "J2a":
            axes[3].axhspan(index - 0.45, index + 0.45, facecolor="#f8e9ec", zorder=0)
        if group == "R1b":
            axes[1].axhspan(index - 0.45, index + 0.45, facecolor="#e8f1f7", zorder=0)
    draw_tree(tree_ax, tree, order, counts)
    stats = summaries(axes[1], order, groups, "dyz3_bp", BLUE, 1e6, 3.05, [0, 1, 2, 3], "DYZ3", "c")
    stats += summaries(axes[2], order, groups, "yq12_bp", GREEN, 1e6, 60, [0, 20, 40, 60], "Yq12", "d")
    stats += summaries(axes[3], order, reference_groups, "unit_bp_array", PURPLE, 1000, 500, [0, 250, 500],
                       "DYZ19", "e", references=True)
    composition_panel(axes[4], order, composition, out)
    fig.suptitle("PRELIMINARY — 1kGP development estimates; assembly references in e and f", color="#914119", fontsize=9.5, y=1.035)
    fig.text(0.39, 0.033, "c–e: interquartile ranges and medians.  f: mean composition; points, individual DYZ1 fractions.",
             ha="center", fontsize=7.5, color="#4d606b")
    fig.text(0.39, 0.012, "Haplogroup branches are schematic.  * Unresolved within NO.", ha="center", fontsize=7.5, color="#4d606b")
    plt.close(fig); print("old Fig. 2 layout not saved: main Fig. 2 comes from manuscript_figures.py", flush=True)
    composition_supplement(composition, order, out)

    map_fig, map_ax = plt.subplots(figsize=(10.0, 4.5))
    world_map(map_ax, locations, land, label=None)
    map_fig.text(0.5, -0.015, "Public cohort population metadata; final sample inclusion and quality filters pending.",
                 ha="center", color="#657780", fontsize=9)
    save(map_fig, ROOT / "results/figures", "cohort_population_coverage")

    pie_fig, pie_ax = plt.subplots(figsize=(12.5, 6.3))
    pie_map(pie_ax, None, pie_samples, land, size_legend=True)
    pie_fig.subplots_adjust(left=0.025, right=0.975, top=0.93, bottom=0.25)
    lineage_legend(pie_fig, pie_samples, 0.23)
    pie_fig.text(0.5, 0.065, "Pie area is proportional to unique individuals per reported location; colours show Y lineages; grey is unassigned.",
                 ha="center", fontsize=9, color="#4d606b")
    pie_fig.text(0.5, 0.035, "Shared individuals are counted once. Leaders connect displaced pies to reported coordinates. HGDP/SGDP final QC is pending.",
                 ha="center", fontsize=9, color="#4d606b")
    save(pie_fig, ROOT / "results/figures", "cohort_world_map")

    with (ROOT / "results/figures/figure2_tree_summary.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["haplogroup", "target", "n", "median_bp", "q25_bp", "q75_bp"])
        writer.writerows(stats)
    with (ROOT / "results/figures/figure2_tree_membership.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["sample", "haplogroup_2019", "display_group", "source"])
        for group in order:
            for row in groups[group]:
                writer.writerow([row["sample"], calls[row["sample"]], group, "development catalog"])
            for row in reference_groups[group]:
                writer.writerow([row["assembly"], calls[row["assembly"]], group, "DYZ19 assembly reference"])
    print("Displayed groups:", ", ".join(f"{name} ({counts[name]})" for name in order))
    print("Mapped cohort populations:", dict(Counter(row["cohort"] for row in locations)))
    print("DYZ19 references displayed:", sum(map(len, reference_groups.values())), "; metadata unavailable:", len(unmapped))
    print("Yq12 composition references:", len(composition), "; assigned to displayed clades:",
          sum(row["display_group"] != "Unassigned" for row in composition))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=ROOT / "data/catalog_d376/out_rule/catalog.tsv")
    parser.add_argument("--calls", type=Path, default=ROOT / "data/method/haplogroups_1kgp_isogg2019_decode1.tsv")
    parser.add_argument("--tree", type=Path, default=ROOT / "data/figure_metadata/isogg2019_tree.tsv")
    parser.add_argument("--locations", type=Path, default=ROOT / "data/figure_metadata/cohort_population_locations.tsv")
    parser.add_argument("--land", type=Path, default=ROOT / "data/figure_metadata/world_land.geojson")
    parser.add_argument("--dyz19", type=Path, default=ROOT / "data/method/dyz19_clean.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "manuscript/figures")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    make_figure(args.catalog, args.calls, args.tree, args.locations, args.land, args.dyz19, args.out)


if __name__ == "__main__":
    main()
