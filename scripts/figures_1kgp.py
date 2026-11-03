#!/usr/bin/env python3
"""Traceable figures and summaries for the dated frozen-model 1kGP snapshot.

Population estimates, calibration references and validation predictions remain
separate. Reading this module never changes the frozen model or stream inputs.
"""

import ast
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import numpy as np
from scipy.stats import spearmanr

from figure2_tree import (
    BLUE, GREEN, PURPLE, LINEAGE_COLORS, read_tsv, save,
    tree_data, pie_map, lineage_legend, align_map_legend_titles, clade_of, pruned_tree,
    load_composition, composition_supplement,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SNAPSHOT = ROOT / "results/estimates_1kgp.tsv"
TREE_SOURCE = ROOT / "results/ytree_v4_2026-11-01/ytree_v4.rooted.nwk"
FROZEN = ROOT / "data/frozen_2026-10-27"
TARGETS = [
    ("dyz3", "DYZ3_bp", "DYZ3", 1e3, "kb"),
    ("yq12", "Yq12_bp", "Yq12", 1e6, "Mb"),
    ("dyz1", "DYZ1_bp", "DYZ1", 1e6, "Mb"),
    ("dyz2", "DYZ2_bp", "DYZ2", 1e6, "Mb"),
    ("dyz18s", "DYZ18S_bp", "DYZ18S", 1e3, "kb"),
    ("dyz19", "DYZ19_bp", "DYZ19", 1e3, "kb"),
    ("dyz17", "DYZ17_bp", "DYZ17", 1e3, "kb"),
]
CUT = {
    "A00", "A0", "A1a", "A1b1", "B", "C", "D", "E", "F1", "F2", "F3", "F4",
    "G", "H", "I", "J1", "J2a1a1a", "J2a1a1b", "J2a1a2", "J2a2", "J2b", "L", "T",
    "N", "O1a", "O1b", "O2a", "O2b", "Q", "R1a", "R1b", "R2", "M", "S",
}
ORDER = [
    "A00", "A0", "A1a", "A1b1", "B", "D", "E", "C", "F1", "F2", "F3", "F4",
    "G", "H", "I", "J1", "J2a1*", "J2a1a1a", "J2a1a1b", "J2a1a2", "J2a2", "J2b",
    "L", "T", "NO*", "N", "O1a", "O1b", "O2a", "O2b", "M", "S", "Q",
    "R1a", "R1b", "R2",
]
from depth_benchmark import SOURCE as DEPTH_SOURCE, load as load_depth_benchmark
DEPTHS = load_depth_benchmark()[2]


def write_rows(path, rows):
    with Path(path).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

def snapshot_tag(path):
    return Path(path).stem.removeprefix("prelim_1kgp_").removeprefix("estimates_")


def load_snapshot(path=DEFAULT_SNAPSHOT):
    raw = read_tsv(path)
    reported_status = {r["sample"]: r for r in read_tsv(ROOT / "data/figure_metadata/y_status_reported.tsv")}
    assert len({r["sample"] for r in raw}) == len(raw)
    parents, _ = tree_data(ROOT / "data/figure_metadata/isogg2019_tree.tsv")
    pedigree = {}
    with (ROOT / "data/catalog_d376/ped.txt").open() as stream:
        header = next(stream).split()
        for line in stream:
            record = dict(zip(header, line.split()))
            pedigree[record["SampleID"]] = record
    valid = []
    for record in raw:
        if record["flag"] not in {"ok", "y_dose_below_0.9"}:
            continue
        assert record["sex_ped"] == "1" and float(record["y_dose"]) > 0
        row = dict(record)
        for _, key, _, _, _ in TARGETS:
            row[key] = float(row[key])
            assert np.isfinite(row[key]) and row[key] > 0, (row["sample"], key)
        row["y_dose"] = float(row["y_dose"])
        row["reference"] = row["reference_male"] == "yes"
        row["family"] = pedigree.get(row["sample"], {}).get("FamilyID", row["sample"])
        row["display_group"] = group_of(row["haplogroup"], parents)
        # Keep the source lineage unchanged; an ancestral NO call is not a resolved N call.
        row["lineage_display"] = "NO" if row["haplogroup"].split()[0] == "NO" else row["lineage"]
        row["ratio_bp"] = row["DYZ1_bp"] / row["DYZ2_bp"]
        supplied = next((record[k] for k in ("y_state", "y_status", "Y_status") if record.get(k)), None)
        status = supplied or reported_status.get(row["sample"], {}).get("y_status", "not_reported")
        status = {"riarrangiamento": "rearrangement", "mosaico": "mosaic"}.get(status.lower(), status)
        row["y_status"] = status
        row["y_status_source"] = (str(Path(path).resolve().relative_to(ROOT)) if Path(path).resolve().is_relative_to(ROOT) else Path(path).name) if supplied else reported_status.get(row["sample"], {}).get("source", "not_reported")
        row["population_included"] = status.lower() not in {"rearrangement", "partial", "yp_deletion", "yq_loss"}
        valid.append(row)
    return raw, valid, pedigree


def in_summary(row, key):
    """Population summary of one target: a PAR2 gain leaves Yq12, DYZ1, DYZ2 and their ratio."""
    return row["population_included"] and not (row["y_status"] == "par2_gain" and key in {"Yq12_bp", "DYZ1_bp", "DYZ2_bp", "ratio_bp", "unit_ratio"})


def population_rows(rows):
    """Keep flagged arm profiles visible, but outside population summaries."""
    return [row for row in rows if row.get("population_included", True)]


def group_of(call, parents):
    node = call.split(" ")[0]
    assert node in parents, call
    original = node
    while node and node not in CUT:
        node = parents[node]
    return node if node in CUT else original + "*"


def legacy_catalog(rows):
    """Adapt eligible frozen rows to the existing pedigree renderer."""
    return [
        {**r, "y_dosage": r["y_dose"], "dyz3_bp": r["DYZ3_bp"], "yq12_bp": r["Yq12_bp"],
         "batch": "698_related" if r["batch"] == "698" else r["batch"]}
        for r in rows
    ]


def parse_newick(path):
    text = Path(path).read_text().strip()
    position = 0

    def read_node():
        nonlocal position
        children = []
        if text[position] == "(":
            position += 1
            while True:
                children.append(read_node())
                if text[position] == ",":
                    position += 1
                else:
                    assert text[position] == ")"
                    position += 1
                    break
        start = position
        while position < len(text) and text[position] not in ":,();":
            position += 1
        name = text[start:position].strip().strip("'\"")
        length = 0.0
        if position < len(text) and text[position] == ":":
            position += 1
            start = position
            while position < len(text) and text[position] not in ",();":
                position += 1
            length = float(text[start:position])
        return {"name": name, "length": length, "children": children}

    root = read_node()
    assert text[position:] in (";", ""), text[position:]
    return root


def prune_tree(node, samples):
    if not node["children"]:
        return {**node} if node["name"] in samples else None
    children = [child for original in node["children"]
                if (child := prune_tree(original, samples)) is not None]
    if not children:
        return None
    if len(children) == 1:
        child = children[0]
        return {**child, "length": child["length"] + node["length"]}
    return {**node, "children": children}


def collapse_tree(tree, samples):
    priority = {group: i for i, group in enumerate(ORDER)}

    def visit(node):
        if not node["children"]:
            record = samples[node["name"]]
            return {**node, "members": [record], "extent": 0.0}
        children = [visit(child) for child in node["children"]]
        members = [row for child in children for row in child["members"]]
        extent = max(child["length"] + child["extent"] for child in children)
        groups = {row["display_group"] for row in members}
        if len(groups) == 1:
            return {**node, "name": next(iter(groups)), "children": [], "members": members, "extent": extent}
        children.sort(key=lambda child: min(priority.get(r["display_group"], 1000) for r in child["members"]))
        return {**node, "children": children, "members": members, "extent": extent}

    collapsed = visit(tree)
    leaves = []

    def place(node, distance):
        node["x"] = distance + node["length"]
        if not node["children"]:
            node["group"] = node["members"][0]["display_group"]
            node["y"] = len(leaves)
            leaves.append(node)
        else:
            for child in node["children"]:
                place(child, node["x"])
            node["y"] = np.mean([c["y"] for c in node["children"]])

    collapsed["length"] = 0.0
    place(collapsed, 0.0)
    assert sum(len(n["members"]) for n in leaves) == len(samples)
    return collapsed, leaves


def validation_summary():
    """Summarize the evaluation of the short-read model."""
    from validation_plots import load
    cv_rows, lineage_rows, _, _ = load()
    cv = {key: {rep: [abs(r["signed"]) for r in rows]
                for rep, rows in groups.items()} for key, groups in cv_rows.items()}
    cv_samples = {key: {r["sample"] for rows in groups.values() for r in rows}
                  for key, groups in cv_rows.items()}
    lolo = {key: [abs(r["signed"]) for rows in groups.values() for r in rows]
            for key, groups in lineage_rows.items()}
    summaries = []
    for key, _, name, _, _ in TARGETS:
        medians = [float(np.median(values)) for values in cv[key].values()]
        assert len(medians) == 10
        summaries.append({
            "target": key, "label": name, "cv_n": len(cv_samples[key]), "cv_repeats": len(medians),
            "cv_mean_median_absolute_pct": float(np.mean(medians)),
            "cv_min_median_absolute_pct": min(medians), "cv_max_median_absolute_pct": max(medians),
            "lolo_n": len(lolo[key]), "lolo_median_absolute_pct": float(np.median(lolo[key])),
        })
    depth_rows = []
    depth_predictions, _, _ = load_depth_benchmark()
    for fraction, depth in DEPTHS.items():
        for key, _, name, _, _ in TARGETS:
            values = [np.median([abs(r["signed"]) for r in rows])
                      for rows in depth_predictions[fraction][key].values()]
            assert len(values) == 10, (fraction, key)
            depth_rows.append({"fraction": fraction, "median_autosomal_depth": depth, "target": key,
                               "label": name, "repeats": len(values),
                               "mean_median_absolute_pct": float(np.mean(values)),
                               "min_median_absolute_pct": float(min(values)),
                               "max_median_absolute_pct": float(max(values))})
    return summaries, depth_rows


def snapshot_summary(raw, valid, pedigree, snapshot_path=DEFAULT_SNAPSHOT):
    population = population_rows(valid)
    index = {r["sample"]: r for r in valid}
    pairs = []
    for line in (ROOT / "data/catalog_d376/fs_pairs.tsv").read_text().splitlines():
        father, son = line.split("\t")
        if father in index and son in index:
            pairs.append((father, son))
    subsets = {}
    for label, selected in (("all_valid", population), ("outside_reference_panel", [r for r in population if not r["reference"]])):
        j = [r for r in selected if r["haplogroup"].startswith("J")]
        short = [r for r in j if r["haplogroup"].startswith("J2a1a1")]
        other = [r for r in j if not r["haplogroup"].startswith("J2a1a1")]
        ratios = [r["ratio_bp"] for r in selected]
        subsets[label] = {
            "n": len(selected), "families": len({r["family"] for r in selected}),
            "j_n": len(j), "j2a1a1_n": len(short), "other_j_n": len(other),
            "j2a1a1_families": len({r["family"] for r in short}),
            "other_j_families": len({r["family"] for r in other}),
            "j2a1a1_dyz19_min_bp": min(r["DYZ19_bp"] for r in short),
            "j2a1a1_dyz19_max_bp": max(r["DYZ19_bp"] for r in short),
            "other_j_dyz19_min_bp": min(r["DYZ19_bp"] for r in other),
            "other_j_dyz19_max_bp": max(r["DYZ19_bp"] for r in other),
            "ratio_bp_median": float(np.median(ratios)),
            "ratio_bp_p1": float(np.percentile(ratios, 1)),
            "ratio_bp_p99": float(np.percentile(ratios, 99)),
            "dyz1_dyz2_spearman": float(spearmanr([r["DYZ1_bp"] for r in selected], [r["DYZ2_bp"] for r in selected]).statistic),
        }
    pair_stats = {}
    for group, selected in (
            ("all", pairs),
            ("cross_batch", [(a,b) for a,b in pairs if index[a]["batch"] != index[b]["batch"]]),
            ("same_batch", [(a,b) for a,b in pairs if index[a]["batch"] == index[b]["batch"]])):
        pair_stats[group] = {"n": len(selected)}
        for key in ("DYZ3_bp", "Yq12_bp"):
            delta = np.array([100 * (index[b][key] / index[a][key] - 1) for a,b in selected])
            pair_stats[group][key] = {
                "median_signed_pct": float(np.median(delta)),
                "median_absolute_pct": float(np.median(abs(delta))),
                "p90_absolute_pct": float(np.percentile(abs(delta), 90))}
    lineages = sorted({r["lineage"] for r in valid})
    return {
        "snapshot": str(Path(snapshot_path).resolve().relative_to(ROOT)), "model_status": "short-read model",
        "completed": len(raw), "stream_target": 1639, "completed_fraction": len(raw) / 1639,
        "valid_males": len(valid), "valid_reference_males": sum(r["reference"] for r in valid),
        "population_summary_n": len(population),
        "population_excluded_samples": [r["sample"] for r in valid if not r["population_included"]],
        "all_reference_males": sum(r["reference_male"] == "yes" for r in raw),
        "flag_counts": dict(Counter(r["flag"] for r in raw)),
        "female_max_y_dose": max(float(r["y_dose"]) for r in raw if r["sex_ped"] == "2"),
        "male_below_0_9_including_no_estimate": sum(r["sex_ped"] == "1" and float(r["y_dose"]) < .9 for r in raw),
        "complete_father_son_pairs": len(pairs),
        "pair_batch_directions": dict(Counter(f"{index[a]['batch']}->{index[b]['batch']}" for a,b in pairs)),
        "pairs_involving_reference": sum(index[f]["reference"] or index[s]["reference"] for f, s in pairs),
        "both_nonreference_pairs": sum(not index[f]["reference"] and not index[s]["reference"] for f, s in pairs),
        "batch_counts": dict(Counter(r["batch"] for r in valid)),
        "superpopulation_counts": dict(Counter(r["superpopulation"] for r in valid)),
        "population_count": len({r["population"] for r in valid}),
        "lineage_counts": dict(Counter(r["lineage"] for r in valid)),
        "display_lineage_counts": dict(Counter(r["lineage_display"] for r in valid)),
        "population_lineage_counts": dict(Counter(r["lineage"] for r in population)),
        "median_bp_by_target": {key: float(np.median([r[key] for r in population])) for _,key,*_ in TARGETS},
        "lineage_medians_bp": {g: {key: float(np.median([r[key] for r in population if r["lineage"] == g]))
                                   for _,key,*_ in TARGETS} for g in lineages},
        "pair_stats": pair_stats,
        "subsets": subsets,
    }


def pedigree_relationships(raw, pedigree):
    """Record documented relationships; a shared family ID alone is not a pair call."""
    relationships = []
    for index, first in enumerate(raw):
        a = pedigree.get(first["sample"], {})
        for second in raw[index + 1:]:
            b = pedigree.get(second["sample"], {})
            parent, child, role, relation = "", "", "", ""
            for one, two, record in ((first, second, b), (second, first, a)):
                if one["sample"] == record.get("FatherID"):
                    parent, child, role, relation = one["sample"], two["sample"], "father", "parent_offspring"
                if one["sample"] == record.get("MotherID"):
                    parent, child, role, relation = one["sample"], two["sample"], "mother", "parent_offspring"
            if not relation:
                known_father = a.get("FatherID") not in (None, "0") and a.get("FatherID") == b.get("FatherID")
                known_mother = a.get("MotherID") not in (None, "0") and a.get("MotherID") == b.get("MotherID")
                if known_father and known_mother:
                    relation = "full_siblings"
                elif known_father or known_mother:
                    relation = "shared_recorded_parent"
            if relation:
                relationships.append({
                    "sample1": first["sample"], "sample2": second["sample"],
                    "relationship": relation, "parent": parent, "child": child, "parent_role": role,
                    "family1": a.get("FamilyID", ""), "family2": b.get("FamilyID", ""),
                    "both_valid_males": first["flag"] == second["flag"] == "ok",
                    "sample1_flag": first["flag"], "sample2_flag": second["flag"],
                    "sample1_reference": first["reference_male"], "sample2_reference": second["reference_male"],
                    "evidence": "documented pedigree; compare separate precomputed public KING estimates",
                })
    return relationships


def autosomal_check(raw, valid, pedigree):
    """Reuse public KING reports; no autosomal kinship is estimated here."""
    metadata = ROOT / "data/figure_metadata"
    sample_ids = {r["sample"] for r in raw}
    valid_index = {r["sample"]: r for r in valid}
    grouped = defaultdict(dict)
    for source in ("a", "b"):
        for row in read_tsv(metadata / f"plink2_deg2_hg38_{source}.kin0"):
            first, second = row["#IID1"], row["IID2"]
            if first in sample_ids and second in sample_ids:
                grouped[tuple(sorted((first, second)))][source] = row
    result = []
    for (first, second), sources in sorted(grouped.items()):
        both_valid = first in valid_index and second in valid_index
        record = {
            "sample1": first, "sample2": second,
            "family1": pedigree[first]["FamilyID"], "family2": pedigree[second]["FamilyID"],
            "same_documented_family": pedigree[first]["FamilyID"] == pedigree[second]["FamilyID"],
            "both_valid_males": both_valid,
            "reported_in_both_tables": len(sources) == 2,
            "king_a": sources.get("a", {}).get("KINSHIP", ""),
            "king_b": sources.get("b", {}).get("KINSHIP", ""),
            "ibs0_a": sources.get("a", {}).get("IBS0", ""),
            "ibs0_b": sources.get("b", {}).get("IBS0", ""),
            "evidence": "precomputed PLINK resource KING-robust report (first/second degree)",
        }
        result.append(record)
    # Conservative sensitivity: one sample per component of reported <=2nd-degree links.
    parents = {sample: sample for sample in valid_index}

    def root(sample):
        while parents[sample] != sample:
            parents[sample] = parents[parents[sample]]
            sample = parents[sample]
        return sample
    for row in result:
        if row["both_valid_males"]:
            parents[root(row["sample1"])] = root(row["sample2"])
    components = defaultdict(list)
    families = defaultdict(list)
    for row in valid:
        components[root(row["sample"])].append(row)
        families[row["family"]].append(row)
    choose = lambda groups: [min(population_rows(group), key=lambda row: row["sample"]) for group in groups.values() if population_rows(group)]
    subsets = {"all_valid": population_rows(valid), "outside_reference_panel": [r for r in population_rows(valid) if not r["reference"]],
               "dosage_at_least_0_9": [r for r in population_rows(valid) if r["y_dose"] >= .9],
               "one_per_documented_family": choose(families),
               "one_per_king_component": choose(components)}
    sensitivities = []
    for name, rows in subsets.items():
        short = [r for r in rows if r["haplogroup"].startswith("J2a1a1")]
        other = [r for r in rows if r["haplogroup"].startswith("J") and not r["haplogroup"].startswith("J2a1a1")]
        sensitivities.append({
            "subset": name, "n": len(rows), "j2a1a1_n": len(short), "other_j_n": len(other),
            "j2a1a1_dyz19_min": min(r["DYZ19_bp"] for r in short),
            "j2a1a1_dyz19_max": max(r["DYZ19_bp"] for r in short),
            "other_j_dyz19_min": min(r["DYZ19_bp"] for r in other),
            "other_j_dyz19_max": max(r["DYZ19_bp"] for r in other),
            "ratio_bp_median": float(np.median([r["ratio_bp"] for r in rows])),
            "dyz1_dyz2_spearman": float(spearmanr([r["DYZ1_bp"] for r in rows], [r["DYZ2_bp"] for r in rows]).statistic),
        })
    extra = [r for r in result if r["both_valid_males"] and not r["same_documented_family"]]
    info = {
        "source": "https://www.cog-genomics.org/plink/2.0/resources",
        "valid_reported_pairs_union": sum(r["both_valid_males"] for r in result),
        "cross_family_pairs_consensus": sum(r["reported_in_both_tables"] for r in extra),
        "cross_family_pairs_union": len(extra),
        "family_representatives": len(subsets["one_per_documented_family"]),
        "king_component_representatives": len(subsets["one_per_king_component"]),
        "selection_rule": "lexicographically first sample per group; no satellite phenotype used",
        "scope": "Published first/second-degree reports; absence does not establish no more-distant relatedness.",
    }
    psam = {r["#IID"]: r for r in read_tsv(metadata / "plink2_hg38_corrected.psam")}
    qc = [{**r, "sex_king_corrected_metadata": psam[r["sample"]]["SEX"],
           "sex_metadata_discordance": r["sex_ped"] != psam[r["sample"]]["SEX"]} for r in raw]
    related = defaultdict(list)
    for row in result:
        if row["both_valid_males"]:
            related[row["sample1"]].append(row)
            related[row["sample2"]].append(row)
    chosen = {row["sample"] for row in choose(components)}
    sample_rows = []
    for row in valid:
        sample = row["sample"]
        component = components[root(sample)]
        links = related[sample]
        sample_rows.append({
            "sample": sample, "family": row["family"],
            "recorded_father": pedigree[sample]["FatherID"],
            "recorded_mother": pedigree[sample]["MotherID"],
            "public_3202_metadata_covered": sample in psam,
            "reported_close_relatives_in_valid_snapshot": len(links),
            "cross_family_links_in_both_tables": sum(
                link["reported_in_both_tables"] and not link["same_documented_family"] for link in links),
            "cross_family_links_either_table": sum(not link["same_documented_family"] for link in links),
            "kinship_component_id": min(r["sample"] for r in component),
            "kinship_component_size": len(component),
            "selected_component_representative": sample in chosen,
            "population_included": row["population_included"],
            "y_status": row["y_status"],
            "reference_male": row["reference_male"],
            "scope": "first/second-degree public reports; no inference of absence of distant kinship",
        })
    return result, sensitivities, info, qc, sample_rows


def validation_panels(fig, summaries, depth_rows):
    # Individual predictions are read from the same archived evaluation as the summaries.
    from figure1b_selected import draw
    draw(fig, compact=True)


def draw_collapsed_tree(ax, root, leaves, individual=False, panel="b"):
    maximum = max(n["x"] + n["extent"] for n in leaves)
    label_x = maximum * 1.04

    def draw(node):
        if node["children"]:
            ys = [n["y"] for n in node["children"]]
            ax.plot([node["x"], node["x"]], [min(ys), max(ys)], color="#60717a", linewidth=.7)
            for child in node["children"]:
                ax.plot([node["x"], child["x"]], [child["y"], child["y"]], color="#60717a", linewidth=.7)
                draw(child)
        else:
            group, y = node["group"], node["y"]
            if individual:
                group = node["members"][0]["haplogroup"]
            color = LINEAGE_COLORS.get("NO" if group.startswith("NO") else group[0], "#65747e")
            end = node["x"] + node["extent"]
            if node["extent"] > 1e-8:
                ax.add_patch(Polygon([[node["x"], y], [end, y - .27], [end, y + .27]],
                                     facecolor=color, edgecolor=color, alpha=.35, linewidth=.45))
            else:
                ax.plot(end, y, "o", color=color, markersize=2)
            ax.plot([end, label_x], [y, y], color="#dce3e7", linewidth=.5, linestyle=":")
            label = (f"{node['members'][0]['sample']}   {group[:10]}" if individual else
                     f"{group} ({len(node['members'])})")
            ax.text(label_x, y, label, va="center", fontsize=7.7,
                    color=color, fontweight="bold" if group.startswith("J2a1a1") or group == "R1b" else "normal")
    draw(root)
    ax.set_xlim(-maximum * .025, maximum * (2.4 if individual else 1.50))
    ax.set_ylim(len(leaves) - .4, -.7)
    ax.axis("off")
    ax.set_title("J SNP subtree; sample / call" if individual else "Y SNP tree; clade (n)",
                 loc="left", fontsize=9.6, pad=11)
    ax.text(-.10, 1.035, panel, transform=ax.transAxes, fontsize=15, fontweight="bold")
    # A scale bar makes the retained sequence-divergence lengths interpretable.
    magnitude = 10 ** np.floor(np.log10(maximum / 3))
    scale = min((magnitude * n for n in (1, 2, 5, 10)), key=lambda v: abs(v - maximum / 3))
    transform = ax.get_xaxis_transform()
    bar_y, label_y = (-.012, -.025) if individual else (-.045, -.085)
    ax.plot([0, scale], [bar_y, bar_y], transform=transform, clip_on=False,
            color="#60717a", linewidth=1.1)
    ax.text(scale / 2, label_y, f"{scale:g} substitutions/site",
            transform=transform, ha="center", va="top", fontsize=6.4, color="#60717a")


def frozen_trait(ax, leaves, key, color, scale, limit, ticks, title, panel, units, highlight, show_panel=True):
    rng = np.random.default_rng(2710 + ord(panel))
    statistics = []
    ax._sample_positions = []
    for y, node in enumerate(leaves):
        group, members = node["group"], node["members"]
        if highlight(group):
            ax.axhspan(y - .44, y + .44, facecolor="#e8f1f7" if panel == "c" else "#f8e9ec", zorder=0)
        values = np.array([row[key] / scale for row in members])
        jitter = y + rng.uniform(-.16, .16, len(members))
        references = np.array([row["reference"] for row in members])
        included = np.array([in_summary(row, key) for row in members])
        flagged = np.array([r["sample"] in {"HG01967", "HG02966"} for r in members]) & (key in {"Yq12_bp", "DYZ19_bp"})
        ax.scatter(values[~flagged], jitter[~flagged], s=4.5, color=color, alpha=.34, linewidths=0, zorder=3)
        if flagged.any():
            ax.scatter(values[flagged], jitter[flagged], s=19, marker="D", facecolors="none",
                       edgecolors="#886333", linewidths=.85, zorder=6, clip_on=False)
        q1, median, q3 = np.percentile(values[included], [25, 50, 75]) if included.any() else [np.nan]*3
        ax.plot([q1, q3], [y, y], linewidth=2.8, color=color, zorder=4)
        ax.plot(median, y, "|", color="#243d49", markersize=7, markeredgewidth=1.1, zorder=5)
        ax._sample_positions.extend((row, value, position) for row, value, position in zip(members, values, jitter))
        statistics.append({"group": group, "target": key, "n": len(members),
                           "n_population_summary": int(included.sum()), "n_excluded_arm_profile": int((~included).sum()),
                           "n_reference": int(references.sum()), "median": float(median * scale),
                           "q25": float(q1 * scale), "q75": float(q3 * scale),
                           "units": "relative dosage" if key == "y_dose" else "bp/bp" if key == "ratio_bp" else "units/units" if key == "unit_ratio" else "bp"})
    ax.set_ylim(len(leaves) - .4, -.7)
    ax.set_xlim(0, limit)
    ax.set_xticks(ticks)
    ax.tick_params(labelsize=7)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.grid(axis="x", color="#e7ecef", linewidth=.55)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=9.5, color=color, pad=11)
    if show_panel:
        ax.text(-.09, 1.035, panel, transform=ax.transAxes, fontsize=15, fontweight="bold", color="#182a34")
    ax.set_xlabel(units, fontsize=7.8)
    return statistics


def annotate_population_extremes(ax, key):
    """Label substantial deviations, with an absolute dosage threshold."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.transforms import Bbox

    points = ax._sample_positions
    q1, q3 = np.percentile([r[key] for r, _, _ in points if in_summary(r, key)], [25, 75])
    iqr = q3 - q1
    low, high = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    by_clade = defaultdict(list)
    for row, _, _ in points:
        if in_summary(row, key):
            by_clade[row["display_group"]].append(row[key])
    medians = {group: float(np.median(values)) for group, values in by_clade.items()}

    def substantial(point):
        row = point[0]
        if row["display_group"] not in medians:
            return False
        return abs(row[key] / medians[row["display_group"]] - 1) >= .50

    selected = {}
    if key == "y_dose":
        for point in sorted(points, key=lambda p: -p[0][key]):
            if point[0][key] > 1.30:
                selected[point[0]["sample"]] = (point, "dosage >1.30; descriptive review threshold")
    else:
        for candidates, reverse, reason in [
            ([p for p in points if p[0][key] < low and substantial(p)], False,
             "lower pooled 1.5-IQR tail and >=50% from clade median; up to three smallest"),
            ([p for p in points if p[0][key] > high and substantial(p)], True,
             "upper pooled 1.5-IQR tail and >=50% from clade median; up to three largest"),
        ]:
            for point in sorted(candidates, key=lambda p: p[0][key], reverse=reverse)[:3]:
                selected[point[0]["sample"]] = (point, reason)
        if key in {"Yq12_bp", "DYZ19_bp"}:
            for point in points:
                if point[0]["sample"] in ({"HG01967", "HG02966", "NA18910", "NA18911"} if key == "DYZ19_bp" else {"HG01967", "HG02966"}) and substantial(point):
                    selected.setdefault(point[0]["sample"], (point, "satellite follow-up candidate; >=50% from clade median"))

    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    bounds = ax.get_window_extent(renderer)
    font = FontProperties(size=5.2)
    anchors = [ax.transData.transform((x, y)) for (r, x, y), _ in selected.values()]
    occupied, placed = [], []
    for (row, value, y), reason in selected.values():
        px, py = ax.transData.transform((value, y))
        width, height, descent = renderer.get_text_width_height_descent(row["sample"], font, False)
        width += 4
        height = max(height + descent + 3, 5.2 * ax.figure.dpi / 72 + 2)
        options = []
        for level in range(9):
            for direction in (1, -1):
                for shift in (0, -.4, .4, -.8, .8):
                    cx = np.clip(px + shift * width, bounds.x0 + width/2 + 1, bounds.x1 - width/2 - 1)
                    cy = np.clip(py + direction * (.8 + level) * height, bounds.y0 + height/2 + 1, bounds.y1 - height/2 - 1)
                    box = Bbox.from_bounds(cx-width/2, cy-height/2, width, height)
                    score = 1e6 * sum(box.overlaps(b) for b in occupied) + 1e4 * sum(box.contains(*p) for p in anchors) + (cx-px)**2 + (cy-py)**2
                    options.append((score, cx, cy, box))
        _, cx, cy, box = min(options, key=lambda option: option[0])
        assert not any(box.overlaps(other) for other in occupied), (key, row["sample"])
        occupied.append(box)
        ax.annotate(row["sample"], (value, y), xytext=ax.transAxes.inverted().transform((cx, cy)),
                    textcoords=ax.transAxes, ha="center", va="center", fontsize=5.2, color="#243d49", zorder=8,
                    bbox=dict(boxstyle="square,pad=0.1", facecolor="white", edgecolor="none"),
                    arrowprops=dict(arrowstyle="-", color="#647780", linewidth=.5, shrinkA=2, shrinkB=1))
        placed.append(dict(target=key, sample=row["sample"], value=row[key], reason=reason,
                          lower_fence=low, upper_fence=high, clade_median=medians[row["display_group"]],
                          relative_deviation=abs(row[key] / medians[row["display_group"]] - 1),
                          dosage_label_threshold=1.30 if key == "y_dose" else ""))
    return placed


def biology_panels(fig, height):
    """Fig. 2b-e: DYZ19 and DYZ3 states with replication, read density, DYZ1 and DYZ2 variability."""
    R = ROOT / "results/model_final_2026-11-02"; GREY, INK, MUTED, ORANGE = "#a9b4ba", "#22313a", "#55666f", "#d9622b"
    C = read_tsv(R / "cohort_1kgp.tsv"); Q = read_tsv(R / "hgsg_calls.tsv"); N = {r["item"]: r["value"] for r in read_tsv(R / "biology_tests.tsv")}
    keep = ("normal", "mosaic_loss", "mosaic_gain", "Y_gain", "par2_gain"); inc = [r for r in C if r["sex"] == "male" and r["estimate_basis"] == "ruler" and r["y_state"] in keep]
    box = lambda x0, x1: fig.add_axes([x0, .72 / height, x1 - x0, 3.33 / height])
    def strip(ax, x, v, color, s=4, a=.55, seed=0, horizontal=False, cap=.3):
        j = np.random.default_rng(seed).uniform(-cap, cap, len(v))
        ax.scatter(*((v, x + j) if horizontal else (x + j, v)), s=s, c=color, lw=0, alpha=a, rasterized=True)
        m = np.median(v); (ax.plot([m, m], [x - .4, x + .4], color=INK, lw=1.3) if horizontal else ax.plot([x - .4, x + .4], [m, m], color=INK, lw=1.3))
    def style(ax, letter):
        ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=7, length=2.5); fig.text({"b": .012, "c": .232, "d": .412, "e": .762}[letter], 4.2 / height, letter, fontsize=15, fontweight="bold")
    # short DYZ19 of J2a1a1, by cohort
    ax = box(.05, .2); own = {r[0]: r[4] for r in csv.reader(open(ROOT / "results/model_checks_2026-11-02/j_calls_hgdp_sgdp.tsv"), delimiter="\t") if r[0].startswith(("HGDP", "SGDP", "HG0", "NA"))}
    Jk = [r for r in inc if r["lineage"] == "J"]; Jq = [r for r in Q if r["call"].startswith("J")]; yes = lambda r: own[r["sample"]].startswith("J2a1a1"); kb = lambda R_: np.array([float(r["DYZ19_bp"]) / 1e3 for r in R_])
    G = [("1kGP", kb([r for r in Jk if r["haplogroup"].startswith("J2a1a1")]), kb([r for r in Jk if not r["haplogroup"].startswith("J2a1a1")]))]
    G += [(c, kb([r for r in Jq if r["cohort"] == c and yes(r)]), kb([r for r in Jq if r["cohort"] == c and not yes(r)])) for c in ("HGDP", "SGDP")]
    assert [(len(a), len(b)) for _, a, b in G] == [(42, 44), (16, 29), (6, 12)], [(len(a), len(b)) for _, a, b in G]
    for i, (c, a, b) in enumerate(G):
        strip(ax, 2.6 * i, a, PURPLE, s=6, a=.8, seed=i, cap=.32); strip(ax, 2.6 * i + 1, b, GREY, s=6, a=.8, seed=i + 5, cap=.32); ax.text(2.6 * i + .5, -.13, c, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7)
    ax.set_xticks([x for i in range(3) for x in (2.6 * i, 2.6 * i + 1)], [f"{len(v)}" for _, a, b in G for v in (a, b)], fontsize=6.1); ax.set_xlim(-.7, 6.9); ax.set_ylim(0, 350); ax.set_ylabel("DYZ19 length (kb)", fontsize=8); style(ax, "b")
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=3.5, color=PURPLE, label="J2a1a1"), Line2D([], [], ls="none", marker="o", ms=3.5, color=GREY, label="Other J")], frameon=False, fontsize=6.2, loc="center", bbox_to_anchor=(.5, .52), handletextpad=.1, ncol=2, columnspacing=.8)
    # unit ratio by lineage
    ax = box(.275, .41); U = read_tsv(R / "dyz1_dyz2_lineages.tsv"); A = {r["measure"].split("_")[1]: float(r["median"]) for r in read_tsv(R / "dyz1_dyz2_truth_vs_yruler.tsv") if r["measure"].startswith("lineage_") and "_truth_" in r["measure"]}
    for i, r in enumerate(U):
        c = GREEN if r["lineage"] == "E" else (BLUE if float(r["ci95_high"]) < 1 else MUTED); ax.plot([float(r["ci95_low"]), float(r["ci95_high"])], [i, i], color=c, lw=1.2); ax.plot(float(r["unit_ratio_median"]), i, "o", ms=3.5, color=c)
        if r["lineage"] in A: ax.plot(A[r["lineage"]], i, "D", ms=3.5, mfc="none", mec=ORANGE, mew=1)
    ax.axvline(1, color=MUTED, lw=.8, ls="--"); ax.set_yticks(range(len(U)), [f'{r["lineage"]} ({r["men"]})' for r in U], fontsize=6.1); ax.set_xlim(.76, 1.2); ax.set_xlabel("DYZ1:DYZ2 unit ratio", fontsize=8); style(ax, "c")
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=3.5, color=MUTED, label="Males"), Line2D([], [], ls="none", marker="D", ms=3.5, mfc="none", mec=ORANGE, label="Assemblies")], frameon=False, fontsize=6.2, loc="center right", bbox_to_anchor=(1, .42), borderaxespad=0, handletextpad=.2)
    # DYZ3 states, by cohort
    ax = box(.455, .745); by = {"1kGP": [(r["haplogroup"], float(r["DYZ3_bp"]) / 1e6) for r in inc], "HGDP": [(r["call"], float(r["DYZ3_bp"]) / 1e6) for r in Q if r["cohort"] == "HGDP"], "SGDP": [(r["call"], float(r["DYZ3_bp"]) / 1e6) for r in Q if r["cohort"] == "SGDP"]}
    ticks, labels = [], []
    for g, (p1, p2, name, x0) in enumerate((("R1b", None, "R1b, other males", 0), ("O2a2b1a1", "O2a2b1a2", "O2a2b1a1, O2a2b1a2", 9.4))):
        for i, c in enumerate(by):
            v1 = np.array([v for h, v in by[c] if h.startswith(p1)]); v2 = np.array([v for h, v in by[c] if (h.startswith(p2) if p2 else not h.startswith(p1))]); x = x0 + 2.9 * i
            strip(ax, x, v1, BLUE, s=4, a=.7, seed=i, cap=.32); strip(ax, x + 1.25, v2, GREY, s=4, a=.45 if len(v2) > 100 else .85, seed=i + 5, cap=.32); ticks += [x, x + 1.25]; labels += [f"{len(v1)}", f"{len(v2):,}"]
            ax.text(x + .6, -.13, c, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7)
        ax.text(x0 + 3.5, .985, name, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7, color=INK)
    ax.set_xticks(ticks, labels, fontsize=5.6); ax.set_xlim(-.7, 17.6); ax.set_ylim(0, 3); ax.axvline(8.6, color="#d8dee2", lw=.8); ax.set_ylabel("DYZ3 length (Mb)", fontsize=8); style(ax, "d")
    # DYZ2 varies more than DYZ1
    ax = box(.805, .968); T = {r["assembly"]: r for r in read_tsv(ROOT / "data/method/yq12_truth_hallast.tsv")}; lin = {r["sample"]: r["lineage"] for r in C}
    kin = {r["sample"]: r["kinship_component_id"] for r in read_tsv(ROOT / "results/estimates_1kgp.relatedness_samples.tsv")}; rep = {}
    for r in sorted((r for r in inc if r["y_state"] != "par2_gain"), key=lambda r: r["sample"]): rep.setdefault(kin.get(r["sample"], r["sample"]), r)
    def centred(P):
        a = np.log2([p[0] for p in P]); b = np.log2([p[1] for p in P]); l = np.array([p[2] for p in P])
        for k in set(l): a[l == k] -= a[l == k].mean(); b[l == k] -= b[l == k].mean()
        return b, a
    men = centred([(float(r["DYZ1_bp"]), float(r["DYZ2_bp"]), r["lineage"]) for r in rep.values() if r["reference_male"] != "yes"]); asm = centred([(float(r["dyz1_units"]), float(r["dyz2_units"]), lin[s]) for s, r in T.items() if s in lin and s != "HG01261"])
    assert len(men[0]) == 1196 and len(asm[0]) == 75, (len(men[0]), len(asm[0]))
    ax.scatter(*men, s=4, c=GREY, lw=0, alpha=.6, label="1kGP males (1,196)", rasterized=True); ax.scatter(*asm, s=9, c=ORANGE, lw=0, alpha=.9, label="Assemblies (75)", rasterized=True)
    L = 1.5; sd = lambda k: float(re.search(r"ratio ([0-9.]+)", N[k]).group(1)); ax.plot([-L, L], [-L, L], color=MUTED, lw=.9, ls="--", label="Equal variation"); ax.plot([-L, L], [-L * sd("composition_sd_1kgp"), L * sd("composition_sd_1kgp")], color=INK, lw=1.2, label="SD ratio 0.762")
    ax.set_xlim(-L, L); ax.set_ylim(-1.2, 1.2); ax.set_xlabel("log2 DYZ2 (centred within lineage)", fontsize=8); ax.set_ylabel("log2 DYZ1 (centred within lineage)", fontsize=8); ax.legend(frameon=False, fontsize=6.2, loc="upper left", borderaxespad=0, handletextpad=.4); style(ax, "e")
    ins = ax.inset_axes([.74, .09, .25, .29]); S = [("1kGP", "1kgp"), ("HGDP", "hgdp"), ("SGDP", "sgdp"), ("Assembly", "assemblies")]
    for i, (n, k) in enumerate(S):
        m_ = re.search(r"ratio ([0-9.]+) \(95% ([0-9.]+) to ([0-9.]+)\)", N[f"composition_sd_{k}"]); v, lo, hi = map(float, m_.groups()); ins.plot([lo, hi], [3 - i, 3 - i], color=INK, lw=1); ins.plot(v, 3 - i, "o", ms=3, color=ORANGE if k == "assemblies" else INK)
    ins.axvline(1, color=MUTED, lw=.8, ls="--"); ins.set_yticks(range(4), [n for n, _ in S[::-1]], fontsize=5.4); ins.set_xlim(.55, 1.05); ins.set_xticks([.6, .8, 1]); ins.tick_params(labelsize=5.4, length=2, pad=1); ins.set_xlabel("SD ratio", fontsize=5.4, labelpad=1); ins.spines[["top", "right"]].set_visible(False); ins.patch.set_alpha(.85)


def population_figure(valid, out, snapshot_path, completed):
    samples = {r["sample"]: r for r in valid}
    tree = prune_tree(parse_newick(TREE_SOURCE), set(samples))
    root, leaves = collapse_tree(tree, samples)
    assert all(n == 1 for n in Counter(leaf["group"] for leaf in leaves).values())
    from cohort_lineages import load_combined_samples
    map_samples = load_combined_samples()
    land = json.loads((ROOT / "data/figure_metadata/world_land.geojson").read_text())
    fig = plt.figure(figsize=(10.6, 9.9))
    height = fig.get_figheight()
    grid = fig.add_gridspec(1, 6, left=.055, right=.965, bottom=5.05 / height, top=9.45 / height,
                            width_ratios=[1.85, .85, .85, .85, 1, .85], wspace=.08)
    axes = [fig.add_subplot(grid[0, i]) for i in range(6)]
    draw_collapsed_tree(axes[0], root, leaves, panel="a")
    stats = frozen_trait(axes[1], leaves, "DYZ3_bp", BLUE, 1e6, 3, [0, 1, 2, 3],
                         "DYZ3", "c", "Length (Mb)", lambda g: g == "R1b", show_panel=False)
    stats += frozen_trait(axes[2], leaves, "Yq12_bp", GREEN, 1e6, 60, [0, 20, 40, 60],
                          "Yq12", "d", "Length (Mb)", lambda g: False, show_panel=False)
    stats += frozen_trait(axes[3], leaves, "DYZ19_bp", PURPLE, 1e3, 550, [0, 250, 500],
                          "DYZ19", "e", "Length (kb)", lambda g: g.startswith("J2a1a1"), show_panel=False)
    for r in valid: r["unit_ratio"] = r["ratio_bp"] * 2404.4 / 3478.9   # bp per Hallast 2026 unit
    ratio_limit = max(2.5, .5 * np.ceil(max(r["unit_ratio"] for r in valid) / .5))
    stats += frozen_trait(axes[4], leaves, "unit_ratio", "#8c7139", 1, ratio_limit, np.arange(0, ratio_limit + .01, 1),
                          "DYZ1:DYZ2", "f", "Unit ratio", lambda g: False, show_panel=False)
    axes[4].axvline(1, color="#8c7139", linestyle="--", linewidth=.8, zorder=1)
    dosage_limit = max(1.5, .25 * np.ceil(max(r["y_dose"] for r in valid) / .25))
    stats += frozen_trait(axes[5], leaves, "y_dose", "#526774", 1, dosage_limit, [.5, 1, 1.5, 2],
                          "Y dosage", "g", "Relative dosage", lambda g: False, show_panel=False)
    axes[5].set_xlim(min(.85, .25 * np.floor(min(r["y_dose"] for r in valid) / .25)), dosage_limit)
    axes[5].axvline(1, color="#526774", linestyle="--", linewidth=.8, zorder=1)
    handles = [
        Line2D([], [], linestyle="none", marker="o", color="#637f8d", markersize=4,
               label="Individual"),
        Line2D([], [], linestyle="none", marker="D", markerfacecolor="none", markeredgecolor="#886333",
               markersize=4, label="Yq satellite depletion"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.8, 4.62 / height), ncol=2,
               frameon=False, fontsize=7.4, handletextpad=.4)
    annotations = []
    for axis, key in zip(axes[1:], ["DYZ3_bp", "Yq12_bp", "DYZ19_bp", "unit_ratio", "y_dose"]):
        annotations.extend(annotate_population_extremes(axis, key))
    write_rows(out / "figure2_labelled_samples.tsv", annotations)
    biology_panels(fig, height)
    save(fig, out, "figure2")
    write_rows(out / "figure2_tree_summary.tsv", stats)
    membership = [{"sample": r["sample"], "haplogroup_2019": r["haplogroup"],
                   "display_group": r["display_group"], "family": r["family"],
                   "reference_male": r["reference_male"],
                   "y_status": r["y_status"], "population_included": r["population_included"],
                   "source": str(Path(snapshot_path).resolve().relative_to(ROOT))}
                  for node in leaves for r in node["members"]]
    write_rows(out / "figure2_tree_membership.tsv", membership)
    counts = Counter(r["lineage"] for r in map_samples)
    write_rows(out / "figure2_lineage_counts.tsv", [
        {"lineage": group, "n_unique_individuals": n, "percent": 100*n/len(map_samples)}
        for group,n in counts.items()])
    pie_fig, pie_ax = plt.subplots(figsize=(12.5, 5.7))
    pie_map(pie_ax, None, map_samples, land, size_legend=True)
    pie_ax.set_title(f"Geographic coverage of 1kGP, HGDP and SGDP  (n={len(map_samples):,} individuals)",
                     fontsize=11, pad=10)
    pie_fig.subplots_adjust(left=.025, right=.975, top=.95, bottom=.16)
    pie_fig.canvas.draw()
    anchor=pie_fig.transFigure.inverted().transform(pie_ax.transData.transform((55.38,-35.63)))
    legend_font=7.6*(pie_ax.get_position().width*pie_fig.get_figwidth())/(.91*10.6)
    map_legend=lineage_legend(pie_fig, map_samples, anchor[1], x=anchor[0], fontsize=legend_font)
    align_map_legend_titles(pie_ax,map_legend)
    pie_fig.text(.5, .047,
                 "Pie area scales with unique individuals per reported location; sectors show broad Y lineages.",
                 ha="center", fontsize=9, color="#4d606b")
    pie_fig.text(.5, .018,
                 "Shared individuals count once. 1kGP contributes length profiles; HGDP and SGDP contribute cohort metadata.",
                 ha="center", fontsize=9, color="#4d606b")
    save(pie_fig, out, "cohort_world_map")
    return root, leaves


def population_supplement(valid, out):
    order = sorted({r["lineage_display"] for r in valid},
                   key=lambda g: (-sum(r["lineage_display"] == g for r in valid), g))
    fig, axes = plt.subplots(3, 3, figsize=(12.6, 9.4), gridspec_kw={"hspace": .44, "wspace": .32})
    rng = np.random.default_rng(271)
    for panel, (ax, (_, key, label, scale, units)) in enumerate(zip(axes.flat, TARGETS)):
        for x, group in enumerate(order):
            selected = [r for r in valid if r["lineage_display"] == group]
            values = np.array([r[key] / scale for r in selected])
            jitter = x + rng.uniform(-.23, .23, len(selected))
            reference = np.array([r["reference"] for r in selected])
            included = np.array([in_summary(r, key) for r in selected])
            ax.scatter(jitter[~reference & included], values[~reference & included], s=5, color=BLUE, alpha=.35, linewidths=0)
            ax.scatter(jitter[reference & included], values[reference & included], s=10, facecolors="none", edgecolors=BLUE, linewidths=.45, alpha=.75)
            ax.scatter(jitter[~included], values[~included], s=19, marker="D", facecolors="none", edgecolors="#886333", linewidths=.8)
            if included.any():
                ax.plot([x - .27, x + .27], [np.median(values[included])]*2, color="#263d48", linewidth=1.6)
        ax.set_xticks(range(len(order)), order, rotation=90, fontsize=7)
        ax.set_ylabel(units, fontsize=8)
        ax.set_title(label, loc="left", fontsize=10)
        ax.text(-.12, 1.045, chr(ord("a") + panel), transform=ax.transAxes,
                fontsize=12, fontweight="bold")
        ax.grid(axis="y", color="#e8edef", linewidth=.5)
    for ax in axes.flat[len(TARGETS) + 1:]:
        ax.set_visible(False)   # unused grid cells
    ax = axes.flat[len(TARGETS)]
    for group in order:
        selected = [r for r in valid if r["lineage_display"] == group]
        for reference in (False, True):
            s = [r for r in population_rows(selected) if r["reference"] == reference]
            color = LINEAGE_COLORS.get(group, LINEAGE_COLORS.get(group[0], BLUE))
            ax.scatter([r["DYZ1_bp"]/1e6 for r in s], [r["DYZ2_bp"]/1e6 for r in s],
                       s=8, facecolors="none" if reference else color,
                       edgecolors=color if reference else "none", linewidths=.5, alpha=.65)
    flagged = [r for r in valid if not r["population_included"]]
    ax.scatter([r["DYZ1_bp"]/1e6 for r in flagged], [r["DYZ2_bp"]/1e6 for r in flagged],
               s=22, marker="D", facecolors="none", edgecolors="#886333", linewidths=.8)
    rho = spearmanr([r["DYZ1_bp"] for r in population_rows(valid)], [r["DYZ2_bp"] for r in population_rows(valid)]).statistic
    ax.set_title(f"DYZ1 versus DYZ2 (Spearman {rho:.2f})", loc="left", fontsize=9.5)
    ax.text(-.12, 1.045, chr(ord("a") + len(TARGETS)), transform=ax.transAxes, fontsize=12, fontweight="bold")
    ax.set_xlabel("DYZ1 (Mb)", fontsize=8)
    ax.set_ylabel("DYZ2 (Mb)", fontsize=8)
    ax.grid(color="#e8edef", linewidth=.5)
    for ax in axes.flat:
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", labelsize=7)
    fig.subplots_adjust(left=.06, right=.98, top=.92, bottom=.135)
    fig.suptitle(f"{ {7: 'Seven', 8: 'Eight'}[len(TARGETS)]} Y satellite measurements in {len(valid)} 1000 Genomes males", fontsize=11, color="#243d49")
    counts = [f"{g}: {sum(r['lineage_display'] == g for r in valid)}" for g in order]
    fig.text(.5, .058, "n per lineage — " + ";  ".join(counts[:10]), ha="center", fontsize=8.8)
    fig.text(.5, .035, ";  ".join(counts[10:]), ha="center", fontsize=8.8)
    fig.text(.5, .010, "Filled points: outside reference panel; open points: reference individuals. Lines: medians.",
             ha="center", fontsize=8)
    save(fig, out, "figureS3_population")


def depth_supplement(depth, out):
    from validation_plots import COLORS
    from validation_ratio import KEY,LABEL,COLOR
    targets=list(TARGETS);colors=list(COLORS)
    if any(r['target']==KEY for r in depth):
        targets.append((KEY,None,LABEL,1,''));colors.append(COLOR)
    if len(targets)==9:fig,axes=plt.subplots(3,3,figsize=(11.8,8.4),sharex=True)
    else:fig,axes=plt.subplots(2,4,figsize=(11.8,5.9),sharex=True)
    for ax, (key, _, label, _, _), color in zip(axes.flat, targets, colors):
        rows = sorted((r for r in depth if r["target"] == key), key=lambda r:r["median_autosomal_depth"])
        ax.plot([r["median_autosomal_depth"] for r in rows], [r["mean_median_absolute_pct"] for r in rows],
                "o-", color=color, markersize=3.5, linewidth=1)
        ax.fill_between([r["median_autosomal_depth"] for r in rows],
                        [r["min_median_absolute_pct"] for r in rows], [r["max_median_absolute_pct"] for r in rows],
                        color=color, alpha=.15, linewidth=0)
        ax.set_xscale("log")
        ax.set_xticks([1, 2, 5, 10, 20, 34], ["1","2","5","10","20","34"])
        ax.set_title(label, loc="left", fontsize=10)
        ax.grid(color="#e8edef", linewidth=.5)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_xlabel("Median autosomal coverage (×)", fontsize=8)
        ax.set_ylabel("Absolute error (%)", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes.flat[len(TARGETS):]: ax.set_visible(False)
    fig.subplots_adjust(wspace=.4, hspace=.45, top=.88, bottom=.13)
    fig.suptitle("Sequencing depth and prediction accuracy", fontsize=11, color="#243d49")
    fig.text(.5, .025, "Points: mean of ten repeat-specific medians. Shading: range across repeats. Training reads remain at full depth.",
             ha="center", fontsize=8)
    save(fig, out, "figureS4_depth")


def j_subtree_figure(valid, out):
    members = [dict(r, display_group=r["sample"]) for r in valid if r["haplogroup"].startswith("J")]
    index = {r["sample"]: r for r in members}
    tree = prune_tree(parse_newick(TREE_SOURCE), set(index))
    root, leaves = collapse_tree(tree, index)
    fig = plt.figure(figsize=(8.2, 10.8))
    grid = fig.add_gridspec(1, 2, left=.08, right=.96, bottom=.09, top=.91, width_ratios=[3.4,1.2], wspace=.1)
    tree_ax, ax = [fig.add_subplot(grid[0, i]) for i in range(2)]
    draw_collapsed_tree(tree_ax, root, leaves, individual=True, panel="a")
    for node in leaves:
        row=node["members"][0];y=node["y"]
        short = row["haplogroup"].startswith("J2a1a1")
        ax.plot(row["DYZ19_bp"]/1e3, y, "o", markersize=4,
                markerfacecolor="none" if row["reference"] else "#bd6078" if short else PURPLE,
                markeredgecolor="#bd6078" if short else PURPLE, markeredgewidth=.7)
    ax.set_ylim(len(leaves)-.4,-.7)
    ax.set_xlim(0,350)
    ax.set_xticks([0,100,200,300])
    ax.set_yticks([])
    ax.set_xlabel("DYZ19 (kb)",fontsize=9)
    ax.set_title("DYZ19 state",fontsize=10,pad=11)
    ax.text(-.10,1.035,"b",transform=ax.transAxes,fontsize=15,fontweight="bold")
    ax.spines[["right","left","top"]].set_visible(False)
    ax.grid(axis="x",color="#e7ecef",linewidth=.6)
    ax.set_axisbelow(True)
    fig.suptitle(f"DYZ19 variation across {len(members)} J chromosomes",fontsize=11,color="#243d49",y=.97)
    fig.text(.5,.04,"Pink: J2a1a1 calls; purple: other J calls. Open points: reference individuals. Labels truncate haplogroup calls.",
             ha="center",fontsize=7.4)
    fig.text(.5,.02,"The tree is inferred from SNPs; labels give haplogroup assignments. Branch lengths represent sequence divergence.",
             ha="center",fontsize=7.2)
    save(fig,out,"figureS5_J_DYZ19")
    write_rows(out/"figureS5_J_membership.tsv",[
        {"sample":r["sample"],"haplogroup":r["haplogroup"],"family":r["family"],
         "reference_male":r["reference_male"],"DYZ19_bp":r["DYZ19_bp"]}
        for n in leaves for r in n["members"]])

def assembly_supplement(out):
    """Preserve the complete-array annotation view, separate from population reads."""
    calls = {r["ID"]: r["Haplogroup"] for r in read_tsv(ROOT / "data/method/haplogroups_1kgp_isogg2019_decode1.tsv")}
    parents, children = tree_data(ROOT / "data/figure_metadata/isogg2019_tree.tsv")
    groups = defaultdict(list)
    for row in load_snapshot()[1]:   # display groups from the current 1kGP estimates (dose >= 0.9)
        if row["y_dose"] >= .9:
            groups[clade_of(calls[row["sample"]], parents)].append(row)
    _, order = pruned_tree(groups, parents, children)
    rows = load_composition(ROOT / "data/method/yq12_truth_hallast.tsv", ROOT / "data/markers_u/truth_regions.tsv",
                            calls, parents, groups)
    composition_supplement(rows, order, out)
    write_rows(out / "figureS2_composition_samples.tsv", rows)


def snapshot_manifest(path, out):
    """Hash the exact dated inputs, figure scripts and generated tables."""
    tag = snapshot_tag(path)
    inputs = [Path(path), TREE_SOURCE]
    inputs += [ROOT / "results/model/lineage_out.tsv"]
    inputs += list(DEPTH_SOURCE.glob("*.tsv"))
    inputs = [p for p in inputs if p.is_file()]
    inputs += [ROOT / name for name in (
        "data/catalog_d376/ped.txt", "data/catalog_d376/fs_pairs.tsv",
        "data/method/dyz19_clean.tsv",
        "data/method/indep_setB.tsv", "data/markers_u/truth_regions.tsv",
        "data/method/haplogroups_1kgp_isogg2019_decode1.tsv",
        "scripts/manuscript_figures.py", "scripts/figure2_tree.py",
        "scripts/figures_1kgp.py", "scripts/prepare_cohort_pies.py",
        "scripts/cohort_lineages.py", "results/unassigned_lineages_calls_2026-10-28.tsv",
        "scripts/figure1b_selected.py", "scripts/validation_plots.py", "scripts/depth_benchmark.py",
        "scripts/pedigree_plots.py", "scripts/validation_ratio.py")]
    inputs += [p for p in (ROOT / "data/figure_metadata").iterdir()
               if p.suffix in {".tsv", ".kin0", ".psam", ".xlsx", ".geojson"} or p.name == "sgdp_metadata.txt"]
    inputs += list((ROOT / "results/cohort_metadata_2026-10-28").glob("*.tsv"))
    outputs = [p for p in out.glob("*.tsv") if not p.name.startswith("snapshot_")]
    outputs += [out / f"snapshot_{tag}_samples.tsv"]
    outputs += list((ROOT / "results").glob(Path(path).stem + ".*.tsv"))
    outputs += [out / f"snapshot_{tag}_summary.json"]
    summary = json.loads(outputs[-1].read_text())
    def record(p):
        return {"path": str(p.absolute().relative_to(ROOT)), "bytes": p.stat().st_size,
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
    manifest = {
        "snapshot": str(Path(path).absolute().relative_to(ROOT)), "date": "2026-10-28",
        "snapshot_version": tag,
        "scope": f"short-read model; {summary['completed']} completed 1kGP records; analysed 1kGP plus HGDP/SGDP metadata map.",
        "validation_arm": "short-read model; per-depth CV predictions and lineage-out predictions in results/model (depth/, lineage_out.tsv)",
        "model_archive_scope": "Final parameter TSV and raw evaluation predictions; binary model is not copied here.",
        "relatedness": "Reused public first/second-degree reports; no new autosomal kinship estimation.",
        "public_king_sources": [
            "https://www.cog-genomics.org/plink/2.0/resources",
            "https://www.dropbox.com/s/xry2angmvmk629j/deg2_hg38.kin0?dl=1",
            "https://www.dropbox.com/s/elf7pk7siurlyak/deg2_hg38.kin0?dl=1",
            "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"],
        "inputs": [record(p) for p in sorted(set(inputs))],
        "derived_outputs": [record(p) for p in sorted(set(outputs))],
    }
    (out / f"snapshot_{tag}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def render_snapshot(path, raw, valid, pedigree, benchmarks, depth, out):
    tag = snapshot_tag(path)
    summary = snapshot_summary(raw, valid, pedigree, path)
    kinship, sensitivities, info, qc, sample_rows = autosomal_check(raw, valid, pedigree)
    summary["relatedness"] = info
    root, leaves = population_figure(valid, out, path, len(raw))
    population_supplement(valid, out)
    from validation_plots import fold_supplement, write_top_samples, write_depth_individuals, lineage_supplement
    fold_supplement(out)
    lineage_supplement(out)
    write_top_samples(out)
    write_depth_individuals(out)
    depth_supplement(depth, out)
    j_subtree_figure(valid, out)
    assembly_supplement(out)
    write_rows(out / "figure1_validation_summary.tsv", benchmarks)
    write_rows(out / "figure1_depth_summary.tsv", depth)
    write_rows(out / f"snapshot_{tag}_samples.tsv", [
        {**{k:v for k,v in r.items() if k not in {"reference"}}, "reference": r["reference"]} for r in valid])
    write_rows(Path(path).with_suffix(".relationships.tsv"), pedigree_relationships(raw, pedigree))
    write_rows(Path(path).with_suffix(".kinship.tsv"), kinship)
    write_rows(Path(path).with_suffix(".sensitivity.tsv"), sensitivities)
    write_rows(Path(path).with_suffix(".sample_qc.tsv"), qc)
    write_rows(Path(path).with_suffix(".relatedness_samples.tsv"), sample_rows)
    write_rows(Path(path).with_suffix(".population_inclusion.tsv"), [
        {k: r[k] for k in ("sample", "y_status", "y_status_source", "population_included")}
        for r in valid])
    follow = []
    for sample, reason in (
        ("HG01967","Short DYZ19 with elevated Y dosage; independent technical/structural follow-up"),
        ("HG02966","Near-absent Yq12/DYZ19 marker signal; estimates below reference range; structural confirmation pending"),
        ("HG01527","Father of HG01529; high DYZ18S is one familial observation"),
        ("HG01529","Son of HG01527; high DYZ18S is one familial observation"),
        ("HG02965","Reference father of HG02966; compare Yq marker depletion in son"),
        ("HG02804","High DYZ1:DYZ2 base-pair ratio; low Yq12 outside reference range; daughter HG02806 (pedigree sex 2); no paternal Y transmission comparison"),
        ("HG02684","Father of HG02686; long DYZ3 beyond reference range and high base-pair ratio require independent checks"),
        ("HG02686","Son of HG02684; familial long DYZ3/high base-pair ratio, not a validated unit-balance exception"),
        ("HG04006","Largest current DYZ3 estimate; beyond calibration range, requires independent measurement")):
        row = next(r for r in valid if r["sample"] == sample)
        follow.append({k:row[k] for k in ("sample","family","haplogroup","reference_male","y_dose","DYZ3_bp","Yq12_bp","DYZ1_bp","DYZ2_bp","DYZ18S_bp","DYZ19_bp","ratio_bp")}|{"follow_up":reason})
    write_rows(Path(path).with_suffix(".follow_up.tsv"),follow)
    (out / f"snapshot_{tag}_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    snapshot_manifest(path, out)
    print("Frozen population:", len(valid), "valid;", summary["valid_reference_males"], "reference;", len(leaves), "SNP-tree rows.")
    return summary


def main():
    raw, valid, pedigree = load_snapshot()
    tree = prune_tree(parse_newick(TREE_SOURCE), {r["sample"] for r in valid})
    collapsed, leaves = collapse_tree(tree, {r["sample"]: r for r in valid})
    print("Collapsed SNP-tree rows:", [(n["group"], len(n["members"])) for n in leaves])
    print("Repeated group labels:", {k: v for k, v in Counter(n["group"] for n in leaves).items() if v > 1})
    print("Summary:", json.dumps(snapshot_summary(raw, valid, pedigree), indent=2))
    print("Validation:", validation_summary()[0])


if __name__ == "__main__":
    main()
