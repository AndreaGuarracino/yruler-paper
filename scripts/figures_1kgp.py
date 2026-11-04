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
from matplotlib.patches import Patch, Polygon
from matplotlib.colors import to_rgba
import numpy as np
from scipy.stats import spearmanr

from palette import RATIO, FS, MM, letter, ORDER, REF, REF_HANDLE, TITLE_PAD
from figure2_tree import (
    BLUE, GREEN, PURPLE, LINEAGE_COLORS, read_tsv, save,
    tree_data, pie_map, lineage_legend, align_map_legend_titles, clade_of, pruned_tree,
    load_composition, composition_supplement, add_reference_calls,
)

ROOT = Path(__file__).resolve().parents[1]
AUX = ROOT / "results/figures"   # figure tables, not figures
DEFAULT_SNAPSHOT = ROOT / "results/estimates_1kgp.tsv"
TREE_SOURCE = ROOT / "results/tree/tree.rooted.nwk"
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
TARGETS.sort(key=lambda t: ORDER.index(t[2]))   # one class order
CUT = {
    "A00", "A0", "A1a", "A1b1", "B", "C", "D", "E", "F1", "F2", "F3", "F4",
    "G", "H", "I", "J1", "J2a1a1a", "J2a1a1b", "J2a1a2", "J2a2", "J2b", "L", "T",
    "N", "O1a", "O1b", "O2a", "O2a1", "O2a2a", "O2a2b1a1", "O2a2b1a2", "O2a2b2", "O2b", "Q", "R1a", "R1b", "R2", "M", "S",
}
ORDER = [
    "A00", "A0", "A1a", "A1b1", "B", "D", "E", "C", "F1", "F2", "F3", "F4",
    "G", "H", "I", "J1", "J2a1*", "J2a1a1a", "J2a1a1b", "J2a1a2", "J2a2", "J2b",
    "L", "T", "NO*", "N", "O1a", "O1b", "O2a", "O2a1", "O2a2a", "O2a2b1a1", "O2a2b1a2", "O2a2b2", "O2b", "M", "S", "Q",
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


def draw_collapsed_tree(ax, root, leaves, individual=False, panel="b", room=1.55, counts=None):
    maximum = max(n["x"] + n["extent"] for n in leaves)
    label_x = maximum * 1.04

    def draw(node):
        if node["children"]:
            ys = [n["y"] for n in node["children"]]
            ax.plot([node["x"], node["x"]], [min(ys), max(ys)], color=LINE, linewidth=.7)
            for child in node["children"]:
                ax.plot([node["x"], child["x"]], [child["y"], child["y"]], color=LINE, linewidth=.7)
                draw(child)
        else:
            group, y = node["group"], node["y"]
            if individual:
                group = node["members"][0]["haplogroup"]
            color = LINEAGE_COLORS.get("NO" if group.startswith("NO") else group[0], LINE)
            end = node["x"] + node["extent"]
            if node["extent"] > 1e-8:
                ax.add_patch(Polygon([[node["x"], y], [end, y - .27], [end, y + .27]],
                                     facecolor=color, edgecolor=color, alpha=.35, linewidth=.3))
            else:
                ax.plot(end, y, "o", color=color, markersize=1.5)
            ax.plot([end, label_x], [y, y], color=GRID, linewidth=.5, linestyle=":")
            label = (f"{node['members'][0]['sample']}   {group[:10]}" if individual else
                     f"{group} ({counts[group] if counts else len(node['members'])})")
            ax.text(label_x, y, label, va="center", fontsize=FS['point'],
                    color=color, fontweight="bold" if group.startswith("J2a1a1") or group in ("R1b", "O2a2b1a1") else "normal")
    draw(root)
    ax.set_xlim(-maximum * .025, maximum * (2.4 if individual else room))
    ax.set_ylim(len(leaves) - .4, -.7)
    ax.axis("off")
    # A scale bar makes the retained sequence-divergence lengths interpretable.
    magnitude = 10 ** np.floor(np.log10(maximum / 3))
    scale = min((magnitude * n for n in (1, 2, 5, 10)), key=lambda v: abs(v - maximum / 3))
    transform = ax.get_xaxis_transform()
    bar_y = -.012 if individual else 1.025   # Fig. 2: header row
    x0 = 0 if individual else .12 * maximum * (room + .025)   # clear of the panel letter
    ax.plot([x0, x0 + scale], [bar_y, bar_y], transform=transform, clip_on=False,
            color=LINE, linewidth=.75)
    ax.text(*((0, -.025) if individual else (x0 + scale * 1.15, bar_y)), f"{scale:g} substitutions/site",
            transform=transform, ha="left", va="top" if individual else "center", fontsize=FS["point"], color=LINE)


def truth_males():
    """Males in each class truth set after QC."""
    out = defaultdict(set); drop = {tuple(l.split()[:2]) for l in open(ROOT / "data/method/truth_qc/x/rule_w4.txt") if l.strip()}
    for r in read_tsv(ROOT / "data/markers_u/truth_regions.tsv"):
        m = r["male"].removesuffix("_R2")
        if (r["class"], m) not in drop: out[r["class"]].add(m)
    for c in ("dyz19", "dyz17"): out[c] = {r["assembly"].removesuffix("_R2") for r in read_tsv(ROOT / f"data/method/{c}_clean.tsv") if r["clean"] == "True"}
    for c in ("dyz1", "dyz2", "dyz18s"): out[c] = {m for m in out["yq12"] if (c, m) not in drop}   # parts of the Yq12 assemblies
    print("truth males per class:", {k: len(v) for k, v in out.items()}, flush=True)
    return out


TRUTH = truth_males()
TRUTH_OF = {"DYZ3_bp": "dyz3", "DYZ17_bp": "dyz17", "DYZ19_bp": "dyz19", "Yq12_bp": "yq12", "DYZ18S_bp": "dyz18s", "DYZ1_bp": "dyz1", "DYZ2_bp": "dyz2", "unit_ratio": "yq12"}   # ratio: Yq12 assemblies


def frozen_trait(ax, leaves, key, color, scale, limit, ticks, title, panel, units, highlight, show_panel=True):
    rng = np.random.default_rng(2710 + ord(panel))
    statistics = []
    ax._sample_positions = []
    for y, node in enumerate(leaves):
        group, members = node["group"], node["members"]
        if highlight(group):
            ax.axhspan(y - .44, y + .44, facecolor=BAND["DYZ3"] if panel == "c" else BAND["DYZ19"], zorder=0)
        values = np.array([row[key] / scale for row in members])
        jitter = y + rng.uniform(-.16, .16, len(members))
        references = np.array([row["reference"] for row in members])
        included = np.array([in_summary(row, key) for row in members])
        truth = np.array([row["sample"] in TRUTH[TRUTH_OF[key]] for row in members])
        ax.scatter(values[~truth], jitter[~truth], s=1.5, color=color, alpha=.34, linewidths=0, zorder=3)
        ax.scatter(values[truth], jitter[truth], edgecolors=color, zorder=3, **REF)   # truth-set males
        ax._triangles = getattr(ax, "_triangles", 0) + int(truth.sum())
        q1, median, q3 = np.percentile(values[included], [25, 50, 75]) if included.any() else [np.nan]*3
        ax.barh(y, q3 - q1, left=q1, height=.3, color=color, linewidth=0, zorder=4)
        ax.plot(median, y, "|", color=INK, markersize=5, markeredgewidth=1, zorder=5)
        ax._sample_positions.extend((row, value, position) for row, value, position in zip(members, values, jitter))
        statistics.append({"group": group, "target": key, "n": len(members),
                           "n_population_summary": int(included.sum()), "n_excluded_arm_profile": int((~included).sum()),
                           "n_reference": int(references.sum()), "median": float(median * scale),
                           "q25": float(q1 * scale), "q75": float(q3 * scale),
                           "units": "relative dosage" if key == "y_dose" else "bp/bp" if key == "ratio_bp" else "units/units" if key == "unit_ratio" else "bp"})
    ax.set_ylim(len(leaves) - .4, -.7)
    column_axis(ax, title, color, units, limit, ticks)
    return statistics


def reads_legend(fig, x=.5, loc="lower center"):
    """Legend for dots and open diamonds, every figure."""
    from palette import REF_HANDLE, top_legend
    top_legend(fig, x, loc, handles=[Line2D([], [], ls="none", marker="o", ms=2.5, color=MUTED, label="From reads"), Line2D([], [], markeredgecolor=MUTED, label="From reads, truth-set male", **REF_HANDLE)], ncol=2)


def column_axis(ax, title, color, units, limit, ticks):
    """One style for every trait column."""
    ax.set_xlim(-.03 * limit, 1.06 * limit)   # end tick labels clear of neighbours
    ax.set_xticks(ticks); ax.set_xticks(np.arange(0, limit + 1e-9, (ticks[1] - ticks[0]) / 2), minor=True)   # minor: half steps
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    from palette import head
    head(ax, title)
    ax.set_xlabel(units)


from palette import GREY, INK, MUTED, ORANGE, LINE, GRID, BAND, J2
TITLE_Y = 5.0   # mm below the axes: axis titles, cohort row


def biology_inputs():
    """1kGP males on the ruler, HGDP/SGDP calls, test results."""
    R = ROOT / "results/model"
    C = read_tsv(R / "cohort_1kgp.tsv"); Q = read_tsv(R / "hgsg_calls.tsv"); N = {r["item"]: r["value"] for r in read_tsv(R / "biology_tests.tsv")}
    keep = ("normal", "mosaic_loss", "mosaic_gain", "Y_gain", "par2_gain"); inc = [r for r in C if r["sex"] == "male" and r["estimate_basis"] == "ruler" and r["y_state"] in keep]
    return R, C, Q, N, inc


def strip(ax, x, v, color, s=4, a=.55, seed=0, cap=.3, ref=None):
    j = np.random.default_rng(seed).uniform(-cap, cap, len(v)); ref = np.zeros(len(v), bool) if ref is None else np.asarray(ref, bool)
    ax.scatter(x + j[~ref], v[~ref], s=s, c=color, lw=0, alpha=a, rasterized=True)
    ax.scatter(x + j[ref], v[ref], edgecolors=color, rasterized=True, **REF)   # truth-set males
    m = np.median(v); ax.plot([x - .4, x + .4], [m, m], color=INK, lw=1)



SUBTREES = {"P": ("Q", "R1a", "R1b", "R2"), "O2a": ("O2a1a", "O2a1b", "O2a2a", "O2a2b1a1", "O2a2b1a2", "O2a2b2")}


def subtree_rows(valid, clade):
    """Males of one clade, grouped by its subclades."""
    sub = SUBTREES[clade]; group = (lambda r: r["display_group"]) if clade == "P" else (lambda r: next((p for p in sub if r["haplogroup"].startswith(p)), None))
    return {r["sample"]: dict(r, display_group=group(r)) for r in valid if group(r) in sub}


COHORT_SHADE = {"1kGP": BLUE, "HGDP": BLUE, "SGDP": BLUE}   # cohorts by position, top to bottom


def subtree_panel(tax, ax, rows, other, highlight, limit, letter):
    """Fig. 2b or c: one subtree, DYZ3 per cohort."""
    from validation_plots import error_box_stats, draw_error_boxes
    from matplotlib.colors import to_rgba
    root, leaves = collapse_tree(prune_tree(parse_newick(TREE_SOURCE), set(rows)), rows)
    assert all(n == 1 for n in Counter(leaf["group"] for leaf in leaves).values()), [leaf["group"] for leaf in leaves]
    data = {n["group"]: {"1kGP": [(r["DYZ3_bp"] / 1e6, r["sample"] in TRUTH["dyz3"]) for r in n["members"] if in_summary(r, "DYZ3_bp")],   # population summary rule
                         **{c: [(float(r["DYZ3_bp"]) / 1e6, False) for r in other if r["cohort"] == c and r["call"].startswith(n["group"])] for c in ("HGDP", "SGDP")}} for n in leaves}
    draw_collapsed_tree(tax, root, leaves, panel=letter, room=3.3, counts={g: "/".join(str(len(d[c])) for c in COHORT_SHADE) for g, d in data.items()})
    stats, rng = [], np.random.default_rng(2710)
    for y, node in enumerate(leaves):
        g = node["group"]
        if g == highlight: ax.axhspan(y - .45, y + .45, facecolor=BAND["DYZ3"], zorder=0)
        for k, (c, col) in enumerate(COHORT_SHADE.items()):
            v = np.array([x for x, _ in data[g][c]]); t = np.array([f for _, f in data[g][c]], bool); yy = y + (k - 1) * .3
            if y == 0: ax.text(limit * .99, yy, c, ha="right", va="center", fontsize=FS["point"], color=INK)   # cohort positions, first row
            if not len(v): continue
            art = draw_error_boxes(ax, [error_box_stats([{"signed": x} for x in v])], [yy], col, .2, showfliers=False)
            for b in art["boxes"]: b.set_alpha(None); b.set_facecolor(to_rgba(col, .25)); b.set_edgecolor(col); b.set_linewidth(.5)
            for w in art["whiskers"]: w.set_alpha(1); w.set_linewidth(.5)
            for m in art["medians"]: m.set_linewidth(1); m.set_zorder(4)
            j = yy + rng.uniform(-.08, .08, len(v))
            ax.scatter(v[~t], j[~t], s=1.2, color=col, alpha=.5, linewidths=0, zorder=3)
            ax.scatter(v[t], j[t], edgecolors=col, zorder=3, **REF)   # truth-set males
            q = np.percentile(v, [10, 25, 50, 75, 90])
            stats.append({"clade": g, "cohort": c, "n": len(v), "n_truth": int(t.sum()), "median_bp": q[2] * 1e6, "q25_bp": q[1] * 1e6, "q75_bp": q[3] * 1e6, "q10_bp": q[0] * 1e6, "q90_bp": q[4] * 1e6})
    ax.set_ylim(len(leaves) - .5, -.5); column_axis(ax, "DYZ3", BLUE, "Length (Mb)", limit, list(np.arange(0, limit + .01, .5)))
    print("Fig. 2b-c n per clade (1kGP/HGDP/SGDP) and 1kGP DYZ3 truth diamonds:", {g: ("/".join(str(len(d[c])) for c in COHORT_SHADE), sum(f for _, f in d["1kGP"])) for g, d in data.items()}, flush=True)
    return stats


def biology_panels(fig, bottom=10, tall=46, sub=(70, 43)):
    """Fig. 2b-f: P and O2a DYZ3, J DYZ19 by cohort, unit ratio, DYZ1 and DYZ2 variability."""
    R, C, Q, N, inc = biology_inputs()
    H = fig.get_figheight() / MM; box = lambda x0, x1, b=bottom, t=tall: fig.add_axes([x0, b / H, x1 - x0, t / H])
    def style(ax, letter, t=tall):
        ax.spines[["top", "right"]].set_visible(False); ax.xaxis.set_label_coords(.5, -TITLE_Y / t)
    def on_top(ax, left=False, ncol=2):   # c-e legends: first entry on one baseline
        p = ax.get_position(); return dict(loc="upper left" if left else "upper center", bbox_to_anchor=(p.x0 if left else (p.x0 + p.x1) / 2, (bottom + tall + 4) / H), bbox_transform=fig.transFigure, ncol=ncol, borderaxespad=0, columnspacing=.8)
    clades = {c: subtree_rows(fig._valid, c) for c in SUBTREES}
    other = [r for r in Q if r["call"].startswith(sum(SUBTREES.values(), ()))]
    limit = max(2, .5 * np.ceil(max([r["DYZ3_bp"] for rows in clades.values() for r in rows.values()] + [float(r["DYZ3_bp"]) for r in other]) / 5e5))
    stats = []
    for (c, highlight), x0 in zip((("P", "R1b"), ("O2a", "O2a2b1a1")), (.03, .52)):
        t, d = box(x0, x0 + .18, *sub), box(x0 + .185, x0 + .465, *sub); stats += subtree_panel(t, d, clades[c], other, highlight, limit, "b" if c == "P" else ""); style(d, "", sub[1])
    # short DYZ19 of J2a1a1, by cohort
    ax = box(.055, .33)
    Jk = [r for r in inc if r["lineage"] == "J"]; Jq = [r for r in Q if r["call"].startswith("J")]; yes = lambda r: r["call"].startswith("J2a1a1"); kb = lambda R_: (np.array([float(r["DYZ19_bp"]) / 1e3 for r in R_]), [r["sample"] in TRUTH["dyz19"] for r in R_])   # HaploGrouper calls, as 1kGP
    G = [("1kGP", kb([r for r in Jk if r["haplogroup"].startswith("J2a1a1")]), kb([r for r in Jk if not r["haplogroup"].startswith("J2a1a1")]))]
    G += [(c, kb([r for r in Jq if r["cohort"] == c and yes(r)]), kb([r for r in Jq if r["cohort"] == c and not yes(r)])) for c in ("HGDP", "SGDP")]
    print("Fig. 2d J2a1a1 vs other J per cohort:", [(c, len(a[0]), len(b[0])) for c, a, b in G], flush=True)
    print("Fig. 2d diamonds (DYZ19 truth):", sum(sum(a[1]) + sum(b[1]) for _, a, b in G), flush=True)
    for i, (c, a, b) in enumerate(G):
        strip(ax, 2.6 * i, a[0], PURPLE, s=3, a=.8, seed=i, cap=.32, ref=a[1]); strip(ax, 2.6 * i + 1, b[0], GREY, s=3, a=.8, seed=i + 5, cap=.32, ref=b[1]); ax.text(2.6 * i + .5, -TITLE_Y / tall, c, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=FS["label"])
    ax.set_xticks([x for i in range(3) for x in (2.6 * i, 2.6 * i + 1)], [f"{len(v[0])}" for _, a, b in G for v in (a, b)]); ax.set_xlim(-.7, 6.9); ax.set_ylim(0, 350); ax.set_ylabel("DYZ19 length (kb)"); style(ax, "c")
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=2.5, color=PURPLE, label="J2a1a1"), Line2D([], [], ls="none", marker="o", ms=2.5, color=GREY, label="Other J")], **on_top(ax))
    # unit ratio by lineage
    ax = box(.40, .625); U = read_tsv(R / "dyz1_dyz2_lineages.tsv"); A = {r["measure"].split("_")[1]: float(r["median"]) for r in read_tsv(R / "dyz1_dyz2_truth_vs_yruler.tsv") if r["measure"].startswith("lineage_") and "_truth_" in r["measure"]}
    for i, r in enumerate(U):
        c = GREEN if r["lineage"] == "E" else (BLUE if float(r["ci95_high"]) < 1 else MUTED); ax.plot([float(r["ci95_low"]), float(r["ci95_high"])], [i, i], color=c, lw=1); ax.plot(float(r["unit_ratio_median"]), i, "o", ms=2.5, color=c)
        if r["lineage"] in A: ax.plot(A[r["lineage"]], i, "x", ms=2.5, color=ORANGE, mew=.6)
    ax.axvline(1, color=MUTED, lw=.6, ls="--"); ax.set_yticks(range(len(U)), [f'{r["lineage"]} ({r["men"]})' for r in U]); ax.set_xlim(.76, 1.2); ax.set_xlabel("DYZ1:DYZ2 unit ratio"); style(ax, "d")
    ax.legend(handles=[Line2D([], [], ls="none", marker="x", ms=2.5, color=ORANGE, mew=.6, label="From assemblies"), Line2D([], [], ls="none", marker="o", ms=2.5, color=MUTED, label="From reads")], **on_top(ax))
    # DYZ2 varies more than DYZ1
    ax = box(.725, .99); T = {r["assembly"]: r for r in read_tsv(ROOT / "data/method/yq12_truth_hallast.tsv")}; lin = {r["sample"]: r["lineage"] for r in C}
    kin = {r["sample"]: r["kinship_component_id"] for r in read_tsv(ROOT / "results/estimates_1kgp.relatedness_samples.tsv")}; rep = {}
    for r in sorted((r for r in inc if r["y_state"] != "par2_gain"), key=lambda r: r["sample"]): rep.setdefault(kin.get(r["sample"], r["sample"]), r)
    def centred(P):
        a = np.log2([p[0] for p in P]); b = np.log2([p[1] for p in P]); l = np.array([p[2] for p in P])
        for k in set(l): a[l == k] -= a[l == k].mean(); b[l == k] -= b[l == k].mean()
        return b, a
    P = [r for r in rep.values() if r["reference_male"] != "yes"]; tri = np.array([r["sample"] in TRUTH["yq12"] for r in P]); print("Fig. 2f diamonds (Yq12 truth):", int(tri.sum()), flush=True)
    men = centred([(float(r["DYZ1_bp"]), float(r["DYZ2_bp"]), r["lineage"]) for r in P]); asm = centred([(float(r["dyz1_units"]), float(r["dyz2_units"]), lin[s]) for s, r in T.items() if s in lin and s != "HG01261"])
    assert len(men[0]) == 1196 and len(asm[0]) == 75, (len(men[0]), len(asm[0]))
    ax.scatter(men[0][~tri], men[1][~tri], s=2, c=GREY, lw=0, alpha=.6, label="From reads (1kGP, 1,196)", rasterized=True); ax.scatter(men[0][tri], men[1][tri], edgecolors=GREY, rasterized=True, **REF); ax.scatter(*asm, s=6, marker="x", c=ORANGE, linewidths=.6, alpha=.9, label="From assemblies (75)", rasterized=True)
    L = .5 * np.ceil(max(np.abs(np.concatenate([*men, *asm]))) / .5); print("Fig. 2f symmetric range", L, flush=True); sd = lambda k: float(re.search(r"ratio ([0-9.]+)", N[k]).group(1)); ax.plot([-L, L], [-L, L], color=MUTED, lw=.6, ls="--", label="Equal variation"); r, a = sd("composition_sd_1kgp"), sd("composition_sd_assemblies"); ax.plot([-L, L], [-L * r, L * r], color=INK, lw=.9); ax.plot([-L, L], [-L * a, L * a], color=ORANGE, lw=.9)   # SD-ratio lines through the origin
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal"); F = {-3: "1/8×", -2: "1/4×", -1: "1/2×", 0: "1×", 1: "2×", 2: "4×", 3: "8×"}; T = [t for t in F if abs(t) <= L]; ax.set_xticks(T, [F[t] for t in T]); ax.set_yticks(T, [F[t] for t in T]); print("Fig. 2f points outside the former +-1.5 x +-1.2 view:", int(sum(((abs(x) > 1.5) | (abs(y) > 1.2)).sum() for x, y in (men, asm))), flush=True); ax.set_xlabel("DYZ2 relative to lineage mean"); ax.set_ylabel("DYZ1 relative to lineage mean"); from matplotlib.legend_handler import HandlerTuple
    dot, cross, eq = Line2D([], [], ls="none", marker="o", ms=2.5, color=GREY), Line2D([], [], ls="none", marker="x", ms=2.5, color=ORANGE, mew=.6), Line2D([], [], color=MUTED, lw=.6, ls="--")
    ax.legend([(Line2D([], [], color=ORANGE, lw=.9), cross), (Line2D([], [], color=INK, lw=.9), dot), eq], [f"From assemblies (75)\n(SD ratio {a:g})", f"From reads (1kGP, 1,196)\n(SD ratio {r:g})", "Equal variation\n(SD ratio 1)"], handler_map={tuple: HandlerTuple(ndivide=1)}, **on_top(ax, left=True, ncol=1)); style(ax, "e")   # points and line share one entry
    B = [.719, .18, .25, .19]; ins = ax.inset_axes(B); x0, y0 = -L + 2 * L * B[0], -L + 2 * L * B[1]; print("Fig. 2f points under the inset:", int(sum(((x >= x0) & (x <= x0 + 2 * L * B[2]) & (y >= y0) & (y <= y0 + 2 * L * B[3])).sum() for x, y in (men, asm))), flush=True); S = [("Assemblies", "assemblies"), ("1kGP", "1kgp"), ("HGDP", "hgdp"), ("SGDP", "sgdp")]
    for i, (n, k) in enumerate(S):
        m_ = re.search(r"ratio ([0-9.]+) \(95% ([0-9.]+) to ([0-9.]+)\)", N[f"composition_sd_{k}"]); v, lo, hi = map(float, m_.groups()); ins.plot([lo, hi], [3 - i, 3 - i], color=INK, lw=.75); ins.plot(v, 3 - i, "x" if k == "assemblies" else "o", ms=2.5 if k == "assemblies" else 2, mew=.6, color=ORANGE if k == "assemblies" else INK)
    ins.axvline(1, color=MUTED, lw=.6, ls="--"); ins.set_yticks(range(4), [n for n, _ in S[::-1]]); ins.set_xlim(.55, 1.05); ins.set_xticks([.7, 1]); ins.set_xticks([.8, .9], minor=True); ins.set_ylim(-.5, 3.5); ins.tick_params(length=1.5, pad=1); ins.set_xlabel("SD ratio", fontsize=FS["tick"], labelpad=3); ins.spines[["top", "right"]].set_visible(False); ins.patch.set_alpha(.85)
    return stats


def population_figure(valid, out, snapshot_path, completed):
    samples = {r["sample"]: r for r in valid}
    tree = prune_tree(parse_newick(TREE_SOURCE), set(samples))
    root, leaves = collapse_tree(tree, samples)
    assert all(n == 1 for n in Counter(leaf["group"] for leaf in leaves).values())
    from cohort_lineages import load_combined_samples
    map_samples = load_combined_samples()
    land = json.loads((ROOT / "data/figure_metadata/world_land.geojson").read_text())
    from palette import MM
    H = 210   # mm: d-f 10-56, b-c 70-113, a 124-197, legend 203-209
    fig = plt.figure(figsize=(180 * MM, H * MM))
    for s, x, y in (("a", 0, H), ("b", 0, 118), ("c", 89.5, 118), ("d", 0, 63), ("e", 64, 63), ("f", 119, 63)): letter(fig, s, (x + 1) / 180, (y - 1) / H)
    grid = fig.add_gridspec(1, 6, left=.02, right=.99, bottom=124 / H, top=197 / H,
                            width_ratios=[1.3, .735, .945, .945, .945, 1.05], wspace=.09)   # tree 35.7 mm as before
    axes = [fig.add_subplot(grid[0, i]) for i in range(6)]; sp = axes.pop(1)   # superpopulation bars next to the tree
    draw_collapsed_tree(axes[0], root, leaves, panel="a", room=1.9)
    from palette import SUPERPOP, SUPERPOP_NAME, NAMES as NM, TITLE as TT
    for y, node in enumerate(leaves):   # 1kGP superpopulation per clade
        c = Counter(r["superpopulation"] for r in node["members"]); left = 0
        for k, col in SUPERPOP.items():
            w = 100 * c[k] / len(node["members"]); sp.barh(y, w, left=left, height=.6, color=col, linewidth=0); left += w
    sp.set_ylim(len(leaves) - .4, -.7); column_axis(sp, "Superpopulation", None, "Males (%)", 100, [0, 50, 100])
    print("Fig. 2a superpopulation counts:", {n["group"]: dict(Counter(r["superpopulation"] for r in n["members"])) for n in leaves}, flush=True)
    p0, p1 = sp.get_position(), axes[-1].get_position()
    fig.legend(handles=[Patch(facecolor=v, label=SUPERPOP_NAME[k]) for k, v in SUPERPOP.items()], loc="upper center", bbox_to_anchor=((p0.x0 + p1.x1) / 2, (H - 1) / H), ncol=5, title=NM["superpop"], title_fontproperties=TT(FS["legend"]), borderaxespad=0)
    fig.legend(handles=[Line2D([], [], markeredgecolor=INK, label="From reads, truth-set male", **REF_HANDLE)], loc="upper left", bbox_to_anchor=(.03, (H - 4.2) / H), fontsize=FS["legend"], borderaxespad=0)   # one diamond key for a-d, f
    stats = frozen_trait(axes[1], leaves, "DYZ3_bp", BLUE, 1e6, 3, [0, 1, 2, 3],
                         "DYZ3", "c", "Length (Mb)", lambda g: g in ("R1b", "O2a2b1a1"), show_panel=False)
    stats += frozen_trait(axes[3], leaves, "Yq12_bp", GREEN, 1e6, 60, [0, 20, 40, 60],
                          "Yq12", "d", "Length (Mb)", lambda g: False, show_panel=False)
    stats += frozen_trait(axes[2], leaves, "DYZ19_bp", PURPLE, 1e3, 550, [0, 200, 400],
                          "DYZ19", "e", "Length (kb)", lambda g: g.startswith("J2a1a1"), show_panel=False)
    for r in valid: r["unit_ratio"] = r["ratio_bp"] * 2404.4 / 3478.9   # bp per Hallast 2026 unit
    ratio_limit = max(2.5, .5 * np.ceil(max(r["unit_ratio"] for r in valid) / .5))
    stats += frozen_trait(axes[4], leaves, "unit_ratio", RATIO, 1, ratio_limit, np.arange(0, ratio_limit + .01, 1),
                          "DYZ1:DYZ2", "f", "Unit ratio", lambda g: False, show_panel=False)
    axes[4].axvline(1, color=RATIO, linestyle="--", linewidth=.6, zorder=1)
    print("Fig. 2a diamonds per column:", {a.get_title(): a._triangles for a in axes[1:]}, flush=True)
    fig._valid = valid; subtrees = biology_panels(fig)
    save(fig, out, "main/fig2_biology")
    write_rows(AUX / "figure2b_subtree_summary.tsv", subtrees)
    write_rows(AUX / "figure2_tree_summary.tsv", stats)
    membership = [{"sample": r["sample"], "haplogroup_2019": r["haplogroup"],
                   "display_group": r["display_group"], "family": r["family"],
                   "reference_male": r["reference_male"],
                   "y_status": r["y_status"], "population_included": r["population_included"],
                   "source": str(Path(snapshot_path).resolve().relative_to(ROOT))}
                  for node in leaves for r in node["members"]]
    write_rows(AUX / "figure2_tree_membership.tsv", membership)
    counts = Counter(r["lineage"] for r in map_samples)
    write_rows(AUX / "figure2_lineage_counts.tsv", [
        {"lineage": group, "n_unique_individuals": n, "percent": 100*n/len(map_samples)}
        for group,n in counts.items()])
    pie_fig, pie_ax = plt.subplots(figsize=(180 * MM, 77 * MM)); pie_ax.set_anchor("N")
    pie_map(pie_ax, None, map_samples, land, size_legend=True, legend_y=-51)   # legends 2 mm (4 degrees) lower
    print(f"Supplementary Fig. 1 map: n={len(map_samples):,} individuals", flush=True)
    pie_fig.subplots_adjust(left=.005, right=.995, top=.995, bottom=.005)
    pie_fig.canvas.draw()
    anchor=pie_fig.transFigure.inverted().transform(pie_ax.transData.transform((55.38,-42.63)))
    map_legend=lineage_legend(pie_fig, map_samples, anchor[1], x=anchor[0], fontsize=FS["legend"])
    align_map_legend_titles(pie_ax,map_legend)
    save(pie_fig, out, "supp/suppfig1_world_map")
    return root, leaves


def population_supplement(valid, out):
    tree = ["A", "B", "D", "E", "C", "F", "G", "H", "I", "J", "L", "T", "NO", "N", "O", "Q", "R1a", "R1b", "R"]   # Fig. 2a tree order
    present = {r["lineage_display"] for r in valid}; order = [g for g in tree if g in present] + sorted(present - set(tree))
    clades = [("I–J", "I", "J"), ("N–O", "NO", "O"), ("P (Q, R)", "Q", "R")]   # brackets over sister lineages
    print("ED7 n per lineage:", ";  ".join(f"{g}: {sum(r['lineage_display'] == g for r in valid)}" for g in order), flush=True)
    panels = TARGETS + [("unit_ratio", "unit_ratio", "DYZ1:DYZ2", 1, "Unit ratio")]   # h as Fig. 2a; Y dosage by lineage: ED5b
    for r in valid: r["unit_ratio"] = r["ratio_bp"] * 2404.4 / 3478.9   # bp per Hallast 2026 unit
    fig, axes = plt.subplots(4, 2, figsize=(180 * MM, 184 * MM), gridspec_kw={"hspace": .5, "wspace": .12})
    fig.subplots_adjust(left=.065, right=.995, top=.94, bottom=.045)
    rng = np.random.default_rng(271)
    for panel, (ax, (_, key, label, scale, units)) in enumerate(zip(axes.flat, panels)):
        from palette import COLOR, head
        color = COLOR[label]
        for x, group in enumerate(order):
            selected = [r for r in valid if r["lineage_display"] == group]
            values = np.array([r[key] / scale for r in selected])
            jitter = x + rng.uniform(-.23, .23, len(selected))
            truth = np.array([r["sample"] in TRUTH.get(TRUTH_OF.get(key), ()) for r in selected])
            included = np.array([in_summary(r, key) for r in selected])   # medians only
            ax.scatter(jitter[~truth], values[~truth], s=2.5, color=color, alpha=.35, linewidths=0)
            ax.scatter(jitter[truth], values[truth], edgecolors=color, **REF)
            if included.any():
                ax.plot([x - .4, x + .4], [np.median(values[included])]*2, color=INK, linewidth=1)   # median, as strip()
        if key == "unit_ratio":   # equal units, as Fig. 2a
            ax.axhline(1, color=color, linestyle="--", linewidth=.6, zorder=1)
        if key == "DYZ3_bp":   # short DYZ3 lineages, as the Fig. 2 bands
            for g in ("R1a", "R1b"):
                if g in order: ax.axvspan(order.index(g) - .45, order.index(g) + .45, color=BAND["DYZ3"], zorder=0, lw=0)
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * .12)
        for name, g0, g1 in clades:
            if g0 in order and g1 in order:
                i0, i1 = order.index(g0), order.index(g1); yb = hi + (hi - lo) * .05
                ax.plot([i0 - .35, i1 + .35], [yb, yb], color=MUTED, lw=.6, clip_on=False)
                ax.text((i0 + i1) / 2, yb + (hi - lo) * .015, name, ha="center", va="bottom", fontsize=FS["point"], color=MUTED)
        last = True   # every panel: name and number of males, vertical (same ticks everywhere)
        ax.set_xticks(range(len(order)), [(f"{g} ({sum(r['lineage_display'] == g for r in valid)})" if last else g) for g in order], rotation=90)   # vertical: 19 lineages per half page
        ax.set_xlim(-.7, len(order) - .3)
        ax.set_ylabel(units)
        head(ax, label)
        letter(fig, chr(ord("a") + panel), .002 if panel % 2 == 0 else .505, ax.get_position().y1 + (5 * MM + (TITLE_PAD - 3) / 72) / fig.get_figheight())
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    for ax in axes.flat[len(panels):]: ax.set_visible(False)
    reads_legend(fig)
    save(fig, out, "ed/ed7_population")


def depth_supplement(valid, out):
    """ED5: Y dosage against coverage, and by lineage."""
    from palette import DOSE, head
    from scipy.stats import kruskal
    R = [dict(r, level="1" if r["level"] == "full" else r["level"]) for r in read_tsv(ROOT / "results/model/depth/y_dose_depth_chr8.tsv")]   # release chr8:91-101 Mb window (um_depthchr8.sh); full depth is level 1
    R = [r for r in R if r["level"] in DEPTHS]
    lv = sorted({r["level"] for r in R}, key=lambda k: DEPTHS[k]); x = [DEPTHS[k] for k in lv]
    D = {k: np.array([abs(float(r["diff_pct"])) for r in R if r["level"] == k]) for k in lv}; q = np.array([np.percentile(D[k], [10, 50, 90]) for k in lv])
    for k, (lo, m, hi) in zip(lv, q): print(f"ED5 Y dosage {k} ({DEPTHS[k]:.1f}x): {len(D[k])} males, |difference from full depth| median {m:.3f}% p10 {lo:.3f}% p90 {hi:.3f}% max {D[k].max():.3f}%", flush=True)
    fig, (ax, b) = plt.subplots(1, 2, figsize=(180 * MM, 66 * MM), gridspec_kw={"width_ratios": [1, 1.7]})
    ax.fill_between(x, q[:, 0], q[:, 2], color=DOSE, alpha=.15, lw=0, label="10th to 90th percentile"); ax.plot(x, q[:, 1], "o-", color=DOSE, ms=2.5, lw=.75, label="Median")
    ax.set_xscale("log"); ax.set_xticks([1, 2, 5, 10, 20, DEPTHS["1"]], ["1", "2", "5", "10", "20", f"{DEPTHS['1']:.0f}"]); ax.minorticks_off(); ax.set_ylim(0, None); ax.grid(True)
    ax.set_xlabel("Median autosomal coverage (×)"); ax.set_ylabel("Absolute difference from full depth (%)"); ax.legend(loc="upper right")
    head(ax, f"Coverage ({len(D[lv[0]])} males)")
    tree = ["A", "B", "D", "E", "C", "F", "G", "H", "I", "J", "L", "T", "NO", "N", "O", "Q", "R1a", "R1b", "R"]   # ED7 and Fig. 2a order
    order = [g for g in tree if any(r["lineage_display"] == g for r in valid)]; lo, hi = .94, 1.06; rng = np.random.default_rng(5)
    V = {g: np.array([float(r["y_dose"]) for r in valid if r["lineage_display"] == g]) for g in order}
    for i, g in enumerate(order):
        v = V[g]; c = np.clip(v, lo, hi); j = i + rng.uniform(-.25, .25, len(v)); inside = (v > lo) & (v < hi)
        b.scatter(j[inside], c[inside], s=2.5, color=DOSE, alpha=.35, lw=0)
        for m, k in (("v", v <= lo), ("^", v >= hi)): b.scatter(j[k], c[k], s=6, marker=m, color=DOSE, lw=0, clip_on=False)
        b.plot([i - .4, i + .4], [np.median(v)] * 2, color=INK, lw=1)
    big = [g for g in order if len(V[g]) > 20]; kw = kruskal(*[V[g] for g in big])
    med = {g: float(np.median(V[g])) for g in big}; g0, g1 = min(med, key=med.get), max(med, key=med.get)
    print(f"ED5b: {len(valid)} males, {len(big)} lineages with > 20 males, medians {g0} {med[g0]:.4f} to {g1} {med[g1]:.4f}, Kruskal-Wallis H {kw.statistic:.2f} P {kw.pvalue:.2e}; outside {lo}-{hi}: {sum(int(((V[g] <= lo) | (V[g] >= hi)).sum()) for g in order)}", flush=True)
    b.text(.01, .97, f"Kruskal–Wallis P = {kw.pvalue:.1e} ({len(big)} lineages with > 20 males)", transform=b.transAxes, va="top", fontsize=FS["point"])
    b.set_xticks(range(len(order)), [f"{g} ({len(V[g])})" for g in order], rotation=90); b.set_xlim(-.7, len(order) - .3)
    b.set_ylim(lo, hi + .016); b.set_yticks(np.arange(lo, hi + .001, .02)); b.set_ylabel("Y dosage (relative)")   # room for the test line; b.grid(axis="y"); b.set_axisbelow(True)
    head(b, f"Lineage (1000 Genomes, {len(valid):,} males)")
    fig.tight_layout(w_pad=2.5)
    for s_, a_ in zip("ab", (ax, b)): fig.text(a_.get_position().x0 - (.075 if s_ == "a" else .06), a_.get_position().y1 + TITLE_PAD / 72 / fig.get_figheight(), s_, fontsize=FS["letter"], fontweight="bold", va="baseline")   # letters on the title baselines
    save(fig, out, "ed/ed5_depth")


def j_subtree_figure(valid, out):
    """ED8: J tree collapsed by haplogroup, DYZ19, bootstrap at J2a1a1."""
    members = [dict(r, display_group=r["haplogroup"][:7]) for r in valid if r["haplogroup"].startswith("J")]   # 7 characters: J2a1a1a, J2a1a1b
    index = {r["sample"]: r for r in members}
    root, leaves = collapse_tree(prune_tree(parse_newick(TREE_SOURCE), set(index)), index)
    B, full = {}, parse_newick(ROOT / "results/tree/RAxML_bipartitions.sup")   # RAxML supports, 100 bootstraps
    def tips(n): return frozenset([n["name"]]) if not n["children"] else frozenset().union(*map(tips, n["children"]))
    every = tips(full)
    def splits(n):
        if n["children"]: B[tips(n)] = B[every - tips(n)] = n["name"]; [splits(c) for c in n["children"]]
    splits(full)
    h = 15 + 3.3 * len(leaves)   # mm: 3.3 per row, title and axis rows
    fig = plt.figure(figsize=(180 * MM, h * MM))
    grid = fig.add_gridspec(1, 3, left=.03, right=.98, bottom=9.5 / h, top=1 - 5.5 / h, width_ratios=[1.6, .62, 1.75], wspace=.05)   # gaps as Fig. 2a
    tree_ax, sp, ax = [fig.add_subplot(grid[0, i]) for i in range(3)]
    from palette import SUPERPOP, SUPERPOP_NAME, NAMES as NM, TITLE as TT
    for y, node in enumerate(leaves):   # 1kGP superpopulation per clade, as Fig. 2a
        c = Counter(r["superpopulation"] for r in node["members"]); left = 0
        for k, col in SUPERPOP.items():
            w = 100 * c[k] / len(node["members"]); sp.barh(y, w, left=left, height=.6, color=col, linewidth=0); left += w
    sp.set_ylim(len(leaves) - .4, -.7); column_axis(sp, "Superpopulation", None, "Males (%)", 100, [0, 50, 100])
    draw_collapsed_tree(tree_ax, root, leaves, panel="a", room=1.3)   # labels end next to the bars
    def boot(n):   # supports of 70 or more at every node
        S = frozenset(r["sample"] for r in n["members"]); v = B.get(S)
        if n["children"]:
            print(f"ED8 node {len(S)} males ({', '.join(sorted({r['haplogroup'][:7] for r in n['members']}))[:60]}): bootstrap {v if v is not None else 'no split in the full tree'}", flush=True)
            if v is not None and int(v) >= 70: tree_ax.annotate(v, (n["x"], n["y"]), xytext=(-2, 1.5), textcoords="offset points", ha="right", va="bottom", fontsize=FS["point"] - 1, color=INK)   # clear of the branches
        for c in n["children"]: boot(c)
    boot(root)
    frozen_trait(ax, leaves, "DYZ19_bp", PURPLE, 1e3, 350, [0, 100, 200, 300], "DYZ19", "b", "Length (kb)", lambda g: g.startswith("J2a1a1"), show_panel=False)
    top = max(r["DYZ19_bp"] for r in members) / 1e3; assert top <= 350, top   # largest J male fits
    ax.set_xlim(-8, 350); ax.set_xticks(range(0, 351, 50)); ax.set_xticks([], minor=True); print(f"ED8 largest DYZ19 {top:.1f} kb", flush=True)
    # one panel as Fig. 2a: no letters
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=2.5, color=MUTED, label="From reads"), Line2D([], [], markeredgecolor=MUTED, label="From reads, truth-set male", **REF_HANDLE)], loc="lower left", borderaxespad=.4)   # empty left side of b
    ax.add_artist(ax.get_legend()); ax.legend(handles=[Patch(facecolor=v, label=SUPERPOP_NAME[k]) for k, v in SUPERPOP.items()], loc="center left", bbox_to_anchor=(0, .5), title=NM["superpop"], title_fontproperties=TT(FS["legend"]), borderaxespad=.4, alignment="left")   # empty left side of b
    save(fig,out,"ed/ed8_J_DYZ19")
    write_rows(AUX / "figureS5_J_membership.tsv",[
        {"sample":r["sample"],"haplogroup":r["haplogroup"],"family":r["family"],
         "reference_male":r["reference_male"],"DYZ19_bp":r["DYZ19_bp"]}
        for n in leaves for r in n["members"]])


def assembly_supplement(out):
    """Preserve the complete-array annotation view, separate from population reads."""
    calls = {r["ID"]: r["Haplogroup"] for r in read_tsv(ROOT / "data/method/haplogroups_1kgp_isogg2019_decode1.tsv")}
    parents, children = tree_data(ROOT / "data/figure_metadata/isogg2019_tree.tsv"); add_reference_calls(calls, parents)
    groups = defaultdict(list)
    for row in load_snapshot()[1]:   # display groups from the current 1kGP estimates (dose >= 0.9)
        if row["y_dose"] >= .9:
            groups[clade_of(calls[row["sample"]], parents)].append(row)
    _, order = pruned_tree(groups, parents, children)
    rows = load_composition(ROOT / "data/method/yq12_truth_hallast.tsv", ROOT / "data/markers_u/truth_regions.tsv",
                            calls, parents, groups)
    pop = [r for r in load_snapshot()[1] if r.get('population_included') in (True, 'True')]   # males of the population analyses
    composition_supplement(rows, order, out, pop)
    write_rows(AUX / "figureS2_composition_samples.tsv", rows)


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
        "scripts/cohort_lineages.py", "results/cohorts/lineage_calls.tsv",
        "scripts/figure1b_selected.py", "scripts/validation_plots.py", "scripts/depth_benchmark.py",
        "scripts/pedigree_plots.py", "scripts/validation_ratio.py")]
    inputs += [p for p in (ROOT / "data/figure_metadata").iterdir()
               if p.suffix in {".tsv", ".kin0", ".psam", ".xlsx", ".geojson"} or p.name == "sgdp_metadata.txt"]
    inputs += list((ROOT / "results/cohorts").glob("*.tsv"))
    outputs = [p for p in AUX.glob("*.tsv") if not p.name.startswith("snapshot_")] + [p for p in out.rglob("*.pdf")]
    outputs += [AUX / f"snapshot_{tag}_samples.tsv"]
    outputs += list((ROOT / "results").glob(Path(path).stem + ".*.tsv"))
    outputs += [AUX / f"snapshot_{tag}_summary.json"]
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
    (AUX / f"snapshot_{tag}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def render_snapshot(path, raw, valid, pedigree, benchmarks, depth, out):
    tag = snapshot_tag(path)
    summary = snapshot_summary(raw, valid, pedigree, path)
    kinship, sensitivities, info, qc, sample_rows = autosomal_check(raw, valid, pedigree)
    summary["relatedness"] = info
    root, leaves = population_figure(valid, out, path, len(raw))
    population_supplement(valid, out)
    from figure1 import ed4   # ED4 in the Fig. 1c style
    ed4(out / 'ed/ed4_lineage_validation.pdf')
    depth_supplement(valid, out)
    j_subtree_figure(valid, out)
    assembly_supplement(out)
    write_rows(AUX / "figure1_validation_summary.tsv", benchmarks)
    write_rows(AUX / "figure1_depth_summary.tsv", depth)
    write_rows(AUX / f"snapshot_{tag}_samples.tsv", [
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
    (AUX / f"snapshot_{tag}_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
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
