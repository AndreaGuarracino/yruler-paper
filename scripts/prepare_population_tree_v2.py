#!/usr/bin/env python3
"""Root the delivered v2 ML tree on its existing A0 split, retaining precision."""
import hashlib
import json
import random
from pathlib import Path

import frozen_snapshot as f
from audit_population_tree import indexed, mrca


def main():
    directory = f.ROOT / "results/ytree_v2_2026-10-27"
    source = directory / "ytree.treefile"
    original = f.parse_newick(source)
    nodes, tips = indexed(original)
    calls = {r["ID"]: r["Haplogroup"] for r in f.read_tsv(
        f.ROOT / "data/method/haplogroups_1kgp_isogg2019_decode1.tsv")}
    a0 = {s for s in tips if calls[s].startswith("A0")}
    assert a0 == {"HG01890", "HG02982", "HG02984"}
    edges = [n for n in nodes if n["parent"] is not None
             and (n["tips"] == a0 or set(tips) - n["tips"] == a0)]
    assert len(edges) == 1, "A0 must already define an edge in the ML tree"
    edge = edges[0]
    graph = [[] for _ in nodes]
    for n in nodes:
        parent = n["parent"]
        if parent is not None:
            graph[parent].append((n["index"], n["length"], n["name"]))
            graph[n["index"]].append((parent, n["length"], n["name"]))

    def orient(i, parent, length, label):
        children = [orient(j, i, distance, support)
                    for j, distance, support in graph[i] if j != parent]
        return {"name": label if children else nodes[i]["name"],
                "length": length, "children": children}

    u, v = edge["index"], edge["parent"]
    rooted = {"name": "", "length": 0., "children": [
        orient(u, v, edge["length"] / 2, edge["name"]),
        orient(v, u, edge["length"] / 2, edge["name"])]}

    def emit(n):
        children = "(" + ",".join(map(emit, n["children"])) + ")" if n["children"] else ""
        return children + n["name"] + ":" + format(n["length"], ".15g")

    f.TREE_SOURCE.write_text(emit(rooted) + ";\n")
    derived = f.parse_newick(f.TREE_SOURCE)
    rn, rt = indexed(derived)
    assert set(rt) == set(tips)
    assert rn[mrca(rn, rt, a0)]["tips"] == a0
    original_splits = {frozenset(n["tips"]) for n in nodes}
    assert all(frozenset(n["tips"]) in original_splits
               or frozenset(set(rt) - n["tips"]) in original_splits for n in rn)

    def distance(ns, ts, a, b):
        anc = mrca(ns, ts, [a, b])
        return ns[ts[a]]["depth"] + ns[ts[b]]["depth"] - 2 * ns[anc]["depth"]

    rng = random.Random(20261027)
    pairs = [rng.sample(sorted(tips), 2) for _ in range(1000)]
    max_error = max(abs(distance(nodes, tips, *p) - distance(rn, rt, *p)) for p in pairs)
    assert max_error < 1e-10
    j = {s for s in tips if calls[s].startswith("J2a1a1")}
    assert rn[mrca(rn, rt, j)]["tips"] == j and len(j) == 42
    h = rt["HG03742"]
    sibling = rn[rn[h]["parent"]]["tips"] - {"HG03742"}
    assert len(sibling) == 280
    assert all(calls[s].startswith(("N", "O")) for s in sibling)
    sites = f.read_tsv(directory / "sites.tsv")
    assert len(sites) == 79941 and sum(int(r["mac"]) == 1 for r in sites) == 40273
    missing = f.read_tsv(directory / "missing_per_sample.tsv")
    assert {r["sample"] for r in missing} == set(tips)
    report = {
        "source": str(source.relative_to(f.ROOT)),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "derived_rooted_tree": str(f.TREE_SOURCE.relative_to(f.ROOT)),
        "derived_sha256": hashlib.sha256(f.TREE_SOURCE.read_bytes()).hexdigest(),
        "root_outgroup": sorted(a0), "rooting": "Midpoint of existing A0 edge; no unrooted split changed",
        "source_tips": len(tips), "sites": len(sites), "singletons": 40273,
        "distance_checks": 1000, "maximum_distance_error": max_error,
        "J2a1a1_carriers_and_MRCA_tips": 42, "HG03742_sister_clade_tips": len(sibling),
        "delivery_issue": "v2val.rooted.nwk is rooted on HG01890 alone, not all three A0 samples; use a derived copy of the full-precision ML tree.",
        "limits": "Provisional source tree; original alignment not delivered, not all fine classification clades agree, no time calibration."
    }
    (f.ROOT / "results/figure2_tree_v2_preparation_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
