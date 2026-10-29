#!/usr/bin/env python3
"""Render main and supplementary figures from the short-read model results."""
import json

import matplotlib.pyplot as plt

import frozen_snapshot as f
import manuscript_figures as m
import validation_plots as v


def main():
    out = f.ROOT / "manuscript/figures"
    cohort = f.ROOT / "results/estimates_1kgp.tsv"; tag = f.snapshot_tag(cohort)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelsize": 10, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "savefig.facecolor": "white", "axes.spines.top": False,
                         "axes.spines.right": False})
    v.CV,v.LO,v.DEP,_ = v.load()
    benchmarks,depth = f.validation_summary()
    from validation_ratio import summary as ratio_summary
    ratio,ratio_depth=ratio_summary(v.CV,v.LO,v.DEP,f.DEPTHS)
    benchmarks.append({key:ratio[key] for key in benchmarks[0]})
    depth.extend(ratio_depth)
    (out/'figure1b_ratio_validation_summary.json').write_text(json.dumps(ratio,indent=2)+'\n')
    raw, valid, pedigree = f.load_snapshot(cohort)
    m.method_figure(f.legacy_catalog(valid), f.ROOT / "data/catalog_d376/fs_pairs.tsv", out,
                    benchmarks=benchmarks, depth=depth)
    f.population_figure(valid, out, cohort, len(raw))
    f.population_supplement(valid, out)
    f.j_subtree_figure(valid, out)
    f.depth_supplement(depth,out)
    from figure1b_selected import main as render_validation
    render_validation()
    from pedigree_plots import main as render_pedigrees
    render_pedigrees()
    f.write_rows(out / "figure1_validation_summary.tsv",benchmarks)
    f.write_rows(out / "figure1_depth_summary.tsv",depth)
    summary = f.snapshot_summary(raw, valid, pedigree, cohort)
    summary['model_status'] = 'short-read model; complete CV, coverage, lineage-out and cohort'
    kinship, sensitivity, info, qc, audit = f.autosomal_audit(raw, valid, pedigree)
    summary["relatedness"] = info
    for suffix, data in [("sensitivity", sensitivity), ("kinship", kinship),
                         ("sample_qc", qc), ("relatedness_samples", audit),
                         ("relationships", f.pedigree_relationships(raw, pedigree))]:
        f.write_rows(cohort.with_suffix("." + suffix + ".tsv"), data)
    f.write_rows(out / f"snapshot_{tag}_samples.tsv", valid)
    (out / f"snapshot_{tag}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    f.snapshot_manifest(cohort,out)
    print(json.dumps({k: summary[k] for k in ["valid_males", "valid_reference_males", "population_summary_n", "complete_father_son_pairs", "median_bp_by_target", "subsets"]}, indent=2))


if __name__ == "__main__":
    main()
