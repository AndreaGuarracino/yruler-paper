# Short-read final model v3 (v3a), 2026-10-28

Complete package to replace the manuscript numbers: CV, depth, lineage-out, holdout, 1kGP cohort, female background, release. It is the same model as the "CV v3a" already used in main.pdf.

## Model

v3 = v2 (grid A + CHM13 + 99 quality-selected HPRC haplotypes at 50%) with two changes to the truth:
- 14 suspicious truth arrays excluded (7 DYZ3 and 7 DYZ19; local duplication or collapse flagged by Flagger/NucFlag; truth_qc/x/rule_final.txt);
- YqHet2.7 and DYZ18sub truth counted over all 2.7 kb copies inside Yq12 (truth_qc/x/truth_all.tsv), not only in the annotated block.

Males evaluated: DYZ3 97, DYZ19 111, the other 6 classes 82.

## Files

| File | Content |
|---|---|
| depth/depth_f<fraction>.tsv | Predictions per repeat, fold, class and male at 7 depth fractions (1 = full CV) |
| depth/summary_mean_of_10_medians.tsv | Mean of the 10 median absolute errors per fraction and class, with the lowest and highest repeat |
| lineage_out.tsv, lineage_out_summary.tsv | One lineage left out at a time (18 lineages): predictions and medians per class and lineage |
| holdout_summary.txt | Fixed holdout of the literature males |
| cohort_1kgp.tsv | Estimates for the 1,639 1kGP samples with the v3 model (flag and y_state columns) |
| female_background_heldout.txt, female_background_filter.txt | Signal in females (no Y) as a percentage of the typical length |
| bed_vs_whole.tsv | BED-mode estimate against the whole file, per sample |
| final_model.tsv, release_model.tsv | Thresholds, number of markers and K per class; release model (PASS check on HG00126) |
| old_vs_new.tsv | Numbers of the frozen 2026-10-27 model against v3 |
| numbers_to_update.md | Lines of the .tex files with the values to update |
| md5.txt | Checksums |

Depth summary (mean of the 10 medians, %):

| Fraction | DYZ3 | Yq12 | DYZ1 | DYZ2 | DYZ18S | YqHet2.7 | DYZ18sub | DYZ19 |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.23 | 1.76 | 2.30 | 2.28 | 3.68 | 5.34 | 4.30 | 2.21 |
| 0.5 | 2.06 | 1.72 | 2.28 | 2.59 | 3.90 | 5.73 | 4.29 | 2.46 |
| 0.167 | 1.93 | 1.68 | 2.25 | 2.39 | 3.80 | 6.18 | 4.40 | 3.03 |
| 0.033 | 2.88 | 1.80 | 2.56 | 2.17 | 4.47 | 6.89 | 4.49 | 3.82 |

## Provenance

Cluster: /path/to/project/pkg/final_model_short_v3_2026-10-28 (markers_u25/cv_v3a, cv_v3a_depth_f*, lolo_v3a; scripts um_pkg.py and um_cvexport.py). Same splits as the main CV (family groups, numpy default_rng(7 + repeat)); thresholds, markers and calibration estimated on the training males at full depth only; in the depth runs the test reads are subsampled.

Note: HG002, HG003 and HG005 (Google NovaSeq reads) are part of this CV. HG002 DYZ2 is -12.6% here and -2.0% with markers of GC >= 38.7% (12 of 31 bases); the later packages correct this.
