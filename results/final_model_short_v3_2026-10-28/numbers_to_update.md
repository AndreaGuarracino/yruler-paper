# Numbers to insert for v3a (2026-10-28)

The v3a CV is already in the manuscript. Depth and lineage exclusion, marked "pending", remain to be inserted. Lines of the .tex files read on 2026-10-28 at 08:50; no manuscript file was modified.

## Pending lines

- online-methods.tex:61 and extended.tex:267: "Leave-one-lineage-out and downsampling evaluations of the current procedure were still incomplete at this draft update"; they are now complete (values below).
- extended.tex:40: draft note "Coverage and excluded-lineage analyses of the current procedure are pending".
- Accuracy table (lineage-out column) and the sentences on depth.

## Lineage exclusion (pooled median absolute error, %; 18 lineages)

| Class | v3a | Old (frozen) |
|---|---|---|
| DYZ3 | 2.20 | 2.39 |
| Yq12 | 1.70 | 1.69 |
| DYZ1 | 2.39 | 2.31 |
| DYZ2 | 2.47 | 2.58 |
| DYZ18S | 3.71 | 3.77 |
| YqHet2.7 | 5.78 | 6.33 |
| DYZ18sub | 4.46 | 4.56 |
| DYZ19 | 2.04 | 2.61 |

Per lineage: lineage_out_summary.tsv.

## Depth (mean of the 10 medians, %; lowest and highest repeat in brackets)

| Fraction | DYZ3 | Yq12 | DYZ19 |
|---|---|---|---|
| 1 | 2.23 (2.12-2.32) | 1.76 (1.62-1.89) | 2.21 (1.99-2.35) |
| 0.667 | 2.23 (2.09-2.56) | 1.76 (1.66-1.83) | 2.37 (2.13-2.55) |
| 0.5 | 2.06 (1.84-2.45) | 1.72 (1.63-1.85) | 2.46 (2.21-2.66) |
| 0.333 | 2.00 (1.88-2.24) | 1.73 (1.62-1.94) | 2.47 (2.24-2.65) |
| 0.167 | 1.93 (1.73-2.13) | 1.68 (1.61-1.74) | 3.03 (2.83-3.28) |
| 0.067 | 2.23 (2.05-2.39) | 1.76 (1.65-1.85) | 3.08 (2.78-3.28) |
| 0.033 | 2.88 (2.62-3.40) | 1.80 (1.73-1.96) | 3.82 (3.52-4.11) |

All 8 classes: depth/summary_mean_of_10_medians.tsv. The actual depth of each subsampled sample is in the chr2_depth_test column of depth/depth_f<fraction>.tsv.

## Holdout (literature males)

DYZ3 2.73% (n 20), Yq12 1.17% (n 12), DYZ19 2.85% (n 20); all classes in holdout_summary.txt.
