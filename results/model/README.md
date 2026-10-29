# Short-read model

The model that yruler uses for Illumina whole-genome reads, and its estimates for the 1000 Genomes Project (1kGP).

| File | Content |
|---|---|
| model.tsv | Settings of each class |
| depth/depth_f<fraction>.tsv | Cross-validation predictions with the test reads thinned to the fraction (1: all reads); columns as in the 10-repeat CV, plus the reference and estimated bp and the chromosome 2 marker depth of the test run |
| depth/summary.tsv | Error statistics per fraction and class |
| lineage_out.tsv, lineage_out_summary.tsv | One prediction per male with his lineage left out of training, and the statistics per class and lineage |
| ../estimates_1kgp.tsv | Estimates, Y dosage, Y state and lineage of the 1,639 1kGP samples |

## Model

For each class:

1. Reference length: the array length in the assemblies. An array is excluded only when the four assembly QC tracks (Flagger HiFi, Flagger ONT, NucFlag HiFi, NucFlag ONT) all flag the same 1 kb or more; this excludes HG01261 for Yq12 and DYZ2. YqHet2.7 is the total amount of 2.7 kb repeat sequence in Yq12, and DYZ18sub the total amount of DYZ18, 3.1 kb and 2.7 kb sequence; they are amounts, not single arrays.
2. Markers: canonical 31-mers found in at least the `presence` fraction of the class assemblies, in at most the `outside` fraction of the assemblies outside the class on Y, not in CHM13 outside Y, and below the female-haplotype filter (`female_t`). DYZ2 markers have at least 16 G+C bases of 31 (GC >= 51.6%), YqHet2.7 markers at least 8 (GC >= 25.8%).
3. Length: K x marker count corrected for G+C / count of the Y single-copy ruler. The G+C correction curve comes from the ruler of the same sample: ruler counts are pooled per G+C bin (counts above 10 x (mean + 1) left out), at any depth; it needs at least 1,000,000 ruler counts, and each bin at least 1,000 markers and 10,000 counts.
4. Yq12: Yq12 - its DYZ2 part (DYZ2 markers without a G+C floor, class `yq12dyz2`) + DYZ2.

`presence`, `outside` and K come from the reference males; cross-validation and lineage-out fit them again in each training set. The G+C floor of each class comes from the same males sequenced in two projects, with leave-one-out.

## model.tsv

| Column | Content |
|---|---|
| class | Marker class |
| presence, outside | Marker selection fractions (step 2) |
| ruler | Y single-copy ruler (y90: Y k-mers in at least 90% of the Y assemblies) |
| female_t | Female-haplotype filter |
| markers, ysinglecopy_markers | Number of class and ruler markers |
| K | Scale factor, bp |
| chr2 | False: the class is normalised to the Y ruler, not to chromosome 2 |
| min_gc | Lowest number of G+C bases of 31 in a marker |
| length | Reported length: the class itself, a signed sum of classes, or `none` (helper class, not reported) |

## Accuracy

Absolute error in %, mean / median / p90 / max. Cross-validation: 5 folds, 10 repeats, family groups; the median is the mean of the 10 repeat medians. Lineage-out: each of 18 lineages left out of training in turn. Fixed holdout: males with published DYZ3 or DYZ1 measurements (Miga 2014, Altemose 2014), with their relatives out of training; not used to choose any setting.

| Class | Males | Cross-validation | Lineage-out | Holdout males | Holdout |
|---|---|---|---|---|---|
| DYZ3 | 104 | 2.74 / 2.31 / 5.20 / 16.04 | 2.74 / 2.33 / 5.14 / 15.60 | 20 | 2.74 / 2.85 / 3.90 / 6.57 |
| Yq12 | 81 | 2.01 / 1.63 / 4.19 / 8.03 | 2.02 / 1.60 / 4.06 / 7.78 | 12 | 1.61 / 1.16 / 3.71 / 5.40 |
| DYZ1 | 82 | 2.67 / 2.28 / 5.42 / 10.34 | 2.70 / 2.28 / 5.14 / 9.94 | 12 | 2.69 / 2.54 / 4.58 / 7.29 |
| DYZ2 | 81 | 2.82 / 2.31 / 5.95 / 13.50 | 2.81 / 2.42 / 5.77 / 13.05 | 12 | 1.90 / 1.06 / 3.44 / 6.21 |
| DYZ18S | 82 | 4.63 / 3.69 / 10.09 / 18.61 | 4.63 / 3.71 / 9.78 / 17.27 | 12 | 5.71 / 4.65 / 11.30 / 13.98 |
| YqHet2.7 | 82 | 6.53 / 5.36 / 13.22 / 26.09 | 6.66 / 5.79 / 13.22 / 21.49 | 12 | 6.91 / 7.64 / 11.80 / 15.31 |
| DYZ18sub | 82 | 4.76 / 4.29 / 10.01 / 16.44 | 4.83 / 4.46 / 10.17 / 15.83 | 12 | 4.56 / 3.51 / 10.28 / 13.97 |
| DYZ19 | 118 | 3.16 / 2.32 / 6.59 / 16.43 | 3.03 / 2.18 / 6.48 / 17.17 | 21 | 3.78 / 2.72 / 8.48 / 9.91 |

## estimates_1kgp.tsv

| Column | Content |
|---|---|
| sample, sex_ped, population, superpopulation | 1kGP sample, pedigree sex (1 male, 2 female), population |
| batch | 2504: original 1kGP set; 698: related samples added later |
| haplogroup, lineage | HaploGrouper call (ISOGG 2019) on the 1kGP chrY genotypes; major lineage (R1a and R1b apart) |
| y_dose | Y dosage: 2 x mean count of the Y ruler markers / mean count of the chromosome 2 markers; 1 is one Y per cell |
| DYZ3_bp ... DYZ19_bp | Estimates in bp |
| reference_male | yes: the male is in the reference set, so his estimate is not independent |
| flag | ok; female; y_dose_below_0.9 (estimates given); no_estimate (fewer than 1,000,000 Y ruler counts: no G+C curve, no Y) |
| yp_ruler_dose, yq_ruler_dose, y_state | Y dosage of the Yp and Yq ruler markers; normal, mosaic_loss, mosaic_gain, partial, yp_deletion, XYY or no_Y |
| yp_deletion_bias_pct | For yp_deletion males: change of the estimates caused by the missing Yp ruler markers |

Flags: ok 1,585, female 40, y_dose_below_0.9 12, no_estimate 2. Medians of the 1,585 ok males, bp: DYZ3 859,326; Yq12 26,984,924; DYZ1 15,341,647; DYZ2 10,781,178; DYZ18S 305,859; YqHet2.7 122,634; DYZ18sub 434,029; DYZ19 330,357.
