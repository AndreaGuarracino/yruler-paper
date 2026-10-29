# Short-read model

The model that yruler uses for Illumina whole-genome reads, and its estimates for the 1000 Genomes Project (1kGP).

| File | Content |
|---|---|
| model.tsv | Settings of each class |
| ../estimates_1kgp.tsv | Estimates, Y dosage, Y state and lineage of the 1,639 1kGP samples |

## Model

For each class:

1. Reference length: the array length in the assemblies. An array is excluded only when the four assembly QC tracks (Flagger HiFi, Flagger ONT, NucFlag HiFi, NucFlag ONT) all flag the same 1 kb or more; this excludes HG01261 for Yq12 and DYZ2. YqHet2.7 is the total amount of 2.7 kb repeat sequence in Yq12, and DYZ18sub the total amount of DYZ18, 3.1 kb and 2.7 kb sequence; they are amounts, not single arrays.
2. Markers: canonical 31-mers found in at least the `presence` fraction of the class assemblies, in at most the `outside` fraction of the assemblies outside the class on Y, not in CHM13 outside Y, and below the female-haplotype filter (`female_t`). DYZ2 markers have at least 16 G+C bases of 31 (GC >= 51.6%), YqHet2.7 markers at least 8 (GC >= 25.8%).
3. Length: K x marker count corrected for G+C / count of the Y single-copy ruler. The G+C correction curve comes from the ruler of the same sample.
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
| DYZ3 | 104 | 2.74 / 2.32 / 5.20 / 16.05 | 2.73 / 2.42 / 5.13 / 15.61 | 20 | 2.74 / 2.86 / 3.91 / 6.56 |
| Yq12 | 81 | 2.01 / 1.64 / 4.18 / 8.01 | 2.03 / 1.61 / 4.15 / 7.77 | 12 | 1.65 / 1.28 / 3.83 / 5.62 |
| DYZ1 | 82 | 2.68 / 2.30 / 5.43 / 10.34 | 2.72 / 2.39 / 5.14 / 9.84 | 12 | 2.69 / 2.54 / 4.59 / 7.29 |
| DYZ2 | 81 | 2.81 / 2.30 / 5.89 / 13.43 | 2.82 / 2.40 / 5.82 / 12.98 | 12 | 2.07 / 1.78 / 3.72 / 5.94 |
| DYZ18S | 82 | 4.63 / 3.68 / 10.09 / 18.60 | 4.62 / 3.71 / 9.78 / 17.26 | 12 | 5.72 / 4.66 / 11.31 / 13.98 |
| YqHet2.7 | 82 | 6.54 / 5.33 / 13.20 / 26.09 | 6.67 / 5.79 / 13.22 / 20.00 | 12 | 6.91 / 7.66 / 11.78 / 15.32 |
| DYZ18sub | 82 | 4.77 / 4.30 / 10.01 / 18.00 | 4.83 / 4.46 / 10.18 / 15.83 | 12 | 4.56 / 3.51 / 10.29 / 13.97 |
| DYZ19 | 118 | 3.16 / 2.32 / 6.60 / 16.45 | 3.02 / 2.18 / 6.48 / 17.18 | 21 | 3.79 / 2.72 / 8.50 / 9.92 |

## estimates_1kgp.tsv

| Column | Content |
|---|---|
| sample, sex_ped, population, superpopulation | 1kGP sample, pedigree sex (1 male, 2 female), population |
| batch | 2504: original 1kGP set; 698: related samples added later |
| haplogroup, lineage | HaploGrouper call (ISOGG 2019) on the 1kGP chrY genotypes; major lineage (R1a and R1b apart) |
| y_dose | Y dosage: 2 x mean count of the Y ruler markers / mean count of the chromosome 2 markers; 1 is one Y per cell |
| DYZ3_bp ... DYZ19_bp | Estimates in bp |
| reference_male | yes: the male is in the reference set, so his estimate is not independent |
| flag | ok; female; y_dose_below_0.9 (estimates given); no_estimate (median count of the Y ruler markers below 10, or no positive ruler baseline for the G+C curve) |
| yp_ruler_dose, yq_ruler_dose, y_state | Y dosage of the Yp and Yq ruler markers; normal, mosaic_loss, mosaic_gain, partial, yp_deletion, XYY or no_Y |
| yp_deletion_bias_pct | For yp_deletion males: change of the estimates caused by the missing Yp ruler markers |

Flags: ok 1,585, female 40, y_dose_below_0.9 8, no_estimate 6. Medians of the 1,585 ok males, bp: DYZ3 858,846; Yq12 26,986,801; DYZ1 15,341,930; DYZ2 10,779,641; DYZ18S 305,846; YqHet2.7 123,192; DYZ18sub 433,913; DYZ19 330,377.
