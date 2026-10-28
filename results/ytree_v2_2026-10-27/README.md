# 1kGP Y tree v2 (2026-10-27): replaces v1 in Fig. 2b

## What changed from v1
| | v1 | v2 |
|---|---|---|
| Sites | PASS biallelic SNPs in chrY:2,781,480-56,887,902, no mask | same NYGC 30x VCF, only inside the Poznik mask (chrY_callable_regions.20130802_b38.bed, 42 intervals, 10,454,819 bp) |
| Singletons | excluded (MAC >= 2) | kept (40,273 of 79,941 sites) |
| Search | IQ-TREE -fast | full search, IQ-TREE 3.1.3, GTR2+FO+ASC, seed 20261027 |
| Support | none | UFBoot 1,000 and SH-aLRT 1,000 |
| Root | common ancestor of all 8 A (A0 + A1a mixed) | on the 3 A0 (HG01890, HG02982, HG02984): separates A0 from A1 + BT |
| Males | 1,597 (HG02300 and NA19226 without Y) | 1,597 (same) |

## Validation (scripts: ytree_validate.py; tables v2val.labels.tsv and v1val.labels.tsv)
A haplogroup label passes if its carriers, subclades included, form exactly one branch of the rooted tree.
| Check | v1 | v2 |
|---|---|---|
| Poznik 2016 labels (ISOGG 2014, 1,234 males) monophyletic | 172 of 309 (55.7%) | 258 of 309 (83.5%) |
| HaploGrouper labels (ISOGG 2019, 1,597 males) monophyletic | 813 of 1,485 (54.7%) | 1,404 of 1,485 (94.5%) |
| HG03742 (Poznik K2a1*) | sister of HG02040 (F) | sister of the NO clade (280 males: O3a 109, O2a 46, N1c 25) |
| J2a1a1 (42 males) | the smallest clade has 47 leaves | exactly 42 |
| A1 (5 males) | monophyletic | monophyletic |
| Label support (1,790 clades) | - | UFBoot >= 95 in 1,655, SH-aLRT >= 80 in 1,674 |

The labels that stay non-monophyletic are mostly deep subclades of rapid expansions (South Asian R1a-Z93, R1b-L11, E1b1a): list in v2val.labels.tsv (column intruders).

## Files
- v2val.rooted.nwk: tree rooted on A0 (for the figure); ytree.treefile: unrooted tree with SH-aLRT/UFBoot support; ytree.contree: UFBoot consensus; ytree.iqtree, ytree.log: IQ-TREE report and log.
- sites.tsv: alignment column, CHROM, POS, build GRCh38, minor allele count (REF and ALT not included: read them from the public VCF 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz).
- missing_per_sample.tsv, ids.txt: missing calls per male and sample order.
- ytree_v2.sh: full command.

## Next
Tree v3 with the Hallast 2023/2026 calls and filters (haploid BCFtools from the CRAMs, depth >= 3, allele >= 85%, sites with < 5% missing, RAxML 8.2.10) and the same validation: replaces v2 if it passes at least the same checks.
