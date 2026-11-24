# Model pipeline

These scripts build the yruler model of the paper (marker 31-mers and `model.tsv`), its cross-validation and the cohort tables. They ran on a Slurm cluster. The `um_` scripts and the download scripts make the inputs; the other scripts go from the inputs to the results.

## Setup

- Copy the scripts to `$P/v1/pipeline/` of a project directory `$P`, put the inputs in `$P/v1/inputs/` (groups below) and run every command from `$P`.
- Edit `config.sh`, the only file with site settings: tool paths and the temporary space shared by the nodes.
- Tools: Python 3.12 with numpy and scipy, samtools 1.24, aria2, KMC 3.2.4, zstd, KmerCamel, clang 18 and cargo (release step), a yruler source checkout in `$P/yruler` (or `YSRC`). The input scripts also use bcftools, agc, impg, MAFFT 7.526 and HMMER 3.4.

## Inputs

| group | holds |
|---|---|
| `markers/` | presence of each universe 31-mer in the 156 reference assemblies (`bits/`), class regions (Hallast et al. 2026 blocks for DYZ1, DYZ2, DYZ18, Yq12), universe keys, pangenome masks, non-Y genome counts, Yp/Yq blocks of the ruler |
| `truth/` | reference lengths in the assemblies (DYZ3, Yq12, DYZ19, DYZ17, repeat units) and the truth exclusions |
| `counts/` | 31-mer counts of the cross-validation males, the nanopore males, thinned reads, females, region slices of the 1kGP, HGDP, SGDP and female CRAMs, the chr8 window, CRAMs counted whole |
| `windows/` | chr8:91-101 Mb window keys, X and PAR rulers, X/PAR window means, chr8 depth per coverage level |
| `lists/` | pedigree, lineages, sex and library metadata, Y states, external and batch-holdout males, nanopore males, mixture pairs |
| `reference/` | GRCh38 with decoys and HLA, the DRAGEN 3.7.6 region BED |

## Steps

Each step submits Slurm jobs:

```shell
O=$PWD/v1/results
OUT=$O bash v1/pipeline/v1.sh model      # model: 14 classes, DYZ2 with the G+C floor 14
OUT=$O bash v1/pipeline/v1.sh cv         # 10 x 5-fold cross-validation
OUT=$O bash v1/pipeline/v1.sh lolo       # lineages left out in turn
OUT=$O bash v1/pipeline/v1.sh holdout    # train on the 2,504 batch, test on the 698 batch
OUT=$O bash v1/pipeline/v1.sh range      # lengths outside the training range
OUT=$O bash v1/pipeline/v1.sh depth      # 9 coverage levels
OUT=$O bash v1/pipeline/v1.sh nanopore   # nanopore cross-validation and model
OUT=$O bash v1/pipeline/v1.sh cohorts    # cohort tables and pairs
OUT=$O bash v1/pipeline/v1.sh release    # markers and model.tsv for yruler, Illumina and nanopore
OUT=$O bash v1/pipeline/v1.sh mixtures   # male and female read mixtures
OUT=$O bash v1/pipeline/v1.sh export     # tables for the figure scripts
```

`model` runs before `cohorts`, `cohorts` before `release`, and `release` before `mixtures`.

## Scripts

| file | does |
|---|---|
| `v1.sh` | settings of the model and submission of each step |
| `config.sh` | site settings |
| `cv_task.sh`, `cv_eval.py` | one class per Slurm task: cross-validation, lineage-out, batch holdout, length range, production markers |
| `model.py` | joins the per-class production markers into `model.npz` and `model.tsv` |
| `summary.py` | error tables of the cross-validation runs, the depth table and the figure tables |
| `cohort_task.sh`, `cohort.py`, `cohort_add.py` | cohort tables from the slice counts, rows of the CRAMs counted whole, pangenome X/PAR doses |
| `pairs.py` | same individual in two projects, parent-offspring pairs |
| `release_task.sh`, `release.py` | markers and `model.tsv` for yruler, masked superstrings, yruler build, checks on HG00126 |
| `mixtures_task.sh`, `mixtures.py` | male and female read mixtures run through yruler |
| `markers_task.sh`, `markers.py` | per-assembly 31-mer sets and their presence over the universe keys (`markers/bits`) |
| `truth.py` | reference lengths read by the other scripts |

## Scripts that make the inputs

They read the Y assemblies of Hallast et al. 2026 and HPRC release 2 (gzipped FASTA in `$P/ks_truth/fa/<assembly>.fa.gz`), the satellite motifs and filter scripts of Hallast et al. 2026 (`$P/hallast_motifs/`), the array QC tables of `data/method/` in this repository, and the 1000 Genomes high-coverage CRAM files.

1. Truth lengths per class and male:
   - `um_hal_nhmmer.sh`: nhmmer with the Hallast et al. motifs and E-values, one task per assembly and motif; the two DYZ17 motifs use `EVAL=1.60E-150`.
   - `um_halfilt.sh`: Hallast et al. filters on the hit tables (`yq12`, `yq12all`, `dyz19`).
   - `um_d17unit.py`: DYZ17 profile HMM from start-to-start repeat units of the `DYZ17_ours` hits (MAFFT, `hmmbuild`, `hmmemit`), written as `hallast_motifs/DYZ17_unit.fasta`; then `um_hal_nhmmer.sh` again with `DYZ17_unit`. It takes the assembly list from the output of `um_newtruth.py`.
   - `um_halregions.py <old marker folder> <new marker folder>`: writes `indep_setB.tsv` (Hallast et al. annotation of the Yq12 classes, the same table as `data/method/yq12_truth_hallast.tsv`), `dyz19_clean.tsv`, `dyz17_truth.tsv` (DYZ17 length: longest block of `DYZ17_unit` hits joined when less than 10 kb apart) and `regions.tsv`.
   - `um_newtruth.py`: DYZ17 arrays from CHM13 31-mers. `indep_truth.py`: Yq12 class amounts from our own profile-HMM hits. `um_prodtruth.py`: all-copies truth.
2. Candidate 31-mers: `markers_task.sh` per assembly, `um_super.py` (class 31-mers, minus non-Y CHM13), `um_v5pos.py` (positions on Y, Yp/Yq blocks of the ruler).
3. Read counts per sample: `pass_cram9.sh` (with `marker_reads.sh`), `keep_dl.sh` and `blockfix.py` for downloads, `um_v5count.sh`, `um_cnpy.py`.
