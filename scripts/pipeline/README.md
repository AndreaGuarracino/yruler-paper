# Model pipeline

These scripts build the released yruler model (marker k-mers and `model.tsv`) and its cross-validation from reference assemblies and short reads. They ran on a Slurm cluster. Copy them to `$P/scripts/` of a project directory and run them from `$P`.

## Environment

- `P`: project directory; `SCRATCH`: local temporary space.
- `YRULER`: yruler binary (default `yruler` in `PATH`); `YRULER_SRC`: yruler repository, for the release build.
- `KMC_BIN`: KMC 3.2.4 binaries; `SPSS_BIN`: masked-superstring tools; `CLANG`: clang 18 prefix for the release build.
- In `PATH`: samtools 1.24, bcftools, zstd, aria2c, agc, impg, Python 3 with numpy.

## Inputs (not in this folder)

- Reference arrays: Y assemblies of Hallast et al. 2026 and HPRC release 2, with the satellite annotation of Hallast et al. 2026 (their filter scripts are in their repository).
- Reads: 1000 Genomes high-coverage CRAM files and GRCh38 (GRCh38_full_analysis_set_plus_decoy_hla.fa).

## Steps

1. Truth lengths per class and male: `truth_cl.py` (library), `indep_truth.py`, `um_newtruth.py` (DYZ17), `um_prodtruth.py`.
2. Candidate k-mers: `um_kmc.sh` and `um_kmc_out.sh` per assembly, `um_super.py` (class k-mers, minus non-Y CHM13), `um_bits.py`, `um_v5pos.py` (positions on Y, autosomal window).
3. Read counts per sample: `pass_cram9.sh` (with `marker_reads.sh`), `keep_dl.sh` and `blockfix.py` for downloads, `um_v5count.sh`, `um_cnpy.py`.
4. Model and accuracy: `um_v9bBx.sh` runs `um_v5run.sh` with the final settings; it calls `um_eval.py` through `um_cvsplit.sh` (repeated 5-fold cross-validation and lineage-out tests), `um_model.py`, `um_apply.py`, `um_holdx.sh` (held-out males), `um_v5sum.py` (accuracy tables) and `um_v9export.py` (tables for the figures).
5. Release: `um_relchain9.sh` builds yruler and runs `um_release.py` (markers and `model.tsv`), `release_pack.sh` (masked superstrings) and `um_relcheck.py` (yruler output against the Python model).
