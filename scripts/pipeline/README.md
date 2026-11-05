# Model pipeline

These scripts build the released yruler model (marker k-mers and `model.tsv`) and its cross-validation from reference assemblies and short reads. They ran on a Slurm cluster. Copy them to `$P/scripts/` of a project directory and run them from `$P`.

## Environment

- `P`: project directory; `SCRATCH`: local temporary space.
- `YRULER`: yruler binary (default `yruler` in `PATH`); `YRULER_SRC`: yruler repository, for the release build.
- `KMC_BIN`: KMC 3.2.4 binaries; `SPSS_BIN`: masked-superstring tools; `CLANG`: clang 18 prefix for the release build.
- `PY`: Python 3 with pandas, for `um_halfilt.sh`.
- In `PATH`: samtools 1.24, bcftools, zstd, aria2c, agc, impg, MAFFT 7.526, HMMER 3.4 (`nhmmer`, `hmmbuild`, `hmmemit`), Python 3 with numpy.

## Inputs (not in this folder)

- Reference arrays: Y assemblies of Hallast et al. 2026 and HPRC release 2, as gzipped FASTA in `$P/ks_truth/fa/<assembly>.fa.gz`.
- Satellite annotation of Hallast et al. 2026: their motifs in `$P/hallast_motifs/<motif>.fasta` (`DYZ1_Yq`, `DYZ2_Con`, `DYZ18_Yq`, `Yqhet_2k7bp`, `Yqhet_3k1bp`, `DYZ19_Yq`) and their filter scripts `filterYq12HMMER.py` and `filterDYZ19HMMER.py` (https://github.com/Markloftus/PanYChromosome), which `um_halfilt.sh` calls as `scripts/filterYq12_hal.py` and `scripts/filterDYZ19_hal.py` with `-i <hit table folder> -o <csv>`.
- DYZ17 seed: `$P/hallast_motifs/DYZ17_ours.fasta`, the consensus (`hmmemit -c`) of a first DYZ17 profile HMM built from repeat units of the assemblies.
- Array QC tables (flanks, gaps, quality-track flags): `dyz19_clean.tsv` and `dyz17_clean.tsv`, as in `data/method/` of this repository.
- Truth exclusions: `$P/truth_qc/x/rule_w4.txt` (class and male per line), as `data/method/truth_qc/x/rule_w4.txt`.
- Reads: 1000 Genomes high-coverage CRAM files and GRCh38 (GRCh38_full_analysis_set_plus_decoy_hla.fa); 1000 Genomes pedigree `$P/lists/1kGP.3202_samples.pedigree_info.txt`.

## Steps

1. Truth lengths per class and male:
   - `um_hal_nhmmer.sh`: nhmmer with the Hallast et al. motifs and E-values, one task per assembly and motif; the two DYZ17 motifs use `EVAL=1.60E-150`.
   - `um_halfilt.sh`: Hallast et al. filters on the hit tables (`yq12`, `yq12all`, `dyz19`).
   - `um_d17unit.py`: DYZ17 profile HMM from start-to-start repeat units of the `DYZ17_ours` hits (MAFFT, `hmmbuild`, `hmmemit`), written as `hallast_motifs/DYZ17_unit.fasta`; then `um_hal_nhmmer.sh` again with `DYZ17_unit`. It takes the assembly list from the output of `um_newtruth.py`.
   - `um_halregions.py <old marker folder> <new marker folder>`: writes `indep_setB.tsv` (Hallast et al. annotation of the Yq12 classes, the same table as `data/method/yq12_truth_hallast.tsv`), `dyz19_clean.tsv`, `dyz17_truth.tsv` (DYZ17 length: longest block of `DYZ17_unit` hits joined when less than 10 kb apart) and `regions.tsv`.
   - `truth_cl.py` (library) reads these tables; set `D17T` to `dyz17_truth.tsv` and `D17C=dyz17_span`, as `um_v9bBx.sh` does.
   - `um_newtruth.py`: DYZ17 arrays from CHM13 31-mers, not the DYZ17 length of the release. `indep_truth.py`: Yq12 class amounts from our own profile-HMM hits. `um_prodtruth.py`: all-copies truth (`TRUTHSET`).
2. Candidate k-mers: `um_kmc.sh` and `um_kmc_out.sh` per assembly, `um_super.py` (class k-mers, minus non-Y CHM13), `um_bits.py`, `um_v5pos.py` (positions on Y, autosomal window).
3. Read counts per sample: `pass_cram9.sh` (with `marker_reads.sh`), `keep_dl.sh` and `blockfix.py` for downloads, `um_v5count.sh`, `um_cnpy.py`.
4. Model and accuracy: `um_v9bBx.sh` runs `um_v5run.sh` with the final settings; it calls `um_eval.py` through `um_cvsplit.sh` (repeated 5-fold cross-validation and lineage-out tests), `um_model.py`, `um_apply.py`, `um_xtech.py` (other sequencing technologies, same-male pairs across projects, parent-offspring pairs), `um_holdx.sh` (held-out males), `um_v5sum.py` (accuracy tables) and `um_v9export.py` (tables for the figures).
5. Release: `um_relchain9.sh` builds yruler and runs `um_release.py` (markers and `model.tsv`), `release_pack.sh` (masked superstrings) and `um_relcheck.py` (yruler output against the Python model).
