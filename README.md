# yruler paper

Manuscript, figures and analysis tables of the yruler paper: lengths of the satellite arrays of the human Y chromosome measured from sequencing reads. The tool is at https://github.com/guarracinolab/yruler.

## Contents

- `manuscript/`: LaTeX source (`main.tex` with `brief-text.tex` and `methods.tex`), bibliography and figure PDFs.
- `results/model/`: per-sample estimates of the final model (`cohort_1kgp.tsv`, `cohort_hgdp_sgdp.tsv`), cross-validation predictions, depth tests and the tables behind the numbers in the text (`text_numbers.tsv`, `methods_numbers.tsv`, `biology_tests.tsv`, `pgls_composition.tsv`).
- `results/checks/`: checks cited in the Methods (relatedness, nanopore cross-validation, sequencing-protocol pairs and others).
- `results/cohorts/`, `results/tree/`, `results/tree_biallelic/`: cohort tables and the Y SNP trees.
- `data/`: reference-array annotations, marker regions and metadata used by the figure scripts.
- `scripts/`: Python scripts that make the figures and tables from the files above; `scripts/pipeline/`: the cluster scripts that build the model and its cross-validation from assemblies and reads.

## Rebuild

Python 3 with numpy, scipy and matplotlib. Run from the repository root.

```shell
R=results/model
python3 scripts/biology_tests.py $R data/method/indep_ours10.tsv data/method/yq12_truth_hallast.tsv $R/biology_tests.tsv $R/dunn_pairs.tsv
python3 scripts/manuscript_figures.py      # Extended Data and Supplementary Fig. 2 into manuscript/figures
python3 scripts/figure1.py                 # Fig. 1 (panels a and b from manuscript/figures/main/fig1/*.svg)
python3 scripts/figure2.py                 # Fig. 2 (layout of manuscript/figures/main/fig2/Fig2_Yruler_v3.svg, data from the tables)
python3 scripts/pedigree_plots.py          # parent-offspring figure
python3 scripts/suppfig1_markerfree.py     # Supplementary Fig. 1
python3 scripts/text_numbers.py $R $R/text_numbers.tsv
python3 scripts/methods_numbers.py $R results/checks/gc_curves.tsv $R/methods_numbers.tsv
python3 scripts/pgls_composition.py $R/cohort_1kgp.tsv $R/pgls_composition.tsv
python3 scripts/relatedness_checks.py $R results/checks/relatedness/gnomad_related_pairs.tsv results/checks/relatedness.tsv results/checks/relatedness/sgdp_king_pairs.tsv
cd manuscript && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

The per-sample estimates come from yruler runs on the sequencing data listed in the Data availability section of the manuscript. Those runs, the marker selection and the cross-validation used a computing cluster; `scripts/pipeline/README.md` lists the steps.
