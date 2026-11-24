# Site settings of the pipeline: the only file to edit on another cluster
# sourced by v1.sh and by every *_task.sh after P (project dir) is set
C=/path/to/tools   # edit: folder with the tools below
export PY=$C/python/bin/python   # Python 3.12 with numpy and scipy, for cv_eval.py; other scripts use python3 in PATH
export SAM=$C/samtools/1.24/bin/samtools ARIA2=$C/aria2/1.37.0/bin/aria2c KMC_BIN=$C/kmc/3.2.4/bin CLANG=$C/clang/18 CARGO=$HOME/.cargo/bin/cargo
export SPSS_BIN=$P/tools/spss/bin BUILDDEPS=$P/tools/builddeps/bin   # KmerCamel and zstd for the masked superstrings, build tools of yruler
export TMPBASE=/path/to/scratch/$USER   # edit: temporary space shared by the nodes
export PATH=$C/zstd/1.5.7/bin:$PATH
