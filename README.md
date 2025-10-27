# eSLAC — Enhanced Single-read Analysis of tRNA Crosstalks

**eSLAC** is a reproducible pipeline that pairs **MSR-seq**–based single-read evidence with robust statistics to quantify:
- modification–modification crosstalk (e.g., **Ψ–Ψ**) and
- **modification–charging** crosstalk
at **isodecoder** and **anticodon** resolutions in human cytosolic and mitochondrial tRNAs.

> Paper: *Interconnected landscape of tRNA pseudouridine modification and charging revealed by enhanced single-read analysis of tRNA crosstalks* (under review).

## Highlights
- Deterministic build of the **tRNA reference** with recorded provenance and checksums.
- Exact **commands** for alignment → BAM → TSV (event tables).
- Transparent **coverage rules**, FDR control, and effect-size filters.
- Outputs at **isodecoder** and **anticodon** levels + publication-ready plots.

## Quick start
```bash
conda env create -f environment.yml
conda activate eslac

# Configure paths/thresholds
cp config/config.yaml.example config/config.yaml
# (edit config/config.yaml, see comments inside)

# Run from FASTQ to tables/plots (Snakemake or Makefile provided)
snakemake -j 8 --use-conda
# or:
make all -j8
