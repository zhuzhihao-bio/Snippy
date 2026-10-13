# Snippy example: clonal structure within host 22M21

Example workflow using [Snippy](https://github.com/tseemann/snippy) (v4.6.0) to assess within-host clonal structure of *Enterococcus faecalis* isolates from host **22M21**, with **X019** as the reference.

## Contents

| Path | Description |
|------|-------------|
| `run_pipeline.sh` | End-to-end pipeline (QC → assemble → Snippy → core → tree/distances) |
| `metadata.txt` | Sample IDs, host, species, role (reference / sample) |
| `Results/core.txt` | snippy-core summary |
| `Results/core.snpdists.tsv` | Pairwise core SNP distances |
| `Results/core.tree` | FastTree Newick tree |
| `Figures/heatmap.pdf` / `.tif` | Core SNP distance heatmap |
| `Figures/tree.pdf` / `.tif` | Core SNP phylogeny |
| `Rscript/plot_heatmap.R` | R script to redraw the heatmap |

## Requirements

- conda / mamba
- `snippy=4.6.0`, `spades=4.0.0`, `fastp`, `snp-dists`, `fasttree` (Bioconda)
- R with `ggplot2`, `scales` (for plotting)

```bash
conda create -y -n snippy -c bioconda -c conda-forge \
  snippy=4.6.0 spades=4.0.0 fastp snp-dists fasttree
conda activate snippy
```

## Data availability

Raw sequencing reads are **not** included. Please contact the author to request access. Place paired FASTQs as:

```text
raw/reads/{id}_1.fq.gz
raw/reads/{id}_2.fq.gz
```

where `{id}` is the first column of `metadata.txt`.

## Quick start

```bash
# 1. Put metadata.txt and run_pipeline.sh in ~/snippy (or edit paths)
# 2. Place reads under raw/reads/
bash run_pipeline.sh
```

Key steps: fastp → SPAdes (`--careful`) → Snippy (`--mincov 10 --minfrac 0.9`, skip self-vs-X019) → `snippy-core` → FastTree → snp-dists.

## Expected result (this example)

Three clones within host 22M21:

- **near-ref**: X019, Y1620, Y1621, Y1624, Y1626 (2–5 SNPs)
- **Cluster A**: Y093, Y097, Y101, Y104 (8–33 SNPs)
- **Cluster B**: Y1464, Y1615 (102 SNPs)

Between-clone distances are ~22,000–24,000 SNPs, supporting **multi-clone coexistence**.

## Plotting

```bash
# from the folder that contains Results/core.snpdists.tsv
Rscript Rscript/plot_heatmap.R Figures
```

Tree figure in `Figures/` was produced from `Results/core.tree` (e.g. with ape / iTOL).

## Citation

- Snippy: Seemann T. https://github.com/tseemann/snippy
- This example: Zhu Zhihao (contact for reads and questions)

## License

Pipeline scripts and example outputs in this repository: MIT. Snippy itself is GPL-2.0.
