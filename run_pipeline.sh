#!/usr/bin/env bash
# Snippy example: E. faecalis isolates
# Author: Zhu Zhihao  Date: 2026.10.12
# Raw reads: Genome Sequence Archive under BioProject PRJCA058049 (run accessions are available in metadata).
# Place paired FASTQs as:
#   raw/reads/{id}_1.fq.gz
#   raw/reads/{id}_2.fq.gz
# (id = first column of metadata.txt)

set -euo pipefail

cd ~/snippy
mkdir -p seq temp/fastp result/assembly ref result/snippy raw/reads

# 0. software
# conda create -y -n snippy -c bioconda -c conda-forge snippy=4.6.0 spades=4.0.0 fastp snp-dists fasttree
# conda activate snippy

# 1. link reads
cut -f1 metadata.txt | tail -n +2 | while read id; do
  ln -sfn "$(pwd)/raw/reads/${id}_1.fq.gz" seq/${id}_1.fq.gz
  ln -sfn "$(pwd)/raw/reads/${id}_2.fq.gz" seq/${id}_2.fq.gz
done

# 2. QC
cut -f1 metadata.txt | tail -n +2 | while read id; do
  mkdir -p temp/fastp/$id
  fastp -w 64 \
    -i seq/${id}_1.fq.gz -I seq/${id}_2.fq.gz \
    -o temp/fastp/$id/${id}_1.fq.gz -O temp/fastp/$id/${id}_2.fq.gz \
    -j temp/fastp/$id/${id}.json -h temp/fastp/$id/${id}.html
done

# 3. assemble
cut -f1 metadata.txt | tail -n +2 | while read id; do
  mkdir -p result/assembly/$id
  spades.py --careful -t 64 \
    -1 temp/fastp/$id/${id}_1.fq.gz -2 temp/fastp/$id/${id}_2.fq.gz \
    -o result/assembly/$id/spades
done

# 4. reference genome for Snippy: X019
ln -sfn "$(pwd)/result/assembly/X019/spades/contigs.fasta" ref/X019.fasta

# 5. Snippy
cut -f1 metadata.txt | tail -n +2 | grep -v '^X019$' | while read id; do
  snippy --cpus 64 --outdir result/snippy/$id \
    --ref "$(realpath ref/X019.fasta)" \
    --R1 temp/fastp/$id/${id}_1.fq.gz --R2 temp/fastp/$id/${id}_2.fq.gz \
    --mincov 10 --minfrac 0.9 --force
done

# 6. core SNPs
snippy-core --ref ref/X019.fasta --prefix result/core result/snippy/*

# 7. tree and pairwise distances
snippy-clean_full_aln result/core.full.aln > result/core.clean.full.aln
FastTree -gtr -nt result/core.aln > result/core.tree
snp-dists result/core.aln > result/core.snpdists.tsv
