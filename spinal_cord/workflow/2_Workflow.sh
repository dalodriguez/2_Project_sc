#!/bin/bash
cd data/

tar -xzvf cellranger-9.0.1.tar.gz
export PATH=$PATH:data/cellranger-9.0.1

cellranger mkref \
  --genome=GRCh38_chr10 \
  --fasta=Reference/chr10.fa \
  --genes=Reference/chr10.gtf

cellranger count \
    --id="SCO-42-2_S21_analysis_ch10" \
    --fastqs=Fastq_subsampled \
    --sample="SCO422" \
    --transcriptome=GRCh38_chr10 \
    --create-bam=false
    --nopreflight \