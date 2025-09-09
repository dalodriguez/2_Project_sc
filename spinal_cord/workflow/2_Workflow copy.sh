#!/bin/bash
cd data/Fastq_subsampled/

gunzip SCO422_S21_L001_R1_001.fastq.gz SCO422_S21_L001_R2_001.fastq.gz

cd data/

tar -xzvf cellranger-9.0.1.tar.gz
export PATH=$PATH:cellranger-9.0.1

cellranger mkref \
  --genome=GRCh38_chr10 \
  --fasta=Reference/chr10.fa \
  --genes=Reference/chr10.gtf

cellranger count \
    --id="SCO-42-2_S21_analysis_ch102" \
    --fastqs=Fastq_subsampled \
    --sample="SCO422" \
    --transcriptome=GRCh38_chr10 \
    --create-bam=false \
    --nopreflight