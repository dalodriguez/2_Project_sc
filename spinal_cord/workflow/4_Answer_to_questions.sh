#!/bin/bash

# QUESTION 1
cd data/Fastq_subsampled
zcat SCO422_S21_L001_R1_001.fastq.gz | echo $((`wc -l`/4))

# QUESTION 2
cd data/Reference
awk '$3=="gene" { 
    n=split($0,a,"gene_id \""); 
    if(n>1){split(a[2],b,"\""); print b[1]}
}' chr10.gtf | sort | uniq | wc -l

# QUESTION 3
cd data/SCO-42-2_S21_analysis_ch10/outs
conda activate python
python3 -c "
import pandas as pd
df = pd.read_csv('metrics_summary.csv')
print(df['Reads Mapped Confidently to Genome'][0])
"
conda deactivate

# QUESTION 4
Rscript -e "
library(Seurat)
seurat <- readRDS('data/SCO_analysis_seurat.rds')
head(VariableFeatures(seurat), 3)
"

# QUESTION 5 
Rscript -e "
library(Seurat)
seurat <- readRDS('data/SCO_analysis_seurat.rds')
markers <- FindMarkers(seurat, ident.1 = 0)
markers[which.max(markers$avg_log2FC), ]
"