#!/bin/bash

# QUESTION 1
cd data/Fastq_subsampled
zcat SCO422_S21_L001_R1_001.fastq.gz | echo $((`wc -l`/4))

# QUESTION 1
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
cd data/SCO-42-2_S21_analysis_ch10/outs/analysis/pca/gene_expression_10_components/
conda activate python
python3 -c "
import pandas as pd
df = pd.read_csv('variance.csv')
print(df['Proportion.Variance.Explained'][0]*100)
"
conda deactivate

# QUESTION 5 
cd data/SCO-42-2_S21_analysis_ch10/outs/filtered_feature_bc_matrix
conda activate python
python3 -c "
import pandas as pd
df = pd.read_csv('barcodes.tsv.gz', sep='\t', header=None)
print(df.shape[0])
"
conda deactivate