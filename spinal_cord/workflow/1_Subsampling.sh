#!/bin/bash

# DO NOT RUN
cd /2_Project_sc/spinal_cord/data/Fastq

# SUBSAMPLE FASTQ FILES
conda activate seqtk
FILES=(*.fastq.gz)
for FILE in "${FILES[@]}"; do
    echo "Subsampling Fastq file: ${FILE}..."
    seqtk sample -s100 "${FILE}" 0.01 | gzip > "${FILE%.fastq.gz}_subsampled.fastq.gz"
done

zcat SCO422_S21_L001_R1_001.fastq.gz | echo $((`wc -l`/4)) "reads"

# DOWNLOAD CELLRANGER
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1756536938&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=HHvEZxFLvSf~rDMB8lYiDM9BYIwegNZOPc6kOvfkCDXkH0oav3QtmgL26Tz-wKMTMDL1csVncAlFGxzScwo0-NuDcup1MPVyur3DMpjSxYW~-kpd6g2Y-K708D1V5w2JBCwhlrnovJZ5sJl0cImJER8GiVnZi2G1Bs6-LcMhrhbNVottFlI-19Wf65kfK~hkDujvrj602aqgIZ9woHnTSKS6-rYRTHS45Kk8UAjj~voIwO0YJjwy7zyoRI2VVMqLTQfjsBIbGAK9dG2G2W6e3ZWJK7WgJcJOKrw2uJuXEHs142DSj~4RYyxV0hYCECNTysOzBF6hPeKXkSj5UGSpJw__"
tar -xzvf cellranger-9.0.1.tar.gz
export PATH=$PATH:/mnt/c/Users/dalod/Desktop/Projects/7_SepalAI/2_Project_sc/sample_submission_boilerplate/data/cellranger-9.0.1
cellranger --version

# DOWNLOAD FULL REFERENCE GRCh38 GENOME
#wget -O refdata-cellranger-GRCh38-2020-A.tar.gz "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
#tar -xzvf refdata-cellranger-GRCh38-2020-A.tar.gz

# DOWNLOAD FULL REFERENCE GRCh38 GENOME AND SUBSET CHROMOSOME 10
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 10 > chr10.fa

wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz
awk '$1 == "10" || $1 ~ /^#/' Homo_sapiens.GRCh38.109.gtf > chr10.gtf