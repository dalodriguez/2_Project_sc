

#  README

## 1. Data sources

This task is based on publicly available sequencing data from the following study:

Title: Spinal Cord Organoids from Human Amniotic Fluid iPSC Recapitulate the Diversity of Cell Phenotypes During Fetal Neural Tube Morphogenesis
DOI: 10.1007/s12035-025-04944-z

The study evaluate different human spinal cord organoids. The dataset a single cell RNA sequencing (scRNAseq) and the analysis in this document is focused on SCO-42-2_S21, originally sequenced  with Illumina NovaSeq 6000, targeting 50K reads/cell. Organoid libraries were prepared with Chromium Next GEM Single Cell 3' Reagent KIits v3.1 (Dual Index).

The subsampled and cleaned FASTQs are stored in `data/Fastq_subsampled` and are used as the inputs for the workflow.

The code to generate the response to the questions are in `workflow/3_Answer to the questions.sh`

## 2. How to download

The following samples were retrieved from Dryad: SCO-42-2_S21_L001_R1_001.fastq.gz and SCO-42-2_S21_L001_R2_001.fastq.gz. 

## 3. Pre-processing / subsampling

Sample were miniaturized using the code in 1_Subsampling, not needed for Taiga submission. 

```bash
conda activate seqtk
FILES=(*.fastq.gz)
for FILE in "${FILES[@]}"; do
    echo "Subsampling Fastq file: ${FILE}..."
    seqtk sample -s100 "${FILE}" 0.01 | gzip > Fastq_subsampled/"${FILE%.fastq.gz}_subsampled.fastq.gz"
done
```

Only 1% of the reads are conserved and the miniaturized Fastq file is stored in data/Fastq_subsampled.

The reference genome (GRCh38) is also miniaturized to contain only chromosome 10, as follows: 

```bash
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 10 > Reference/chr10.fa

wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz
awk '$1 == "10" || $1 ~ /^#/' Homo_sapiens.GRCh38.109.gtf > Reference/chr10.gtf
```
The output is sored in data/Reference. 

## 4. How the workflow works

The workflow files is stored in `workflow/`.

### Step 1 – Compile Reference Genome

**Purpose:** Prepare reference genome
**Tools:** `cellranger`
**Inputs:**  gtf and fa files of the human GRCh38 genome
**Outputs:** GRCh38_chr10 folder containing the compiled reference genome
**Command:**

```bash
cellranger mkref \
  --genome=GRCh38_chr10 \
  --fasta=Reference/chr10.fa \
  --genes=Reference/chr10.gtf
```

---

### Step 2 Aligned, quality control and quantification

**Purpose:** Perform the cellranger pipeline
**Tools:** `cellranger`
**Inputs:** Miniaturized Fastq files and compiled reference genome
**Outputs:** Processed data
**Command:**

```bash
cellranger count \
    --id="SCO-42-2_S21_analysis_ch10" \
    --fastqs=Fastq_subsampled \
    --sample="SCO422" \
    --transcriptome=GRCh38_chr10 \
    --create-bam=false
```