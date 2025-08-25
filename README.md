# TbGOF

A pipeline for the identification of differentially represented ORFs arising from the Gain-of-Function library following forward genetic screening, including drug resistance screens. 

---

Publication: 

---

## Overview

**TbGOF** is a set of scripts and tools designed to analyze gain-of-function (GoF) screens in _T. brucei_. It identifies open reading frames (ORFs) that are significantly overrepresented in treated samples compared to untreated controls, helping to discover genes associated different stresses like exposure to drug.

This pipeline was initially developed for analyzing screens in kinetoplastid parasites but can be adapted to other organisms. This pipeline revamps, add options, and automatises the data analysis from Carter M, et al. [A Trypanosoma brucei ORFeome-Based Gain-of-Function Library Identifies Genes That Promote Survival during Melarsoprol Treatment](https://journals.asm.org/doi/full/10.1128/msphere.00769-20?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org). mSphere. 2020 Oct 7;5(5):e00769-20. doi: 10.1128/mSphere.00769-20. PMID: 33028684; PMCID: PMC7568655.

---

## Key Features

- Handles multiple replicates per condition  
- Uses DESeq2 for statistical testing
- Computes differential representation of ORFs, set your own Fold Change threshold
- Test the diversity of the control culture
- Highly customizable for different datasets and experimental designs
- support macOS and Linux

---

## Dependencies
The pipeline requires the following software and R packages. All are installed automatically when using the provided conda lock files, but the list is provided here for users installing manually.

### Core Software:
bash (≥4.0) – Required to run the ORF-enrich.sh shell pipeline.
R (≥4.0.0) – Required for statistical analysis with DESeq2.

### Bioinformatics Tools:
Bowtie (1.3.1) – Aligns reads to the reference genome/ORFeome.
samtools (≥1.10) – Converts, sorts, and indexes alignment files.
SRA Toolkit (≥2.11.0) – Downloads sequencing data from NCBI when using -R.
Trimmomatic (≥0.39) – Optional, used for read trimming/quality filtering if included in your workflow.

### R Packages:
DESeq2 (≥1.30.0) – Performs differential representation analysis.
tidyverse (≥1.3.0) – Data wrangling and plotting.
optparse – Command-line argument parsing in R scripts.

### System Libraries:
libstdc++ – C++ standard library (required for Bowtie on Linux; fixed via ldlib-links.sh).

---

## Conda Environment Setup
### Clone the git repository and enter the the ORFeome directory
```
git clone https://github.com/Franck-Dumetz/TbGOF.git
cd TbGOF
```
### Linux 
```
conda create -n TbGOF --file conda-linux-64.lock
```
### MacOS
```
conda create -n TbGOF --file conda-osx-64.lock
```
### Activate the conda environment
```
conda activate TbGOF
```
```
# [Linux only] Fix potential C++ library compatibility issues (required by bowtie)
./ldlib-links.sh
```
---

## Running the Pipeline

### Usage

```
./ORF-enrich.sh -A <annotation.gff> -G <genome.fasta> -T <treatments.csv> -C <fold_change> [-R <sra_list.txt> | -F <fastq_directory>] [-u] [-m]
```

### Argument Descriptions

| Flag | Required? | Description |
|------|-----------|-------------|
| `-A` | ✅ | **Genome annotation file** in GFF format. |
| `-G` | ✅ | **Genome file** in FASTA format. |
| `-T` | ✅ | **Treatments CSV file** for grouping samples.<br><br>Must be a **two-column, headerless CSV** with:<br>• **Column 1 – Sample names**<br>&nbsp;&nbsp;• If using **SRA IDs** (`-R`): names must exactly match the SRA accession IDs.<br>&nbsp;&nbsp;• If using **local FASTQ files** (`-F`): names must exactly match the FASTQ filenames, without the `.fastq` extension (e.g., `sample1.fastq` → `sample1`).<br>• **Column 2 – Condition/library names**, starting with `treated_` or `untreated_` (e.g., `untreated_NewMP`).<br><br>An example is provided in the `test-data/` folder. |
| `-C` | ✅ | **Fold change** as an integer. |
| `-R` | One of `-R` or `-F` is required | Path to a text file containing SRA accession numbers, one per line:<br>```SRR10846669\nSRR10846670\n...``` |
| `-F` | One of `-R` or `-F` is required | Path to a directory containing `.fastq` files. |
| `-u` | Optional | Process **uniquely aligned** reads only. |
| `-m` | Optional | Process **multi-mapped** reads only.<br>If neither `-u` nor `-m` is used, both types will be processed. |

> **Note:** The order of arguments does not matter, but each flag must come before its corresponding input.
---

## Output


### Running with Test Data
Use the following command to run with the provided test data set:
```
./ORF-enrich.sh -A test-data/test-annot.gff -R test-data/test-sra.txt -G test-data/test-genome.fasta -u -T test-data/test-treatments.csv -C 4
```
### Rerunning with Different Fold Changes
Once you have successfully run the full pipeline once, you can rerun the differential expression step only using different fold-change thresholds without reprocessing all the data:
```
Rscript Deseq2.R <fold-change>
```
