# TbGOF

A pipeline for the identification of differentially represented ORFs arising from the Gain-of-Function library following forward genetic screening, including drug resistance screens. 

---

Publication: 

---

## Overview

**TbGOF** is a set of scripts and tools designed to analyze gain-of-function (GoF) screens in _T. brucei_. It identifies open reading frames (ORFs) that are significantly overrepresented in treated samples compared to untreated controls, helping to discover genes associated different stresses like exposure to drug.

The pipeline was validated using sequencing data from Carter et al. (2020) [A Trypanosoma brucei ORFeome-Based Gain-of-Function Library Identifies Genes That Promote Survival during Melarsoprol Treatment](https://journals.asm.org/doi/full/10.1128/msphere.00769-20?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org) (mSphere. 2020 Oct 7;5(5):e00769-20. doi: 10.1128/mSphere.00769-20. PMID: 33028684; PMCID: PMC7568655), with both unique and multiple alignments assessed. Using the TriTrypDB _T. brucei_ reference genome (v52, approximate) and a fold change threshold of 4, the pipeline reproduced the original analysis, identifying 56 overrepresented genes, of which 54 overlapped with the 57 reported in the original study.

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

**Core Software**:
- bash (≥4.0) – Required to run the ORF-enrich.sh shell pipeline.
- R (≥4.0.0) – Required for statistical analysis with DESeq2.

**Bioinformatics Tools**:
- Bowtie (≥2.5.5) – Aligns reads to the reference genome/ORFeome.
- samtools (≥1.10) – Converts, sorts, and indexes alignment files.
- SRA Toolkit (≥2.11.0) – Downloads sequencing data from NCBI when using -R.
- trim_galore (≥2.3.0) – Used for read trimming/quality filtering.

**R Packages**:
- DESeq2 (≥1.30.0) – Performs differential representation analysis.
- tidyverse (≥1.3.0) – Data wrangling and plotting.
- optparse – Command-line argument parsing in R scripts.

---

## Conda Environment Setup
### Clone the git repository and enter the the ORFeome directory
```bash
git clone https://github.com/Franck-Dumetz/TbGOF.git
cd TbGOF
```
### Set up the conda environment
```
conda create --file environment.yml
```
### Activate the conda environment
```bash
conda activate TbGOF
```

---

## Running the Pipeline

### Usage

```bash
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

### Running with Test Data
Use the following command to run with the provided test data set:
```bash
./ORF-enrich.sh -A test-data/test-annot.gff -R test-data/test-sra.txt -G test-data/test-genome.fasta -u -T test-data/test-treatments.csv -C 4
```
### Rerunning with Different Fold Changes
Once you have successfully run the full pipeline once, you can rerun the differential expression step only using different fold-change thresholds without reprocessing all the data:
```bash
Rscript Deseq2.R <fold-change>
```
---

## Output

### Low-count gene reports
Three CSV files listing genes with zero counts:
- `no_counts_all.csv` – Genes with zero counts in all samples.
- `no_counts_treated.csv` – Genes with zero counts in all treated samples.
- `no_counts_untreated.csv` – Genes with zero counts in all untreated samples.

### Differential expression results
For each treated vs. untreated comparison, the pipeline generates a CSV containing:
- Normalized counts for each replicate
- Average normalized counts for each condition
- Fold change
- Adjusted p-value

Only genes with an adjusted p-value < 0.05 are included in these files.

For example, if there are two treated samples and two untreated samples, the pipeline produces:
- Four CSVs identifying genes overexpressed in the treated sample.
- Four CSVs identifying genes overexpressed in the untreated sample.

### Common overexpressed genes
Two additional CSV files identify genes that are consistently overexpressed across all comparisons:
- Genes overexpressed in every treated-vs-untreated comparison.
- Genes overexpressed in every untreated-vs-treated comparison.

Example output files are provided in the `example-results/` folder. These demonstrate the expected CSV formats generated by the pipeline. These are the results from running the test command given above. 

#### ⚠️ Troubleshooting
If the `results/` folder is empty, missing, or only partially populated, the pipeline may have failed silently. Check `output.log` for details on any errors. Common causes include insufficient disk space, which can prevent intermediate files from being created and cause DESeq2 analysis to fail.

---

## Terminal Output

When you run the pipeline, progress messages will be printed to the terminal. These indicate each major stage of the workflow.  

Example run with the test dataset:  

```bash
(TbGOF) computer1 TbGOF % ./ORF-enrich.sh -A test-data/test-annot.gff -R test-data/test-sra.txt -G test-data/test-genome.fasta -u -T test-data/test-treatments.csv -C 4
<<SRAs converted to FASTQs>>
<<FASTQs trimmed>>
<<Bowtie complete>>
<<BAM files made>>
<<GFF trimmed>>
<<Read count complete>>
<<DESeq2 analysis complete.>>
<<Results saved in: results/foldchange_4_<comparison_name>.csv>>
<<Genes with no counts are listed in: results/no-counts.csv>>
```

---

## Intermediate Files

The pipeline also saves files generated during intermediate steps:

- **FASTQ files from SRA** (if `-R` is used)  
- **Trimmed FASTQ files**  
- **BAM files**  
- **Raw read counts** in `counts.csv`  

To save space, **SAM files are removed** after conversion, but all their information is preserved in the corresponding BAM files.

---

## Supported GFF3 annotation formats

This pipeline currently supports two GFF3 annotation formats:

- Ribosomal profiling annotation
    - Gene identifier is taken from the ID attributes.
    - Gene description is taken from the Name attribute.
- TriTrypDB annotation
    - Gene identifier is taken from the gene_id attribute (or ID when gene_id is unavailable).
    - Gene description is taken from the description attribute of the corresponding mRNA feature.

Other GFF3 formats may require modifications to the attribute-parsing step.

---

## Contact & Support
For any issues or questions, please open an issue in this repository or contact us at fdumetz@som.umaryland.edu.
