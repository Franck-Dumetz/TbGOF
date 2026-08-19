#!/bin/bash

# TbGOF – Analyzing T. brucei Gain-of-Function drug-screening data 
# Copyright (C) 2025 Anushka Shome and Franck Dumetz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/.

# This script automates the analysis of ORFeome screening data. It accepts
# either SRA accession IDs or FASTQ files (gzipped or not) as input, performs quality trimming,
# alignment to a reference genome using Bowtie, BAM file generation and sorting,
# gene quantification via SeqMonk, and differential expression analysis with DESeq2.
# Both uniquely and/or multi-mapped reads can be analyzed depending on user flags.
# Differential expression results are output as Excel files with fold change filtering.


set -euo pipefail

gff=""
sras=""
fasta=""
unique=0
multiple=0
treatments=""
fastqs=""
fold=""
threads=4
force=0

usage="
Usage: ./ORF-enrich.sh -A <annotation.gff> -G <genome.fasta> -T <treatments.csv> -C <fold_change> [-R <sra_list.txt> | -F <fastq_directory>] [-u] [-m] [-p <threads>] [-f]

Required arguments:
  -A  Path to genome annotation file in GFF format.
  -G  Path to genome sequence file in FASTA format.
  -T  Path to treatments CSV file (no header, two columns: sample name, condition).
  -C  Fold change for differential expression analysis.

Input source (choose one):
  -R  Path to text file with list of SRA accession IDs.
  -F  Path to directory containing FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz).

Optional flags:
  -u  Only process uniquely aligned reads.
  -m  Only process multi-mapped reads.
      (If neither -u nor -m is used, both types will be processed.)
  -p  Number of threads to use for trimming/alignment/sorting (default: 4).
  -f  Force re-trimming and re-alignment even if a complete set of BAM
      files already exists in bam/ (default: skip if complete — useful
      when re-running after a downstream failure, e.g. DESeq2/count-reads
      erroring after alignment already finished).

Note:
  - Sample names in the treatments file must match the FASTQ filenames
    (without extension) or SRA IDs.
  - Condition names in the treatments file must start with 'treated_' or
    'untreated_'.
"

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -A) gff="$2"; shift 2 ;;
    -R) sras="$2"; shift 2 ;;
    -G) fasta="$2"; shift 2 ;;
    -u) unique=1; shift 1 ;;
    -m) multiple=1; shift 1 ;;
    -T) treatments="$2"; shift 2 ;;
    -F) fastqs="$2"; shift 2 ;;
    -C) fold="$2"; shift 2 ;;
    -p) threads="$2"; shift 2 ;;
    -f) force=1; shift 1 ;;
    *)
      echo "<<Invalid flag: $1>>"
      echo "$usage"
      exit 1
      ;;
  esac
done

# Argument checking
if [[ -z "$sras" && -z "$fastqs" ]]; then
  echo "<<No SRA ID file or FASTQ directory provided>>"
  echo "$usage"
  exit 1
fi

if [[ ! -z "$fastqs" ]]; then
  fastqs="${fastqs%/}"
  fastqs_abs="$(realpath "$fastqs")"
  mkdir -p fastqs
  # Symlink rather than copy: keeps the rest of the pipeline working off a
  # standardized local fastqs/ path without duplicating (potentially large)
  # raw read files on disk. trim_galore/bowtie2 follow symlinks transparently.
  shopt -s nullglob
  src_fastqs=("$fastqs_abs"/*.fastq "$fastqs_abs"/*.fq "$fastqs_abs"/*.fastq.gz "$fastqs_abs"/*.fq.gz)
  shopt -u nullglob
  if [[ ${#src_fastqs[@]} -eq 0 ]]; then
    echo "<<No .fastq/.fq/.fastq.gz/.fq.gz files found in $fastqs>>"
    exit 1
  fi
  ln -sf "${src_fastqs[@]}" fastqs/
fi

if [[ -z "$gff" ]]; then
  echo "<<No gff file provided>>"; echo "$usage"; exit 1
fi
if [[ -z "$fasta" ]]; then
  echo "<<No fasta file provided>>"; echo "$usage"; exit 1
fi
if [[ "$unique" -eq 0 && "$multiple" -eq 0 ]]; then
  echo "<<No alignment flags provided. Defaulting to unique and multiple alignments>>"
  unique=1
  multiple=1
fi
if [[ -z "$treatments" ]]; then
  echo "<<No treatment file provided>>"; echo "$usage"; exit 1
fi

# Preflight: confirm every external tool the pipeline calls is actually on
# PATH before doing any work. This is what would have caught the missing
# trim_galore immediately (e.g. because a conda env wasn't activated in a
# SLURM job) instead of letting the pipeline silently limp through empty
# trimmed/bam directories and print misleading "complete" messages.
missing=0
for tool in trim_galore bowtie2 bowtie2-build samtools python Rscript; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "<<Required tool not found on PATH: $tool>>"
    missing=1
  fi
done
if [[ "$missing" -eq 1 ]]; then
  echo "<<Preflight check failed. Is the correct conda/environment active? Aborting.>>"
  exit 1
fi

echo "<<Using $threads thread(s)>>"

# SRA -> FASTQ (populates & creates the fastq directory)
if [[ ! -d fastqs ]]; then
  ./fastq.sh $sras >> output.log
  echo "<<SRAs converted to FASTQs>>"
fi

# FASTQ -> Trimmed FASTQ -> sorted, indexed BAM. Trimming and alignment are
# each skipped independently if their outputs already look complete, unless
# -f forces a re-run. This lets you re-run after a downstream failure (e.g.
# count-reads.py/DESeq2 erroring, or a crash mid-alignment) without
# repeating whichever earlier step(s) already finished.
# NOTE: neither check detects changed input files — use -f if fastqs/
# contents changed since trimming/alignment were last run.
shopt -s nullglob
raw_fastqs=(fastqs/*.fastq fastqs/*.fq fastqs/*.fastq.gz fastqs/*.fq.gz)
shopt -u nullglob
if [[ ${#raw_fastqs[@]} -eq 0 ]]; then
  echo "<<No fastq files found in fastqs/>>"
  exit 1
fi
n_expected=${#raw_fastqs[@]}

mkdir -p bam
align_complete=1
if [[ "$force" -eq 1 ]]; then
  align_complete=0
else
  if [[ "$unique" -eq 1 ]]; then
    shopt -s nullglob; have_m1=(bam/*_m1.bam); shopt -u nullglob
    [[ ${#have_m1[@]} -eq "$n_expected" ]] || align_complete=0
  fi
  if [[ "$multiple" -eq 1 ]]; then
    shopt -s nullglob; have_m10=(bam/*_m10.bam); shopt -u nullglob
    [[ ${#have_m10[@]} -eq "$n_expected" ]] || align_complete=0
  fi
fi

if [[ "$align_complete" -eq 1 ]]; then
  echo "<<Complete BAM set already present in bam/ — skipping trimming and alignment (use -f to force)>>"
else

mkdir -p trimmed
shopt -s nullglob
trimmed_fastqs=(trimmed/*.fq trimmed/*.fastq trimmed/*.fq.gz trimmed/*.fastq.gz)
shopt -u nullglob

trim_complete=1
if [[ "$force" -eq 1 ]]; then
  trim_complete=0
else
  [[ ${#trimmed_fastqs[@]} -eq "$n_expected" ]] || trim_complete=0
fi

if [[ "$trim_complete" -eq 1 ]]; then
  echo "<<Complete set of trimmed FASTQs already present in trimmed/ — skipping trimming (use -f to force)>>"
else
  # trim_galore auto-detects gzip input and gzip-compresses output by default.
  trim_galore --cores "$threads" --output_dir trimmed "${raw_fastqs[@]}" >> output.log 2>&1
  echo "<<FASTQs trimmed>>"
  shopt -s nullglob
  trimmed_fastqs=(trimmed/*.fq trimmed/*.fastq trimmed/*.fq.gz trimmed/*.fastq.gz)
  shopt -u nullglob
fi

if [[ ${#trimmed_fastqs[@]} -eq 0 ]]; then
  echo "<<No trimmed fastq files found in trimmed/ to align>>"
  exit 1
fi

# Trimmed FASTQ -> sorted, indexed BAM directly (bowtie2 piped into samtools
# sort; no SAM files are ever written to disk). bowtie2 reads gzipped FASTQ
# natively, so trimmed .fq.gz files are used as-is.
bowtie2-build --threads "$threads" "$fasta" index >> output.log 2>&1

for file in "${trimmed_fastqs[@]}"; do
  base="${file##*/}"
  base="${base%_*}"   # strip trailing _trimmed(.fq|.fastq)(.gz)

  if [ "$unique" -eq 1 ]; then
    bowtie2 --local -p "$threads" -x index -U "$file" 2>> output.log \
      | samtools sort -@ "$threads" -o bam/"${base}_m1.bam" -
    samtools index -@ "$threads" bam/"${base}_m1.bam"
  fi

  if [ "$multiple" -eq 1 ]; then
    bowtie2 --local -k 10 -p "$threads" -x index -U "$file" 2>> output.log \
      | samtools sort -@ "$threads" -o bam/"${base}_m10.bam" -
    samtools index -@ "$threads" bam/"${base}_m10.bam"
  fi
done
echo "<<Bowtie2 alignment + sorted/indexed BAM files complete>>"

fi

# Trim the gff
python trim-gff.py $gff
echo "<<GFF trimmed>>"

# Count reads per gene per BAM (pysam-based, no per-gene subprocess spawning)
python count-reads.py
echo "<<Read count complete>>"

# Create the treatments.csv
if [[ "$unique" -eq 1 && "$multiple" -eq 1 ]]; then
  python make-csv.py -um $treatments
elif [ "$unique" -eq 1 ]; then
  python make-csv.py -u $treatments
else
  python make-csv.py -m $treatments
fi

# Run Deseq2
mkdir -p results
Rscript Deseq2.R $fold 1 >> output.log 2>&1
python switch-treatments.py
Rscript Deseq2.R $fold >> output.log 2>&1
python switch-treatments.py
echo "<<DESeq2 analysis complete.>>"
echo "<<Results saved in: results/foldchange_${fold}_<comparison_name>.csv>>"
echo "<<Genes with no counts are listed in: results/no-counts.csv>>"

# Cleaning all created directories
# rm -r trimmed
