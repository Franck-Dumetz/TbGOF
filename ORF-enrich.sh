#!/bin/bash

# ORFeome – Analyzing ORFeome screening data
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
# either SRA accession IDs or FASTQ files as input, performs quality trimming,
# alignment to a reference genome using Bowtie, BAM file generation and sorting,
# gene quantification via SeqMonk, and differential expression analysis with DESeq2.
# Both uniquely and/or multi-mapped reads can be analyzed depending on user flags.
# Differential expression results are output as Excel files with fold change filtering.

gff=""
sras=""
fasta=""
unique=0
multiple=0
treatments=""
fastqs=""
fold=""

usage="
Usage: ./ORF.sh -A <annotation.gff> -G <genome.fasta> -T <treatments.csv> [-R <sra_list.txt> | -F <fastq_directory>] [-u] [-m]

Required arguments:
  -A  Path to genome annotation file in GFF format.
  -G  Path to genome sequence file in FASTA format.
  -T  Path to treatments CSV file (no header, two columns: sample name, condition).
  -C  Fold change for differential expression analysis.

Input source (choose one):
  -R  Path to text file with list of SRA accession IDs.
  -F  Path to directory containing FASTQ files.

Optional flags:
  -u  Only process uniquely aligned reads.
  -m  Only process multi-mapped reads.
      (If neither -u nor -m is used, both types will be processed.)

Note:
  • Sample names in the treatments file must match the FASTQ filenames or SRA IDs.
  • Condition names in the treatments file must start with 'treated_' or 'untreated_'.
"

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -A)
      gff="$2"
      shift 2
      ;;
    -R)
      sras="$2"
      shift 2
      ;;
    -G)
      fasta="$2"
      shift 2
      ;;
    -u)
      unique=1
      shift 1
      ;;
    -m)
      multiple=1
      shift 1
      ;;
    -T)
      treatments="$2"
      shift 2
      ;;
    -F)
      fastqs="$2"
      shift 2
      ;;
    -C)
      fold="$2"
      shift 2
      ;;
    *)
      echo "<<Invalid flag: $1>>"
      echo "$usage"
      exit 1
      ;;
  esac
done

#Argument checking
if [[ -z "$sras" && -z "$fastqs" ]]; then
  echo "<<No SRA ID file or FASTQ directory provided>>"
  echo "$usage"
  exit 1
fi

if [[ ! -z "$fastqs" ]]; then
  fastqs="${fastqs%/}"
  mkdir fastqs
  cp "${fastqs}"/* fastqs/
fi

if [[ -z "$gff" ]]; then
  echo "<<No gff file provided>>"
  echo "$usage"
  exit 1
fi

if [[ -z "$fasta" ]]; then
  echo "<<No fasta file provided>>"
  echo "$usage"
  exit 1
fi

if [[ "$unique" -eq 0 && "$multiple" -eq 0 ]]; then
  echo "<<No alignment flags provided. Defaulting to unique and multiple alignments>>"
  unique=1
  multiple=1
fi

if [[ -z "$treatments" ]]; then
  echo "<<No treatment file provided>>"
  echo "$usage"
  exit 1
fi

#SRA -> FASTQ (populates & creates the fastq directory)
if [[ ! -d fastqs ]]; then
  ./fastq.sh $sras >> output.log
  echo "<<SRAs converted to FASTQs>>"
fi

#FASTQ -> Trimmed FASTQ (populates & creates the trimmed directory)
mkdir trimmed
trim_galore --output_dir trimmed fastqs/* >> output.log 2>&1
echo "<<FASTQs trimmed>>"

#Trimmed FASTQ -> SAM (populates & creates the sam directory)
mkdir sam

#bowtie-build $fasta index >> output.log
bowtie2-build $fasta index >> output.log 2>&1

for file in trimmed/*.fq; do
  basename="${file##*/}"
  basename="${basename%_*}"
  if [ "$unique" -eq 1 ]; then
    bowtie2 --local -x index -U $file -S sam/"${basename}_m1.sam" 2>> output.log
  fi
  if [ "$multiple" -eq 1 ]; then
    bowtie2 --local -k 10 -x index -U $file -S sam/"${basename}_m10.sam" 2>> output.log
  fi
  #bowtie2 --local -x index -U $file -S sam/"${basename}_m10.sam" 2>> output.log
done
#bowtie --best --strata -t -v 2 -a -m 10 -S index trimmed/
echo "<<Bowtie complete>>"

# SAM -> BAM
mkdir bam

for file in sam/*; do
  basename="${file##*/}"
  basename=${basename%.*}
  samtools view -S -b $file > bam/"${basename}_unsorted.bam" 2>> output.log
  samtools sort -o bam/"${basename}.bam" bam/"${basename}_unsorted.bam" 2>> output.log
  samtools index bam/"${basename}.bam" 2>> output.log
done

rm bam/*_unsorted.bam
echo "<<BAM files made>>"

#Trim the gff
python trim-gff.py $gff
echo "<<GFF trimmed>>"

#Run seq-monk (generates counts.csv)
python count-reads.py
echo "<<Read count complete>>"

#Create the treatments.csv
if [[ "$unique" -eq 1 && "$multiple" -eq 1 ]]; then
  python make-csv.py -um $treatments
elif [ "$unique" -eq 1 ]; then
  python make-csv.py -u $treatments
else
  python make-csv.py -m $treatments
fi

#echo "Treatments file created"

#Run Deseq2
mkdir results
Rscript Deseq2.R $fold 1 >> output.log 2>&1
python switch-treatments.py
Rscript Deseq2.R $fold >> output.log 2>&1
python switch-treatments.py
echo "<<DESeq2 analysis complete.>>"
echo "<<Results saved in: results/foldchange_${fold}_<comparison_name>.csv>>"
echo "<<Genes with no counts are listed in: results/no-counts.csv>>"

#Cleaning all created directories
#rm -r fastqs
#rm -r trimmed
rm -r sam
#rm -r bam
