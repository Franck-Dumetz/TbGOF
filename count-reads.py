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

# This script counts the number of raw reads aligning to each gene in a genome annotation (GFF) file, using one or more input BAM files.
# It outputs a CSV file where:
# - Rows represent genes
# - Columns represent samples (each sample corresponds to a BAM file)
# - Each cell contains the raw read count for that gene in that sample
# Requires samtools to be installed in the same environment, as well as the packages in the imports above.
# Input bam files should also have index files in the same directory.

import sys
import subprocess
import csv
from collections import defaultdict
import os

def read_gff(fn, bam):
    file_counts[bam] = 0
    fields = []
    with open(fn) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if "gene_id" in fields[8]:
                gene = fields[8].split("gene_id=")[1].split(";")[0]
                count_reads(bam, fields, gene)
                des[fields[8].split("gene_id=")[1].split(";")[0]] = fields[8].split("Description=")[1].split(";")[0]
            else:
                gene = fields[8].split("ID=")[1].split(";")[0]
                count_reads(bam, fields, gene)
                des[fields[8].split("ID=")[1].split(";")[0]] = fields[8].split("Description=")[1].split(";")[0]


def count_reads(bam, fields, name):
    chrom = fields[0]
    start = int(fields[3])
    end = int(fields[4])
    #name = fields[8].split(":")[0]
    #name = fields[8].split("gene_id=")[1].split(";")[0]
    cmd = ["samtools", "view", f"bam/{bam}", f"{chrom}:{start}-{end}"]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True, check=True)
    output = result.stdout.strip()
    if not output:
        counts[bam][name] = 0
    else:
        lines = result.stdout.strip().split("\n")
        counts[bam][name] = len(lines)
        file_counts[bam] += len(lines)
    #print(f"{bam} {counts[bam][name]}")
    return counts

counts = defaultdict(dict)
des = defaultdict(dict)
file_counts = defaultdict(dict)

#bam_files = ["SRR10846669.bam", "SRR10846670.bam", "SRR10846671.bam", "SRR10846672.bam", "SRR10846673.bam", "SRR10846674.bam", "SRR10846675.bam", "SRR10846676.bam", "SRR10846677.bam", "SRR10846678.bam", "SRR10846679.bam", "SRR10846680.bam"] # change the name of the bam files here
bam_files = [f for f in os.listdir("bam") if f.endswith(".bam") and os.path.isfile(os.path.join("bam", f))]

for bam in bam_files:
    read_gff("trimmed_genome.gff", bam) # change the genome name here

#Counts is full
row_headers = sorted({row for col in counts.values() for row in col})
col_headers = sorted(counts.keys())
with open('counts.csv', 'w', newline='') as file: # change the output csv file name here
    writer = csv.writer(file)
    writer.writerow(["Gene"] + ["Description"] + [col.split('.')[0] for col in col_headers])
    for row in row_headers:
        writer.writerow([row] + [des[row]] + [counts[col].get(row, 0) for col in col_headers])
