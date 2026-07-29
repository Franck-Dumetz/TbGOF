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

# This script processes a GFF genome annotation file to retain only CDS (coding sequence) entries.
# For each CDS, it trims the coordinates to keep only the first 100 bp of the region.
# - For '+' strand genes, it keeps positions from start to start+99 (or to the original end if shorter)
# - For '-' strand genes, it keeps positions from end-99 to end (or from start if shorter)
# The resulting trimmed GFF is written to a new file.

import sys
from collections import defaultdict

# Change these file names
inp = sys.argv[1]
out = "trimmed_genome.gff"
des = defaultdict(dict)

with open(inp, 'r') as infile, open(out, 'w') as outfile:
    for line in infile:
        if line.startswith("#") or not line.strip(): # if the line is a comment or an empty line, write it to the new gff file
            outfile.write(line)
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9: #if the line is incorrectly formatted, skip over it
            continue

        # GFF files have 9 fields: chromosome, source, type, start amino acid, end amino acid, alignment score, strand, phase, attributes
        if fields[2] == "mRNA":
            gene_id = fields[8].split("ID=")[1].split(";")[0]
            description = fields[8].split("description=")[1].split(";")[0]
            des[gene_id] = description
            
        if fields[2] != "CDS":
            continue

        start = int(fields[3])
        end = int(fields[4])

        # + are the forward sense strands (add 99 to start) and - are backwards (antisense) strands (subtract 99 from start)
        if fields[6] == "+":
            e = start + 99
            if e > end:
                e = end
            s = start
        else:
            s = end - 99
            if s < start:
                s = start
            e = end
        fields[3] = str(s)
        fields[4] = str(e)
        parent = fields[8].split("ID=")[1].split(";")[0]
        if "Parent=" in fields[8]:
            parent = fields[8].split("Parent=")[1].split(";")[0]
        if parent in des:
            fields[8] = f"{fields[8]};Description={des[parent]}"
        else:
            name = fields[8].split("Name=")[1].split(";")[0]
            fields[8] = f"{fields[8]};Description={name}"

        outfile.write("\t".join(fields) + "\n")
