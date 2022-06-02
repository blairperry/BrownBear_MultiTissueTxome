#!/usr/bin/env python2.7

# This file takes as input an Orthogroups.tsv file (output from OrthoFinder run), a species name, and the column in Orthogroups.tsv pertaining to that species, and outputs a TSV that can be used by agat_sq_add_attributes_from_tsv.pl to add orthogroup annotations to an existing GFF. This script assumes that the description field linking orthogroups to GFF entries is "Parent" but feel free to change this below if needed.

# Usage:
# python2 parseOrthogroups.py Orthogroups.tsv hSapiens 3

import sys

orthogroups = sys.argv[1]
species_name = sys.argv[2]
column = int(sys.argv[3]) - 1

outfile = orthogroups[:-3] + species_name + '.tsv'

print ""
print "Input file: " + orthogroups
print "Column with info for " + species_name + ": " + sys.argv[3]
print "Saving output to: " + outfile
print ""

gene_to_ortho = {}

with open(orthogroups) as a:
    for line in a.readlines()[1:]:
        line = line.rstrip().split('\t')
        if len(line) > column:
            if len(line[column]) > 0:
                og_id = line[0]
                for entry in line[column].split(","):
                    entry = entry.lstrip().rstrip()
                    if entry not in gene_to_ortho:
                        gene_to_ortho[entry] = og_id


with open(outfile,'w') as out:
    print >> out, '\t'.join(['Parent','Orthogroup'])
    for entry in gene_to_ortho:
        output_entry = [entry,gene_to_ortho[entry]]
        print >> out, '\t'.join(output_entry)

