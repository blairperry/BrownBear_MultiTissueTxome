#!/usr/bin/env python2.7

# This script takes as input a GFF and the output from parseOrthogroups.py. For all exon entries in the GFF, it will add orthogroup information into the attribute field. All non-exon entries will be printed as-is.

# Usage:
# python2 annotateGFF.py [path/to/GFF] [path/to/annotation/tsv] [path/to/output/file]

import sys

in_gff = sys.argv[1]
in_tsv = sys.argv[2]
outfile = sys.argv[3]

att_dict = {}
key_field = ''
new_att = ''

with open(in_tsv) as a:
    for i,line in enumerate(a.readlines()):
        if i == 0:
            key_field = line.rstrip().split('\t')[0]
            new_att = line.rstrip().split('\t')[1]
        else:
            key_id = line.rstrip().split('\t')[0]
            og_id = line.rstrip().split('\t')[1]
            att_dict[key_id] = og_id

with open(outfile,'w') as out:
    with open(in_gff) as a:
        for line in a.readlines():
            if '\texon\t' in line:
                line = line.rstrip().split('\t')
                full_att = line[8].split(';')
                key_orig = [x for x in full_att if key_field in x][0].split('=')[1]
                if key_orig in att_dict:
                    og_id = att_dict[key_orig]
                else:
                    og_id = 'NA'
                new_att_entry = new_att + '=' + og_id
                full_att.append(new_att_entry)
                full_att = ';'.join(full_att)
                new_line = '\t'.join(line[:8]) + '\t' + full_att
                print >> out, new_line
            else:
                print >> out, line.rstrip()

