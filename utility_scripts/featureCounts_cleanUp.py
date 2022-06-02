#!/usr/bin/env python2.7

## Cleans up some of the clunky and redundant info from the featureCounts output file.

## Usage: python featureCounts_cleanUp.py <path/to/featureCounts/output/file>

import sys

infile = sys.argv[1]
outfile = infile[:-4] + '.cleaned.txt'

star_bam_suffix = 'Aligned.sortedByCoord.out.bam' # suffix of star result files (i.e. everything after the sample ID, shared by all files). Change if needed

with open(outfile,'w') as out:
    with open(infile) as a:
        for line in a.readlines():
            if line[0] == '#':
                print >> out, line.rstrip()

            else:
                line = line.rstrip().split('\t')

                # Clean column headers
                if line[0] == 'Geneid':
                    star_dir = '/'.join(line[6].split('/')[:-1]) + '/' # get path to STAR bam output directory
                    new_header = []
                    for item in line:
                        if star_bam_suffix in item:
                            item = item.replace(star_dir,'').replace('Aligned.sortedByCoord.out.bam','')
                            new_header.append(item)
                        else:
                            new_header.append(item)
                    print >> out, '\t'.join(new_header)

                # Clean redundant/duplicated chr entries and simplify start, end, strand
                else:
                    if len(line[1].split(';')) >= 1:
                        line[1] = line[1].split(';')[0]
                        line[2] = line[2].split(';')[0]
                        line[3] = line[3].split(';')[-1]
                        line[4] = line[4].split(';')[0]
                        print >> out, '\t'.join(line)

