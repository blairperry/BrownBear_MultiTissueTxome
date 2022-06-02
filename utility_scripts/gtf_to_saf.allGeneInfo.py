#!/usr/bin/env python2.7

## Converts a GTF file to a simplified annotation file (SAF) for use with featureCounts.

## Usage: python gtf_to_saf.py <path/to/input/GTF> <feature>
## Assumes that file conforms to standard GTF format with description fields in the 9th column.
## feature : feature from 3rd column of GTF that featureCounts will use to count read overlap (ex., exon)

## Output example:
## GeneID Chr Start End Strand
## gene12345_497097_SymbolA chr1 3204563 3207049 -
## gene12345_497097_SymbolA chr1 3411783 3411982 -
## gene12345_497097_SymbolA chr1 3660633 3661579 -
## gene5678_100503874_SymbolB chr1 3637390 3640590 -
## gene5678_100503874_SymbolB chr1 3648928 3648985 -
## gene5678_100503874_SymbolB chr1 3670236 3671869 -
## *** In this example, each line represents an exon entry (feature) and is grouped by the gene_id, geneID, and locus (meta-feature)
##      that you want featureCounts to use to measure expression.

## Will output a .SAF file with the same basename as the input GTF.

import sys

infile = sys.argv[1]
feature = sys.argv[2]

header = ['GeneID','Chr','Start','End','Strand']

outfile = infile[:-4] + '.allGeneInfo.saf'

with open(outfile,'w') as out:
    print >> out, '\t'.join(header)

    with open(infile) as a:
        for line in a.readlines():
            if line[0] != '#':
                line=line.rstrip().split('\t')
                if line[2] == feature:
                    #print line
                    chr = line[0]
                    start = line[3]
                    end = line[4]
                    strand = line[6]
                    descrip = line[8].split('; ')
                    meta_entry1 = [entry for entry in descrip if 'GeneID' in entry]
                    if len(meta_entry1) > 0:
                        meta_entry1 = meta_entry1[0].replace('GeneID:','')
                        meta_id1 = meta_entry1.split(' ')[1].replace('"', '')
                    else:
                        meta_id1 = 'NA'

                    meta_entry2 = [entry for entry in descrip if 'gene_id' in entry][0]
                    meta_id2 = meta_entry2.split(' ')[1].replace('"', '')

                    meta_entry3 = [entry for entry in descrip if 'gene "' in entry]
                    if len(meta_entry3) > 0:
                        meta_entry3=meta_entry3[0]
                    else:
                        meta_entry3 = [entry for entry in descrip if 'product "' in entry][0]
                    meta_id3 = meta_entry3.split(' ')[1].replace('"', '')

                    meta_output = ':'.join([meta_id1,meta_id2,meta_id3])
                    outlist = [meta_output,chr,start,end,strand]
                    print >> out, '\t'.join(outlist)

