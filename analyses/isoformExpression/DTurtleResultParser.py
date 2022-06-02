
## Script to parse through the *.err file produced by a Kamiak run of
##  dturtle on all pairwise combinations of tissues. Want to pull out the
##  number of retained genes for each analysis and number of significant
##  DTU transcripts.


results = [['comparison','num_retained','num_sig_genes','num_sig_tx']]

with open('7_DTUanalyses/R_39966666.err') as a:
    comparison = ''
    num_retained = -1
    num_sig_genes = -1
    num_sig_tx = -1


    for line in a.readlines():
        if line[0:3] == 'Com':
            comparison = line.split(': ')[1].replace('\'','').replace(' ','_').rstrip()
            # print comparison
        elif line[0:3] == 'Ret':
            num_retained = int(line.split(' ')[1])
            # print num_retained
        elif line[0:3] == 'Fou':
            num_sig_genes = int(line.split(' ')[1])
            # print num_sig_genes
            num_sig_tx = int(line.split(' ')[5])
            # print num_sig_tx
            res_list = [comparison,str(num_retained),str(num_sig_genes),str(num_sig_tx)]
            results.append(res_list)

            comparison = ''
            num_retained = -1
            num_sig_genes = -1
            num_sig_tx = -1

with open('7_DTUanalyses/parsedResultSummary_02.09.22.tsv','w') as out:
    for line in results:
        print >> out, '\t'.join(line)
