
library(DTUrtle)
library(tidyverse)

#the BiocParallel framework is used to parallelize the computations.
#Using 12 cores:
biocpar <- BiocParallel::MulticoreParam(12)
#or standard serial computation (only 1 core)
# biocpar <- BiocParallel::SerialParam()
#multiple other options available for computational clusters.

load('dTurtleInputData_forKamiak_02.04.22')

for (row in 1:nrow(all.combos)) {
  tissue1 <- all.combos[row, 1]
  tissue2  <- all.combos[row, 2]
  
  analysis <- paste(str_replace_all(tissue1,' ','_'),str_replace_all(tissue2,' ','_'),sep = '_vs_')
  
  dturtle.res <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd.necropsy, id_col = "id",
                                         cond_col = 'group', cond_levels = c(tissue1,tissue2), 
                                         filtering_strategy = "bulk", 
                                         BPPARAM = biocpar,
                                         filter_only=F)
  
  dturtle.res <- posthoc_and_stager(dturtle = dturtle.res, ofdr = 0.05)
  
  dturtle.res <- create_dtu_table(dturtle = dturtle.res)
  
  write_tsv(dturtle.res$dtu_table,paste('dturtleRes_',analysis,'_DTUTable.tsv',sep=''))
  write_tsv(as.data.frame(dturtle.res$sig_gene),paste('dturtleRes_',analysis,'_SigGenes.tsv',sep=''))
  write_tsv(as.data.frame(dturtle.res$sig_tx),paste('dturtleRes_',analysis,'_SigTx.tsv',sep=''))
}

