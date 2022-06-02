
library(tidyverse)
library(igraph)
library(pheatmap)

res.summary <- read_tsv('analyses/isoformExpression/7_DTUanalyses/parsedResultSummary_02.09.22.tsv') %>% 
  mutate(tissue1 = str_split_fixed(comparison,'_vs_',2)[,1],
         tissue2 = str_split_fixed(comparison,'_vs_',2)[,2])

### Heatmap of number of DTU transcripts between each tissue

res.sigTx <- res.summary %>% 
  select(tissue1,tissue2,num_sig_tx)

g <- graph.data.frame(res.sigTx, directed=FALSE)

# add value as a weight attribute
res.sigTx.matrix <- as.data.frame(get.adjacency(g, attr='num_sig_tx', sparse=FALSE))
res.sigTx.matrix[res.sigTx.matrix == 0] <- NA

row.names(res.sigTx.matrix) <- str_replace_all(row.names(res.sigTx.matrix),'_',' ')
names(res.sigTx.matrix) <- str_replace_all(names(res.sigTx.matrix),'_',' ')


pheatmap(res.sigTx.matrix,
         cellwidth = 20,
         cellheight = 20,
         border_color = 'black',
         cluster_rows = T,
         cluster_cols = T,
         cutree_rows = 4,
         cutree_cols = 4,
         main='Total num sig DTU tx',
         color = viridis::viridis(20),
         na_col = 'white')


### Read in and merge all result tables

files <- list.files(pattern = '*.DTUTable.tsv', full.names = TRUE,path = 'analyses/isoformExpression/7_DTUanalyses/')
all_data <- map_df(files, ~read.delim(.x,sep = '\t') %>% mutate(comparison = basename(.x))) %>% 
  select(1:9) %>% 
  filter(!is.na(gene_ID)) %>% 
  mutate(comparison = str_remove_all(comparison,'dturtleRes_'),
         comparison = str_remove_all(comparison,'_DTUTable.tsv'),
         tissue1 = str_split_fixed(comparison,'_vs_',2)[,1],
         tissue2 = str_split_fixed(comparison,'_vs_',2)[,2],
         genesWOneOrMorSigTx = ifelse(number_significant_tx > 0, T, F))

# Get number of genes with at least one significant DTU transcript (gene q < 0.05 and one or more tx q < 0.05)
all_data.numNonZeroTx <- all_data %>% 
  group_by(tissue1, tissue2, genesWOneOrMorSigTx) %>% 
  tally() %>% 
  filter(genesWOneOrMorSigTx == T) %>% 
  select(tissue1,tissue2,n)

g2 <- graph.data.frame(all_data.numNonZeroTx, directed=FALSE)

# add value as a weight attribute
res.numNonZeroTx.matrix <- as.data.frame(get.adjacency(g2, attr='n', sparse=FALSE))
res.numNonZeroTx.matrix[res.numNonZeroTx.matrix == 0] <- NA

row.names(res.numNonZeroTx.matrix) <- str_replace_all(row.names(res.numNonZeroTx.matrix),'_',' ')
names(res.numNonZeroTx.matrix) <- str_replace_all(names(res.numNonZeroTx.matrix),'_',' ')

pheatmap(res.numNonZeroTx.matrix,
         cellwidth = 20,
         cellheight = 20,
         border_color = 'black',
         color = viridis::viridis(20),
         cluster_rows = T,
         cluster_cols = T,
         cutree_rows = 4,
         cutree_cols = 4,
         main='Genes w/ one or more sig DTU tx',
         na_col = 'white')


