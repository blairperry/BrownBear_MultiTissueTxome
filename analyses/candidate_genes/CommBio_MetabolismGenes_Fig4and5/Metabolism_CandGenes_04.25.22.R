
library(ggpmisc)
library(tidyverse)
library(pheatmap)
library(viridis)
library(colorspace)


# Read in bear gene-level counts ------------------------------------------

bear.gene_counts <- read_tsv('data/raw_counts/a_geneLevelQuants/bearNecropsy_geneLevelCounts_10.12.21_bwp.cleaned.txt',comment = '#') %>% 
  janitor::clean_names() %>% 
  dplyr::select(-start,-end,-chr,-strand) %>% 
  mutate(geneid = str_split_fixed(geneid,'[:]',3)[,2]) %>% 
  mutate(rowSum = rowSums(as.matrix(.[,-c(1,2)]))) %>%
  filter(rowSum > 0) %>%
  dplyr::select(-rowSum) %>% 
  dplyr::select(-cfa_s9,-ittfl_s44)

bear.gene_counts.mat <- as.matrix(bear.gene_counts[,-c(1,2)])
row.names(bear.gene_counts.mat) <- bear.gene_counts$geneid

bear.x <- bear.gene_counts.mat / bear.gene_counts$length
tpm.mat <- t( t(bear.x) * 1e6 / colSums(bear.x) )
bear.tpm.mat <- tpm.mat

bear.tpm.mat.df <- bear.tpm.mat %>% as.data.frame() %>% rownames_to_column(var='gene_id') %>% 
  filter(rowSums(.[,-1])>0)


# Read in sample info -----------------------------------------------------

sample_info_necrop <- readxl::read_xlsx('data/sample_info/MultiTissue_mapping_stats.xlsx') %>% 
  janitor::clean_names() %>% 
  dplyr::select(sample_id,tissue)

sample_info_prevRNA <- readxl::read_xlsx('data/sample_info/prevRNA_SampleInfo.xlsx') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = paste(tissue,phys_state,sep='_')) %>% 
  mutate(sample_id = str_to_upper(str_split_fixed(sample,'[_]',2)[,1])) %>% 
  dplyr::select(sample_id, tissue)

all.sample_info <- sample_info_necrop %>% 
  bind_rows(sample_info_prevRNA) %>% 
  filter(sample_id != 'CFA')


# Calculate within-tissue z-scores ----------------------------------------

bear.tpm.mat.df.zscores <- bear.tpm.mat.df %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'tpm') %>% 
  filter(tpm > 0) %>% 
  mutate(logTpm = log10(tpm+1)) %>% 
  group_by(sample) %>% 
  mutate(zscore = (logTpm - mean(logTpm)) / sd(logTpm))


# Read in and subset Jansen 2019 supp table to relevant genes -------------

metabGenes <- readxl::read_xlsx('analyses/candidate_genes/CommBio_MetabolismGenes_Fig4and5/Fig4and5GeneCuration_04.25.22.xlsx',col_names = T,na = 'NA') %>% 
  janitor::clean_names()

head(metabGenes)

# Read in and join StringDB mapping and KEGG characterization 

string.mapping <- read_tsv('analyses/candidate_genes/CommBio_MetabolismGenes_Fig4and5/string_mapping.tsv') %>% 
  select(queryItem, preferredName)

string.kegg <- readxl::read_xlsx('analyses/candidate_genes/CommBio_MetabolismGenes_Fig4and5/stringdb.enrichment.KEGG.xlsx') %>% 
  janitor::clean_names() %>% 
  filter(focal_pathway=='yes') %>% 
  separate_rows(matching_proteins_in_your_network_labels,sep = ',') %>% 
  left_join(string.mapping,by=c('matching_proteins_in_your_network_labels'='preferredName')) %>% 
  select(gene_id = queryItem, string_id = matching_proteins_in_your_network_labels, pathway=term_description)

metabGenes.wPathways <- metabGenes %>% 
  left_join(string.kegg,by=c('alternate_symbol'='gene_id')) %>% 
  mutate(pathway = ifelse(is.na(pathway),'Other',pathway))

# Candidate gene within sample z-scores

cand.zscores <- bear.tpm.mat.df.zscores %>% 
  filter(gene_id %in% metabGenes$bear_gene_id) %>% 
  left_join(metabGenes,by=c('gene_id'='bear_gene_id')) %>% 
  mutate(sample = str_replace_all(sample,'_hy','hy')) %>% 
  mutate(sample = str_replace_all(sample,'_hi','hi')) %>% 
  mutate(sample = str_to_upper(str_split_fixed(sample,'[_]',2)[,1])) %>% 
  left_join(all.sample_info,by=c('sample'='sample_id'))

cand.zscores.heatdata <- cand.zscores %>% 
  mutate(sample = paste(sample, tissue, sep=' ')) %>% 
  select(gene_symbol,sample,zscore) %>% 
  pivot_wider(names_from = sample,values_from = zscore,values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames('gene_symbol')

cand.zscores.annotdata <- bear.tpm.mat.df.zscores %>%
  ungroup() %>% 
  filter(gene_id %in% metabGenes$bear_gene_id) %>% 
  left_join(metabGenes.wPathways,by=c('gene_id'='bear_gene_id')) %>% 
  filter(pathway != 'Other') %>% 
  select(gene_symbol,pathway) %>% 
  unique() %>% 
  mutate(present = 1) %>% 
  pivot_wider(names_from = 'pathway',values_from = present,values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames('gene_symbol')


paletteLength <- 50
lim <-  max(abs(min(cand.zscores.heatdata,na.rm = T)),abs(max(cand.zscores.heatdata,na.rm = T))) 
myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

pheatmap(cand.zscores.heatdata,
         scale = 'none',
         color=diverge_hcl(50,palette = 'Blue-Red 3'),
         breaks=myBreaks,
         cutree_cols = 6,
         cutree_rows = 6,
         border_color = NA,
         annotation_row = cand.zscores.annotdata)

