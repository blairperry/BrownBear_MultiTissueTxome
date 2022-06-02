# BiocManager::install("biomaRt",force = T)
# install.packages('tidyverse')
# install.packages('ggbeeswarm')

library(biomaRt)
library(patchwork)
library(ggbeeswarm)
library(tidyverse)
library(cols4all)

############################################################################################################

## Housekeeping genes downloaded from https://www.tau.ac.il/~elieis/HKG/
## Read in and use BioMart to get ensembl IDs

# hks <- read_tsv('analyses/candidate_genes/housekeeping/HK_genes.txt',col_names = c('symbol','refseq'))
# ensembl <- useMart('ensembl',dataset='hsapiens_gene_ensembl')
# 
# hks.ensembl <- getBM(attributes = c('hgnc_symbol','refseq_mrna','ensembl_gene_id'),filters = 'refseq_mrna',values = hks$refseq,mart=ensembl)
# 
# write_tsv(hks.ensembl,'analyses/candidate_genes/housekeeping/HK_genes.ensemblIDs.txt')


############################################################################################################

## Read in housekeeping ids, ortholog info, and combine

hks <- read_tsv('analyses/candidate_genes/housekeeping/HK_genes.ensemblIDs.txt')
og.info <- read_tsv('analyses/orthoFinder/id_conversion/HumanBearOrthologs_wIdConvert_10.20.21.tsv') %>% 
  mutate(humanGeneID = str_split_fixed(humanGeneID,'[.]',2)[,1])

hks.og <- hks %>% 
  left_join(og.info,by=c('ensembl_gene_id'='humanGeneID')) %>% 
  filter(!is.na(orthogroup1))

length(unique(hks.og$ensembl_gene_id))
sum(unique(hks.og$ensembl_gene_id) %in% hks$ensembl_gene_id) / length(hks$ensembl_gene_id) # ~88% of HKs have OG info


sample_info_necrop <- readxl::read_xlsx('data/sample_info/MultiTissue_mapping_stats.xlsx') %>% 
  janitor::clean_names() %>% 
  dplyr::select(sample_id,tissue)

sample_info_prevRNA <- readxl::read_xlsx('data/sample_info/prevRNA_SampleInfo.xlsx') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = paste(tissue,phys_state,sep='_')) %>% 
  mutate(sample_id = str_to_upper(str_split_fixed(sample,'[_]',2)[,1])) %>% 
  dplyr::select(sample_id, tissue)

all.sample_info <- sample_info_necrop %>% 
  bind_rows(sample_info_prevRNA)


############################################################################################################

## Read in gene expression data - raw, TPM, and TMM normalized

### RAW

bear.raw.mat <- read_tsv('analyses/rnaseq_normalization/bear_oglevel_rawCounts.txt')
human.raw.mat <- read_tsv('analyses/rnaseq_normalization/human_oglevel_rawCounts.txt')

bear.raw.df <- bear.raw.mat %>% 
  dplyr::rename(og_id=orthogroup1) %>% 
  filter(str_detect(og_id,'OG')) %>% 
  dplyr::select(-og_length,-ittfl_s44,-cfa_s9)

human.raw.df <- human.raw.mat %>% 
  dplyr::rename(og_id=orthogroup1) %>% 
  filter(str_detect(og_id,'OG')) %>% 
  dplyr::select(-og_length)

sum(bear.raw.df$og_id %in% human.raw.df$og_id)

combined.raw.df <- bear.raw.df %>% 
  left_join(human.raw.df,by='og_id') %>% 
  drop_na()

hk.raw <- combined.raw.df %>% 
  mutate(is_hk = ifelse(og_id %in% hks.og$orthogroup1,T,F))

combined.raw.forPlot <- hk.raw %>% 
  pivot_longer(-c(1,286),names_to = 'sample',values_to = 'raw_count') %>% 
  mutate(species = ifelse(str_detect(sample,'GTEX'),'human','bear')) %>% 
  mutate(sample = str_replace_all(sample,'_hy','hy')) %>% 
  mutate(sample = str_replace_all(sample,'_hi','hi')) %>% 
  mutate(sample = str_replace_all(sample,'_fa','fa')) %>% 
  mutate(sample_id = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  mutate(sample_id = ifelse(species=='human',sample,sample_id)) %>% 
  left_join(all.sample_info,by='sample_id') %>% 
  mutate(tissue = ifelse(species=='human',str_split_fixed(sample,'[_]',2)[,1],tissue)) 


### TPM

bear.tpm.mat <- read_tsv('analyses/rnaseq_normalization/bear_oglevel_tpmCounts.txt')
human.tpm.mat <- read_tsv('analyses/rnaseq_normalization/human_oglevel_tpmCounts.txt')


bear.tpm.df <- bear.tpm.mat %>% 
  dplyr::rename(og_id=rowname) %>% 
  filter(str_detect(og_id,'OG')) %>% 
  dplyr::select(-ittfl_s44,-cfa_s9)


human.tpm.df <- human.tpm.mat %>% 
  dplyr::rename(og_id=rowname) %>% 
  filter(str_detect(og_id,'OG'))

combined.tpm.df <- bear.tpm.df %>% 
  left_join(human.tpm.df) %>% 
  drop_na()

hk.tpms <- combined.tpm.df %>% 
  mutate(is_hk = ifelse(og_id %in% hks.og$orthogroup1,T,F))

combined.tpm.forPlot <- hk.tpms %>% 
  pivot_longer(-c(1,286),names_to = 'sample',values_to = 'tpm_count') %>% 
  mutate(species = ifelse(str_detect(sample,'GTEX'),'human','bear')) %>% 
  mutate(sample = str_replace_all(sample,'_hy','hy')) %>% 
  mutate(sample = str_replace_all(sample,'_hi','hi')) %>% 
  mutate(sample = str_replace_all(sample,'_fa','fa')) %>% 
  mutate(sample_id = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  mutate(sample_id = ifelse(species=='human',sample,sample_id)) %>% 
  left_join(all.sample_info,by='sample_id') %>% 
  mutate(tissue = ifelse(species=='human',str_split_fixed(sample,'[_]',2)[,1],tissue)) 


### TMM

tmms <- read_tsv('analyses/rnaseq_normalization/combined_oglevel_tmmCounts.txt') %>% 
  dplyr::select(-ittfl_s44)

hk.tmms <- tmms %>% 
  mutate(is_hk = ifelse(og_id %in% hks.og$orthogroup1,T,F))

combined.tmm.forPlot <- hk.tmms %>% 
  pivot_longer(-c(1,286),names_to = 'sample',values_to = 'tmm_count') %>% 
  mutate(species = ifelse(str_detect(sample,'GTEX'),'human','bear')) %>% 
  mutate(sample = str_replace_all(sample,'_hy','hy')) %>% 
  mutate(sample = str_replace_all(sample,'_hi','hi')) %>% 
  mutate(sample = str_replace_all(sample,'_fa','fa')) %>% 
  mutate(sample_id = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  mutate(sample_id = ifelse(species=='human',sample,sample_id)) %>% 
  left_join(all.sample_info,by='sample_id') %>% 
  mutate(tissue = ifelse(species=='human',str_split_fixed(sample,'[_]',2)[,1],tissue)) 

combined.tmm.forPlot


### Combine all

combined.allCounts.forPlot <- combined.tmm.forPlot %>% 
  left_join(combined.raw.forPlot) %>% 
  left_join(combined.tpm.forPlot)


###

combined.tmm.forPlot.active <- combined.tmm.forPlot %>% 
  filter((species == 'bear' & str_detect(tissue,'Active')) | (species == 'human')) %>% 
  group_by(og_id,is_hk,species,tissue) %>% 
  summarise(median_count = median(tmm_count),
            avg_count = mean(tmm_count)) %>% 
  group_by(og_id) %>% 
  mutate(tissue = str_split_fixed(tissue,'[_]',2)[,1]) %>% 
  filter(tissue %in% c('Adipose','Liver','Muscle'))

combined.tmm.forPlot.active.v2 <- combined.tmm.forPlot.active %>% 
  dplyr::select(og_id,is_hk,species,tissue,median_count) %>% 
  pivot_wider(names_from = species,values_from = median_count)

combined.tmm.forPlot.active.goodGenes <- combined.tmm.forPlot.active.v2 %>% 
  filter(bear > 0 & human > 0)


p1 <- ggplot(combined.tmm.forPlot.active.goodGenes,aes(x=human,y=bear,color=is_hk)) +
  facet_wrap(~tissue) +
  geom_point(data=subset(combined.tmm.forPlot.active.goodGenes,is_hk=='FALSE'),alpha=0.2) + 
  geom_point(data=subset(combined.tmm.forPlot.active.goodGenes,is_hk=='TRUE'),alpha=0.3) + 
  geom_abline(slope=1,intercept = 0,lty=2) +
  scale_color_manual(values = c('FALSE'='#E69F00','TRUE'='#0072B2')) +
  scale_y_log10() +
  labs(x='Human - Median Expression per Orthogroup\n(log10 scaled axis)',y='Bear - Median Expression per Orthogroup\n(log10 scaled axis)') +
  scale_x_log10() +
  theme_linedraw()
p1

p2 <- ggplot(combined.tmm.forPlot.active.v2,aes(x=log2(human / bear),color=is_hk)) +
  facet_wrap(~tissue) +
  geom_density(lwd=1) +
  geom_vline(xintercept = 0,lty=2) +
  scale_color_manual(values = c('FALSE'='#E69F00','TRUE'='#0072B2')) +
  xlab('Log2(Human Expression / Bear Expression)\n (Median Orthogroup-level Counts)') +
  theme_linedraw()

p1 / p2


## Comparing standard deviation across tissues in both species

combined.tmm.stDevs.bySpecies <- combined.tmm.forPlot.active.v2 %>% 
  group_by(og_id) %>% 
  mutate(bear_isZero = sum(bear==0),human_isZero = sum(human==0)) %>% 
  filter(bear_isZero < 3 & human_isZero < 3) %>% 
  dplyr::select(-bear_isZero,-human_isZero) %>% 
  mutate(bear = log10(bear+1),human=log10(human+1)) %>% 
  mutate(bear_stdev = sd(bear),human_stdev = sd(human)) %>% 
  dplyr::select(og_id,is_hk,bear=bear_stdev,human=human_stdev) %>% 
  unique() %>% 
  pivot_longer(-c(1,2),names_to = 'species',values_to = 'stdev')

hk.percUnder1.human <- nrow(combined.tmm.stDevs.bySpecies[which(combined.tmm.stDevs.bySpecies$species=='human' & combined.tmm.stDevs.bySpecies$is_hk==T & combined.tmm.stDevs.bySpecies$stdev< 1),])  / nrow(combined.tmm.stDevs.bySpecies[which(combined.tmm.stDevs.bySpecies$species=='human' & combined.tmm.stDevs.bySpecies$is_hk==T),1]) * 100
hk.percUnder1.bear <- nrow(combined.tmm.stDevs.bySpecies[which(combined.tmm.stDevs.bySpecies$species=='bear' & combined.tmm.stDevs.bySpecies$is_hk==T & combined.tmm.stDevs.bySpecies$stdev< 1),])  / nrow(combined.tmm.stDevs.bySpecies[which(combined.tmm.stDevs.bySpecies$species=='bear' & combined.tmm.stDevs.bySpecies$is_hk==T),1]) * 100

p.separate <- ggplot(combined.tmm.stDevs.bySpecies,aes(x=is_hk,y=stdev,fill=is_hk)) +
  # geom_jitter(alpha=0.1) +
  geom_boxplot(outlier.alpha = 0.5,show.legend = F) +
  # ggrepel::geom_text_repel(data=subset(combined.tmm.stDevs.bySpecies,(species=='bear' & is_hk & stdev > 1)),aes(label=og_id)) +
  geom_hline(yintercept = 1,lty=2) +
  facet_wrap(~species,nrow = 1) +
  scale_y_continuous(limits = c(0,3.2)) +
  labs(x='Is Housekeeping Orthogroup?',y='standard deviation (log10(tmm count + 0.1)',title = 'Standard deviation across tissues (human and bear separate)') +
  theme_linedraw() + theme(panel.grid = element_blank())

p.human <-  ggplot(subset(combined.tmm.stDevs.bySpecies,species=='human'),aes(x=is_hk,y=stdev,fill=is_hk)) +
  # geom_jitter(alpha=0.1) +
  geom_boxplot(outlier.alpha = 0.5,show.legend = F) +
  # ggrepel::geom_text_repel(data=subset(combined.tmm.stDevs.bySpecies,(species=='bear' & is_hk & stdev > 1)),aes(label=og_id)) +
  scale_fill_manual(values = c('FALSE'='#E69F00','TRUE'='#0072B2')) +
  annotate('text',label=paste('Percent HKOs\nwith SD < 1: ',format(round(hk.percUnder1.human, 2), nsmall = 2),'%',sep=''),x=T,y=2.5) +
  geom_hline(yintercept = 1,lty=2) +
  scale_y_continuous(limits = c(0,3.2)) +
  labs(x='Is Housekeeping Orthogroup?',y='standard deviation (log10(tmm count + 0.1)',title = 'Human') +
  theme_linedraw() + theme(panel.grid = element_blank())


  
p.bear <-  ggplot(subset(combined.tmm.stDevs.bySpecies,species=='bear'),aes(x=is_hk,y=stdev,fill=is_hk)) +
  # geom_jitter(alpha=0.1) +
  geom_boxplot(outlier.alpha = 0.5,show.legend = F) +
  # ggrepel::geom_text_repel(data=subset(combined.tmm.stDevs.bySpecies,(species=='bear' & is_hk & stdev > 1)),aes(label=og_id)) +
  geom_hline(yintercept = 1,lty=2) +
  scale_fill_manual(values = c('FALSE'='#E69F00','TRUE'='#0072B2')) +
  annotate('text',label=paste('Percent HKOs\nwith SD < 1: ',format(round(hk.percUnder1.bear, 2), nsmall = 2),'%',sep=''),x=T,y=2.5) +
  scale_y_continuous(limits = c(0,3.2)) +
  labs(x='Is Housekeeping Orthogroup?',y='standard deviation (log10(tmm count + 0.1)',title = 'Bear') +
  theme_linedraw() + theme(panel.grid = element_blank())

combined.tmm.stDevs.combined <- combined.tmm.forPlot.active.v2 %>% 
  group_by(og_id) %>% 
  mutate(bear_isZero = sum(bear==0),human_isZero = sum(human==0)) %>% 
  filter(bear_isZero < 3 & human_isZero < 3) %>% 
  dplyr::select(-bear_isZero,-human_isZero) %>% 
  pivot_longer(-c(1,2,3),names_to = 'species',values_to = 'tmmcount') %>% 
  mutate(stdev = sd(log10(tmmcount+0.1))) %>% 
  dplyr::select(og_id,is_hk,stdev) %>% 
  unique() 

hk.percUnder1.combined <- nrow(combined.tmm.stDevs.combined[which(combined.tmm.stDevs.combined$is_hk==T & combined.tmm.stDevs.combined$stdev< 1),])  / nrow(combined.tmm.stDevs.combined[which(combined.tmm.stDevs.combined$is_hk==T),1]) * 100


p.combined <- ggplot(combined.tmm.stDevs.combined,aes(x=is_hk,y=stdev,fill=is_hk)) +
  geom_boxplot(outlier.alpha = 0.5,show.legend = F) +
  # geom_violin(trim = T,) +
  geom_hline(yintercept = 1,lty=2) +
  scale_fill_manual(values = c('FALSE'='#E69F00','TRUE'='#0072B2')) +
  annotate('text',label=paste('Percent HKOs\nwith SD < 1: ',format(round(hk.percUnder1.combined, 2), nsmall = 2),'%',sep=''),x=T,y=2.5) +
  labs(x='Is Housekeeping Orthogroup?',y='',title = 'Human and Bear') +
  scale_y_continuous(limits = c(0,3.2)) +
  theme_linedraw() + theme(panel.grid = element_blank())


combined.tmm.stDevs.combined %>% 
  filter(is_hk==T) %>% 
  mutate(lessThan1 = ifelse(stdev < 1,T,F)) %>% 
  group_by(lessThan1) %>% 
  tally() %>% 
  mutate(percent = n / sum(n)) # 99.3% of housekeeping gene orthogroups have st dev < 1 across all species, tissues


# p.separate + p.combined

p.human + p.bear + p.combined




