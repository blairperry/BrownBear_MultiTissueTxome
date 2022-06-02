
# BiocManager::install("TissueEnrich")

library(TissueEnrich)
library(SummarizedExperiment)
library(patchwork)
library(tidyverse)

load('analyses/isoformExpression/7_DTUanalyses/dTurtleInputData_forKamiak_02.04.22.Rdata')

dir <- '/Volumes/WorkingDrive_BWP/_bearMultiTissueRNAseq_Sept2021/data/isoforms/kallisto_new_Feb2022/kallisto_quants_Feb2022/'

filenames <- list.files(dir)

files <- file.path(dir, filenames, "abundance.h5")
names(files) <- str_split_fixed(filenames,'[_]',2)[,1]

txi.kallisto <- tximport::tximport(files, type = "kallisto", txOut = TRUE)

id.convert <- read_csv('data/isoforms/new_merge.combined.renamed.mapping.txt')

allTissue.counts <- as.data.frame(txi.kallisto$abundance) %>% # Reading in the TPM estimates from Kallisto
  rownames_to_column(var='old_id') %>% 
  left_join(id.convert) %>% 
  dplyr::select(-old_id) %>% 
  dplyr::select(tx_id=new_id,everything())


bear.tx_level.median_counts <- allTissue.counts %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'counts') %>%
  mutate(sample = str_to_upper(sample)) %>% 
  filter(sample != 'CFA') %>% 
  left_join(pd.multiTissue,by=c('sample'='id')) %>%
  filter(!is.na(group)) %>%
  mutate(group = ifelse(group %in% c('Ventricle','Atrium'),'Heart (Ventricle + Atrium)',group)) %>%
  mutate(group = ifelse(str_detect(group,'Kidney'),'Kidney (Medulla + Cortex)',group)) %>%
  group_by(tx_id,group) %>%
  summarise(median_count = median(counts)) %>%
  pivot_wider(names_from = group,values_from = median_count) %>%
  column_to_rownames(var='tx_id') %>% 
  filter(rowSums(as.matrix(.)>1)>1)

bear.tx_level.se<-SummarizedExperiment(assays = SimpleList(as.matrix(bear.tx_level.median_counts)),rowData = row.names(bear.tx_level.median_counts),colData = colnames(bear.tx_level.median_counts))

bear.tx_level.TisseEnrichTest <- teGeneRetrieval(bear.tx_level.se,expressedGeneThreshold = 5)

bear.tx_level.tissueSpecRes <- assay(bear.tx_level.TisseEnrichTest) %>% as.data.frame()

bear.tx_level.tissueSpecRes.filtered <- bear.tx_level.tissueSpecRes %>%
  filter(str_detect(Group,'Enriched|Enhanced'))

p1 <- bear.tx_level.tissueSpecRes.filtered %>% filter(str_detect(Group,'Tissue-')) %>% 
  group_by(Tissue) %>% 
  tally() %>% 
  arrange(-n) %>% 
  ggplot(aes(x=reorder(Tissue,-n),y=n)) +
  geom_bar(stat='identity',fill='dodgerblue4') + 
  labs(x='Tissue',y='Number of Tissue-Specific Transcripts',title = 'Bear - Tissue Specific Transcripts') +
  # scale_y_continuous(expand = c(0,0),limits =c(0,3000)) +
  theme_linedraw() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1,vjust=1))
p1

median_counts.forPlot <- as.data.frame(bear.tx_level.median_counts) %>% rownames_to_column(var='tx_id') %>% 
  pivot_longer(-1,names_to='tissue',values_to = 'median_count')

# ggplot(subset(median_counts.forPlot,tx_id =='PB.10074.7'),aes(x=reorder(tissue,-median_count),y=median_count)) +
#   geom_bar(stat='identity',fill='dodgerblue4') +
#   labs(x='Tissue',y='Median Count') +
#   theme_linedraw() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1,vjust=1))


# tissueSpecRes.filtered.annot <- tissueSpecRes.filtered %>% 
#   left_join(tx2gene,by=c('Gene'='transcript_id')) %>% 
#   dplyr::select(1:4)


# Bear gene-level counts --------------------------------------------------

bear.gene_counts <- read_tsv('data/raw_counts/a_geneLevelQuants/bearNecropsy_geneLevelCounts_10.12.21_bwp.cleaned.txt',comment = '#') %>% 
  janitor::clean_names() %>% 
  dplyr::select(-start,-end,-chr,-strand) %>% 
  mutate(geneid = str_split_fixed(geneid,'[:]',3)[,2]) %>% 
  mutate(rowSum = rowSums(as.matrix(.[,-c(1,2)]))) %>%
  filter(rowSum > 0) %>%
  dplyr::select(-rowSum)

bear.gene_counts.mat <- as.matrix(bear.gene_counts[,-c(1,2)])
row.names(bear.gene_counts.mat) <- bear.gene_counts$geneid

bear.x <- bear.gene_counts.mat / bear.gene_counts$length
tpm.mat <- t( t(bear.x) * 1e6 / colSums(bear.x) )
bear.tpm.mat <- tpm.mat

bear.tpm.mat.df <- bear.tpm.mat %>% as.data.frame() %>% rownames_to_column()

# write_csv(bear.tpm.mat.df,'analyses/rnaseq_normalization/bear_TPMnorm_geneLevelCounts.csv')

bear.gene.median_counts <- bear.tpm.mat.df %>% 
  dplyr::select('gene_id'=rowname,everything()) %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'counts') %>% 
  mutate(sample = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  left_join(pd.multiTissue,by=c('sample'='id')) %>% 
  filter(!is.na(group)) %>% 
  mutate(group = ifelse(group %in% c('Ventricle','Atrium'),'Heart (Ventricle + Atrium)',group)) %>% 
  mutate(group = ifelse(str_detect(group,'Kidney'),'Kidney (Medulla + Cortex)',group)) %>% 
  group_by(gene_id,group) %>% 
  summarise(median_count = median(counts)) %>% 
  pivot_wider(names_from = group,values_from = median_count) %>% 
  column_to_rownames(var='gene_id')


bear.gene_level.se<-SummarizedExperiment(assays = SimpleList(as.matrix(bear.gene.median_counts)),rowData = row.names(bear.gene.median_counts),colData = colnames(bear.gene.median_counts))

bear.gene_level.TissueEnrichTest <- teGeneRetrieval(bear.gene_level.se,expressedGeneThreshold = 5)

bear.gene_level.tissueSpecRes <- assay(bear.gene_level.TissueEnrichTest) %>% as.data.frame()

bear.gene_level.tissueSpecRes.filtered <- bear.gene_level.tissueSpecRes %>%
  filter(str_detect(Group,'Enriched|Enhanced'))

bear.gene_level.tissueSpecRes.filtered %>% 
  bind_rows(bear.tx_level.tissueSpecRes.filtered) %>% 
  mutate(type = ifelse(str_detect(Gene,'gene'),'gene','transcript')) %>% 
  filter(str_detect(Group,'Tissue-')) %>% 
  group_by(Tissue,type) %>% 
  tally() %>% 
  arrange(desc(type),n) %>% 
  mutate(Tissue = ifelse(str_detect(Tissue,'Kidney'),'Kidney',Tissue)) %>% 
  pivot_wider(names_from = type,values_from = n) %>% 
  ggplot(aes(x=gene,y=transcript,color=Tissue)) +
  geom_abline(intercept = 0,slope = 1,lty=2,color='grey') +
  geom_point(aes(fill=Tissue),pch=21,color='black',size=3) +
  ggrepel::geom_text_repel(aes(label=Tissue)) +
  scale_x_continuous(limits=c(0,700)) +
  scale_y_continuous(limits=c(0,700)) +
  labs(x='Number of tissue-specific genes',y='Number of tissue-specific transcripts') +
  theme_linedraw() + theme(aspect.ratio = 1)


