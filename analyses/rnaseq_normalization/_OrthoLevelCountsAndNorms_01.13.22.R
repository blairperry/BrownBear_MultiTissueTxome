
library(tidyverse)
library(janitor)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(viridis)
library(patchwork)

###
### Converting gene-level raw counts to orthogroup-level counts - Human RNA-seq data
###


# Read in gff and calculate gene length (sum of exons) for all genes

human.gff <- read_tsv('data/reference/hsapiens/homSap.exons.gff',col_names = c('chr','type','feature','start','stop','ignore','strand','ignore2','description')) %>% 
  select(-chr,-type,-feature,-ignore,-strand,-ignore2) %>% 
  mutate(geneID = str_split_fixed(description,';',4)[,3] %>% str_remove('gene_id=')) %>% 
  mutate(length = abs(stop-start)) %>% 
  select(geneID,length) %>% 
  group_by(geneID) %>% 
  summarise(exonLength = sum(length))
  
head(human.gff)


# Read in Orthogroup information ------------------------------------------

human.ogs <- read_tsv('analyses/orthoFinder/id_conversion/HumanBearOrthologs_wIdConvert_10.20.21.tsv') %>% 
  select(humanGeneID,orthogroup1) 
  
head(human.ogs)

# Read in raw counts and sample info ----------------------------------

gtex.sampleInfo <- read_tsv('data/gtex_human_rnaseq/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt') %>% 
  select(sample=SAMPID,tissue=SMTS) %>% 
  mutate(tissue = str_replace_all(tissue,' ','_')) %>% 
  mutate(tissue = str_remove_all(tissue,'_Tissue'))

# Note: May need to unzip gzipped file before reading into R
all.raw_counts <- read_tsv('data/gtex_human_rnaseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.subsample.txt') %>% 
  pivot_longer(-c(1,2),names_to = 'sample',values_to = 'rawcounts') %>% 
  left_join(gtex.sampleInfo,by='sample') %>% 
  mutate(sample = paste(tissue,sample,sep='_')) %>% 
  select(-tissue) %>% 
  pivot_wider(names_from=sample,values_from=rawcounts,values_fill = 0) %>%
  mutate(rowSum = rowSums(as.matrix(.[,-c(1,2)]))) %>%
  filter(rowSum > 0) %>%
  select(-rowSum) %>%
  left_join(human.ogs,by=c('Name'='humanGeneID')) %>%
  left_join(human.gff,by=c('Name'='geneID')) %>%
  mutate(orthogroup1 = ifelse(is.na(orthogroup1),Name,orthogroup1)) # Currently, any gene without an orthogroup with bear is being considered it's own orthogroup

head(all.raw_counts$orthogroup1,20)
all.raw_counts


# Sum gene-level counts by orthogroup to get raw orthogroup-level  --------

all.raw_counts.oglevel <- all.raw_counts %>% 
  pivot_longer(-c(1,2,203,204),names_to = 'sample',values_to = 'counts') %>% 
  group_by(orthogroup1,sample) %>% 
  summarise(og_counts = sum(counts),
            og_length = sum(exonLength)) %>% 
  pivot_wider(names_from = sample,values_from = og_counts)

# write_tsv(all.raw_counts.oglevel,'analyses/rnaseq_normalization/human_oglevel_rawCounts.txt')
# all.raw_counts.oglevel <- read_tsv('analyses/rnaseq_normalization/human_oglevel_rawCounts.txt')

head(all.raw_counts.oglevel)

### Normalizing data
### Following a tpm -> TMM approach based on: https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-021-04414-y.pdf

counts.mat <- as.matrix(all.raw_counts.oglevel[,-c(1,2)])
row.names(counts.mat) <- all.raw_counts.oglevel$orthogroup1

x <- counts.mat / all.raw_counts.oglevel$og_length
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
human.tpm.mat <- tpm.mat

human.tpm.mat.df <- human.tpm.mat %>% as.data.frame() %>% rownames_to_column()
# write_tsv(human.tpm.mat.df,'analyses/rnaseq_normalization/human_oglevel_tpmCounts.txt')

edgeR.delist <- DGEList(counts = tpm.mat)
edgeR.delist <- calcNormFactors(edgeR.delist)

human.tmm.counts <- as.data.frame(cpm(edgeR.delist))

head(human.tmm.counts)


# Reformat table for easier plotting --------------------------------------

human.tmm.forPlot <- human.tmm.counts %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(str_detect(rowname,'OG')) %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'norm_count') %>% 
  mutate(treatment = str_split_fixed(sample,'[_]',2)[,1]) %>% 
  filter(norm_count > 0)

ggplot(human.tmm.forPlot,aes(x=norm_count,group=sample,color=treatment)) +
  geom_density(lwd=0.1,alpha=0.5) +
  scale_x_log10() +
  theme_linedraw()+
  ggtitle('Human - TPM+TMM Normalized Counts') +
  xlab('Orthogroup-level Counts (log10 scaled axis)') +
  ylab('Density') +
  theme_linedraw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'),aspect.ratio = 1)

human.tpm.forPlot <- human.tpm.mat %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(str_detect(rowname,'OG')) %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'norm_count') %>% 
  mutate(treatment = str_split_fixed(sample,'[_]',2)[,1]) %>% 
  filter(norm_count > 0)

ggplot(human.tpm.forPlot,aes(x=norm_count,group=sample,color=treatment)) +
  geom_density(lwd=0.1,alpha=0.5) +
  scale_x_log10() +
  theme_linedraw()+
  ggtitle('Human - TPM Normalized Counts') +
  xlab('Orthogroup-level Counts (log10 scaled axis)') +
  ylab('Density') +
  theme_linedraw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'),aspect.ratio = 1)

human.raw.forPlot <- all.raw_counts.oglevel %>% 
  as.data.frame() %>% 
  filter(str_detect(orthogroup1,'OG')) %>% 
  pivot_longer(-c(1,2),names_to = 'sample',values_to = 'norm_count') %>% 
  mutate(treatment = str_split_fixed(sample,'[_]',2)[,1]) %>% 
  filter(norm_count > 0)

ggplot(human.raw.forPlot,aes(x=norm_count,group=sample,color=treatment)) +
  geom_density(lwd=0.1,alpha=0.5) +
  scale_x_log10() +
  ggtitle('Human - Raw Counts') +
  xlab('Orthogroup-level Counts (log10 scaled axis)') +
  ylab('Density') +
  theme_linedraw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'),aspect.ratio = 1)



###
### Converting gene-level raw counts to orthogroup-level counts - Bear RNA-seq data
###

# Read in gff and calculate gene length (sum of exons) for all genes

bear.gff <- read_tsv('data/reference/uarctos/ursArc.longCDS.exons.gff',col_names = c('chr','type','feature','start','stop','ignore','strand','ignore2','description'),comment = "#") %>% 
  select(-chr,-type,-feature,-ignore,-strand,-ignore2) %>% 
  mutate(geneID = str_split_fixed(description,';',3)[,2] %>% str_remove('Parent=')) %>% 
  mutate(length = abs(stop-start)) %>% 
  select(geneID,length) %>% 
  group_by(geneID) %>% 
  summarise(exonLength = sum(length))

head(bear.gff)

# Read in Orthogroup information ------------------------------------------

bear.ogs <- read_tsv('analyses/orthoFinder/id_conversion/HumanBearOrthologs_wIdConvert_10.20.21.tsv') %>% 
  select(bearID, bearGeneID,orthogroup1) 

head(bear.ogs)


# Read in raw counts and sample info ----------------------------------

all.raw_counts <- read_tsv('data/raw_counts/a_geneLevelQuants/bearNecropsy_geneLevelCounts_10.12.21_bwp.cleaned.txt',comment = '#') %>% 
  clean_names() %>% 
  select(-start,-end,-chr,-strand) %>% 
  mutate(geneid = str_split_fixed(geneid,'[:]',3)[,2]) %>% 
  mutate(rowSum = rowSums(as.matrix(.[,-c(1,2)]))) %>%
  filter(rowSum > 0) %>%
  select(-rowSum) %>%
  left_join(bear.ogs,by=c('geneid'='bearGeneID')) %>%
  left_join(bear.gff,by=c('bearID'='geneID')) %>%
  mutate(orthogroup1 = ifelse(is.na(orthogroup1),geneid,orthogroup1)) %>%  # Currently, any gene without an orthogroup with bear is being considered it's own orthogroup
  mutate(exonLength = ifelse(is.na(exonLength),length,exonLength))
head(all.raw_counts$orthogroup1,20)

all.raw_counts.oglevel <- all.raw_counts %>% 
  pivot_longer(-c(1,2,89,90,91),names_to = 'sample',values_to = 'counts') %>% 
  group_by(orthogroup1,sample) %>% 
  summarise(og_counts = sum(counts),
            og_length = sum(exonLength)) %>% 
  pivot_wider(names_from = sample,values_from = og_counts)

# write_tsv(all.raw_counts.oglevel,'analyses/rnaseq_normalization/bear_oglevel_rawCounts.txt')
# all.raw_counts.oglevel <- read_tsv('analyses/rnaseq_normalization/bear_oglevel_rawCounts.txt')

head(all.raw_counts.oglevel,20)

### Normalizing data
### Following a tpm -> TMM approach based on: https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-021-04414-y.pdf

counts.mat <- as.matrix(all.raw_counts.oglevel[,-c(1,2)])
row.names(counts.mat) <- all.raw_counts.oglevel$orthogroup1

x <- counts.mat / all.raw_counts.oglevel$og_length
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
bear.tpm.mat <- tpm.mat

bear.tpm.mat.df <- bear.tpm.mat %>% as.data.frame() %>% rownames_to_column()

# write_tsv(bear.tpm.mat.df,'analyses/rnaseq_normalization/bear_oglevel_tpmCounts.txt')


edgeR.delist <- DGEList(counts = tpm.mat)
edgeR.delist <- calcNormFactors(edgeR.delist)

bear.tmm.counts <- as.data.frame(cpm(edgeR.delist))

head(bear.tmm.counts)

# write_csv(tmm.counts,'analyses/rnaseq_normalization/bear_OrthoLevelTPMs_11.11.21.csv')


# Read in and combine sample info from this and previous study ------------

sample_info_necrop <- readxl::read_xlsx('data/sample_info/MultiTissue_mapping_stats.xlsx') %>% 
  clean_names() %>% 
  select(sample_id,tissue)

sample_info_prevRNA <- readxl::read_xlsx('data/sample_info/prevRNA_SampleInfo.xlsx') %>% 
  clean_names() %>% 
  mutate(tissue = paste(tissue,phys_state,sep='_')) %>% 
  mutate(sample_id = str_to_upper(str_split_fixed(sample,'[_]',2)[,1])) %>% 
  select(sample_id, tissue)

all.sample_info <- sample_info_necrop %>% 
  bind_rows(sample_info_prevRNA)


# Reformat table for easier plotting --------------------------------------

bear.tmm.forPlot <- bear.tmm.counts %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(str_detect(rowname,'OG')) %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'norm_count') %>% 
  mutate(sample = str_replace_all(sample,'_hy','hy')) %>% 
  mutate(sample = str_replace_all(sample,'_hi','hi')) %>% 
  mutate(sample = str_replace_all(sample,'_fa','fa')) %>% 
  mutate(sample_id = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  left_join(all.sample_info,by='sample_id') %>% 
  filter(norm_count > 0 ) %>% 
  filter(str_detect(tissue,'Active')) 

ggplot(bear.tmm.forPlot,aes(x=norm_count,group=sample,color=tissue)) +
  geom_density(lwd=0.3,alpha=0.85) +
  scale_x_log10() +
  theme_linedraw()+
  ggtitle('Bear - TPM+TMM Normalized Counts') +
  xlab('Orthogroup-level Counts (log10 scaled axis)') +
  ylab('Density') +
  theme_linedraw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'),aspect.ratio = 1)

bear.tpm.forPlot <- bear.tpm.mat %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(str_detect(rowname,'OG')) %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'norm_count') %>% 
  filter(sample %in% bear.tmm.forPlot$sample) %>% 
  filter(norm_count > 0) %>% 
  mutate(sample_id = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  left_join(all.sample_info,by='sample_id') 

head(bear.tpm.forPlot)

ggplot(bear.tpm.forPlot,aes(x=norm_count,group=sample,color=tissue)) +
  geom_density(lwd=0.3,alpha=0.85) +
  scale_x_log10() +
  theme_linedraw()+
  ggtitle('Bear - TPM Normalized Counts') +
  xlab('Orthogroup-level Counts (log10 scaled axis)') +
  ylab('Density') +
  theme_linedraw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'),aspect.ratio = 1)


bear.raw.forPlot <- all.raw_counts.oglevel %>% 
  as.data.frame() %>% 
  filter(str_detect(orthogroup1,'OG')) %>% 
  pivot_longer(-c(1,2),names_to = 'sample',values_to = 'norm_count') %>% 
  filter(sample %in% bear.tmm.forPlot$sample) %>% 
  mutate(sample_id = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  left_join(all.sample_info,by='sample_id') %>% 
  filter(norm_count > 0)

ggplot(bear.raw.forPlot,aes(x=norm_count,group=sample,color=tissue)) +
  geom_density(lwd=0.3,alpha=0.85) +
  scale_x_log10() +
  ggtitle('Bear - Raw Counts') +
  xlab('Orthogroup-level Counts (log10 scaled axis)') +
  ylab('Density') +
  theme_linedraw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'),aspect.ratio = 1)



##########################################################################
# Combined TMM normalization of bear and human samples simultaneously ####
##########################################################################


# Read in TPM data (generated above) --------------------------------------

bear.tpm.mat <- read_tsv('analyses/rnaseq_normalization/bear_oglevel_tpmCounts.txt')
human.tpm.mat <- read_tsv('analyses/rnaseq_normalization/human_oglevel_tpmCounts.txt')

bear.tpm.df <- bear.tpm.mat %>% 
  dplyr::rename(og_id=rowname) %>% 
  filter(str_detect(og_id,'OG')) %>% 
  select(-cfa_s9) # Remove contaminated adipose sample (see Jansen et al 2019)

human.tpm.df <- human.tpm.mat %>% 
  dplyr::rename(og_id=rowname) %>% 
  filter(str_detect(og_id,'OG'))

# Combine human and bear tpm tables and normalize with TMM method in edgeR
combined.tpm.df <- bear.tpm.df %>% 
  left_join(human.tpm.df) %>% 
  drop_na()

combined.tpm.mat <- as.matrix(combined.tpm.df[,-1])
row.names(combined.tpm.mat) <- combined.tpm.df$og_id

edgeR.delist <- DGEList(counts = combined.tpm.mat)

edgeR.delist <- calcNormFactors(edgeR.delist)

combined.tmm.counts <- cpm(edgeR.delist) %>% as.data.frame() %>% rownames_to_column(var='og_id')

head(combined.tmm.counts)

# write_tsv(combined.tmm.counts,'analyses/rnaseq_normalization/combined_oglevel_tmmCounts.txt')
# combined.tmm.counts <- read_tsv('analyses/rnaseq_normalization/combined_oglevel_tmmCounts.txt')



# Reformat for easier plotting --------------------------------------------

combined.tmm.forPlot <- combined.tmm.counts %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'tmm_count') %>% 
  mutate(species = ifelse(str_detect(sample,'GTEX'),'human','bear')) %>% 
  mutate(sample = str_replace_all(sample,'_hy','hy')) %>% 
  mutate(sample = str_replace_all(sample,'_hi','hi')) %>% 
  mutate(sample = str_replace_all(sample,'_fa','fa')) %>% 
  mutate(sample_id = str_split_fixed(sample,'[_]',2)[,1] %>% str_to_upper()) %>% 
  left_join(all.sample_info,by='sample_id') %>% 
  mutate(sample_id = ifelse(species=='human',sample,sample_id)) %>% 
  mutate(tissue = ifelse(species=='human',str_split_fixed(sample,'[_]',2)[,1],tissue)) 

combined.tmm.forPlot.active <- combined.tmm.forPlot %>% 
  filter((species == 'bear' & str_detect(tissue,'Active')) | (species == 'human')) %>% 
  group_by(og_id,species,tissue) %>% 
  summarise(median_count = median(tmm_count),
            avg_count = mean(tmm_count)) %>% 
  group_by(og_id) %>% 
  mutate(tissue = str_split_fixed(tissue,'[_]',2)[,1]) %>% 
  filter(tissue %in% c('Adipose','Liver','Muscle'))

combined.tmm.forPlot.active.v2 <- combined.tmm.forPlot.active %>% 
  select(og_id,species,tissue,median_count) %>% 
  pivot_wider(names_from = species,values_from = median_count)

combined.tmm.forPlot.active.goodGenes <- combined.tmm.forPlot.active.v2 %>% 
  filter(bear > 0 & human > 0)


p1 <- ggplot(combined.tmm.forPlot.active.goodGenes,aes(x=human,y=bear,color=tissue)) +
  facet_wrap(~tissue) +
  geom_point(alpha=0.2,pch=1) + 
  geom_abline(slope=1,intercept = 0,lty=2) +
  # geom_smooth() +
  scale_y_log10() +
  labs(x='Human - Median Expression per Orthogroup\n(log10 scaled axis)',y='Bear - Median Expression per Orthogroup\n(log10 scaled axis)') +
  scale_x_log10() +
  theme_linedraw()
p1

p2 <- ggplot(combined.tmm.forPlot.active.v2,aes(x=log2(human / bear),color=tissue)) +
  facet_wrap(~tissue) +
  geom_density(lwd=1) +
  geom_vline(xintercept = 0,lty=2) +
  xlab('Log2(Human Expression / Bear Expression)\n (Median Orthogroup-level Counts)') +
  theme_linedraw()

p1 / p2



# Calculate correlation between tissues and species -----------------------

corr.data <- combined.tmm.forPlot.active %>% 
  mutate(tissue = paste(species, tissue, sep = '_')) %>% 
  select(-avg_count,-species) %>% 
  pivot_wider(names_from = tissue, values_from = median_count,values_fill = 0) 
  
corr.res <- cor(corr.data[,-1],method = 'spearman') %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(str_detect(rowname,'human')) %>% 
  select(1,contains('bear')) %>% 
  column_to_rownames(var='rowname')

pheatmap(corr.res,
        display_numbers = T,
        border_color = 'white',
        number_color = 'black',
        color=viridis(50,begin = .45))

cor.test(x=corr.data$bear_Muscle,y=corr.data$human_Muscle,method = 'spearman')


# Reformat tables for PCA calculation and plotting ------------------------

pca.data.bear <- combined.tmm.forPlot %>% 
  # filter(str_detect(tissue,'dip') | str_detect(tissue,'us') | str_detect(tissue,'iv')) %>% 
  select(ogid = 1,sample,tmm_count,species,tissue) %>% 
  mutate(treatment = ifelse(species=='human','human',
                            ifelse(str_detect(tissue,'Hibern'),'Hibernation',
                                   ifelse(str_detect(tissue,'Hyper'),'Hyperphagia',
                                          ifelse(str_detect(tissue,'Active'),'Active','Necropsy'))))) %>% 
  mutate(tissue = ifelse(tissue == 'Adipose (Subcutaneous)','Adipose',tissue)) %>% 
  mutate(tissue_simple = ifelse(str_detect(tissue,'dip'),'Adipose',
                                ifelse(str_detect(tissue,'us'),'Muscle','Liver'))) %>% 
  filter(species != 'human') %>% 
  mutate(sample_detail = paste(sample,species,treatment,tissue,sep = ':')) %>% 
  select(ogid,sample_detail,tmm_count) %>% 
  mutate(tmm_count = log10(tmm_count+0.001)) %>% 
  pivot_wider(names_from = sample_detail,values_from = tmm_count) %>% 
  filter(ogid %in% combined.tmm.forPlot.active.goodGenes$og_id) %>% 
  rowwise() %>% 
  mutate(median = median(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = 'ogid') %>% 
  ungroup() %>% 
  arrange(-median) %>% 
  top_n(n=10000,wt=median) %>% # Select top 10000 genes
  select(-median) %>% 
  as.matrix()

pca.data.all.sds <- combined.tmm.forPlot %>% 
  select(ogid = 1,sample,tmm_count,species,tissue) %>% 
  mutate(treatment = ifelse(species=='human','human',
                            ifelse(str_detect(tissue,'Hibern'),'Hibernation',
                                   ifelse(str_detect(tissue,'Hyper'),'Hyperphagia',
                                          ifelse(str_detect(tissue,'Active'),'Active','Necropsy'))))) %>% 
  mutate(tissue_simple = str_split_fixed(tissue,'[_ ]',2)[,1]) %>% 
  mutate(sample_detail = paste(sample,species,treatment,tissue,sep = ':')) %>% 
  select(ogid,sample_detail,tmm_count) %>% 
  mutate(tmm_count = log10(tmm_count+1)) %>% 
  pivot_wider(names_from = sample_detail,values_from = tmm_count) %>% 
  filter(ogid %in% combined.tmm.forPlot.active.goodGenes$og_id) %>% 
  mutate(stdev = rowSds(as.matrix(.[,-1]))) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = 'ogid') %>% 
  ungroup() %>% 
  arrange(-stdev) %>% 
  top_n(n=250,wt=stdev) %>% # select top 250 genes based on stdev
  select(-stdev) %>% 
  as.matrix()

pca.res.bear <- prcomp(t(pca.data.bear))
pca.res.all.sds <- prcomp(t(pca.data.all.sds))

# Bear only
pca.scores.bear <- pca.res.bear$x %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  separate(rowname,sep='[:]',into = c('sample','species','treatment','tissue')) %>%
  mutate(tissue = str_split_fixed(tissue,'[_]',2)[,1]) %>%
  filter(tissue != 'NA')

##
## For main figure
##

p.bears.all <- pca.scores.bear %>% 
  mutate(tissue_simple = str_split_fixed(tissue,'[ (]',2)[,1]) %>% 
  mutate(tissue_simple = ifelse(str_detect(tissue,'Kidney'),tissue,tissue_simple)) %>% 
  ggplot(aes(x=PC1,y=PC2,fill=tissue_simple,color=tissue_simple,shape=treatment)) +
  geom_point(size=4,alpha=0.8,show.legend = T,color='black') +
  scale_shape_manual(values=c('Active'=21,'Hibernation'=24,'Hyperphagia'=22,'Necropsy'=23))+
  ggrepel::geom_text_repel(aes(label=tissue_simple),show.legend = T) +
  cols4all::scale_fill_discrete_c4a_cat(palette = 'rainbow') +
  cols4all::scale_color_discrete_c4a_cat(palette = 'rainbow') +
  xlab(paste('PC1 (Variance explained: ',summary(pca.res.bear)$importance[2,1]*100,'%)',sep='')) +
  ylab(paste('PC2 (Variance explained: ',summary(pca.res.bear)$importance[2,2]*100,'%)',sep='')) +
  theme_linedraw() + theme(aspect.ratio = 1)


# All - by st dev

pca.scores.all.sds <- pca.res.all.sds$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname,sep='[:]',into = c('sample','species','treatment','tissue')) %>% 
  mutate(tissue = str_split_fixed(tissue,'[_ ]',2)[,1]) %>% 
  filter(tissue != 'NA') %>% 
  mutate(tissue = str_replace_all(tissue,'Atrium','Heart'),
         tissue = str_replace_all(tissue,'Ventricle','Heart'),
         tissue = str_replace_all(tissue,'Small','Small Intestine'),
         tissue = str_replace_all(tissue,'Gall','Gall Bladder'))


h.pca1 <- ggplot(pca.scores.all.sds,aes(x=PC1,y=PC2,color=tissue,fill=tissue,shape=species)) +
  geom_point(size=4,alpha=0.7,show.legend = T,color='black') +
  # scale_shape_manual(values=c('Active'=1,'Hibernation'=2,'human'=0,'Hyperphagia'=4,'Necropsy'=5))+ ### Use when shape=treatment
  scale_shape_manual(values=c('human'=24,'bear'=21))+### use when shape=species
  cols4all::scale_fill_discrete_c4a_cat(palette = 'rainbow') +
  cols4all::scale_color_discrete_c4a_cat(palette = 'rainbow') +
  ggrepel::geom_text_repel(aes(label=tissue)) +
  scale_x_reverse() +
  xlab(paste('PC1 (Variance explained: ',summary(pca.res.all.sds)$importance[2,1]*100,'%)',sep='')) +
  ylab(paste('PC2 (Variance explained: ',summary(pca.res.all.sds)$importance[2,2]*100,'%)',sep='')) +
  theme_linedraw() + theme(panel.grid = element_line(color = 'grey60'),aspect.ratio = 1)

h.pca1

h.pca2 <- ggplot(pca.scores.all.sds,aes(x=PC1,y=PC3,color=tissue,fill=tissue,shape=species)) +
  geom_point(size=4,alpha=0.8,show.legend = T,color='black') +
  # scale_shape_manual(values=c('Active'=1,'Hibernation'=2,'human'=0,'Hyperphagia'=4,'Necropsy'=5))+  ### Use when shape=treatment
  scale_shape_manual(values=c('human'=24,'bear'=21))+### use when shape=species
  cols4all::scale_fill_discrete_c4a_cat(palette = 'rainbow') +
  cols4all::scale_color_discrete_c4a_cat(palette = 'rainbow') +
  ggrepel::geom_text_repel(aes(label=tissue)) +
  scale_x_reverse() +  xlab(paste('PC1 (Variance explained: ',summary(pca.res.all.sds)$importance[2,1]*100,'%)',sep='')) +
  ylab(paste('PC3 (Variance explained: ',summary(pca.res.all.sds)$importance[2,3]*100,'%)',sep='')) +
  theme_linedraw() + theme(panel.grid = element_line(color = 'grey60'),aspect.ratio = 1)

h.pca3 <- ggplot(pca.scores.all.sds,aes(x=PC2,y=PC3,color=tissue,fill=tissue,shape=species)) +
  geom_point(size=4,alpha=0.8,show.legend = T,color='black') +
  # scale_shape_manual(values=c('Active'=1,'Hibernation'=2,'human'=0,'Hyperphagia'=4,'Necropsy'=5))+  ### Use when shape=treatment
  scale_shape_manual(values=c('human'=24,'bear'=21))+### use when shape=species
  cols4all::scale_fill_discrete_c4a_cat(palette = 'rainbow') +
  cols4all::scale_color_discrete_c4a_cat(palette = 'rainbow') +
  ggrepel::geom_text_repel(aes(label=tissue)) +
  scale_x_reverse() +  xlab(paste('PC2 (Variance explained: ',summary(pca.res.all.sds)$importance[2,2]*100,'%)',sep='')) +
  ylab(paste('PC3 (Variance explained: ',summary(pca.res.all.sds)$importance[2,3]*100,'%)',sep='')) +
  theme_linedraw() + theme(panel.grid = element_line(color = 'grey60'),aspect.ratio = 1)


# Fig 1 panels

p.bears.all + theme(panel.grid = element_blank()) | h.pca1 + theme(panel.grid = element_blank())


