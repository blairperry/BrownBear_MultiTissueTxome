
# BiocManager::install(c("BiocParallel", "GenomicRanges", "Gviz", "rtracklayer", "stageR", "tximport", "DESeq2"))
# remotes::install_github("TobiTekath/DTUrtle",upgrade = 'never')

library(DTUrtle)
library(tidyverse)

#the BiocParallel framework is used to parallelize the computations.
#Using 4 cores:
# biocpar <- BiocParallel::MulticoreParam(2)
#or standard serial computation (only 1 core)
biocpar <- BiocParallel::SerialParam()
#multiple other options available for computational clusters.

# Import and format data --------------------------------------------------

##### import gtf Annotation to get transcript to gene mapping
tx2gene <- import_gtf(gtf_file = "data/isoforms/new_merge.combined.renamed_tappAS_annot_from_SQANTI3.gff3",feature_type = 'gene')

tx2gene <- move_columns_to_front(df = tx2gene, 
                                 columns = c("seqnames", "ID"))

tx2gene <- tx2gene %>% 
  mutate(ID = ifelse(ID == '', 'None', ID)) %>% 
  select(transcript_id=1,gene_id=2,everything())


id.convert <- read_csv('data/isoforms/new_merge.combined.renamed.mapping.txt')

tx2gene_tempForImport <- tx2gene %>% 
  left_join(id.convert,by=c('transcript_id'='new_id')) %>% 
  select(old_id,gene_id)


##### import transcript-level quantification data, for example from Salmon, and perform appropriate scaling
files <- Sys.glob("/Volumes/WorkingDrive_BWP/_bearMultiTissueRNAseq_Sept2021/data/isoforms/kallisto_new_Feb2022/kallisto_quants_Feb2022/*/abundance.h5")
names(files) <- gsub(".*/","",gsub("/abundance.h5","",files))

cts_temp <- import_counts(files = files, type = "kallisto",
                          countsFromAbundance='dtuScaledTPM',
                          countsCol='est_counts',
                          tx2gene=tx2gene_tempForImport)

##### rename transcripts using ID conversion table

cts_temp2 <- cts_temp %>% 
  as.data.frame() %>% 
  rownames_to_column(var='old_id') %>% 
  left_join(id.convert) %>% 
  select(-old_id) %>% 
  select(new_id,everything()) %>% 
  column_to_rownames(var='new_id')

names(cts_temp2) <- str_to_upper(str_split_fixed(names(cts_temp2),'[_]',2)[,1])
names(cts_temp2)


cts <- as.matrix(cts_temp2)
head(cts)



##### create a sample data sheet, specifying which sample belongs to which group

sample_info_necrop <- readxl::read_xlsx('data/sample_info/MultiTissue_mapping_stats.xlsx') %>% 
  janitor::clean_names() %>% 
  select(sample_id,tissue)

sample_info_prevRNA <- readxl::read_xlsx('data/sample_info/prevRNA_SampleInfo.xlsx') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = paste(tissue,phys_state,sep='_')) %>% 
  mutate(sample_id = str_to_upper(str_split_fixed(sample,'[_]',2)[,1])) %>% 
  select(sample_id, tissue)

all.sample_info <- sample_info_necrop %>% 
  bind_rows(sample_info_prevRNA)

pd <- data.frame("id"=colnames(cts), "group"="your_grouping_variable", 
                 stringsAsFactors = FALSE) %>% 
  mutate(id = str_split_fixed(id,'[_]',2)[,1]) %>% 
  mutate(id = str_to_upper(id)) %>% 
  left_join(all.sample_info,by=c('id'='sample_id')) %>% 
  dplyr::select(id,group=tissue)

# write_tsv(pd,'analyses/isoformExpression/fullSampleInfo.tsv')



# Format and save Rdata object for working on Kamiak ----------------------



pd.multiTissue <- pd %>% 
  filter(str_detect(group,'Active|Hibernation|Hyperphagia',negate = T)) %>% 
  mutate(group = ifelse(str_detect(group,'Cortex'),'Kidney_cortex',ifelse(str_detect(group,'Medulla'),'Kidney_medulla',group))) %>% 
  mutate(group = ifelse(str_detect(group,'[(]'),str_split_fixed(group,'[ ]',2)[,1],group)) 

temp <- pd.multiTissue$group %>% unique()
all.combos <- as.data.frame(t(utils::combn(temp,2)))

save(cts,pd.multiTissue,all.combos,tx2gene,file = 'analyses/isoformExpression/7_DTUanalyses/dTurtleInputData_forKamiak_02.04.22.Rdata')



######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### Write script to conduct all pairwise comparisons (will run on Kamiak)
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 

# BiocManager::install(c("BiocParallel", "GenomicRanges", "Gviz", "rtracklayer", "stageR", "tximport", "DESeq2"))
# remotes::install_github("TobiTekath/DTUrtle",upgrade = 'never')

library(DTUrtle)
library(tidyverse)

#the BiocParallel framework is used to parallelize the computations.
#Using 12 cores:
biocpar <- BiocParallel::MulticoreParam(12)
#or standard serial computation (only 1 core)
# biocpar <- BiocParallel::SerialParam()
#multiple other options available for computational clusters.

load('dTurtleInputData_forKamiak.Rdata')

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
