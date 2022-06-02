
library(tidyverse)


# Read in data ------------------------------------------------------------

sample.info <- read_tsv('data/gtex_human_rnaseq/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
head(sample.info)

# NOTE: Relevant tissue information is in columns SMTS and SMTSD, with latter being more detailed.
#       Want to include stomach, spleen, small intestine, skin, lung, liver, kidney (medulla and cortex if available), heart (ventricle and atrium), gall bladder, adipose, and muscle types

sample.info %>% 
  filter(str_detect(SMTS,'Liver|Adipose|Muscle|Stomach|Spleen|Small Intestine|Skin|Lung|Kidney|Heart|Gall')) %>% 
  filter(str_detect(SMAFRZE,'RNASEQ')) %>% 
  select(SMTSD) %>% 
  unique() 

# Subcategories to exclude: Heart Atrial Appendage, Adipose Visceral, Skin - Sun Exposed (Lower Leg), Cells Cultured fibroblasts

sample.info.filtered <- sample.info %>% 
  filter(str_detect(SMTS,'Liver|Adipose|Muscle|Stomach|Spleen|Small Intestine|Skin|Lung|Kidney|Heart|Gall')) %>% 
  filter(str_detect(SMTSD,'Appendage|Visceral|Lower|Cultured',negate = T)) %>% 
  filter(str_detect(SMAFRZE,'RNASEQ')) %>% 
  arrange(SMTS,SAMPID)

sample.info.filtered %>% select(SMTS,SMTSD) %>% unique() 

sample.info.filtered %>% 
  group_by(SMTS) %>% 
  tally() 

## Randomly subsample to 20 samples per tissue
set.seed(100)

sample.info.filtered.sample <- sample.info.filtered %>% group_by(SMTS) %>% sample_n(size = 20)

## Read in raw gene count headers and extract relevant column indices
## Because the raw count table is so large, it's slow to subset columns in R. Going to get indices for desired columns and subset in bash w/ awk

gtex.header <- read_tsv('data/gtex_human_rnaseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct',skip=2,n_max = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  select(1) %>% 
  rowid_to_column()

gtex.header.indices <- gtex.header %>% 
  filter(rowname == 'Name' | rowname == 'Description' | rowname %in% sample.info.filtered.sample$SAMPID) %>% 
  mutate(rowid = paste('$',rowid,sep = ''))

cat(paste(gtex.header.indices$rowid,collapse = ','))

## Run in terminal
# tail -n +3 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct | awk 'BEGIN {OFS="\t"} {print $1,$2,$27,$120,$805,$843,$877,$934,$947,$1085,$1184,$1204,$1219,$1449,$1545,$1693,$1694,$1892,$2007,$2043,$2314,$2344,$2636,$2728,$2740,$2742,$2756,$2783,$2786,$3156,$3610,$3615,$3625,$3797,$3885,$3978,$4021,$4065,$4067,$4095,$4110,$4530,$4661,$4680,$4731,$4818,$4902,$4924,$5052,$5063,$5334,$5335,$5499,$5533,$5558,$5576,$5903,$5994,$6026,$6028,$6256,$6464,$6574,$6706,$6731,$6792,$7024,$7043,$7081,$7126,$7160,$7206,$7273,$7395,$7555,$7717,$7900,$8106,$8222,$8434,$8435,$8516,$8639,$8682,$8827,$8829,$8935,$9072,$9086,$9164,$9259,$9367,$9385,$9603,$9628,$9647,$9699,$9701,$9810,$9868,$9928,$9937,$9959,$10035,$10055,$10317,$10585,$10633,$10717,$10732,$10743,$10787,$11045,$11182,$11328,$11422,$11437,$11510,$11565,$11663,$11773,$11798,$11816,$11825,$11833,$11936,$11981,$12003,$12027,$12101,$12117,$12184,$12234,$12288,$12353,$12421,$12427,$12502,$12561,$12620,$12656,$12681,$12796,$13028,$13133,$13170,$13313,$13330,$13488,$13524,$13574,$13581,$13617,$13618,$13638,$13688,$13708,$13719,$13799,$13894,$14072,$14076,$14219,$14343,$14401,$14575,$14677,$14709,$14871,$14918,$14952,$14967,$14969,$15097,$15247,$15311,$15394,$15514,$15516,$15763,$15814,$15917,$15971,$16040,$16049,$16052,$16079,$16137,$16241,$16424,$16446,$16544,$16621,$16656,$16700,$16888,$16967,$17133,$17217,$17262,$17266,$17345}' > GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.subsample.txt

## Read in subsampled data and verify that all samples are present (need to unzip gzipped file first)

gtex.sub <- read_tsv('data/gtex_human_rnaseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.subsample.txt')

gtex.sub.samps <- gtex.sub %>% 
  pivot_longer(-c(1,2),names_to = 'sample',values_to = 'rawCounts') %>% 
  select(sample) %>% 
  unique()

sum(gtex.sub.samps$sample %in% sample.info.filtered.sample$SAMPID) # All Good!


