
library(tidyverse)

human.raw <- read_tsv('analyses/orthoFinder/id_convserion/homSap_proteinHeaders.txt',col_names = F)
bear.raw <- read_tsv('analyses/orthoFinder/id_convserion/ursArc_proteinHeaders.txt',col_names = F)

human.convertTable <- human.raw %>% 
  separate(X1,sep=' ',into = c('txid','geneid','chr','type')) %>% 
  select(-type) %>% 
  mutate(txid = str_remove(txid,'>'),
         geneid = str_remove(geneid,'gene='),
         chr = str_remove(chr,'seq_id=')
         ) %>% 
  arrange(geneid)

head(human.convertTable)

bear.convertTable <- bear.raw %>% 
  separate(X1,sep=' ',into = c('txid','geneid','name','chr','type')) %>% 
  select(-type) %>% 
  mutate(txid = str_remove(txid,'>'),
         geneid = str_remove(geneid,'gene='),
         chr = str_remove(chr,'seq_id='),
         name = str_remove(name,'name=')
  ) %>% 
  arrange(geneid)

head(bear.convertTable)

# write_tsv(bear.convertTable,'analyses/orthoFinder/id_convserion/bear_conversionTable_10.20.21.txt')
# write_tsv(human.convertTable,'analyses/orthoFinder/id_convserion/human_conversionTable_10.20.21.txt')


### Make full human-bear ortholog table with converted ids 

human.bear.orthos <- read_tsv('analyses/orthoFinder/Orthologues/Orthologues_homSap.longCDS.prots/homSap.longCDS.prots__v__ursArc.longCDS.prots.tsv') %>% 
  select(orthogroup1=1,humanID=2,bearID=3) %>% 
  mutate(oneToOne = ifelse(str_detect(humanID,', '),F,
                           ifelse(str_detect(bearID,', '),F,
                                  T))) %>% 
  mutate(humanID = str_split(humanID, ', '),
         bearID = str_split(bearID,', ')) %>% 
  unnest(humanID) %>% 
  unnest(bearID)

orthos_idConvert <- human.bear.orthos %>% 
  left_join(human.convertTable,by=c('humanID'='txid')) %>% 
  mutate(humanGeneID = geneid) %>% 
  select(-geneid,-chr) %>% 
  left_join(bear.convertTable,by=c('bearID'='txid')) %>% 
  mutate(bearGeneID = geneid) %>% 
  select(-geneid,-name,-chr)
  
head(orthos_idConvert)

write_tsv(orthos_idConvert,'analyses/orthoFinder/id_convserion/HumanBearOrthologs_wIdConvert_10.20.21.tsv')
