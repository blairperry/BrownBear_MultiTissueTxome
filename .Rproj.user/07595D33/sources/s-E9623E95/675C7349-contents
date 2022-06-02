# BiocManager::install("tximport")
#BiocManager::install("rhdf5")

library(tidyverse)
library(tximport)
library(rhdf5)


# Read in Kallisto results and combine ------------------------------------

dir <- '/Volumes/WorkingDrive_BWP/_bearNecropsyRNAseq_Sept2021/data/isoforms/kallisto_quants'

filenames <- list.files(dir)

files <- file.path(dir, filenames, "abundance.h5")
names(files) <- str_split_fixed(filenames,'[_]',2)[,1]

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

head(txi.kallisto$counts)
head(txi.kallisto$abundance)

# Read in ID conversion table (BBEAR to PB IDs) ---------------------------

id.convert <- read_csv('data/isoforms/new_merge.combined.renamed.mapping.txt')

allTissue.counts <- as.data.frame(txi.kallisto$counts) %>% 
  rownames_to_column(var='old_id') %>% 
  left_join(id.convert) %>% 
  select(old_id,new_id,everything())

sum(is.na(allTissue.counts$new_id)) # Check that everything has a new ID

#