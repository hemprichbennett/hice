library(tidyverse)
library(seqinr)


library(tidyverse)
field_data <- read_csv('data/Edited_all_files_with_lat_long_VKedits.csv')

hice_faeces_ids <- field_data %>%
  filter(Species == 'Hice') %>%
  pull(Faeces_no1) %>%
  gsub('T', '', .)

post_qc_fa <- read.fasta('data/processed_dna_data/25_april_strict_lengths/for_r/91/all_post_QC.txt',
                         whole.header = T,
                         as.string = T)

seqname_tib <- tibble(seqnames = names(post_qc_fa),
                      seqs = NA) %>%
  separate(seqnames, into =c('seqnames', 'ncopies'), sep = '-') %>%
  mutate(ncopies = as.numeric(ncopies))

for(i in 1:length(post_qc_fa)){
  seqname_tib$seqs[i] <- as.character(post_qc_fa[[i]])
}


seqname_tib$seqnames <- gsub('\\.DAVE.+', '', seqname_tib$seqnames)
seqname_tib$seqnames
seqname_tib$seqnames <- gsub('GC.EC.7330.Dave.2.[1-9].', '', seqname_tib$seqnames)
seqname_tib$seqnames
seqname_tib$seqnames <- gsub('GC.EC.7330.Dave.2017.[1-9].', '', seqname_tib$seqnames)
seqname_tib$seqnames
seqname_tib$seqnames <- gsub('\\.Dave.+', '', seqname_tib$seqnames)
seqname_tib$seqnames
seqname_tib$seqnames <- gsub('_.+', '', seqname_tib$seqnames)
seqname_tib$seqnames
seqname_tib$seqnames <- gsub('.+\\.', '', seqname_tib$seqnames)
seqname_tib$seqnames
seqname_tib$seqnames <- gsub('T', '', seqname_tib$seqnames)


hice_tib <- filter(seqname_tib, seqnames %in% hice_faeces_ids)
sum(hice_tib$ncopies)
