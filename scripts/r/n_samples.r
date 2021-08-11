if(interactive()==TRUE){
  library('here')
  library(ggplot2)
  library(reshape2)
  library(tidyverse)
  library(ggridges)
  library(gridExtra)
  library(forcats)
  library(corrplot)
  library(DataExplorer)
  
}else{
  library(here, lib.loc = '/data/home/btw863/r_packages/')
}

setwd(here())
source('scripts/r/r_network_gen.r')

field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
field_data_2 <- field_data
#The bit below deals with the inevitable merging issues with some samples being from faeces_no1, some from faeces_no2
field_data$Faeces_no2 <- NULL

field_data_2$Faeces_no1 <- field_data_2$Faeces_no2
field_data_2$Faeces_no2 <- NULL

field_data <- rbind(field_data, field_data_2)
all_interactions <- r_network_gen(collapse_species = F, desired_species = 'Hice', include_malua = T, lulu = T)

# get the colnames from all_interactions, as they contain the site and year
bats <- colnames(all_interactions) %>%
  as_tibble() %>%
  separate(col = value, sep = ', ', into = c('site', 'year')) %>%
  table(.)

write.csv(bats, file = 'results/captures.csv')
