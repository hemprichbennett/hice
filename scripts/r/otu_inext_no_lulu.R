#### Header ####
## Project:
## Script purpose:
## Date:
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

#Make a nested list, where the main list has a value for each clustering level, and each item within it is a distinct network generated at that level #####

library(here)
library(bipartite)
library(ggplot2)
library(iNEXT)
library(DataExplorer)
library(magrittr)

source('scripts/r/r_network_gen.r')



setwd(here())
# getwd()

inpath <- 'data/processed_dna_data/lulu/'
filenames <- list.files(pattern = '*binary*', path = inpath, recursive = T)
#filenames <- filenames
filenames
filenames <- paste(inpath, filenames, sep = '')
#filenames <- filenames[grep('lulu', filenames)]

rawnets <- lapply(filenames, read.table, header = F, stringsAsFactors = F, row.names=1)
names(rawnets) <- gsub('data/processed_dna_data/lulu/', '', filenames)
names(rawnets) <- gsub('/.+', '', names(rawnets))
for(i in 1:length(rawnets)){
  all_interactions <- rawnets[[i]]
  all_interactions[1,1] <- 'MOTU'
  all_interactions[1,] <- gsub('X', '', all_interactions[1,])
  all_interactions[2:nrow(all_interactions),2:ncol(all_interactions)] <- ifelse(all_interactions[2:nrow(all_interactions), 2:ncol(all_interactions)]==0,0,1)
  rawnets[[i]] <- all_interactions
}


netlists <- lapply(rawnets, function(x) r_network_gen(input = x, collapse_species = F, desired_species='Hice', include_malua = T))

names(netlists) <- names(rawnets)

#Make them by site, not site and year
for(i in 1:length(netlists)){
  #colnames(netlists[[i]]) <- gsub(',.+', '', colnames(netlists[[i]]))
  #Get rid of the superfluous first two rows
  netlists[[i]] <- netlists[[i]][-c(1,2),]}

#Get rid of the redundant first two rows
#for(i in 1:length(netlists)){netlists[[i]] <- netlists[[i]][-c(1,2),]}


#Each item of netlists is a matrix of interactions, where the colnames are where each individual cervinus
#was caught. Make a nested list so we can compare sites
new_netlists <- list()
for(i in 1:length(netlists)){
  clust <- names(netlists)[i]
  new_netlists[[clust]] <- list()
  print(clust)
  for(a in 1:length(unique(colnames(netlists[[i]])))){
    site <- unique(colnames(netlists[[i]]))[a]
    desired_cols <- which(colnames(netlists[[i]])==site)
    #print(desired_cols)
    new_netlists[[clust]][[site]] <- netlists[[i]][,desired_cols]
    new_netlists[[clust]][[site]] <- ifelse(new_netlists[[clust]][[site]]==0, 0, 1)
    #print(site)
  }
}
#plot_str(new_netlists)

inext_in <- list()
for(z in 1:length(new_netlists)){
  clustering <- names(new_netlists)[z]
  hice_list <- list()
  for(i in 1:length(new_netlists[[z]])){
    nam <- names(new_netlists[[z]])[i]
    hice_list[[nam]] <- as.numeric(c(ncol(new_netlists[[z]][[i]]), rowSums(new_netlists[[z]][[i]])))
    
  }
  inext_in[[clustering]] <- hice_list
}

#temp <- hice_list$SAFE

all_inext <- lapply(inext_in, function(x) iNEXT(x, datatype = 'incidence_freq'))



asy_list <- list()

for(i in 1:length(all_inext)){
  n <- names(all_inext)[i]
  asy_list[[n]] <- all_inext[[i]]$AsyEst
  asy_list[[n]] <- asy_list[[n]][which(asy_list[[n]]$Diversity=='Species richness'),]
  asy_list[[n]]$clustering <- n
  asy_list[[n]]$n_samples <- all_inext[[i]]$DataInfo$T
}

asymptote_ests <- do.call(rbind, asy_list)


####Rearrange our stats ####
#asymptote_ests <- inc_all$AsyEst
#asymptote_ests <- asymptote_ests[asymptote_ests$Diversity=='Species richness',]
asymptote_ests$Site <- as.character(asymptote_ests$Site)
#asymptote_ests <- cbind(asymptote_ests, inc_all$DataInfo$T)
#asymptote_ests <- asymptote_ests[order(asymptote_ests$Site),]
#colnames(asymptote_ests)[8] <- 'number of samples'
asymptote_ests$percent_completeness <- (asymptote_ests$Observed*100)/asymptote_ests$Estimator
asymptote_ests$N_samples_reqd <- (asymptote_ests$`n_samples`/asymptote_ests$percent_completeness)*100
#asymptote_ests <- asymptote_ests[,c(1,3,4,9,8,10,5,6,7)]

write.csv(asymptote_ests, 'results/non_lulu_otu_inext_non_melted.csv')

palette <- c("#8960b3",
             "#56ae6c",
             "#ba495b",
             "#b0923b")

asymptote_ests$Site %<>%
  gsub('DANUM', 'Danum', .)%<>%
  gsub('MALIAU', 'Maliau', .)

n_samples_plot <- ggplot(asymptote_ests, aes(x = clustering, y = N_samples_reqd, colour = Site)) + geom_point(size = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values =palette, name = 'Site')+
  labs(x= 'Clustering %', y = 'Number of samples required')

pdf('plots/inext_otus/n_samples.pdf')
n_samples_plot
dev.off()


melted_asymptote <- melt(asymptote_ests[,-c(2,5,6,7)], id.vars = c('Site', 'clustering','n_samples'))


melted_asymptote$variable %<>%
  gsub('Observed', 'Observed MOTU richness', .)%<>%
  gsub('Estimator', 'Estimated MOTU richness', .)%<>%
  gsub('percent_completeness', 'Estimated % completeness', .) %<>%
  gsub('N_samples_reqd', 'Number of samples required\nfor complete community', .)

melted_asymptote$Treatment <- 'non_LULU'

write.csv(melted_asymptote, 'results/non_lulu_otu_inext.csv')

facetted <- ggplot(melted_asymptote, aes(x = clustering, y = value, colour = Site)) + 
  geom_point(size = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values =palette, name = 'Site')+
  labs(x= 'Clustering %', y = NULL)+ facet_wrap( ~ variable, scales = 'free')

facetted

# pdf('plots/inext_otus/facets.pdf')
# facetted
# dev.off()
