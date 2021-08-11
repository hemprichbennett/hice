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
filenames <- list.files(pattern = '.csv', path = inpath)
#filenames <- filenames
filenames
filenames <- paste(inpath, filenames, sep = '')
#filenames <- filenames[grep('lulu', filenames)]

rawnets <- lapply(filenames, read.csv, header = F, stringsAsFactors = F, row.names=1)
names(rawnets) <- gsub('.*\\/', '', filenames)
names(rawnets) <- gsub('_.+', '', names(rawnets))
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
    newnam <- paste(nam, ' (', ncol(new_netlists[[z]][[i]]), ' bats)', sep = '')
    hice_list[[newnam]] <- as.numeric(c(ncol(new_netlists[[z]][[i]]), rowSums(new_netlists[[z]][[i]])))
    
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
  
# pdf('plots/inext_otus/n_samples.pdf')
# n_samples_plot
# dev.off()


melted_asymptote <- melt(asymptote_ests[,-c(2,5,6,7)], id.vars = c('Site', 'clustering','n_samples'))


melted_asymptote$variable %<>%
  gsub('Observed', 'Observed MOTU richness', .)%<>%
  gsub('Estimator', 'Estimated MOTU richness', .)%<>%
    gsub('percent_completeness', 'Estimated % completeness', .) %<>%
    gsub('N_samples_reqd', 'Number of samples required\nfor complete community', .)

melted_asymptote$Treatment <- 'LULU'

#write.csv(melted_asymptote, 'results/lulu_otu_inext.csv')

facetted <- ggplot(melted_asymptote, aes(x = clustering, y = value, colour = Site)) + 
  geom_point(size = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values =palette, name = 'Site')+
  labs(x= 'Clustering %', y = NULL)+ facet_wrap( ~ variable, scales = 'free')

facetted

# pdf('plots/inext_otus/facets.pdf')
# facetted
# dev.off()
  


palette <- c("#698fca",
             "#c29d31",
             "#7968d7",
             "#739f48",
             "#b963b7",
             "#4aaa86",
             "#c75980",
             "#c76347")


div_list <- list()

for(i in 1:length(all_inext)){
  n <- names(all_inext)[i]
  div_list[[n]] <- all_inext[[i]]$AsyEst
  div_list[[n]] <- div_list[[n]][which(div_list[[n]]$Diversity=='Shannon diversity'),]
  div_list[[n]]$clustering <- n
  div_list[[n]]$n_samples <- all_inext[[i]]$DataInfo$T
}

div_ests <- do.call(rbind, div_list)

ggplot(div_ests, aes(x=n_samples, y = Estimator, colour = Site, size=clustering))+ geom_point()


ggplot(div_ests, aes(x=clustering, y = Estimator, colour = Site))+ geom_point()+ 
  geom_smooth(method = lm, se = T)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_manual(values = palette)
  

div_ests$clustering <- as.numeric(div_ests$clustering)
library(lme4)
summary(lm(Estimator ~ clustering + n_samples + Site, data = div_ests))


sites <- c()
mins <- c()
maxs <- c()
for(i in 1:length(unique(div_ests$Site))){
  s <- as.character(unique(div_ests$Site)[i])
  minimum <- min(div_ests[which(div_ests$Site==s),'LCL'])
  maximum <- max(div_ests[which(div_ests$Site==s),'UCL'])
  sites <- c(sites, s)
  mins <- c(mins, minimum)
  maxs <- c(maxs, maximum)
}

df <- data.frame(sites, mins, maxs)
df_melt <- melt(df)
ggplot(df_melt, aes(x = variable, y = value, colour = sites)) + geom_point() + geom_abline()



combs <- expand.grid(df$sites, df$sites) # Make a df of all the pairwise combinations
colnames(combs) <- c('Site1', 'Site2')
combs <- combs[-which(combs$Site1==combs$Site2),] #Get rid of the self-pairs
combs$in_range <- 0
for(i in 1:nrow(combs)){
  s1 <- combs$Site1[i]
  s2 <- combs$Site2[i]
  s1min <- df$mins[which(df$sites==s1)]
  s1max <- df$maxs[which(df$sites==s1)]
  s2min <- df$mins[which(df$sites==s2)]
  s2max <- df$maxs[which(df$sites==s2)]
  
  if(s2max > s1min && s2max < s1max){
    combs$in_range[i] <- 1
  }
  if(s1max > s2min && s1max < s2max){
    combs$in_range[i] <- 1
  }
}

#Reorder the factors by sample size
ord <- gsub('.+ \\(', '', levels(combs$Site1))
ord <- gsub(' .+', '', ord)

ord2 <- gsub('.+ \\(', '', levels(combs$Site2))
ord2 <- gsub(' .+', '', ord2)

combs$Site1 <- factor(combs$Site1, levels(combs$Site1)[order(as.numeric(ord))])
combs$Site2 <- factor(combs$Site2, levels(combs$Site2)[order(as.numeric(ord2))])

ggplot(combs, aes(x=Site1, y=Site2, fill = in_range))+ geom_tile()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  

library(igraph)
library(qgraph)


#This bit removes the duplicate pairs, courtesy of https://stackoverflow.com/questions/29170099/remove-duplicate-column-pairs-sort-rows-based-on-2-columns?rq=1
combs[1:2] <- t( apply(combs[1:2], 1, sort) )
combs <- combs[!duplicated(combs[1:2]),]

#igraph_obj <- get.adjacency(g,sparse=FALSE)
#Trying to make a column of bat counts
sizes <- gsub('.+\\(', '', combs$Site1)
sizes <- gsub(' .+', '', sizes)

a <- qgraph(combs, layout = "spring", label.cex = 2, vsize = as.numeric(sizes))
