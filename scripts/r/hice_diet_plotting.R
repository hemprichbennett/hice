####This script makes dietary barplots, and runs a chi-squared test to look at the difference between Hice diet
####Between networks

library(stringr)
library('here')
library('bipartite')
library(iNEXT)
library(tidyverse)
library(reshape2)


setwd(here())
getwd()
#####Load in the data, reformat it to make it more useful ####

#all_interactions <- read.table(here('data/processed_dna_data/galaxy_r_workflow/95/all_post_QC_otus.txt.table_binary.out'), sep = '\t', header = F, stringsAsFactors = F, row.names = 1)
all_interactions <- read.csv('data/processed_dna_data/lulu/95/lulu_95.csv', header = F, stringsAsFactors = F)#
all_interactions[1,1] <- 'MOTU'
all_interactions[1,] <- gsub('X', '', all_interactions[1,])
all_interactions[2:nrow(all_interactions),2:ncol(all_interactions)] <- ifelse(all_interactions[2:nrow(all_interactions), 2:ncol(all_interactions)]==0,0,1)


#This finds all samples with 'GC' in the name and gives them a useful name
gc <- grep('GC',all_interactions[1,])
for(i in 1:length(gc)){
  #print(str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7])
  temp <- str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7]
  all_interactions[1,gc[i]] <- str_split(temp, pattern='_')[[1]][1]
}

#This finds all samples without 'GC' in the name and gives them a useful name
non_gc <- seq(1, ncol(all_interactions))[-gc]
str_split(all_interactions[1,non_gc[2]], pattern = '\\.')[[1]][1] #Finish this
for(i in 1:length(non_gc)){
  all_interactions[1,non_gc[i]] <- str_split(all_interactions[1,non_gc[i]], pattern = '\\.')[[1]][1]
}

badcols <- c('1774','4437', '2070', '2275')#Sadly these columns match two different samples, so must be removed for now until checked against the field data

all_interactions <- all_interactions[,-which(all_interactions[1,] %in% badcols)]


field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
field_data$Site <- gsub('DVCA', 'Danum', field_data$Site)
field_data$Site <- gsub('DANUM', 'Danum', field_data$Site)
field_data$Site <- gsub('MALIAU', 'Maliau', field_data$Site)
field_data$Site <- gsub('MALUA', 'SBE', field_data$Site)
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)

#####This loop filters for hice and makes a list entry for each Site ####
batnames <- c()
siteslist <- list()
for(i in 1:ncol(all_interactions)){
  batname <- all_interactions[1,i]
  if(!batname %in%field_data$Faeces_no1 && !batname %in% field_data$Faeces_no2){
    next()
  }
  if(batname %in% field_data$Faeces_no1){
    site <- field_data[which(field_data$Faeces_no1==batname), 'SiteAndYear']
    if(length(which(field_data$Faeces_no1==batname))>1){next()}#If theres more than one match
    if(field_data[which(field_data$Faeces_no1==batname), 'Species']!= 'Hice'){#Only allow hice into the upcoming list
      next()
    }
    
  }else if(batname %in% field_data$Faeces_no2){
    site <- field_data[which(field_data$Faeces_no2==batname), 'SiteAndYear']
    if(length(which(field_data$Faeces_no2==batname))>1){next()}
    if(field_data[which(field_data$Faeces_no2==batname), 'Species']!= 'Hice'){#Only allow hice into the upcoming list
      next()
    }
  
    }
  if(site %in% names(siteslist)){
    pos <- which(names(siteslist) == site)
    siteslist[[pos]] <- cbind(siteslist[[pos]], as.numeric(all_interactions[,i]))
  }else{
    siteslist[[site]] <- matrix(nrow = nrow(all_interactions), ncol = 1, as.numeric(all_interactions[,i]))
    rownames(siteslist[[site]]) <- all_interactions[,1]
  }
}


#The bat names are still in row 1
for(i in 1:length(siteslist)){
  colnames(siteslist[[i]]) <- siteslist[[i]][1,]
  siteslist[[i]] <- siteslist[[i]][-1,]
  print(names(siteslist)[i])
  print(ncol(siteslist[[i]]))
}
#####Add the taxonomic information, format the data #####
orderdata <- read.csv('data/taxonomy/order.csv', stringsAsFactors = F)
familydata <- read.csv('data/taxonomy/family.csv', stringsAsFactors = F)

iNEXT_making <- function(prey_data, sites_list){
  colnames(prey_data) <- c('MOTU', 'Taxa')
  taxa_list <- list()
  
  #Put the taxonomic information onto the sites data
  
  for(b in 1: length(siteslist)){
    all_interactions <- siteslist[[b]]
    taxa_mat <- matrix(nrow=0, ncol=ncol(all_interactions))
    
    z <- 1
    for(i in 1: nrow(all_interactions)){ #Make a matrix of simplified interactions
      rowname = rownames(all_interactions)[i]
      if(!rowname %in% prey_data$MOTU){
        next()
      }
      tax = as.character(prey_data[which(prey_data$MOTU == rowname),'Taxa'])
      if(is.null(nrow(taxa_mat))){ #If its the first iteration there won't be any rownames yet, so the next if statement will fail
        taxa_mat <- rbind(taxa_mat, all_interactions[i,])
        rownames(taxa_mat)[z] <- tax
        z <- z+1
      }
      if(tax %in% rownames(taxa_mat)){
        to_merge = which(rownames(taxa_mat)==tax)
        taxa_mat[to_merge,] <- taxa_mat[to_merge,]+ all_interactions[i,]
      }else{
        taxa_mat <- rbind(taxa_mat, all_interactions[i,])
        rownames(taxa_mat)[z] <- tax
        z <- z+1
      }
    }
    taxa_list[[b]] <- taxa_mat
    names(taxa_list)[b] <- names(siteslist)[b]
  }
  return(taxa_list)

}

orderlist <- iNEXT_making(prey_data = orderdata, sites_list = siteslist)
prop_list <- list()
for(i in 1:length(orderlist)){
  prop_list[[i]] <- data.frame(rowSums(orderlist[[i]])/sum(orderlist[[i]]), rownames(orderlist[[i]]))
}
names(prop_list) <- names(orderlist)



#####Format the data a bit for the chi-squared etc ####
molten_proportions <- melt(prop_list )
molten_proportions <- molten_proportions[,-2] #This column sucks, get rid of it
colnames(molten_proportions) <- c('prey_taxa', 'value', 'Site')

#Annoyingly having melted the data we must now cast it to form our contingency table
hice_wide <- acast(molten_proportions, Site ~ prey_taxa)
hice_wide[is.na(hice_wide)] <- 0
hice_wide <- t(hice_wide) #Has to be the other way round for the chi-squared test
chisq <- chisq.test(hice_wide)
chisq
chisq$expected




#PCA

for_pca <- t(sapply(orderlist, function(x) rowSums(x)))

for_pca <- for_pca[,-which(apply(for_pca, 2, var)==0)] #Get rid of columns with zero variance, they add nothing and screw
#up the scaling

my_prcomp <- prcomp(for_pca, scale = T)

library(factoextra)

powers <- get_eigenvalue(my_prcomp)

write.csv(powers, 'results/PCA.csv')

fviz_eig(my_prcomp, addlabels = TRUE, ylim = c(0, 70))

z <- fviz_pca_ind(my_prcomp,
             repel = TRUE
)

z <- z+ ggtitle(NULL)+ xlab(paste('PC1 (', round(powers$variance.percent[1], 1), '%)', sep = ''))+
  ylab(paste('PC2 (', round(powers$variance.percent[2], 1), '%)', sep = ''))

pdf('plots/pca/sites.pdf')
z
dev.off()

fviz_pca_var(my_prcomp,
             col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_contrib(my_prcomp, choice = "var", axes = 1, top = 10)
fviz_contrib(my_prcomp, choice = "var", axes = 2, top = 10)



var <- get_pca_var(my_prcomp)
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)
corrplot(var$contrib, is.corr=FALSE)    

ind <- get_pca_ind(my_prcomp)

#Make the dataframe for NMDS and simper etc
taxa <- as.data.frame(matrix(ncol = 1+ nrow(orderlist[[1]]),nrow = 0))
names(taxa) <- c('site', rownames(orderlist$SAFE))

for(i in 1:length(orderlist)){
  taxa <- rbind(taxa, cbind(names(orderlist)[i], t(orderlist[[i]])))
}
colnames(taxa)[1] <- 'site'
# 
taxa[,2:ncol(taxa)] <- apply(taxa[,2:ncol(taxa)], 2, as.numeric)

#remove columns with no interactions
badcols <- names(which(colSums(taxa[,2:ncol(taxa)])==0))
taxa[badcols] <- NULL

temp_for_MDS <- taxa[,2:ncol(taxa)]

taxa <- taxa[-which(rowSums(temp_for_MDS)==0),]

#removing the bat which consumes only a mite and so skews everything hugely
taxa <- taxa[-which(rownames(taxa)=='1473'),]

ECHOtest2=adonis(taxa[,2:ncol(taxa)]~taxa$site)
ECHOtest2
summary(ECHOtest2)

ECHOtest3=simper(taxa[,2:ncol(taxa)], taxa$site, permutations = 1000)

summary(ECHOtest3, ordered = TRUE, digits = max(3,getOption("digits") - 3))

ECHOtest3

#Get all the significant simper values into one place
simper_outputs <- summary(ECHOtest3, ordered = TRUE, digits = max(3,getOption("digits") - 3))
library(reshape2)
simper_simplified <- simper_outputs[[1]][FALSE,]
namesvec <- c()
for(i in 1:length(simper_outputs)){
  temp <- simper_outputs[[i]]
  temp_2 <- temp[which(temp$p < 0.05),]
  simper_simplified <- rbind(simper_simplified, temp_2)
  if(nrow(temp_2)> 0){
    namesvec <- c(namesvec, rep(names(simper_outputs)[i], nrow(temp_2)))
  }
}

#Add site names and taxa names
simper_simplified$sites <- namesvec
simper_simplified$taxa <- gsub('[1-9]', '', rownames(simper_simplified))

#Rearrange columns
simper_simplified <- simper_simplified[,c(8, 9, seq(1,7))]
simper_simplified$sites <- gsub('_', '-', simper_simplified$sites)

#simper_simplified <- simper_simplified[,c(1,2,8, 9)]
simper_simplified[,c(3:9)] <- round(simper_simplified[,c(3:9)], digits = 3)

write.csv(simper_simplified, 'results/significant_simper.csv')


library(vegan)

###### Intra-site NMDS #####
#Make a list of communities by species matrixes, where each list item is a site (exclude SBE as it was only sampled in one year)

sites <- unique(gsub(',.+', '', taxa$site))
sites <- sites[-which(sites=='SBE')]

sites_for_MDS <- list()
intrasite_MDS <- list()

set.seed(2)


not_SBE <- taxa[-grep('SBE', taxa$site),]

not_SBE_MDS <- metaMDS(not_SBE[,2:ncol(not_SBE)], # Our community-by-species matrix
                   trymax = 200,
        k=2) # The number of reduced dimensions

data.scores <- as.data.frame(scores(not_SBE_MDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$SiteAndYear <- not_SBE$site  #  add the Site variable created earlier
data.scores$Site <- gsub(',.+', '', data.scores$SiteAndYear)
data.scores$Year <- gsub('.+,', '', data.scores$SiteAndYear)

intrasite_plot <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2, col=Year)) +
  geom_point(aes(size = 1.8)) +
  stat_ellipse() +
  theme_bw()+ facet_wrap(. ~ Site)+
  scale_color_viridis_d() + 
  theme(text = element_text(size = 20),
        legend.position = 'bottom')+
  # remove point size from legend
  guides(size = FALSE)

intrasite_plot

pdf('plots/NMDS/intrasite.pdf', width = 14)
intrasite_plot
dev.off()

jpeg('plots/NMDS/intrasite.jpeg', width = 14, height = 7, units = 'in', res = 300)
intrasite_plot
dev.off()

##### 2016-only NMDS #####

#get rid of the non-2016 stuff
sixteen_sites <- taxa[grep('2016', taxa$site),]

sixteen_MDS <- metaMDS(sixteen_sites[,2:ncol(sixteen_sites)], # Our community-by-species matrix
               k=2) # The number of reduced dimensions

#stuff copied and modified from https://chrischizinski.github.io/rstats/vegan-ggplot2/
sixteen_scores <- as.data.frame(scores(sixteen_MDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame

sixteen_scores$Site <- sixteen_sites$site  #  add the Site variable created earlier


sixteen_MDS_plot <- ggplot(sixteen_scores, aes(x=NMDS1, y=NMDS2, col=Site)) +
  geom_point(aes(size = 1.8)) +
  stat_ellipse() +
  theme_bw()+ 
  scale_color_viridis_d() +
  theme(text = element_text(size = 20), 
        legend.position = 'bottom')+
  # remove point size from legend
  guides(size = FALSE)

sixteen_MDS_plot
  
pdf('plots/NMDS/2016.pdf')
sixteen_MDS_plot
dev.off()

jpeg('plots/NMDS/2016.jpeg', width = 9, height = 7, units = 'in', res = 300)
sixteen_MDS_plot
dev.off()



# Clump all years together,  so its just site -----------------------------


#get rid of the year identifier

clumped_sites <- taxa

#clumped_sites$site <- gsub(',.+', '',  clumped_sites$site)

clumped_MDS <- metaMDS(clumped_sites[,2:ncol(clumped_sites)], # Our community-by-species matrix
                       k=2, # The number of reduced dimensions
                       trymax = 100)

#stuff copied and modified from https://chrischizinski.github.io/rstats/vegan-ggplot2/
clumped_scores <- as.data.frame(scores(clumped_MDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame

clumped_scores$Site <- clumped_sites$site  #  add the Site variable created earlier


clumped_MDS_plot <- ggplot(clumped_scores, aes(x=NMDS1, y=NMDS2, col=Site)) +
  geom_point() +
  stat_ellipse() +
  theme_bw()+ 
  scale_color_viridis_d() +
  theme(text = element_text(size = 20), 
        legend.position = 'bottom')+
  # remove point size from legend
  guides(size = FALSE)

clumped_MDS_plot

