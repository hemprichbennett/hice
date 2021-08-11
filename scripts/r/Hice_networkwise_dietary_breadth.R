##################################################
## Project: Bat-diet
## Script purpose: calculating diversity etc in cervinus diet for ALL INDIVIDUALS AT SITE
## Date: 12/06/18
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################
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
  library(ggpubr)
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

#Currently the networks are named 'site, year'. For this version we just want 'site'

unique(colnames(all_interactions))

prey_data <- read.csv('data/taxonomy/order.csv')
colnames(prey_data) <- c('MOTU', 'Taxa')
prey_data$MOTU <- as.character(prey_data$MOTU)




###Add the taxonomic information to the interactions for everything
taxa_mat <- all_interactions[1,]
z <- 2
for(i in 1: nrow(all_interactions)){
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
#The colnames currently have the sample numbers, we want them to be a useable value in the data frame
colnames(taxa_mat) <- taxa_mat[1,]
taxa_mat <- taxa_mat[-1,]

#####Make a list with a network for each site####

sites_list <- list()
sites <- unique(colnames(all_interactions))




for(i in 1:length(sites)){
  m =  all_interactions[,which(colnames(all_interactions )==sites[i])]
  colnames(m) = m[1,]
  m = m[-c(1,2),]
  m <- m[-which(rowSums(m)==0),]
  sites_list[[i]] <- m
}
print(sites)
sites <- gsub('DANUM', 'Danum', sites)
sites <- gsub('MALIAU', 'Maliau', sites)


print(sites)
names(sites_list) <- sites




#####Do some ecology ####
hice_ecology <- matrix(nrow=0, ncol=1+ncol(specieslevel(matrix(sample(0: 1, size =100, replace = T), nrow = 10, ncol = 10), level = 'higher')))

for(i in 1:length(sites_list)){
  starttime <- Sys.time()
  if(length(which(duplicated(colnames(sites_list[[i]]))))>0){
    sites_list[[i]] <- sites_list[[i]][,-which(duplicated(colnames(sites_list[[i]])))]
  }
  splevel = specieslevel(sites_list[[i]], level = 'higher')
  endtime <- Sys.time()
  cat(names(sites_list)[i],'took', endtime-starttime,'\n')
  hice_ecology <- rbind(hice_ecology, cbind(rep(names(sites_list)[i], nrow(splevel)), splevel))
}
hice_ecology <- cbind(rownames(hice_ecology), hice_ecology)
# write.csv(hice_ecology, 'hice_ecology_1.csv')
# hice_ecology <- read.csv('hice_ecology_1.csv')
#####Add some taxonomic information
colnames(hice_ecology)[c(1:2)] <- c('Sample','Site')
t_taxa <- t(taxa_mat)
t_taxa <- cbind(t_taxa, rownames(t_taxa))
t_taxa <- as.data.frame(t_taxa)
colnames(t_taxa)[ncol(t_taxa)] <- 'Sample_no'

hice_ecology <- merge(hice_ecology, t_taxa, by.x='Sample' , by.y= 'Sample_no')

hice_ecology <- merge(x = hice_ecology, y = field_data, by.x = 'Sample', by.y = 'Faeces_no1')

write.csv(hice_ecology, 'hice_ecology_3.csv')
#####Work out the nestedness of each bat, then add it to the df#####



#####Make some graphs comparing sites and testing them with an anova ####

hice_ecology$habitat_type <- rep(NA, nrow(hice_ecology))
colnames(hice_ecology)[which(colnames(hice_ecology)=='Site.y')] <- 'Site'
hice_ecology$Site <- gsub('MALIAU', 'Maliau', hice_ecology$Site)
hice_ecology$Site <- gsub('DANUM', 'Danum', hice_ecology$Site)
hice_ecology$Site <- gsub('DVCA', 'Danum', hice_ecology$Site)
hice_ecology$Site <- gsub('MALUA', 'SBE', hice_ecology$Site)

hice_ecology$SiteAndYear <- gsub('MALIAU', 'Maliau', hice_ecology$SiteAndYear)
hice_ecology$SiteAndYear <- gsub('DANUM', 'Danum', hice_ecology$SiteAndYear)
hice_ecology$SiteAndYear <- gsub('DVCA', 'Danum', hice_ecology$SiteAndYear)
hice_ecology$SiteAndYear <- gsub('MALUA', 'SBE', hice_ecology$SiteAndYear)

hice_ecology[grep('SAFE', hice_ecology$Site), 'habitat_type'] <- 'Logged'
hice_ecology[grep('SBE', hice_ecology$Site), 'habitat_type'] <- 'Logged'
hice_ecology[grep('Danum', hice_ecology$Site), 'habitat_type'] <- 'Primary'
hice_ecology[grep('Maliau', hice_ecology$Site), 'habitat_type'] <- 'Primary'

#We need the number of nodes per network to be stored, for plotting later
hice_ecology$n_nodes <- rep(NA, nrow(hice_ecology))
for(i in 1:length(unique(hice_ecology$SiteAndYear))){
  var= unique(hice_ecology$SiteAndYear)[i]
  rows = which(hice_ecology$SiteAndYear==var)
  val = ncol(sites_list[[which(names(sites_list)==var)]])
  hice_ecology$n_nodes[rows] <- val
}


hice_ecology$Site <- fct_rev(hice_ecology$Site) #ggridges plots the factors in an annoying order, this rectifies it
#hice_ecology$SiteYearAndSex <- paste(hice_ecology$SiteAndYear, hice_ecology$Sex)
colnames(hice_ecology) <- gsub('PDI', 'Resource range', colnames(hice_ecology)) #As the interactions are binary what we output was resource range, not PDI

# #The docs for ridgeplots are here https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
# ridge
# pdf('plots/Site comparisons/normalised_degree_ridgeplot.pdf')
# ridge
# dev.off()

melted_hice <- melt(hice_ecology[,c('SiteAndYear', 'habitat_type', 'degree', 'normalised.degree', 'proportional.similarity', 'Resource range', 'n_nodes', 'Age', 'Reproductive_condition', 'Sex')],
                    id.vars = c( 'habitat_type', 'n_nodes', 'Age', 'Reproductive_condition', 'Sex', 'SiteAndYear'))
#The orders of the metrics aren't alphabetical, this alphabetises them and makes their names prettier
melted_hice$variable <- gsub('normalised.degree', 'Normalised degree', melted_hice$variable)
melted_hice$variable <- gsub('partner.diversity','Partner diversity', melted_hice$variable)
melted_hice$variable <- gsub('^degree$','Degree', melted_hice$variable)
melted_hice$variable <- gsub('proportional.similarity','Proportional similarity', melted_hice$variable)
melted_hice$variable <- gsub('PDI','Resource range', melted_hice$variable)
melted_hice$variable <- gsub('closeness','Closeness centrality', melted_hice$variable)
melted_hice$variable <- ordered(melted_hice$variable, levels=unique(melted_hice$variable)[order(as.character(unique(melted_hice$variable)))])



melted_hice$for_labels <- paste(melted_hice$Site, ' (', melted_hice$n_nodes, ' bats)', sep = '')

facet_ridge <- ggplot(melted_hice, aes (y=fct_rev(for_labels), x =value, fill=habitat_type)) + 
  geom_density_ridges(scale= 0.5)+ #The scale determines the space between the rows
  theme_ridges()+ #This changes the theme to make it more aesthetically pleasing
  scale_fill_cyclical(values = c("#85d7da","#d0ca9f"), guide = 'legend', name = 'Habitat type')+
  scale_x_continuous(expand = c(0.01, 0)) + #Make the space between the labels and plot smaller
  scale_y_discrete(expand = c(0.01, 0))+ #Make it so the top series actually fits in the plot
  ylab(NULL)+ xlab(NULL)+
  facet_wrap( ~ variable, ncol=2, scales = 'free_x', strip.position = 'bottom')+ #free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"),#strip stuff sorts the facet labels, spacing adjusts the space between facets
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        text = element_text(size=12))+
  theme(legend.text =element_text(size = 10), 
        legend.position="bottom") #Trying to standardise the sizes across both plots
facet_ridge

pdf('plots/Hice/networks_facet_ridgeplot.pdf', width = 7, height = 9)
facet_ridge
dev.off()


#####Now comes the significance testing ####
anova_list <- list()
for(i in 1:length(unique(melted_hice$variable))){ #Could probably just use lapply now that a load of the loop has been killed
  v <- as.character(unique(melted_hice$variable)[i])
  s <- melted_hice[which(melted_hice$variable==v),]
  my_anova <- aov(s$value ~ s$SiteAndYear)
  summary(my_anova)
  anova_list[[i]] <- TukeyHSD(my_anova)
  names(anova_list)[i] <- v
}

mets <- c()
net1s <- c()
net2s <- c()
signif <- c()
for(l in 1:length(anova_list)){
  lis <- anova_list[[l]]$`s$SiteAndYear`
  for(i in 1:nrow(lis)){
    net1 <- strsplit(rownames(lis),split = '-')[[i]][1]
    net2 <- strsplit(rownames(lis),split = '-')[[i]][2]
    
    if(lis[i,4]<=0.05){
      sig <- 's'
    }else{sig <- 'n-s'}
    
    net1s <- c(net1s, net1)
    net2s <- c(net2s, net2)
    mets <- c(mets, names(anova_list)[l])
    signif <- c(signif, sig)
  }
}
sig_df <- data.frame(as.factor(mets), as.factor(net1s), as.factor(net2s), as.factor(signif))
colnames(sig_df) <- c('metrics', 'network_1', 'network_2', 'significance')

for(i in 1:length(unique(sig_df$metrics))){
  var <- unique(sig_df$metrics)[i]
  temp <-  sig_df[which(sig_df$metrics==var),]
  temp <- temp[,-1]
  mat <- dcast( temp, network_1 ~ network_2, value.var = 'significance')
  write.csv(mat, paste('results/tukey/', var, '.csv', sep = ''), row.names = FALSE)
}

temp <-  sig_df[which(sig_df$metrics=='Degree'),]
temp <- temp[,-1]
dcast( temp, network_1 ~ network_2, value.var = 'significance')

mat <- matrix(0, 4, 4)
mat[temp[,1:2]] <- temp[,3]
mat

sig_df$network_2 <- fct_rev(sig_df$network_2) #ggridges plots the factors in an annoying order, this rectifies it

write.csv(sig_df, 'results/full_size_sig.csv')


melted_hice[which(!melted_hice$Age %in% c('A', 'J')), 'Age'] <- NA
melted_hice$Age <- as.factor(melted_hice$Age)
melted_hice[which(!melted_hice$Sex %in% c('M', 'F')), 'Sex'] <- NA
melted_hice$Sex <- factor(melted_hice$Sex)
melted_hice$Reproductive_condition <- gsub("LA", "L", melted_hice$Reproductive_condition)
melted_hice[which(!melted_hice$Reproductive_condition %in% c('L', 'PL', 'PR')), 'Reproductive_condition'] <- 'NR'
melted_hice$Site <- gsub('DVCA', 'Danum', melted_hice$Site)
melted_hice$Site <- gsub('DANUM', 'Danum', melted_hice$Site)
melted_hice$Site <- gsub('MALIAU', 'Maliau', melted_hice$Site)
melted_hice$Site <- gsub('MALUA', 'SBE', melted_hice$Site)
melted_hice$Reproductive_condition <- factor(melted_hice$Reproductive_condition)
melted_hice$Site <- factor(melted_hice$Site)
melted_hice$habitat_type <- factor(melted_hice$habitat_type)
mets <- c('Degree', "Normalised degree", 'Resource range', 'Proportional similarity')
# 
# 
#Do a multiple regression per network
 for(i in 1:length(mets)){
   
   met = mets[i]
   mod <- lm(value ~ Sex + Age + Reproductive_condition + habitat_type,  data = melted_hice[which(melted_hice$variable==mets[i]),])
   temp_df <- hice_ecology[which(!is.na(hice_ecology$Sex)),]
   sink(paste('data/output_data/hice_stats/', met, '.txt', sep = ''))
   print(summary(mod))
   sink()
 }
#Now do stepwise model selection
library(MASS)
for(i in 1:length(mets)){
  
  met = mets[i]
  
  
  for_step <- melted_hice[which(melted_hice$variable==met),]
  badrows <- c()
  for(r in 1:nrow(for_step)){
    if(NA %in% unlist(for_step[r,])){
      badrows <- c(badrows, r)
    }
  }
  for_step <- for_step[-badrows,]
  temp_degree <- lm(value ~ Sex + Age + Reproductive_condition + habitat_type,  data = for_step)
  
  
  step.model <- stepAIC(temp_degree, direction = "backward", 
                        trace = FALSE)

  mod <- lm(value ~ Sex + Age + Reproductive_condition + habitat_type,  data = melted_hice[which(melted_hice$variable==mets[i]),])
  temp_df <- hice_ecology[which(!is.na(hice_ecology$Sex)),]
  sink(paste('data/output_data/hice_stats/', met, '_step.txt', sep = ''))
  print(summary(step.model))
  sink()
}




melted_degree <- melted_hice[which(melted_hice$variable=='Degree'),]
melted_prop_sim <- melted_hice[which(melted_hice$variable=='Proportional similarity'),]
#1-off ANOVA on degree
# fit_both_sexes <- lm(log(value) ~ Sex + Age + Year +  Site, data = melted_degree)
# fitrep <- lm(log(value) ~Year + Site + Reproductive_condition,  data = melted_degree[which(melted_degree$Sex=='F'),])
# 

degree_df <- melted_hice[which(melted_hice$variable=='Degree'),]
badrows <- c()
for(r in 1:nrow(degree_df)){
  if(NA %in% unlist(degree_df[r,])){
    badrows <- c(badrows, r)
  }
}
degree_df <- degree_df[-badrows,]
temp_degree <- lm(value ~ Sex + Age + Reproductive_condition + habitat_type,  data = degree_df)

library(MASS)
step.model <- stepAIC(temp_degree, direction = "backward", 
                            trace = FALSE)
summary(step.model)

####Look at correlations ####

taxa_for_cor <- t(taxa_mat)
taxa_for_cor <- taxa_for_cor[which(rowSums(taxa_for_cor)>2),]
taxa_for_cor <- taxa_for_cor[,which(colSums(taxa_for_cor)>20)]
cormat <- round(cor(taxa_for_cor),2)
res1 <- cor.mtest(taxa_mat, conf.level = .99)



###Here we use pearsons similarity to plot only the significantly correlated points, 
#ordered along the angular order of the eigenvectors
pdf('plots/Hice/network_taxa_correlations.pdf')
taxa_corr <- corrplot(cormat, method = "circle", p.mat = res1$p, sig.level = .05, type = 'upper', order = 'AOE',
                      tl.col = "black", tl.srt = 45, insig = 'blank', col = c('black', 'white'),
                      bg = "lightblue",
                      cl.pos = "b")

dev.off()



hice_stats <- hice_ecology[,seq(4,22)]
#I get rid of these values here as they never change and so are worthless for this
hice_stats <- hice_stats[,-which(colnames(hice_stats) %in% c("node.specialisation.index.NSI", "betweenness", "weighted.betweenness"))]

hice_cor <- round(cor(hice_stats),2)
stat_sig <- cor.mtest(hice_stats, conf.level = .95)


#pdf('plots/Hice/metric_correlations.pdf')
#corrplot(hice_cor, method = "circle", p.mat = stat_sig$p, sig.level = .05, type = 'upper', order = 'AOE',
#         tl.col = "black", tl.srt = 45, insig = 'blank')

#dev.off()




tax_df <- hice_ecology[,c(1,2 , seq(23,47))]
tax_df <- melt(tax_df, id.vars = c('Sample', 'Site.x'))
colnames(tax_df)[c(3,4)] <- c('Order', 'Present/absent')
tax_df$`Present/absent` <- as.integer(as.character(tax_df$`Present/absent`))
tax_df$`Present/absent` <- ifelse(tax_df$`Present/absent`== 0, 0, 1)

prop_present <- sapply(unique(tax_df[,c('Site.x', 'Order')]),  function(x) as.character(x))
prop <- c()
nbats <- c()
#for(i in 1: 1){
for(i in 1: nrow(prop_present)){
  tem <- tax_df[which(tax_df$Site== as.character(prop_present[i,1])
                      & tax_df$Order==as.character(prop_present[i,2])),]
  prop <- c(prop, sum(tem$`Present/absent`)/nrow(tem)) #The number of bats that consumed the order, divided by total bats
  nbats <- c(nbats, nrow(tem))
}


prop_present <- cbind(prop_present, prop, nbats)
prop_present <- as.data.frame(prop_present)
prop_present$prop <- as.numeric(as.character(prop_present$prop))
prop_present$nbats <- as.integer(as.character(prop_present$nbats))

prop_present$treatment_and_site <- as.character(prop_present$Site.x) %>%
  gsub('Danum', 'Old-growth: Danum', .) %>%
  gsub('Maliau', 'Old-growth: Maliau', .) %>%
  gsub('SAFE', 'Logged: SAFE', .) %>%
  gsub('SBE', 'Logged: SBE', .)

wide_tiles <- ggplot(data = prop_present[which(prop_present$prop!=0),], aes(x = Order, y =treatment_and_site)) +
  geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", name ='Proportion of\nindividuals\ncontaining\norder', 
                       midpoint = 0.5, 
                       breaks= c(0.001,0.5,1),
                       labels = c(0.001,0.5,1),
                       limits = c(0.001,1)
                       ) +
  labs(y ="Site", x = 'Prey order')+
  theme_dark(base_size = 14)+
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(colour=, size=15),
        legend.title = element_text(size=15, face="bold"),
        legend.position="bottom")#+
  # guides(fill=guide_legend(title="Proportion of\nindividuals containing\norder"),
  #       strip.background =element_rect(fill="black", size = 3))
wide_tiles

pdf('plots/network_proportion_of_individuals_containing_wide.pdf', width = 17)
wide_tiles
dev.off()

jpeg('plots/network_proportion_of_individuals_containing_wide.jpg', width = 17, height = 7, units = 'in', res = 300)
wide_tiles
dev.off()