library(stringr)
library('here')
library('bipartite')
library(iNEXT)
library(ggplot2)
source(here('scripts/r/iNEXT_prep.R'))
# source('scripts/r/hernani_comparisons.R')

setwd(here())
getwd()
all_interactions <- read.csv('data/processed_dna_data/lulu/95/lulu_95.csv', header = F, stringsAsFactors = F)#
all_interactions[1, 1] <- "MOTU"
all_interactions[1, ] <- gsub("X", "", all_interactions[1, 
                                                        ])
all_interactions[2:nrow(all_interactions), 2:ncol(all_interactions)] <- ifelse(all_interactions[2:nrow(all_interactions), 
                                                                                                2:ncol(all_interactions)] == 0, 0, 1)


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

badcols <- c('1774','4437', '2070', '2275', '4260', '4531')#Sadly these columns match two different samples, so must be removed for now until checked against the field data

all_interactions <- all_interactions[,-which(all_interactions[1,] %in% badcols)]



#all_interactions_with_extra <- rbind(all_interactions[1,], all_interactions)
#I'm putting in an extra row here with the colnames so we can retain the sample numbers without having to rewrite a load of the.matrix.reloaders code
field_data <- read.csv(here('data/Edited_all_files_with_lat_long_VKedits.csv'))
field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)
field_data$SiteAndYear <- paste(field_data$Site, field_data$Year)

all_hice <- iNEXT.prep(master.data = field_data, ID.column.1 = "Faeces_no1", ID.column.2 = "Faeces_no2",species.column = "SiteAndYear", split.by.column = "Species", split.by.var = "Hice", OTU.matrix = all_interactions)

all_hice <- all_hice[,which(all_hice[1,]>7)]
hice_list <- list()
for (i in 1: ncol(all_hice)){
  hice_list[[i]] <- all_hice[which(rowSums(all_hice)>0),i]
  
}
names(hice_list) <- colnames(all_hice)
names(hice_list) <- gsub(' ', ', ', names(hice_list))
names(hice_list) <- gsub('DANUM', 'Danum', names(hice_list))
names(hice_list) <- gsub('DVCA', 'Danum', names(hice_list))
names(hice_list) <- gsub('MALUA', 'SBE', names(hice_list))
names(hice_list) <- gsub('MALIAU', 'Maliau', names(hice_list))


a <- iNEXT(hice_list, datatype = 'incidence_freq')

#####Here I replicate the gginext command, but tweak it a bit because it doesn't natively allow it ####

z <- fortify.iNEXT(a)
z$site_type <- rep(NA, nrow(z))
z[grep('SAFE', z$site),'site_type'] <- 'Logged'
z[grep('SBE', z$site),'site_type'] <- 'Logged, replanted'
z[grep('Danum', z$site),'site_type'] <- 'Primary'
z[grep('Maliau', z$site),'site_type'] <- 'Primary'
z$site_type <- as.factor(z$site_type)
type = 1
se = TRUE
facet.var = "site"
color.var = "site"
grey = T

if(ncol(z) ==7) {se <- FALSE}
datatype <- unique(z$datatype)
if(color.var=="none"){
  if(levels(factor(z$order))>1 & "site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
    color.var <- "both"
    z$col <- z$shape <- paste(z$site, z$order, sep="-")
    
  }else if("site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
    color.var <- "site"
    z$col <- z$shape <- z$site
  }else if(levels(factor(z$order))>1){
    warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
    color.var <- "order"
    z$col <- z$shape <- factor(z$order)
  }else{
    z$col <- z$shape <- rep(1, nrow(z))
  }
}else if(color.var=="order"){     
  z$col <- z$shape <- factor(z$order)
}else if(color.var=="site"){
  if(!"site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
    z$col <- z$shape <- factor(z$order)
  }
  z$col <- z$shape <- z$site
}else if(color.var=="both"){
  if(!"site"%in%names(z)){
    warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
    z$col <- z$shape <- factor(z$order)
  }
  z$col <- z$shape <- paste(z$site, z$order, sep="-")
}

z$lty <- factor(z$method, c("interpolated", "observed", "extrapolated"), c("interpolation", "observation", "extrapolation"))
z$col <- factor(z$col)
data.sub <- z[which(z$method=="observed"),]

palette <- c("#E69F00", "#56B4E9", "#009E73")

g <- ggplot(z, aes_string(x="x", y="y", color = 'site_type')) + 
  geom_point(size=3, data=data.sub)+
  ylab('OTU diversity') + xlab('Number of bats sampled')+
  geom_line(aes_string(linetype="lty"), lwd=1.5)+
  scale_color_manual(values=palette)



g <- g +
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #Get rid of a load of the default crap
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x = element_text(size=10),
      axis.text.y = element_text(size=8))


g <- g + facet_wrap( ~ site, nrow=2, scales = 'free_x', strip.position = 'top')+ #free_x is required so that the x-axes aren't all constrained to showing the same thing
  theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"))#strip stuff sorts the facet labels, spacing adjusts the space between facets
g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr"), alpha=0.2)
g

pdf('plots/iNEXT/hice_inext.pdf', width = 10, height = 7)
g
dev.off()


####Rearrange our stats ####
asymptote_ests <- a$AsyEst
asymptote_ests <- asymptote_ests[asymptote_ests$Diversity=='Species richness',]
asymptote_ests$Site <- as.character(asymptote_ests$Site)
asymptote_ests <- asymptote_ests[order(asymptote_ests$Site),]
asymptote_ests <- cbind(asymptote_ests, a$DataInfo$T)
colnames(asymptote_ests)[8] <- 'number of samples'
asymptote_ests$percent_completeness <- asymptote_ests$Estimator/asymptote_ests$Observed
write.csv(asymptote_ests, 'data/output_data/hice_stats/hice_diversity_ests.csv')
