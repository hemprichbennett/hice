  if(interactive()==TRUE){
    library('here')
    library(ggplot2)
    library(tidyverse)
    library(ggridges)
    library(gridExtra)
    library(forcats)
    library(reshape2)
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
  
  field_data$Faeces_no2 <- NULL
  
  field_data_2$Faeces_no1 <- field_data_2$Faeces_no2
  field_data_2$Faeces_no2 <- NULL
  
  field_data <- rbind(field_data, field_data_2)
  
  all_interactions <- r_network_gen(collapse_species = F, desired_species = 'Hice', include_malua = T, lulu = T)
  
  
  prey_data <- read.csv('data/taxonomy/family.csv')
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
  
  t_taxa_mat <- data.frame(t(taxa_mat))
  
  t_taxa_mat$Sample <- rownames(t_taxa_mat)
  
  t_taxa_mat <- merge(x=t_taxa_mat, y =field_data, by.x = 'Sample', by.y = 'Faeces_no1')
  


tax_df <- t_taxa_mat[,c(seq(2,146), 162)]
tax_df <- melt(tax_df, id.vars = c('Site'))
colnames(tax_df)[c(2,3)] <- c('Family', 'Present/absent')
tax_df$`Present/absent` <- as.integer(as.character(tax_df$`Present/absent`))
tax_df$`Present/absent` <- ifelse(tax_df$`Present/absent`== 0, 0, 1)

prop_present <- sapply(unique(tax_df[,c('Site', 'Family')]),  function(x) as.character(x))
prop <- c()
nbats <- c()
#for(i in 1: 1){
for(i in 1: nrow(prop_present)){
  tem <- tax_df[which(tax_df$Site== as.character(prop_present[i,1])
                      & tax_df$Family==as.character(prop_present[i,2])),]
  prop <- c(prop, sum(tem$`Present/absent`)/nrow(tem)) #The number of bats that consumed the Family, divided by total bats
  nbats <- c(nbats, nrow(tem))
}


prop_present <- cbind(prop_present, prop, nbats)
prop_present <- as.data.frame(prop_present)
prop_present$prop <- as.numeric(as.character(prop_present$prop))
prop_present$nbats <- as.integer(as.character(prop_present$nbats))

prop_present$Site <- gsub('DANUM', 'Danum', prop_present$Site)
prop_present$Site <- gsub('DVCA', 'Danum', prop_present$Site)
prop_present$Site <- gsub('MALIAU', 'Maliau', prop_present$Site)
prop_present$Site <- gsub('MALUA', 'SBE', prop_present$Site)


prop_present$order <- rep(NA, nrow(prop_present)) #Make a column of the order which each family belongs to
prop_present$order[prop_present$Family %in% c('Sparassidae', 'Salticidae', 'Araneidae', 'Theridiidae', 'Clubionidae', 'Linyphiidae', 'Pholcidae', 
                                              'Uloboridae')] <- 'Araneae (Spiders)'
prop_present$order[prop_present$Family %in% c('Ectobiidae', 'Blaberidae', 'Termitidae', 'Blattidae', 'Rhinotermitidae', 'Cryptocercidae')] <- 'Blattodea (Termites\nand cockroaches)'
prop_present$order[prop_present$Family %in% c('Chrysomelidae', 'Cerambycidae', 'Curculionidae', 'Elateridae', 'Carabidae', 'Ptilodactylidae', 'Mordellidae', 'Eucnemidae',
                                              'Latridiidae', 'Dermestidae', 'Mycetophagidae', 'Dytiscidae', 'Nitidulidae', 'Cleridae', 'Leiodidae', 'Byrrhidae',
                                              'Throscidae', 'Scarabaeidae', 'Psephenidae', 'Tenebrionidae', 'Zopheridae', 'Cantharidae')] <- 'Coleoptera (Beetles)'
prop_present$order[prop_present$Family %in% c('Culicidae', 'Chironomidae', 'Tachinidae', 'Sciaridae', 'Psychodidae', 'Phoridae', 'Ceratopogonidae', 
                                              'Muscidae', 'Mycetophilidae', 'Chloropidae', 'Calliphoridae', 'Tabanidae', 'Stratiomyidae', 'Cecidomyiidae', 'Milichiidae',
                                              'Syrphidae', 'Asilidae', 'Ephydridae', 'Dolichopodidae', 'Pediciidae', 'Drosophilidae', 'Keroplatidae', 'Tipulidae',
                                              'Platypezidae')] <- 'Diptera (Flies)'
prop_present$order[prop_present$Family %in% c('Entomobryidae', 'Isotomidae')] <- 'Entomobryomorpha (Springtails)'
prop_present$order[prop_present$Family %in% c('Baetidae', 'Heptageniidae', 'Ephemerellidae')] <- 'Ephemeroptera (Mayflies)'
prop_present$order[prop_present$Family %in% c('Euphausiidae')] <- 'Euphausiacea (krill)'
prop_present$order[prop_present$Family %in% c('Armadillidiidae')] <- 'Isopoda (Woodlice)'
prop_present$order[prop_present$Family %in% c('Cicadidae', 'Cicadellidae', 'Pentatomidae', 'Derbidae', 'Aphididae', 'Miridae', 'Fulgoridae', 'Flatidae', 'Cydnidae',
                                              'Pemphigidae', 'Kinnaridae', 'Rhyparochromidae', 'Dictyopharidae', 'Delphacidae', 'Nogodinidae', 'Cixiidae', 'Hormaphididae',
                                              'Lygaeidae')] <- 'Hemiptera (True bugs)'
prop_present$order[prop_present$Family %in% c('Ichneumonidae', 'Mutillidae', 'Perilampidae', 'Formicidae', 'Braconidae', 'Dryinidae', 'Agaonidae', 'Vespidae',
                                              'Tenthredinidae')] <- 'Hymenoptera (Ants, wasps, etc)'
prop_present$order[prop_present$Family %in% c('Geometridae', 'Lymantriidae', 'Tortricidae', 'Noctuidae', 'Lycaenidae', 'Sphingidae', 'Crambidae', 'Tineidae', 'Erebidae',
                                              'Nolidae', 'Adelidae', 'Pterophoridae', 'Limacodidae', 'Pyralidae', 'Blastobasidae', 'Hesperiidae', 'Lecithoceridae')] <- 'Lepidoptera (Butterflies and moths)'
prop_present$order[prop_present$Family %in% c('Mantidae', 'Tarachodidae')] <- 'Mantodea (Mantises)'
prop_present$order[prop_present$Family %in% c('Laelapidae')] <- 'Mesostigmata (Mites)'
prop_present$order[prop_present$Family %in% c('Chrysopidae', 'Hemerobiidae', 'Coniopterygidae')] <- 'Neuroptera (Net-winged insects)'
prop_present$order[prop_present$Family %in% c('Brachychthoniidae')] <- 'Oribatida (Moss mites)'
prop_present$order[prop_present$Family %in% c('Acrididae', 'Gryllidae', 'Tettigoniidae', 'Acrididae', 'Chorotypidae')] <- 'Orthoptera (Grasshoppers)'
prop_present$order[prop_present$Family %in% c('Diapheromeridae', 'Phasmatidae')] <- 'Phasmatodea (Stick insects)'
prop_present$order[prop_present$Family %in% c('Paradoxosomatidae')] <- 'Polydesmida (Millipedes)'
prop_present$order[prop_present$Family %in% c('Caeciliusidae', 'Archipsocidae', 'Stenopsocidae', 'Lepidopsocidae', 'Psocidae', 'Epipsocidae')] <- 'Psocoptera (Barkflies)'
prop_present$order[prop_present$Family %in% c('Eremaeidae', 'Brachychthoniida', 'Scheloribatidae', 'Ceratozetidae', 'Oppiidae')] <- 'Sarcoptiformes (Mites)'
prop_present$order[prop_present$Family %in% c('Philopotamidae', 'Psychomyiidae', 'Leptoceridae')] <- 'Trichoptera (Caddisflies)'
prop_present$order[prop_present$Family %in% c('Phlaeothripidae')] <- 'Thysanoptera (Thrips)'
prop_present$order[prop_present$Family %in% c('Tarsonemidae', 'Eupodidae', 'Bdellidae', 'Cheyletidae', 'Ereynetidae', 'Teneriffiidae', 'Cunaxidae')] <- 'Trombidiformes (Mites)'


unique(prop_present$Family[which(is.na(prop_present$order))])


tiles <- ggplot(data = prop_present[which(prop_present$prop!=0),], aes(y = fct_rev(Family), x =Site)) + geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 0.175,  limits = c(0,0.35)) +
  labs(fill='Proportion of\nsamples\ncontaining',
       x ="Site", y = 'Prey')+
  theme(panel.grid.major = element_blank(),
        panel.background= element_rect(fill = "darkgray", colour = "darkgray",size = 0.5, linetype = NULL),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(order ~ ., scales = 'free', space="free", switch = 'both')+
  theme(strip.text.y = element_text(angle = 180, colour = 'black', size = 8), strip.placement = 'outside',
        strip.background =element_rect(fill="white", size = 3))
tiles

pdf('plots/proportion_of_individuals_containing_family.pdf', height = 13)
tiles
dev.off()

widetiles <- ggplot(data = prop_present[which(prop_present$prop!=0),], aes(x = Family, y = fct_rev(Site))) + 
  geom_tile(aes(fill=prop), colour = 'white')+
  scale_fill_gradient2(low = "white", mid = "blue",
                       high = "black", midpoint = 0.175, name ='Proportion of\nsamples\ncontaining', limits = c(0,0.35)) +
  labs(y ="Bat species", x = 'Prey')+
  theme_linedraw(base_size = 14)+
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text.x = element_text(colour = 'white', size = 8, angle = 90), strip.placement = "outside",
        panel.background = element_rect(fill = "darkgray",
                                        colour = "darkgray",
                                        size = 0.5, linetype = NULL),
        strip.background =element_rect(fill="black", size = 3))+
  facet_grid(. ~order, scales = 'free', space="free", switch = 'both')#+

widetiles

pdf('plots/proportion_of_individuals_containing_family_order.pdf', width = 17)
widetiles
dev.off()
