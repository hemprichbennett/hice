#### Header ####
## Project: cervinus ecology
## Script purpose: Combining the outputs from otu_inext_lulu.R and otu_inext_no_lulu.R
## Date:
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(ggplot2)


lulu <- read.csv('results/lulu_otu_inext.csv')
non_lulu <- read.csv('results/non_lulu_otu_inext.csv')


both <- rbind(lulu, non_lulu)

colnames(both) <- gsub('Site', 'Network', colnames(both))

both$Treatment <- gsub('non_', 'Non-', both$Treatment)

both$variable <- gsub('Estimated MOTU richness', 'Estimated MOTU\nrichness', both$variable)
both$variable <- gsub('Observed MOTU richness', 'Observed MOTU\nrichness', both$variable)
both$variable <- gsub('Estimated % completeness', 'Estimated %\ncompleteness', both$variable)
both$variable <- gsub("Number of samples required\nfor complete community", "Number of samples\nrequired for\ncomplete community", both$variable)

both$Site <- gsub(',.+','', both$Network)
both$Year <- gsub('.+ ', '', both$Network)

palette <- c("#8960b3",
             "#56ae6c",
             "#ba495b",
             "#b0923b")

facetted <- ggplot(both, aes(x = clustering, y = value, colour = Site, shape = Year)) + 
  geom_point(size = 2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values =palette, name = 'Site')+
  labs(x= 'Clustering %', y = NULL)+ 
  facet_grid(variable ~ Treatment, scales = 'free_y', switch = 'y')+
  theme(strip.placement = "outside", legend.position="bottom",
        text = element_text(size = 18),
        # split the legend over two lines, reduce the margin sizes of it
        legend.box="vertical", legend.margin=margin()) +
        # increase the point size in legend
        guides(shape = guide_legend(override.aes = list(size = 5)),
               colour = guide_legend(override.aes = list(size = 5)))
  

facetted


pdf('plots/inext_otus/both_treatments.pdf', height = 10)
facetted
dev.off()

jpeg('plots/inext_otus/both_treatments.jpeg', width = 7, height = 10, units = 'in', res = 300)
facetted
dev.off()
