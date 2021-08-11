#### Header ####
## Project: bat-diet
## Script purpose: Finding the effect of including different numbers of bats on degree
## Date: 02/07/21
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes: this is a different analysis to that which I did in my thesis, which 
## was to identify the effects of removing specific species. Here I look at the
## effect of the number of bat species present in a network
##################################################

library(tidyverse)
library(bipartite)
library(iNEXT)


source('scripts/r/r_network_gen.r')

big_net_matrix <- r_network_gen(collapse_species = F, desired_species = 'Hice', 
                      include_malua = T, lulu = T, split_by = 'site')

# split the matrix into multiple matrixes in a list, one for each site

site_names <- unique(gsub('\\,.+', '', colnames(big_net_matrix)))

networks <- list()
for(i in site_names){
  networks[[i]] <- big_net_matrix[, grep(i, colnames(big_net_matrix))]
  
  # get rid of empty rows
  networks[[i]] <- networks[[i]][-which(rowSums(networks[[i]]) == 0),]
  
  # get rid of the two weirdly formatted first rows
  networks[[i]] <- networks[[i]][-c(1:2),]
}

# give the networks nicer names
names(networks) <- gsub('DANUM', 'Danum', names(networks))
names(networks) <- gsub('MALIAU', 'Maliau', names(networks))

inext_plot <- ggiNEXT(iNEXT(networks, datatype = 'incidence_raw', 
                            endpoint = 300))

inext_plot

# Analysis ----------------------------------------------------------------

# now that we have a nicely formatted dataset, lets retain only n bats, at 
# random. This is currently failing as the number of possible subsets is huge


# the number of iterations per network and n_bats
n_iterations <- 100000 # (100 thousand)


# the numbers of bats we are going to keep per iteration
nbat_numbers <- seq(10, 100, by = 10)


out_list <- list() # a list to store the output objects in while the loop is going
z <- 1 # an iterator for use in specifying where in out_list to store things
set.seed(12345) # the random seed, saved for reproducibility

for(n in 1:length(networks)){
  # select a network to work on
  full_net <- networks[[n]]
  net_name <- names(networks)[n]
  print(net_name)
  
  for(n_bats in nbat_numbers){
    if(n_bats > ncol(full_net)){ next()}
    print(n_bats)
    
    for(i in 1:n_iterations){
      # randomly choose a subset of bats
      cols_to_use <- sample(1:ncol(full_net), n_bats, replace = F)
      
      subnet <- full_net[,cols_to_use]
      
      # get the number of OTUs consumed
      n_OTUs <- length(which(rowSums(subnet) > 0))
      
      # store the output
      out_list[[z]] <- data.frame(net = net_name, 
                                  number_of_bats = n_bats, 
                                  number_of_OTUs = n_OTUs)
      # increase the iterator
      z <- z + 1
    }
  }
}


# Combine and plot output -------------------------------------------------


out_df <- bind_rows(out_list) %>%
# add a column for if the site is OG or logged
  mutate(logging = ifelse(net %in% c('Danum', 'Maliau'), 'Old growth', 'Logged')) %>%
  # add the word 'bats to the number_of_bats variable, for prettier formatting
  mutate(number_of_bats = paste(number_of_bats, 'bats')) %>%
  mutate(number_of_bats = fct_relevel(number_of_bats, 
                                      levels = paste(seq(10,100, by = 10), 'bats'))) %>%
  # put the network factor into alphabetical order
  mutate(net = fct_relevel(net, sort))


write_csv(out_df, 'results/bat_rarefying.csv')

nbat_plot <- ggplot(out_df, aes(x = net, y = number_of_OTUs, fill = logging)) +
  geom_violin() +
  facet_wrap(. ~ number_of_bats, ncol = 5) +
  scale_fill_viridis_d() + 
  theme_classic() +
  theme(legend.position = 'bottom', 
        # rotate the x-axis labels
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        # increase the text size
        text = element_text(size = 17)) +
  labs(x = 'Site', y = 'Number of MOTUs detected',
       fill = 'Logging treatment')

nbat_plot  
ggsave('plots/nbat_plot.jpeg', nbat_plot)
ggsave('plots/nbat_plot.pdf', nbat_plot)

