

rm(list = ls())
library(tidyverse)
library(data.table)

####################
####################
## Cell type markers
load("./data/markers_human.rda")

## Extract the top markers of each cluster
CellMarkers <- lapply(markers_human, function(x, i){
  i <- 1:length(x)
  all <- x[i][, c("human", 'avg_log2FC')]
  # Up-regulated
  up <- all %>%
    filter(avg_log2FC > 0) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::select('human') %>%
    dplyr::rename(up = human)
  # Down-regulated
  down <- all %>%
    filter(avg_log2FC < 0) %>%
    dplyr::arrange(avg_log2FC) %>%
    dplyr::select('human') %>%
    dplyr::rename(down = human)
  return(c(up, down))
}) 

CellMarkers_c5_c6_c7 <- CellMarkers[c('c5', 'c6', 'c7')]
  
CellMarkers_c5_c6_c7 <- lapply(rapply(CellMarkers_c5_c6_c7, enquote, how="unlist"), eval)

CellMarkers_c5_c6_c7_top100 <- lapply(CellMarkers_c5_c6_c7, function(x){
  x <- x[1:100]
})

#################################
## make pairs by pairing up- with down-regulated genes

c5_pairs <- expand_grid(CellMarkers_c5_c6_c7$c5.up, CellMarkers_c5_c6_c7$c5.down) 
colnames(c5_pairs) <- c('up', 'down')

c6_pairs <- expand_grid(CellMarkers_c5_c6_c7$c6.up, CellMarkers_c5_c6_c7$c6.down) 
colnames(c6_pairs) <- c('up', 'down')

c7_pairs <- expand_grid(CellMarkers_c5_c6_c7$c7.up, CellMarkers_c5_c6_c7$c7.down) 
colnames(c7_pairs) <- c('up', 'down')

#########
## make pairs by pairing up- with down-regulated genes: using only top 100 markers

c5_pairs_top <- expand_grid(CellMarkers_c5_c6_c7_top100$c5.up, CellMarkers_c5_c6_c7_top100$c5.down) 
colnames(c5_pairs_top) <- c('up', 'down')

c6_pairs_top <- expand_grid(CellMarkers_c5_c6_c7_top100$c6.up, CellMarkers_c5_c6_c7_top100$c6.down) 
colnames(c6_pairs_top) <- c('up', 'down')

c7_pairs_top <- expand_grid(CellMarkers_c5_c6_c7_top100$c7.up, CellMarkers_c5_c6_c7_top100$c7.down) 
colnames(c7_pairs_top) <- c('up', 'down')

#########
# bind together

PRN_pairs <- rbind(c5_pairs, c6_pairs, c7_pairs)

PRN_pairs_top <- rbind(c5_pairs_top, c6_pairs_top, c7_pairs_top)

##########
# save
save(PRN_pairs, file = './data/PRN_pairs.rda')


save(PRN_pairs_top, file = './data/PRN_pairs_top.rda')

