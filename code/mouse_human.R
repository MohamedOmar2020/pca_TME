



rm(list = ls())

library(superheat)
library(switchBox)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(tidyverse)



###########################
## Load the markers
temp <- list.files("./data/mast2", full.names = T)

markers <- lapply(temp, read.csv)
names(markers) <- temp
names(markers) <- gsub(".+\\/", "", names(markers))
names(markers) <- gsub(".csv", "", names(markers))
names(markers) <- gsub("mast", "", names(markers))

################################
## Convert to human gene symbols

# function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 <- getLDS(attributes = c("hgnc_symbol"),
                    filters = "hgnc_symbol",
                    values = x,
                    mart = human,
                    attributesL = c("mgi_symbol"),
                    martL = mouse,
                    uniqueRows=TRUE)
  
  #humanx <- unique(genesV2[, 2])
  humanx <- genesV2[, 1:2]
  
  print(head(humanx))
  return(humanx)
}

markers_human <- lapply(markers, function(x, y, z, i) {
  i <- 1:length(x)
  y <- convertMouseGeneList(x[, 'X'])
  colnames(y) <- c("human", "X") 
  z <- merge(x, y, by = "X")
  z <- z[order(z[,"avg_log2FC"], decreasing = T), ]
  return(z)
})


save(markers_human, file = "./objs/markers_human.rda")


