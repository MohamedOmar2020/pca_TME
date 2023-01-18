
rm(list = ls())

# load libraries
library(data.table)
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(reticulate)
library(corrplot)

pkgs <- c(
  "Seurat", "SeuratWrappers", "ggplot2", "batchelor", "circlize",
  "dplyr", "optparse", "reshape2", "data.table", "magrittr",
  "patchwork", "scales", "GSVA", "RColorBrewer", "ggridges", "ggridges",
  "clusterProfiler", "survminer", "survminer", "monocle", "tidyverse",
  "psych", "ggrepel", "pheatmap", "escape", "multcomp", "agricolae"
)

color_for_use <- c(
  "#8F2A47", "#B83A4D", "#D25C52", "#E27E56",
  "#ECA86B", "#F4CB85", "#F8E8A2", "#FAF8C7", "#EBF0AF",
  "#CEE2A2", "#ABD3A6", "#82C3A5", "#609EB0", "#4C78B1",
  "#5C519B"
)

lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) 



#######################
# set the active conda environment
use_condaenv(condaenv = "scutils", required = TRUE)
use_condaenv(condaenv = "pathml2", required = TRUE)


# load the anndata module
ad <- import("anndata", convert = FALSE)
reticulate::py_version()

#############################
# load the combined CAFs and mouse h5ad object
adata_all <- ad$read_h5ad("data/adata_all.h5ad")
# access normalized data matrix
adata_all_dataInput <- t(py_to_r(adata_all$X))
range(adata_all_dataInput)
rownames(adata_all_dataInput) <- rownames(py_to_r(adata_all$var))
colnames(adata_all_dataInput) <- rownames(py_to_r(adata_all$obs))

## access meta data
adata_all_metaData <- py_to_r(adata_all$obs)
adata_all_meta <- adata_all_metaData
# some cleaning
#adata_all_meta$cluster <- paste0("c", adata_all_meta$cluster)
#adata_all_meta <- adata_all_meta[!(adata_all_meta$cluster == 'cNA'), ]
adata_all_dataInput <- adata_all_dataInput[, rownames(adata_all_meta)]

adata_all <- CreateSeuratObject(adata_all_dataInput)
adata_all <- AddMetaData(adata_all, adata_all_meta)

table(adata_all$cluster)

# change the active idents
Idents(adata_all) <- adata_all$cluster

## Similarity analysis
col1 <- c(
  "#8F2A47", "#B83A4D", "#D25C52", "#E27E56",
  "#ECA86B", "#F4CB85", "#F8E8A2", "#FAF8C7", "#EBF0AF",
  "#CEE2A2", "#ABD3A6", "#82C3A5", "#609EB0", "#4C78B1",
  "#5C519B"
)
simi <- function(x) {
  Object_Select <- subset(adata_all, ident = x)
  Object_Select <- FindVariableFeatures(
    Object_Select,
    selection.method = "vst", nfeatures = 1000
  )
  Expr_Matrix <- as.data.frame(t(
    as.data.frame(
      Object_Select@assays$RNA@data[VariableFeatures(Object_Select), ]
    )
  ))
  Expr_Matrix$cluster <- as.character(Object_Select@active.ident)
  Expr_Matrix_Mean <- aggregate(
    Expr_Matrix[, 1:(ncol(Expr_Matrix) - 1)],
    by = list(Expr_Matrix$cluster), FUN = mean
  )
  rownames(Expr_Matrix_Mean) <- Expr_Matrix_Mean$Group.1
  Expr_Matrix_Mean <- Expr_Matrix_Mean[, 2:ncol(Expr_Matrix_Mean)]
  Expr_Matrix_Mean <- as.data.frame(t(Expr_Matrix_Mean))
  Corr_Result <- cor(Expr_Matrix_Mean)
  pdf("MouseStromal_humanCAFs_Cor.pdf", width = 7, height = 7)
  corrplot(Corr_Result,
           method = "color",
           col = rev(brewer.pal(11, "RdYlBu")),
           type = "full",
           tl.col = "black",
           order = "hclust", is.corr = FALSE, addCoef.col = "grey"
  )
  dev.off()
}
clusters <- unique(adata_all$cluster)
simi(clusters)


