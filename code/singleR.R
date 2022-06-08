

rm(list = ls())

Sys.setenv(RETICULATE_PYTHON = "/Users/mohamedomar/opt/anaconda3/envs/scutils/bin/python3.10")
library(reticulate)

library(data.table)
library(readxl)
library(ggplot2)
library(patchwork)
library(SingleR)
library(Seurat)
library(SeuratDisk)
library(tools)
library(stringr)

######################################
# convert the human h5ad to h5seurat
#Convert("./human/erg_fibroblasts_scvi_v6_regulons2.h5ad", dest = "h5seurat", assay = 'RNA', overwrite = TRUE)
#human_obj <- LoadH5Seurat('./human/erg_fibroblasts_scvi_v6_regulons2.h5seurat', assay = 'RNA')
#human_obj

# load the mouse data
#mouse_obj <- readRDS('./outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_4.rds')

# conda env and scanpy
use_condaenv(condaenv = "scutils", required = TRUE)
sc <- import("scanpy")

## convert the human h5ad to seurat
mouse_adata <- sc$read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons2_5.h5ad')

counts <- t(mouse_adata$X)
colnames(counts) <-  mouse_adata$obs_names$to_list()
rownames(counts) <-  mouse_adata$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

mouse_obj <- CreateSeuratObject(counts)
mouse_obj <- AddMetaData(mouse_obj,  mouse_adata$obs)

# Remove the unlabelled libraries here.
mouse_obj <- mouse_obj[,!is.na(mouse_obj$cluster)]

############
## convert the human h5ad to seurat
human_adata <- sc$read_h5ad('./human/erg_fibroblasts_scvi_v6_regulons2.h5ad')

counts <- t(human_adata$X)
colnames(counts) <-  human_adata$obs_names$to_list()
rownames(counts) <-  human_adata$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

human_obj <- CreateSeuratObject(counts)
human_obj <- AddMetaData(human_obj,  human_adata$obs)

###############################
# convert both to singlecellexperiment object
mouse_sce <- as.SingleCellExperiment(mouse_obj)
human_sce <- as.SingleCellExperiment(human_obj)

rownames(human_sce) <- str_to_title(rownames(human_sce))


# divide the mouse data into erg and non-erg
table(mouse_sce$key)
mouse_erg_sce <- mouse_sce[, mouse_sce$key == 'terg']
mouse_nonerg_sce <- mouse_sce[, mouse_sce$key != 'terg']

# divide the human data into erg and non-erg
table(human_sce$erg)
human_erg_sce <- human_sce[, human_sce$erg == 'positive']
human_nonerg_sce <- human_sce[, human_sce$erg == 'negative']

###############################
# mapping on all data

# SingleR() expects reference datasets to be normalized and log-transformed.
library(scuttle)
mouse_sce <- logNormCounts(mouse_sce)

# the same for human_sce
human_sce <- logNormCounts(human_sce)

by.t <- scran::pairwiseTTests(assay(mouse_sce, 2), mouse_sce$cluster, direction="up")
all_markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=50)
pred_human <- SingleR(test=human_sce, ref=mouse_sce, labels=mouse_sce$cluster, genes = markers)

#pred_human <- SingleR(test=human_sce, ref=mouse_sce, labels=mouse_sce$cluster, de.method="wilcox")
table(pred_human$labels)
table(mouse_sce$cluster)

plotScoreHeatmap(pred_human)
plotDeltaDistribution(pred_human, ncol = 8)

all.markers <- metadata(pred_human)$de.genes
human_sce$labels <- pred_human$labels

# Beta cell-related markers
library(scater)
plotHeatmap(human_sce, order_columns_by="labels",
            features=unique(unlist(all.markers))) 

################################################
## mapping based on erg

# erg positive

# normalization and log-transformation.
mouse_erg_sce <- logNormCounts(mouse_erg_sce)

# the same for human_sce
human_erg_sce <- logNormCounts(human_erg_sce)

by.t <- scran::pairwiseTTests(assay(mouse_erg_sce, 2), mouse_erg_sce$cluster, direction="up")
markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=50)
pred_human_erg <- SingleR(test=human_erg_sce, ref=mouse_erg_sce, labels=mouse_erg_sce$cluster, genes = markers)

#pred_human_erg <- SingleR(test=human_erg_sce, ref=mouse_erg_sce, labels=mouse_erg_sce$cluster, de.method="classic", de.n = 500, genes = 'de')

table(pred_human_erg$labels)
table(mouse_erg_sce$cluster)

plotScoreHeatmap(pred_human_erg)
plotDeltaDistribution(pred_human_erg, ncol = 4)

human_erg_markers <- metadata(pred_human_erg)$de.genes
human_erg_sce$labels <- pred_human_erg$labels

plotHeatmap(human_erg_sce, order_columns_by="labels",
            features=unique(unlist(human_erg_markers))) 

################
# erg negative

# normalization and log-transformation.
mouse_nonerg_sce <- logNormCounts(mouse_nonerg_sce)

# the same for human_sce
human_nonerg_sce <- logNormCounts(human_nonerg_sce)

by.t <- scran::pairwiseTTests(assay(mouse_nonerg_sce, 2), mouse_nonerg_sce$cluster, direction="up")
markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=10)
pred_human_nonerg <- SingleR(test=human_nonerg_sce, ref=mouse_nonerg_sce, labels=mouse_nonerg_sce$cluster, genes = markers)

#pred_human_nonerg <- SingleR(test=human_nonerg_sce, ref=mouse_nonerg_sce, labels=mouse_nonerg_sce$cluster, de.method="wilcox")
table(pred_human_nonerg$labels)
table(mouse_nonerg_sce$cluster)

plotScoreHeatmap(pred_human_nonerg)

plotDeltaDistribution(pred_human_nonerg, ncol = 4)

human_nonerg_sce$labels <- pred_human_nonerg$labels

##################
# combine sce
human_sce2 <- cbind(human_erg_sce, human_nonerg_sce)
table(human_sce2$labels)
table(human_sce$labels)
plotHeatmap(human_sce, order_columns_by="labels",
            features=unique(unlist(all.markers))) 

# combine predictions
pred_human_erg_df <- as.data.frame(pred_human_erg)
pred_human_nonerg_df <- as.data.frame(pred_human_nonerg)

pred_human2 <- dplyr::bind_rows(pred_human_erg_df, pred_human_nonerg_df)
pred_human2 <- DataFrame(pred_human2, row.names = rownames(pred_human2))

plotScoreHeatmap(pred_human2)




