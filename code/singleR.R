

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
library(scater)
library(scran)

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
pred_human <- SingleR(test=human_sce, ref=mouse_sce, labels=mouse_sce$cluster, genes = all_markers)

#pred_human <- SingleR(test=human_sce, ref=mouse_sce, labels=mouse_sce$cluster, de.method="wilcox")
table(pred_human$labels)
table(mouse_sce$cluster)

# png('./singleR/scoreHeatmaps_all.png', width = 2000, height = 2000, res = 200)
# plotScoreHeatmap(pred_human)
# dev.off()


plotDeltaDistribution(pred_human, ncol = 8)

pred_markers_all <- metadata(pred_human)$de.genes
#pred_markers_all <- unlist(pred_markers_all)

pred_markers_all <- lapply(pred_markers_all, function(x){
  x = c(unique(unlist(x)))
})

pred_markers_all_df <- do.call('cbind', pred_markers_all)

human_sce$labels <- pred_human$labels

# Beta cell-related markers
#plotHeatmap(human_sce, order_columns_by="labels",
#            features=unique(unlist(all.markers))) 

################################################
## mapping based on erg

# erg positive

# normalization and log-transformation.
mouse_erg_sce <- logNormCounts(mouse_erg_sce)

# the same for human_sce
human_erg_sce <- logNormCounts(human_erg_sce)

by.t_erg <- scran::pairwiseTTests(assay(mouse_erg_sce, 2), mouse_erg_sce$cluster, direction="up")
markers_erg <- scran::getTopMarkers(by.t_erg[[1]], by.t_erg[[2]], n=50)
pred_human_erg <- SingleR(test=human_erg_sce, ref=mouse_erg_sce, labels=mouse_erg_sce$cluster, genes = markers_erg)

#pred_human_erg <- SingleR(test=human_erg_sce, ref=mouse_erg_sce, labels=mouse_erg_sce$cluster, de.method="classic", de.n = 500, genes = 'de')

table(pred_human_erg$labels)
table(mouse_erg_sce$cluster)

#plotScoreHeatmap(pred_human_erg)
#plotDeltaDistribution(pred_human_erg, ncol = 4)

human_erg_markers <- metadata(pred_human_erg)$de.genes
human_erg_sce$labels <- pred_human_erg$labels

#plotHeatmap(human_erg_sce, order_columns_by="labels",
#            features=unique(unlist(human_erg_markers))) 

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

#plotScoreHeatmap(pred_human_nonerg)

#plotDeltaDistribution(pred_human_nonerg, ncol = 4)
human_nonerg_markers <- metadata(pred_human_nonerg)$de.genes

human_nonerg_sce$labels <- pred_human_nonerg$labels

##################
# combine sce
human_sce2 <- cbind(human_erg_sce, human_nonerg_sce)
table(human_sce2$labels)
table(human_sce$labels)
#plotHeatmap(human_sce, order_columns_by="labels",

### combine the predictions
## pred_human_erg doesn't have a column for cluster 5 in the scores matrix
pred_human_erg$scores <- pred_human_erg$scores %>% 
  as.data.frame() %>% 
  dplyr::mutate('5'=NA) %>%
  dplyr::relocate('5', .before = 6)


pred_human2 <- rbind(pred_human_erg, pred_human_nonerg)

human_pred_markers <- metadata(pred_human2)$de.genes

human_pred_markers <- lapply(human_pred_markers, function(x){
  x = c(unique(unlist(x)))
})

human_pred_markers_df2 <- do.call('cbind', human_pred_markers)

################################
## markers human predictions
pred_markers_all <- findMarkers(human_sce2, groups=human_sce2$labels, test.type = 'binom', pval.type = 'any')
markers_c3 <- pred_markers_all[[4]]
markers_c3 <- as.data.frame(markers_c3)
markers_c3['Wnt2', ]

by.t <- scran::pairwiseTTests(assay(human_sce2, 2), human_sce2$labels, direction="up")
pred_markers_all <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=100, pairwise = TRUE)
markers_c3 <- pred_markers_all[[4]]
markers_c3 <- as.data.frame(markers_c3)


png('./singleR/scoreHeatmaps_all.png', width = 2000, height = 2000, res = 200)
plotScoreHeatmap(pred_human2)
dev.off()


pred_markers <- metadata(pred_human2)$de.genes
pred_markers_c3 <- pred_markers$`3` 
pred_markers_c3 <- unique(unlist(pred_markers_c3))
human_sce2$labels <- pred_human2$labels

# Beta cell-related markers
library(scater)
plotHeatmap(human_sce2, order_columns_by="labels",
            features=c('Wnt2', 'Wnt6')) 


##############################
## convert to seurat object

human_seurat <- as.Seurat(human_sce, counts = "counts", data = "logcounts")

## convert to loom object
SaveH5Seurat(human_seurat, filename = "./human/singleR/human_pred.h5seurat", overwrite = T)
Convert("./human/singleR/human_pred.h5seurat", dest = "h5ad")


human_seurat2 <- as.Seurat(human_sce2, counts = "counts", data = "logcounts")
## convert to loom object
SaveH5Seurat(human_seurat2, filename = "./human/singleR/human_pred2.h5seurat", overwrite = T)
Convert("./human/singleR/human_pred2.h5seurat", dest = "h5ad")



###################################################################################
human_seurat2[["percent.mt"]] <- PercentageFeatureSet(human_seurat2, pattern = "^Mt-")

# Visualize QC metrics as a violin plot
VlnPlot(human_seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(human_seurat2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_seurat2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#human_seurat2 <- subset(human_seurat2, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 25)

human_seurat2 <- NormalizeData(human_seurat2, normalization.method = "LogNormalize", scale.factor = 10000)

human_seurat2 <- FindVariableFeatures(human_seurat2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(human_seurat2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(human_seurat2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(human_seurat2)
human_seurat2 <- ScaleData(human_seurat2, features = all.genes)


human_seurat2 <- RunPCA(human_seurat2, features = VariableFeatures(object = human_seurat2))


human_seurat2 <- FindNeighbors(human_seurat2, dims = 1:10)
human_seurat2 <- FindClusters(human_seurat2, resolution = 0.5)

# 'umap-learn')
human_seurat2 <- RunUMAP(human_seurat2, dims = 1:10)

Idents(human_seurat2) <- human_seurat2$labels
# individual clusters
DimPlot(human_seurat2, reduction = "umap")


VlnPlot(human_seurat2, features = c("Ar"), slot = "counts", log = TRUE)

















