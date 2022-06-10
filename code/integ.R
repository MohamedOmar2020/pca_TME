rm(list = ls())

Sys.setenv(RETICULATE_PYTHON = "/Users/mohamedomar/opt/anaconda3/envs/scutils/bin/python3.10")
library(reticulate)
library(data.table)
library(readxl)
library(ggplot2)
library(patchwork)
library(projectR)
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
#mouse_obj <- readRDS('./outs/h5ads/prostate_mouse_loda.rds')

# conda env and scanpy
use_condaenv(condaenv = "scutils", required = TRUE)
sc <- import("scanpy")
scp <- import('scipy')

## convert the human h5ad to seurat
mouse_adata <- sc$read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad')
counts <- scp$sparse$csr_matrix$todense(mouse_adata$layers['counts'])
counts <- t(counts)
#counts <- t(mouse_adata$layers['counts'])
colnames(counts) <-  mouse_adata$obs_names$to_list()
rownames(counts) <-  mouse_adata$var_names$to_list()
rownames(counts) <- str_to_title(rownames(counts))
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

# counts <- scp$sparse$csr_matrix$todense(mouse_adata$raw$X)
# counts <- t(counts)
# colnames(counts) <-  mouse_adata$raw$obs_names$to_list()
# rownames(counts) <-  mouse_adata$raw$var_names$to_list()
# counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

mouse_obj <- CreateSeuratObject(counts)
mouse_obj <- AddMetaData(mouse_obj,  mouse_adata$obs)

# Remove the unlabelled libraries here.
mouse_obj <- mouse_obj[,!is.na(mouse_obj$cluster)]

############
## convert the human h5ad to seurat
human_adata <- sc$read_h5ad('./human/erg_fibroblasts_scvi_v6_regulons2.h5ad')

# counts <- scp$sparse$csr_matrix$todense(human_adata$raw$X)
# counts <- t(counts)
# colnames(counts) <-  human_adata$raw$obs_names$to_list()
# rownames(counts) <-  human_adata$raw$var_names$to_list()
# rownames(counts) <- str_to_title(rownames(counts))
# counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

counts <- scp$sparse$csr_matrix$todense(human_adata$layers['counts'])
counts <- t(counts)
colnames(counts) <-  human_adata$obs_names$to_list()
rownames(counts) <-  human_adata$var_names$to_list()
rownames(counts) <- str_to_title(rownames(counts))
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

human_obj <- CreateSeuratObject(counts)
human_obj <- AddMetaData(human_obj,  human_adata$obs)

# find common genes
commonGns <- intersect(rownames(mouse_obj), rownames(human_obj))

# subset
mouse_obj <- mouse_obj[commonGns, ]
human_obj <- human_obj[commonGns, ]

###############################
# seurat standard workflow
seurat_list <- list(mouse = mouse_obj, human = human_obj)

# normalization and identification of highly variable features
#for (i in 1:length(seurat_list)) {
  #seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
#  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", 
#                                             nfeatures = 2000, verbose = FALSE)
#}

for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- SCTransform(seurat_list[[i]], verbose = FALSE)
}

## integration:
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features, 
                                  verbose = FALSE)


query <- seurat_list[["human"]]


anchors <- FindTransferAnchors(reference = seurat_list[["mouse"]], query = seurat_list[["human"]], 
                               dims = 1:30, features = commonGns, normalization.method = "SCT")

predictions <- TransferData(anchorset = anchors, refdata = mouse_obj$cluster, 
                            dims = 1:30)

query <- AddMetaData(query, metadata = predictions)
table(query$predicted.id)

VlnPlot(query, c("Wnt2"), group.by = "predicted.id")



#query <- ScaleData(query, verbose = FALSE)
query <- RunPCA(query, npcs = 30, verbose = FALSE)
query <- RunUMAP(query, reduction = "pca", dims = 1:30)
Idents(query) <- query$predicted.id
DimPlot(query, reduction = "umap", group.by = "predicted.id")

# find all markers of cluster 2
cluster0.markers <- FindMarkers(query, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 20)


cluster1.markers <- FindMarkers(query, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 20)

cluster2.markers <- FindMarkers(query, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)

cluster3.markers <- FindMarkers(query, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 20)


FeaturePlot(query, c("Wnt2", "Wnt6"))
















