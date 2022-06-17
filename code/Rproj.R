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
library(projectR)

############################################################
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

# commonGns
commonGns <- intersect(rownames(mouse_mat), rownames(human_mat))

######################################
## Obtaining PCs to project
mouse_mat <- as.matrix(mouse_obj@assays$RNA@counts)
mouse_metadata <- mouse_obj@meta.data

# do PCA on mouse expression data
pc_mouse<-prcomp(t(mouse_mat))
pcVAR <- round(((pc_mouse$sdev)^2/sum(pc_mouse$sdev^2))*100,2)
dPCA <- data.frame(cbind(pc_mouse$x, mouse_metadata))

#plot pca
library(ggplot2)
#setCOL <- scale_colour_manual(values = c("blue","black","red"), name="Condition:")
#setFILL <- scale_fill_manual(values = c("blue","black","red"),guide = FALSE)
#setPCH <- scale_shape_manual(values=c(23,22,25,25,21,24),name="Cell Line:")
pPCA <- ggplot(dPCA, aes(x=PC1, y=PC2, colour=cluster,
                         fill=cluster)) +
  geom_point(alpha=.6)+
  #setCOL + setPCH + setFILL +
  #scale_size_area(breaks = c(2,4,6), name="Day") +
  theme(legend.position=c(0,0), legend.justification=c(0,0),
        legend.direction = "horizontal",
        panel.background = element_rect(fill = "white",colour=NA),
        legend.background = element_rect(fill = "transparent",colour=NA),
        plot.title = element_text(vjust = 0,hjust=0,face="bold")) +
  labs(title = "PCA of mouse clusters",
       x=paste("PC1 (", pcVAR[1],"% of varience)",sep=""),
       y=paste("PC2 (", pcVAR[2],"% of varience)",sep=""))

pPCA

#############################################
## Projection
# data to project into PCs from human expression data
human_mat <- as.matrix(human_obj@assays$RNA@counts)
mouse_metadata <- human_obj@meta.data

library(projectR)
PCA2human <- projectR(data = human_mat,loadings=pc_mouse,
                      full = TRUE, dataNames = rownames(human_mat))

pd.ESepiGen4c1l<-data.frame(Condition=sapply(colnames(p.ESepiGen4c1l$mRNA.Seq),
                                             function(x) unlist(strsplit(x,'_'))[1]),stringsAsFactors=FALSE)
pd.ESepiGen4c1l$color<-c(rep("red",2),rep("green",3),rep("blue",2),rep("black",2))

names(pd.ESepiGen4c1l$color)<-pd.ESepiGen4c1l$Cond
dPCA2ESepi<- data.frame(cbind(t(PCA2ESepi[[1]]),pd.ESepiGen4c1l))

#plot pca
setEpiCOL <- scale_colour_manual(values = c("red","green","blue","black"),
                                 guide = guide_legend(title="Lineage"))
pPC2ESepiGen4c1l <- ggplot(dPCA2ESepi, aes(x=PC1, y=PC2, colour=Condition)) +
  geom_point(size=5) + setEpiCOL +
  theme(legend.position=c(0,0), legend.justification=c(0,0),
        panel.background = element_rect(fill = "white"),
        legend.direction = "horizontal",
        plot.title = element_text(vjust = 0,hjust=0,face="bold")) +
  labs(title = "Encode RNAseq in target PC1 & PC2",
       x=paste("Projected PC1 (",round(PCA2ESepi[[2]][1],2),"% of varience)",sep=""),
       y=paste("Projected PC2 (",round(PCA2ESepi[[2]][2],2),"% of varience)",sep=""))





