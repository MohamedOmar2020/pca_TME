
##############Libraries####
library(Matrix)
library(reshape)
library(readxl)
library(ggplot2)
library(devtools)
library(circlize)
library(ggpubr)
library(plyr)
library(dplyr)
library(scales)
library(viridis)
library(grid)
library(Seurat)
library(mgcv)
library(MASS)
library(colorRamps)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpmisc)
library(cowplot)
library(M3C)
library(ConsensusClusterPlus)
library(clusterProfiler)
library(rstatix)
library(fgsea)
library(org.Mm.eg.db)
library(GO.db)
library(limma)
library(reticulate)
library(Seurat)
library(SeuratDisk)

############################
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

##############################
# set the active conda environment
use_condaenv(condaenv = "scutils", required = TRUE)

# load the anndata module
ad <- import("anndata", convert = FALSE)

# load the mouse h5ad object
adata_mouse <- ad$read_h5ad("outs/h5ads/adata_mouse.h5ad")
# access normalized data matrix
mouse_dataInput <- t(py_to_r(adata_mouse$X))
range(mouse_dataInput)
rownames(mouse_dataInput) <- rownames(py_to_r(adata_mouse$var))
colnames(mouse_dataInput) <- rownames(py_to_r(adata_mouse$obs))

## access meta data
mouse_metaData <- py_to_r(adata_mouse$obs)
mouse_meta <- mouse_metaData
# some cleaning
mouse_meta$cluster <- paste0("c", mouse_meta$cluster)
mouse_dataInput <- mouse_dataInput[, rownames(mouse_meta)]

#######################
## create seurat object
mouse_object <- CreateSeuratObject(counts = mouse_dataInput, meta.data = mouse_meta)
Idents(mouse_object) <- mouse_object$cluster

# add the umap
umapCoord <- py_to_r(adata_mouse$obsm['X_umap'])
rownames(umapCoord) <- rownames(mouse_meta)
umapCoord <- as(umapCoord, "matrix")
mouse_object@reductions$umap <- CreateDimReducObject(embeddings = umapCoord, key = "UMAP_", global = T, assay = "RNA")

# get the colors
clusterColors <- py_to_r(adata_mouse$uns['cluster_colors'])

# re-order the levels
mouse_object@active.ident <- factor(mouse_object@active.ident, levels = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'))
table(mouse_object@active.ident)



DimPlot(mouse_object, reduction = "umap", cols = clusterColors, pt.size = 2)

####################################################
stress_genes= c("G0s2", "Jun", "Junb", "Jund", "Fos", "Dusp1", "Cdkn1a", "Fosb", "Btg2", "Klf6", "Klf4")
cc_genes =  unique(c(paste0(toupper(substr(tolower(as.character(unlist(cc.genes.updated.2019))), 1, 1)), substr(tolower(as.character(unlist(cc.genes.updated.2019))), 2, nchar(tolower(as.character(unlist(cc.genes.updated.2019)))))), "1700020L24Rik" ,"5730416F02Rik" ,"Agpat4","Asf1b", "Aspm","Ccdc18","Ccr9","Clspn","Cyp4f17","Dek","Dnmt3b" ,"Dtl","Fancm","Fignl1","Gm14214","Gm14730","Gm15428","Gm15448","Gm21992","Gm23248","Gm26465","Gm5145" ,"Mcm4","Mcm8","Mki67","Oip5","Pcna","Pcna-ps2","Pigv","Rnd1","Snrpa","Ube2c"))
IFN_genes = unique(c(grep("Irf", rownames(mouse_object@assays$RNA@counts), v = T), grep("Ifi", rownames(rownames(mouse_object@assays$RNA@counts)), v = T), "Cd2","Cd3d", "Cmpk2","Cxcl10", "Isg15","Isg20","Oasl1" ,"Phf11b" ,"Plac8", "Rsad2","Rtp4","Sdf4", "Slfn4","Tnfrsf18" ,"Tnfrsf4","Tnfrsf9","Trex1","Usp18","Xaf1"))
ccl_genes = grep("Ccl", rownames(mouse_object@assays$RNA@counts), v = T)
MHC_genes =  c(grep("^H2-", rownames(mouse_object@assays$RNA@counts), v = T), "Cd74", "B2m")
comp_genes =  c(grep("^C1q", rownames(mouse_object@assays$RNA@counts), v = T))
ig_genes = c(grep("^Igj|^Igh|^Igk|^Igl", rownames(mouse_object@assays$RNA@counts), v = T))
hb_genes = c(grep("^Hba|^Hbb", rownames(mouse_object@assays$RNA@counts), v = T))
#markers = read.csv(paste0(workdir, "/", assay, "_markers.csv"))[, 1:2]
#LR_pairs = read.table(paste0(projdir, "/mouse_lr_pair.txt"), header = T)

########################################################################################################
########################################################################################################
##############fig 1h-k####
dot_clusters = c('c3', 'c4')
dot_mm = unique(mouse_object$key)
num_cells = setdiff(names(mouse_object$orig.ident[as.character(mouse_object@active.ident) %in% dot_clusters]), mouse_object$key %in% dot_mm); length(num_cells)


MouseModel = as.vector(mouse_object$key[num_cells])
ct = as.vector(mouse_object@active.ident[num_cells])
dot_data = data.frame(Mouse = MouseModel, ct = ct)
dot_factors = c("Mouse")

dot_sum1 = table(dot_data[,c("Mouse", "ct")])/rowSums(table(dot_data[,c("Mouse", "ct")]))*100
attributes(dot_sum1)$class <- "matrix"
dot_sum2 = uniquecombs(dot_data[, dot_factors], ordered = F); rownames(dot_sum2) = as.vector(dot_sum2$x);dot_sum2 = dot_sum2[rownames(dot_sum1),]
rownames(dot_sum2) = rownames(dot_sum1) = c()
dot_sum = cbind(dot_sum1, dot_sum2)
colnames(dot_sum) <- c('c3', 'c4', 'model')
# colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)] = paste0("ct_", colnames(dot_sum)[0:(length(unique(seurat@active.ident))-1)])
#dot_sum$Treatment = mapvalues(factor(dot_sum$Treatment, levels = c("UT","m1", "m2a")), from = c("UT","m1", "m2a"), to = c("Ctrl","m1", "m2a"))
dot_plot = melt(dot_sum)
dot_sum_lymph = dot_sum
dot_show = c("c3", "c4")
ggplot(dot_plot[dot_plot$ct %in% dot_show,], aes(x = Mouse, y = value)) +
  facet_wrap(. ~ ct,  strip.position = 'top', nrow = 1, scales = "free") +
  stat_summary(aes(fill = Treatment), color = "black", fun = mean, geom = "bar", size = 0.5, width = 0.6) +
  stat_summary(fun.data=mean_sdl, geom = "errorbar", width = 0.2, size = 0.3, color = "black") +
  scale_fill_manual(values = c("grey90", "grey60", "grey20")) + #c("cornflowerblue", "orange", "red")
  geom_point(size = 1.5, position = position_jitter(seed = 1, width = 0.2, height = 0), shape = 21, fill = "white") +
  geom_signif(comparisons = list(c("m1", "m2a"), c("Ctrl", "m2a")), test = "t.test", map_signif_level = F, family = "serif", textsize = 2, step_increase = 0.15) +
  theme_classic() +
  theme(text = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, angle = 0, vjust = 0.7),
        strip.background = element_blank(),
        axis.title.x =  element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm")) +
  ylab("% of lymphoid") 
ggsave(file = paste0("fig1_makeup_l.png"), width = 50, height = 60, units = "mm", dpi = 150, device = "png", bg = 'transparent')
ggsave(file = paste0("fig1_makeup_l.pdf"), width = 50, height = 60, units = "mm", dpi = 150, device = "pdf", bg = 'transparent')





