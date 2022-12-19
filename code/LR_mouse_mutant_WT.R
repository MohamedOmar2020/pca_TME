# clean workspace
rm(list = ls())

# load libraries
library(data.table)
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(reticulate)
library(ComplexHeatmap)

# set the active conda environment
use_condaenv(condaenv = "scutils", required = TRUE)

# load the anndata module
ad <- import("anndata", convert = FALSE)

# load the mouse h5ad object
adata_mouse_all_wt <- ad$read_h5ad("./forCellChat/mouse_all_raw_wt.h5ad")
adata_mouse_all_mutant <- ad$read_h5ad("./forCellChat/mouse_all_raw_mutant.h5ad")

##########################
# access normalized data matrix
mouse_all_dataInput_wt <- t(py_to_r(adata_mouse_all_wt$X))
range(mouse_all_dataInput_wt)
rownames(mouse_all_dataInput_wt) <- rownames(py_to_r(adata_mouse_all_wt$var))
colnames(mouse_all_dataInput_wt) <- rownames(py_to_r(adata_mouse_all_wt$obs))

#####
mouse_all_dataInput_mutant <- t(py_to_r(adata_mouse_all_mutant$X))
range(mouse_all_dataInput_mutant)
rownames(mouse_all_dataInput_mutant) <- rownames(py_to_r(adata_mouse_all_mutant$var))
colnames(mouse_all_dataInput_mutant) <- rownames(py_to_r(adata_mouse_all_mutant$obs))

########################
## access meta data
mouse_all_metaData_wt <- py_to_r(adata_mouse_all_wt$obs)
mouse_all_meta_wt <- mouse_all_metaData_wt
# some cleaning
table(mouse_all_meta_wt$cluster)
mouse_all_meta_wt <- mouse_all_meta_wt[!is.na(mouse_all_meta_wt$cluster), ]
#levels(mouse_all_meta$cluster) <- c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages")
mouse_all_dataInput_wt <- mouse_all_dataInput_wt[, rownames(mouse_all_meta_wt)]

#########
mouse_all_metaData_mutant <- py_to_r(adata_mouse_all_mutant$obs)
mouse_all_meta_mutant <- mouse_all_metaData_mutant
# some cleaning
table(mouse_all_meta_mutant$cluster)
mouse_all_meta_mutant <- mouse_all_meta_mutant[!is.na(mouse_all_meta_mutant$cluster), ]
#levels(mouse_all_meta$cluster) <- c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages")
mouse_all_dataInput_mutant <- mouse_all_dataInput_mutant[, rownames(mouse_all_meta_mutant)]

#############################################
# creat a cellchat object
cellchat_mouse_all_wt <- createCellChat(object = mouse_all_dataInput_wt, meta = mouse_all_meta_wt, group.by = "compartment")
cellchat_mouse_all_mutant <- createCellChat(object = mouse_all_dataInput_mutant, meta = mouse_all_meta_mutant, group.by = "compartment")

# Add cell information into meta slot of the object (Optional)
cellchat_mouse_all_wt <- addMeta(cellchat_mouse_all_wt, meta = mouse_all_meta_wt)
cellchat_mouse_all_wt <- setIdent(cellchat_mouse_all_wt, ident.use = "compartment") # set "labels" as default cell identity
levels(cellchat_mouse_all_wt@idents) # show factor levels of the cell labels
groupSize_wt <- as.numeric(table(cellchat_mouse_all_wt@idents)) # number of cells in each cell group

####
cellchat_mouse_all_mutant <- addMeta(cellchat_mouse_all_mutant, meta = mouse_all_meta_mutant)
cellchat_mouse_all_mutant <- setIdent(cellchat_mouse_all_mutant, ident.use = "compartment") # set "labels" as default cell identity
levels(cellchat_mouse_all_mutant@idents) # show factor levels of the cell labels
groupSize_mutant <- as.numeric(table(cellchat_mouse_all_mutant@idents)) # number of cells in each cell group

######################################################
# clean
rm(adata_mouse_all_mutant, adata_mouse_all_wt, mouse_all_dataInput_mutant, mouse_all_dataInput_wt, mouse_all_meta_mutant, mouse_all_meta_wt, mouse_all_metaData_mutant, mouse_all_metaData_wt)

#######################################################################################################################################
# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB
cellchat_mouse_all_wt@DB <- CellChatDB.use
cellchat_mouse_all_mutant@DB <- CellChatDB.use


#################################################
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_mouse_all_wt <- subsetData(cellchat_mouse_all_wt) # This step is necessary even if using the whole database
future::plan("multicore", workers = 8) # do parallel
cellchat_mouse_all_wt <- identifyOverExpressedGenes(cellchat_mouse_all_wt)
cellchat_mouse_all_wt <- identifyOverExpressedInteractions(cellchat_mouse_all_wt)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_wt <- projectData(cellchat_mouse_all_wt, PPI.mouse)

#####
cellchat_mouse_all_mutant <- subsetData(cellchat_mouse_all_mutant) # This step is necessary even if using the whole database
cellchat_mouse_all_mutant <- identifyOverExpressedGenes(cellchat_mouse_all_mutant)
cellchat_mouse_all_mutant <- identifyOverExpressedInteractions(cellchat_mouse_all_mutant)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_mutant <- projectData(cellchat_mouse_all_mutant, PPI.mouse)

#####################################################
## Inference of cell-cell communication network
options(future.globals.maxSize=20*1024^3)

# Compute the communication probability and infer cellular communication network
cellchat_mouse_all_wt <- computeCommunProb(cellchat_mouse_all_wt)
cellchat_mouse_all_mutant <- computeCommunProb(cellchat_mouse_all_mutant)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_mouse_all_wt <- filterCommunication(cellchat_mouse_all_wt, min.cells = 10)
cellchat_mouse_all_mutant <- filterCommunication(cellchat_mouse_all_mutant, min.cells = 10)

######################
# Infer the cell-cell communication at a signaling pathway level
cellchat_mouse_all_wt <- computeCommunProbPathway(cellchat_mouse_all_wt)
cellchat_mouse_all_mutant <- computeCommunProbPathway(cellchat_mouse_all_mutant)

# Calculate the aggregated cell-cell communication network
cellchat_mouse_all_wt <- aggregateNet(cellchat_mouse_all_wt)
cellchat_mouse_all_mutant <- aggregateNet(cellchat_mouse_all_mutant)

## Get the dataframe of communications at the L/R level
df.net_mouse_wt <- subsetCommunication(cellchat_mouse_all_wt)
df.net_mouse_mutant <- subsetCommunication(cellchat_mouse_all_mutant)

## Get the dataframe of communications at the pathway level
df.net_mouse_pathway_wt <- subsetCommunication(cellchat_mouse_all_wt, slot.name = "netP")
df.net_mouse_pathway_mutant <- subsetCommunication(cellchat_mouse_all_mutant, slot.name = "netP")

# all the pathways
table(df.net_mouse_pathway_wt$pathway_name)
table(df.net_mouse_pathway_mutant$pathway_name)


#####################################3
# compute centrality
cellchat_mouse_all_wt <- netAnalysis_computeCentrality(cellchat_mouse_all_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_mutant <- netAnalysis_computeCentrality(cellchat_mouse_all_mutant, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

######################################
## Merge the 2 objects
object.list <- list(WT = cellchat_mouse_all_wt, mutant = cellchat_mouse_all_mutant)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

##########################################################################################################
# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

########################
# Compare the number of interactions and interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


#############################################
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

# Number of interactions between wt and mutant
png('./figures/cellchat_mutant_wt/Diff_N_interactions.png', width = 2000, height = 1500, res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# Compute and visualize the pathway distance in the learned joint manifold
#rankSimilarity(cellchat, type = "functional")

# Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('cyan3', 'indianred1'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c('cyan3', 'indianred1'))

png('./figures/cellchat_mutant_wt/Diff_informationFlow.png', width = 1500, height = 2000, res = 300)
gg1 
dev.off()

#####################################################
# Identify the upgulated and down-regulated signaling ligand-receptor pairs
gg1 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in mutants", angle.x = 45, remove.isolate = T, n.colors = 5, min.quantile = 0.25, line.on= T)
gg2 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in mutants", angle.x = 45, remove.isolate = T, n.colors = 5)

gg1 + gg2

#####################################################
# Identify dysfunctional signaling by using differential expression analysis

# define mutant as the positive dataset
pos.dataset = "mutant"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in mutant
net.up <- subsetCommunication(cellchat, net = net, datasets = "mutant",ligand.logFC = 0.2, receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in wildtypes, i.e.,downregulated in mutants
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.1, receptor.logFC = -0.1)

# do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# Chord diagram
png('./figures/cellchat_mutant_wt/up_mutants.png', width = 4500, height = 4000, res = 300)
netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = net.up, 
                     thresh = 0.01, lab.cex = 0.6, small.gap = 15, big.gap = 20,
                     legend.pos.y = 275, 
                     legend.pos.x = 70, 
                     title.name = "")
dev.off()

png('./figures/cellchat_mutant_wt/up_WT.png', width = 4500, height = 4000, res = 300)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = net.down, 
                     thresh = 0.01, lab.cex = 0.8, small.gap = 15, big.gap = 20, 
                     legend.pos.y = 275, 
                     legend.pos.x = 70, 
                     title.name = "")
dev.off()

####################3
# pathway level
png('./figures/cellchat_mutant_wt/up_mutants_pathways.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', net = net.up, 
                     thresh = 0.05, lab.cex = 0.5, small.gap = 1,
                     legend.pos.y = 180, 
                     legend.pos.x = 20, 
                     title.name = "")
dev.off()

png('./figures/cellchat_mutant_wt/up_WT_pathways.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list[[1]], slot.name = 'netP', net = net.down, 
                     thresh = 0.05, lab.cex = 0.8, small.gap = 3, 
                     legend.pos.y = 180, 
                     legend.pos.x = 20, 
                     title.name = "")
dev.off()

####################################################
# Compare outgoing (or incoming) signaling associated with each cell population
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = 'Wild types', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 12)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = 'Mutants', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 12)

png('./figures/cellchat_mutant_wt/Diff_heatmap2.png', width = 4000, height = 3000, res = 300)
draw(ht2 + ht1, ht_gap = unit(2, "cm"))
dev.off()

############
# save
saveRDS(cellchat, file = "./forCellChat/cellchat_mouse_WT_vs_mutants.rds")

# load
cellchat = readRDS("./forCellChat/cellchat_mouse_WT_vs_mutants.rds")




