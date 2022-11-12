
# clean workspace
rm(list = ls())

# load libraries
library(data.table)
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(reticulate)

# set the active conda environment
use_condaenv(condaenv = "scutils", required = TRUE)

# load the anndata module
ad <- import("anndata", convert = FALSE)

# load the mouse h5ad object
adata_mouse_stroma <- ad$read_h5ad("./forCellChat/adata_mouse_stroma_raw.h5ad")
# access normalized data matrix
mouse_stroma_dataInput <- t(py_to_r(adata_mouse_stroma$X))
range(mouse_stroma_dataInput)
rownames(mouse_stroma_dataInput) <- rownames(py_to_r(adata_mouse_stroma$var))
colnames(mouse_stroma_dataInput) <- rownames(py_to_r(adata_mouse_stroma$obs))

## access meta data
mouse_stroma_metaData <- py_to_r(adata_mouse_stroma$obs)
mouse_stroma_meta <- mouse_stroma_metaData
# some cleaning
mouse_stroma_meta$cluster <- paste0("c", mouse_stroma_meta$cluster)
mouse_stroma_meta <- mouse_stroma_meta[!(mouse_stroma_meta$cluster == 'cNA'), ]
mouse_stroma_dataInput <- mouse_stroma_dataInput[, rownames(mouse_stroma_meta)]

#############################################
# creat a cellchat object
cellchat_mouse_stroma <- createCellChat(object = mouse_stroma_dataInput, meta = mouse_stroma_meta, group.by = "cluster")

# Add cell information into meta slot of the object (Optional)
cellchat_mouse_stroma <- addMeta(cellchat_mouse_stroma, meta = mouse_meta)
cellchat_mouse_stroma <- setIdent(cellchat_mouse_stroma, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_stroma@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_mouse_stroma@idents)) # number of cells in each cell group

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB
cellchat_mouse_stroma@DB <- CellChatDB.use

#################################################
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_mouse_stroma <- subsetData(cellchat_mouse_stroma) # This step is necessary even if using the whole database
future::plan("multicore", workers = 5) # do parallel
cellchat_mouse_stroma <- identifyOverExpressedGenes(cellchat_mouse_stroma)
cellchat_mouse_stroma <- identifyOverExpressedInteractions(cellchat_mouse_stroma)
# project gene expression data onto PPI network (optional)
cellchat_mouse_stroma <- projectData(cellchat_mouse_stroma, PPI.mouse)

#####################################################
## Inference of cell-cell communication network

# Compute the communication probability and infer cellular communication network
cellchat_mouse_stroma <- computeCommunProb(cellchat_mouse_stroma)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_mouse_stroma <- filterCommunication(cellchat_mouse_stroma, min.cells = 10)

######################
# Infer the cell-cell communication at a signaling pathway level
cellchat_mouse_stroma <- computeCommunProbPathway(cellchat_mouse_stroma)

# Calculate the aggregated cell-cell communication network
cellchat_mouse_stroma <- aggregateNet(cellchat_mouse_stroma)

## Get the dataframe of communications at the L/R level
df.net_mouse <- subsetCommunication(cellchat_mouse_stroma)

## Get the dataframe of communications at the pathway level
df.net_mouse_pathway <- subsetCommunication(cellchat_mouse_stroma, slot.name = "netP")

# all the pathways
table(df.net_mouse_pathway$pathway_name)

#####################################################
#  visualize the aggregated cell-cell communication network:  the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat_mouse_stroma@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_mouse_stroma@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_mouse_stroma@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

###########
# examine the signaling sent from each cell group.
mat_mouse <- cellchat_mouse_stroma@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat_mouse)) {
  mat2 <- matrix(0, nrow = nrow(mat_mouse), ncol = ncol(mat_mouse), dimnames = dimnames(mat_mouse))
  mat2[i, ] <- mat_mouse[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_mouse), title.name = rownames(mat_mouse)[i])
}

##########################################################
## c5,6,7: NEPC, High expression of TGFB, PERIOSTIN, and WNT signaling 
## c3 and c4: High AR -- WNT signaling

# Visualization of cell-cell communication network
pathways.show <- c("WNT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = c(6,7,8) # a numeric vector. 
netVisual_aggregate(cellchat_mouse_stroma, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = 'hierarchy')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_mouse_stroma, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_mouse_stroma, signaling = pathways.show, layout = "chord")

# Heatmap
png('./figures/cellchat/PERIOSTIN_heatmap.png', width = 2000, height = 2000, res = 300)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_mouse_stroma, signaling = 'PERIOSTIN', color.heatmap = "Reds")
dev.off()


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
png('./figures/cellchat/c3.png', width = 2000, height = 2000, res = 300)
netVisual_chord_gene(cellchat_mouse_stroma, sources.use = 4, lab.cex = 0.2,legend.pos.y = 20)
dev.off()

# Plot the signaling gene expression distribution using violin/dot plot
png('./figures/cellchat/Wnt_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_stroma, signaling = "WNT")
dev.off()

png('./figures/cellchat/PERIOSTIN_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_stroma, signaling = "PERIOSTIN")
dev.off()

png('./figures/cellchat/PDGF_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_stroma, signaling = "PDGF")
dev.off()

plotGeneExpression(cellchat_mouse_stroma, signaling = "THY1")

##############
# Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat_mouse_stroma <- netAnalysis_computeCentrality(cellchat_mouse_stroma, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
png('./figures/cellchat/WNT_signaling_role.png', width = 2000, height = 2000, res = 300)
netAnalysis_signalingRole_network(cellchat_mouse_stroma, signaling = 'WNT', width = 8, height = 2.5, font.size = 10)
dev.off()

png('./figures/cellchat/PERIOSTIN_signaling_role.png', width = 2000, height = 2000, res = 300)
netAnalysis_signalingRole_network(cellchat_mouse_stroma, signaling = 'PERIOSTIN', width = 8, height = 2.5, font.size = 10)
dev.off()

png('./figures/cellchat/PDGF_signaling_role.png', width = 2000, height = 2000, res = 300)
netAnalysis_signalingRole_network(cellchat_mouse_stroma, signaling = 'PDGF', width = 8, height = 2.5, font.size = 10)
dev.off()

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_mouse_stroma)
gg2 <- netAnalysis_signalingRole_scatter(cellchat_mouse_stroma, signaling = c("WNT", "ncWNT"))
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_mouse_stroma, pattern = "outgoing", font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_mouse_stroma, pattern = "incoming", font.size = 5)

png('./figures/cellchat/outgoing_signaling_role.png', width = 4000, height = 2000, res = 300)
ht1
dev.off()

#####################################
## Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
selectK(cellchat_mouse_stroma, pattern = "outgoing")

nPatterns = 4
cellchat_mouse_stroma <- identifyCommunicationPatterns(cellchat_mouse_stroma, pattern = "outgoing", k = nPatterns)

# river plot
png('./figures/cellchat/outgoing_signaling_pattern.png', width = 4000, height = 2000, res = 300)
netAnalysis_river(cellchat_mouse_stroma, pattern = "outgoing")
dev.off()

## Identify and visualize incoming communication pattern of secreting cells
selectK(cellchat_mouse_stroma, pattern = "incoming")

nPatterns = 4
cellchat_mouse_stroma <- identifyCommunicationPatterns(cellchat_mouse_stroma, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat_mouse_stroma, pattern = "incoming")

###############################################
## Manifold and classification learning analysis of signaling networks
# Identify signaling groups based on their functional similarity
cellchat_mouse_stroma <- computeNetSimilarity(cellchat_mouse_stroma, type = "functional")
cellchat_mouse_stroma <- netEmbedding(cellchat_mouse_stroma, type = "functional")
cellchat_mouse_stroma <- netClustering(cellchat_mouse_stroma, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat_mouse_stroma, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
cellchat_mouse_stroma <- computeNetSimilarity(cellchat_mouse_stroma, type = "structural")
cellchat_mouse_stroma <- netEmbedding(cellchat_mouse_stroma, type = "structural")
cellchat_mouse_stroma <- netClustering(cellchat_mouse_stroma, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat_mouse_stroma, type = "structural", label.size = 3.5)


netVisual_embeddingZoomIn(cellchat_mouse_stroma, type = "structural", nCol = 2)



############
# save
saveRDS(cellchat_mouse_stroma, file = "./forCellChat/cellchat_mouse_stroma_LS.rds")







