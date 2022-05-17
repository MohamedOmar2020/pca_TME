
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
adata_mouse <- ad$read_h5ad("./forCellChat/adata_mouse_norm.h5ad")
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
mouse_meta <- mouse_meta[!(mouse_meta$cluster == 'cNA'), ]
mouse_dataInput <- mouse_dataInput[, rownames(mouse_meta)]

#############################################
# creat a cellchat object
cellchat_mouse <- createCellChat(object = mouse_dataInput, meta = mouse_meta, group.by = "cluster")

# Add cell information into meta slot of the object (Optional)
cellchat_mouse <- addMeta(cellchat_mouse, meta = mouse_meta)
cellchat_mouse <- setIdent(cellchat_mouse, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_mouse@idents)) # number of cells in each cell group

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB
cellchat_mouse@DB <- CellChatDB.use

#################################################
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_mouse <- subsetData(cellchat_mouse) # This step is necessary even if using the whole database
future::plan("multicore", workers = 5) # do parallel
cellchat_mouse <- identifyOverExpressedGenes(cellchat_mouse)
cellchat_mouse <- identifyOverExpressedInteractions(cellchat_mouse)
# project gene expression data onto PPI network (optional)
cellchat_mouse <- projectData(cellchat_mouse, PPI.mouse)

#####################################################
## Inference of cell-cell communication network

# Compute the communication probability and infer cellular communication network
cellchat_mouse <- computeCommunProb(cellchat_mouse)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_mouse <- filterCommunication(cellchat_mouse, min.cells = 10)

######################
# Infer the cell-cell communication at a signaling pathway level
cellchat_mouse <- computeCommunProbPathway(cellchat_mouse)

# Calculate the aggregated cell-cell communication network
cellchat_mouse <- aggregateNet(cellchat_mouse)

## Get the dataframe of communications at the L/R level
df.net_mouse <- subsetCommunication(cellchat_mouse)

## Get the dataframe of communications at the pathway level
df.net_mouse_pathway <- subsetCommunication(cellchat_mouse, slot.name = "netP")

# all the pathways
table(df.net_mouse_pathway$pathway_name)

#####################################################
#  visualize the aggregated cell-cell communication network:  the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat_mouse@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_mouse@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_mouse@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

###########
# examine the signaling sent from each cell group.
mat_mouse <- cellchat_mouse@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat_mouse)) {
  mat2 <- matrix(0, nrow = nrow(mat_mouse), ncol = ncol(mat_mouse), dimnames = dimnames(mat_mouse))
  mat2[i, ] <- mat_mouse[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_mouse), title.name = rownames(mat_mouse)[i])
}

##########################################################
## c5,6,7: NEPC, High expression of TGFB and WNT signaling 
## c3 and c4: High AR -- WNT signaling

# Visualization of cell-cell communication network
pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(6,7) # a numeric vector. 
netVisual_aggregate(cellchat_mouse, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = 'hierarchy')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_mouse, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_mouse, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_mouse, signaling = pathways.show, color.heatmap = "Reds")


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat_mouse, sources.use = 4, lab.cex = 0.2,legend.pos.y = 20)

# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat_mouse, signaling = "WNT")

##############
# Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat_mouse <- netAnalysis_computeCentrality(cellchat_mouse, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_mouse, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_mouse)
gg2 <- netAnalysis_signalingRole_scatter(cellchat_mouse, signaling = c("WNT", "ncWNT"))
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_mouse, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_mouse, pattern = "incoming")
ht1 + ht2

#####################################
## Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
selectK(cellchat_mouse, pattern = "outgoing")

nPatterns = 4
cellchat_mouse <- identifyCommunicationPatterns(cellchat_mouse, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat_mouse, pattern = "outgoing")

## Identify and visualize incoming communication pattern of secreting cells
selectK(cellchat_mouse, pattern = "incoming")

nPatterns = 4
cellchat_mouse <- identifyCommunicationPatterns(cellchat_mouse, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat_mouse, pattern = "incoming")

###############################################
## Manifold and classification learning analysis of signaling networks
# Identify signaling groups based on their functional similarity
cellchat_mouse <- computeNetSimilarity(cellchat_mouse, type = "functional")
cellchat_mouse <- netEmbedding(cellchat_mouse, type = "functional")
cellchat_mouse <- netClustering(cellchat_mouse, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat_mouse, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
cellchat_mouse <- computeNetSimilarity(cellchat_mouse, type = "structural")
cellchat_mouse <- netEmbedding(cellchat_mouse, type = "structural")
cellchat_mouse <- netClustering(cellchat_mouse, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat_mouse, type = "structural", label.size = 3.5)


netVisual_embeddingZoomIn(cellchat_mouse, type = "structural", nCol = 2)



############
# save
saveRDS(cellchat_mouse, file = "./forCellChat/cellchat_mouse_LS.rds")







