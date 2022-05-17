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

# load the human h5ad object
adata_human <- ad$read_h5ad("./forCellChat/adata_human_norm.h5ad")
# access normalized data matrix
human_dataInput <- t(py_to_r(adata_human$X))
range(human_dataInput)
rownames(human_dataInput) <- rownames(py_to_r(adata_human$var))
colnames(human_dataInput) <- rownames(py_to_r(adata_human$obs))

## access meta data
human_metaData <- py_to_r(adata_human$obs)
human_meta <- human_metaData
# some cleaning
human_meta$cluster <- paste0("c", human_meta$cluster)
table(human_meta$cluster)
human_meta <- human_meta[!(human_meta$cluster == 'cnan'), ]
human_dataInput <- human_dataInput[, rownames(human_meta)]

#############################################
# creat a cellchat object
cellchat_human <- createCellChat(object = human_dataInput, meta = human_meta, group.by = "cluster")

# Add cell information into meta slot of the object (Optional)
cellchat_human <- addMeta(cellchat_human, meta = human_meta)
cellchat_human <- setIdent(cellchat_human, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_human@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_human@idents)) # number of cells in each cell group

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB
cellchat_human@DB <- CellChatDB.use

#################################################
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_human <- subsetData(cellchat_human) # This step is necessary even if using the whole database
future::plan("multicore", workers = 5) # do parallel
cellchat_human <- identifyOverExpressedGenes(cellchat_human)
cellchat_human <- identifyOverExpressedInteractions(cellchat_human)
# project gene expression data onto PPI network (optional)
cellchat_human <- projectData(cellchat_human, PPI.mouse)

#####################################################
## Inference of cell-cell communication network

# Compute the communication probability and infer cellular communication network
cellchat_human <- computeCommunProb(cellchat_human)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_human <- filterCommunication(cellchat_human, min.cells = 10)

######################
# Infer the cell-cell communication at a signaling pathway level
cellchat_human <- computeCommunProbPathway(cellchat_human)

# Calculate the aggregated cell-cell communication network
cellchat_human <- aggregateNet(cellchat_human)

## Get the dataframe of communications at the L/R level
df.net_human <- subsetCommunication(cellchat_human)

## Get the dataframe of communications at the pathway level
df.net_human_pathway <- subsetCommunication(cellchat_human, slot.name = "netP")

# all the pathways
table(df.net_human_pathway$pathway_name)

#####################################################
#  visualize the aggregated cell-cell communication network:  the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat_human@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_human@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_human@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

###########
# examine the signaling sent from each cell group.
mat_human <- cellchat_human@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat_human)) {
  mat2 <- matrix(0, nrow = nrow(mat_human), ncol = ncol(mat_human), dimnames = dimnames(mat_human))
  mat2[i, ] <- mat_human[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_human), title.name = rownames(mat_human)[i])
}

##########################################################
## c5,6,7: NEPC, High expression of TGFB and WNT signaling 
## c3 and c4: High AR -- WNT signaling

# Visualization of cell-cell communication network
pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(6,7) # a numeric vector. 
netVisual_aggregate(cellchat_human, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = 'hierarchy')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_human, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_human, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_human, signaling = pathways.show, color.heatmap = "Reds")


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat_human, sources.use = 4, lab.cex = 0.2,legend.pos.y = 20)

# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat_human, signaling = "WNT")

##############
# Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat_human <- netAnalysis_computeCentrality(cellchat_human, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_human, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_human)
gg2 <- netAnalysis_signalingRole_scatter(cellchat_human, signaling = c("WNT", "ncWNT"))
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_human, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_human, pattern = "incoming")
ht1 + ht2

#####################################
## Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
selectK(cellchat_human, pattern = "outgoing")

nPatterns = 4
cellchat_human <- identifyCommunicationPatterns(cellchat_human, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat_human, pattern = "outgoing")

## Identify and visualize incoming communication pattern of secreting cells
selectK(cellchat_human, pattern = "incoming")

nPatterns = 3
cellchat_human <- identifyCommunicationPatterns(cellchat_human, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat_human, pattern = "incoming")

###############################################
## Manifold and classification learning analysis of signaling networks
# Identify signaling groups based on their functional similarity
cellchat_human <- computeNetSimilarity(cellchat_human, type = "functional")
cellchat_human <- netEmbedding(cellchat_human, type = "functional")
cellchat_human <- netClustering(cellchat_human, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat_human, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
cellchat_human <- computeNetSimilarity(cellchat_human, type = "structural")
cellchat_human <- netEmbedding(cellchat_human, type = "structural")
cellchat_human <- netClustering(cellchat_human, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat_human, type = "structural", label.size = 3.5)


netVisual_embeddingZoomIn(cellchat_human, type = "structural", nCol = 2)



############
# save
saveRDS(cellchat_human, file = "./forCellChat/cellchat_human_LS.rds")






