library(zellkonverter)
library(CellChat)
library(patchwork)
library(SeuratDisk)
library(Seurat)

# from scanpy via SeuratDisk [WORKING]
Convert("/athena/lodalab/scratch/ryc4001/scproj/human/erg/outs/h5ads/human_erg_lr_mut.h5ad", dest="h5seurat", overwrite=TRUE)
mut <- LoadH5Seurat("/athena/lodalab/scratch/ryc4001/scproj/human/erg/outs/h5ads/human_erg_lr_mut.h5seurat", assays="RNA")
counts <- GetAssayData(object=mut, slot="counts")
scaled <- normalizeData(counts, scale.factor=10000, do.log=TRUE)
mut <- SetAssayData(object=mut, slot="scale.data", new.data=scaled)
cellchat <- createCellChat(object=mut, group.by="phenotypes")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "human_erg_interactions.csv")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
pathways.show <- cellchat@netP$pathways

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

for (pathway in pathways.show){
  par(mfrow=c(1,1))
  try(netVisual_aggregate(cellchat, signaling = pathway, layout = "chord"))
}

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
try(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10))

cellchat.mut <- cellchat

Convert("/athena/lodalab/scratch/ryc4001/scproj/human/erg/outs/h5ads/human_erg_lr_wt.h5ad", dest="h5seurat", overwrite=TRUE)
wt <- LoadH5Seurat("/athena/lodalab/scratch/ryc4001/scproj/human/erg/outs/h5ads/human_erg_lr_wt.h5seurat", assays="RNA")
counts <- GetAssayData(object=wt, slot="counts")
scaled <- normalizeData(counts, scale.factor=10000, do.log=TRUE)
wt <- SetAssayData(object=wt, slot="scale.data", new.data=scaled)
cellchat <- createCellChat(object=wt, group.by="phenotypes")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "untreated_interactions.csv")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
pathways.show <- cellchat@netP$pathways

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

for (pathway in pathways.show){
  par(mfrow=c(1,1))
  try(netVisual_aggregate(cellchat, signaling = pathway, layout = "chord"))
}

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


cellchat.wt <- cellchat

object.list <- list(WT=cellchat.wt, MUT=cellchat.mut)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

try(gg1 <- netVisual_heatmap(cellchat))
try(gg2 <- netVisual_heatmap(cellchat, measure = "weight"))
try(gg1 + gg2)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T)

#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use="phenotypes")
#gg1

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
par(mfrow = c(1,2), xpd=TRUE)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
ht1 + ht2

#rankSimilarity(cellchat, type = "functional")

pos.dataset = "MUT"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "MUT",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
