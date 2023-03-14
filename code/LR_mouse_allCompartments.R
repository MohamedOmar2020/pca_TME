# clean workspace
rm(list = ls())

# load libraries
library(data.table)
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(reticulate)

########################
# 
netVisual_chord_gene <- function (object, slot.name = "net", color.use = NULL, signaling = NULL, 
                                  pairLR.use = NULL, net = NULL, sources.use = NULL, targets.use = NULL, 
                                  lab.cex = 0.8, small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03), 
                                  link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, 
                                  reduce = -1, transparency = 0.4, link.border = NA, title.name = NULL, 
                                  legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE, 
                                  thresh = 0.05, ...) 
{
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }
  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }
  if (is.null(net)) {
    prob <- slot(object, "net")$prob
    pval <- slot(object, "net")$pval
    prob[pval > thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source", "target", "interaction_name")
    pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", 
                                              "pathway_name", "ligand", "receptor", "annotation", 
                                              "evidence"))
    idx <- match(net$interaction_name, rownames(pairLR))
    temp <- pairLR[idx, ]
    net <- cbind(net, temp)
  }
  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], 
                                 pairLR.use = object@LR$LRsig, key = "pathway_name", 
                                 matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }
  if (!is.null(pairLR.use)) {
    if ("interaction_name" %in% colnames(pairLR.use)) {
      net <- subset(net, interaction_name %in% pairLR.use$interaction_name)
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
    }
  }
  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source", "target", "pathway_name", 
                                "prob"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net <- net %>% dplyr::group_by(source_target, pathway_name) %>% 
      dplyr::summarize(prob = sum(prob))
    a <- stringr::str_split(net$source_target, "sourceTotarget", 
                            simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net$ligand <- net$pathway_name
    net$receptor <- " "
  }
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      sources.use <- levels(object@idents)[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  else {
    sources.use <- levels(object@idents)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      targets.use <- levels(object@idents)[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  else {
    targets.use <- levels(object@idents)
  }
  df <- subset(net, prob > 0)
  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }
  if (length(unique(net$ligand)) == 1) {
    message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  }
  df$id <- 1:nrow(df)
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(" ", 
                                                       j - 1), collapse = ""))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(" ", 
                                                           j - 1), collapse = ""))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }
  cell.order.sources <- levels(object@idents)[levels(object@idents) %in% 
                                                sources.use]
  cell.order.targets <- levels(object@idents)[levels(object@idents) %in% 
                                                targets.use]
  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]
  order.source <- unique(df.ordered.source[, c("ligand", "source")])
  order.target <- unique(df.ordered.target[, c("receptor", 
                                               "target")])
  order.sector <- c(order.source$ligand, order.target$receptor)
  if (is.null(color.use)) {
    color.use = scPalette(nlevels(object@idents))
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source, 
                                                                         df$target))]
  }
  else if (is.null(names(color.use))) {
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source, 
                                                                         df$target))]
  }
  edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector
  df.plot <- df.ordered.source[, c("ligand", "receptor", "prob")]
  if (directional == 2) {
    link.arr.type = "triangle"
  }
  else {
    link.arr.type = "big.arrow"
  }
  circlize::circos.clear()
  circlize::chordDiagram(df.plot, order = order.sector, col = edge.color, 
                         grid.col = grid.col, transparency = transparency, link.border = link.border, 
                         directional = directional, direction.type = c("diffHeight", "arrows"), link.arr.type = link.arr.type, link.arr.length = 0.1, annotationTrack = "grid", 
                         annotationTrackHeight = annotationTrackHeight, preAllocateTracks = list(track.height = max(strwidth(order.sector))), 
                         small.gap = small.gap, big.gap = big.gap, link.visible = link.visible, 
                         scale = scale, link.target.prop = link.target.prop, reduce = reduce, 
                         ...)
  circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = circlize::get.cell.meta.data("xlim")
    xplot = circlize::get.cell.meta.data("xplot")
    ylim = circlize::get.cell.meta.data("ylim")
    sector.name = circlize::get.cell.meta.data("sector.index")
    circlize::circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                          niceFacing = TRUE, adj = c(0, 0.5), cex = lab.cex)
  }, bg.border = NA)
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), 
                                  type = "grid", legend_gp = grid::gpar(fill = color.use), 
                                  title = "immune cell types", legend_height = unit(12, "cm"), labels_gp = grid::gpar(fontsize=12), title_gp = grid::gpar(fontsize = 14))
    ComplexHeatmap::draw(lgd, x = unit(1, "npc") - unit(legend.pos.x, 
                                                        "mm"), y = unit(legend.pos.y, "mm"), just = c("right", 
                                                                                                      "bottom"))
  }
  circlize::circos.clear()
  if (!is.null(title.name)) {
    text(-0, 1.02, title.name, cex = 1)
  }
  gg <- recordPlot()
  return(gg)
}

###################################
# set the active conda environment
use_condaenv(condaenv = "scutils", required = TRUE)

# load the anndata module
ad <- import("anndata", convert = FALSE)

# load the mouse h5ad object
adata_mouse_all <- ad$read_h5ad("./forCellChat/mouse_all_raw.h5ad")
# access normalized data matrix
mouse_all_dataInput <- t(py_to_r(adata_mouse_all$X))
range(mouse_all_dataInput)
rownames(mouse_all_dataInput) <- rownames(py_to_r(adata_mouse_all$var))
colnames(mouse_all_dataInput) <- rownames(py_to_r(adata_mouse_all$obs))

## access meta data
mouse_all_metaData <- py_to_r(adata_mouse_all$obs)
mouse_all_meta <- mouse_all_metaData
# some cleaning
table(mouse_all_meta$cluster)
mouse_all_meta <- mouse_all_meta[!is.na(mouse_all_meta$cluster), ]
#levels(mouse_all_meta$cluster) <- c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages")
mouse_all_dataInput <- mouse_all_dataInput[, rownames(mouse_all_meta)]

#############################################
# creat a cellchat object
cellchat_mouse_all <- createCellChat(object = mouse_all_dataInput, meta = mouse_all_meta, group.by = "cluster")

# Add cell information into meta slot of the object (Optional)
cellchat_mouse_all <- addMeta(cellchat_mouse_all, meta = mouse_all_meta)
cellchat_mouse_all <- setIdent(cellchat_mouse_all, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_mouse_all@idents)) # number of cells in each cell group

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB
cellchat_mouse_all@DB <- CellChatDB.use

######################################################
# clean
rm(mouse_all_dataInput, mouse_all_meta)

#################################################
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_mouse_all <- subsetData(cellchat_mouse_all) # This step is necessary even if using the whole database
future::plan("multicore", workers = 8) # do parallel
cellchat_mouse_all <- identifyOverExpressedGenes(cellchat_mouse_all)
cellchat_mouse_all <- identifyOverExpressedInteractions(cellchat_mouse_all)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all <- projectData(cellchat_mouse_all, PPI.mouse)

#####################################################
## Inference of cell-cell communication network
options(future.globals.maxSize=20*1024^3)

# Compute the communication probability and infer cellular communication network
cellchat_mouse_all <- computeCommunProb(cellchat_mouse_all, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_mouse_all <- filterCommunication(cellchat_mouse_all, min.cells = 10)

######################
# Infer the cell-cell communication at a signaling pathway level
cellchat_mouse_all <- computeCommunProbPathway(cellchat_mouse_all)

# Calculate the aggregated cell-cell communication network
cellchat_mouse_all <- aggregateNet(cellchat_mouse_all)

## Get the dataframe of communications at the L/R level
df.net_mouse <- subsetCommunication(cellchat_mouse_all)

## Get the dataframe of communications at the pathway level
df.net_mouse_pathway <- subsetCommunication(cellchat_mouse_all, slot.name = "netP")

# all the pathways
table(df.net_mouse_pathway$pathway_name)


Periostin_int <- df.net_mouse %>%
  dplyr::filter(pathway_name == "PERIOSTIN")

#####################################################
#  visualize the aggregated cell-cell communication network:  the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat_mouse_all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_mouse_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_mouse_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

###########
# examine the signaling sent from each cell group.
mat_mouse <- cellchat_mouse_all@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat_mouse)) {
  mat2 <- matrix(0, nrow = nrow(mat_mouse), ncol = ncol(mat_mouse), dimnames = dimnames(mat_mouse))
  mat2[i, ] <- mat_mouse[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_mouse), title.name = rownames(mat_mouse)[i])
}

##########################################################
levels(cellchat_mouse_all@idents)


################################################################
# Chord diagram
# png('./figures/cellchat_allCompartments/WNT_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'WNT', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/PERIOSTIN_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'PERIOSTIN', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/TGFb_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'TGFb', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/ncWNT_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'ncWNT', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/PDGF_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'PDGF', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/PTN_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'PTN', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/COLLAGEN_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'COLLAGEN', sources.use= c('c5', 'c6', 'c7'), targets.use = 'epithelium', title.space	= 1, layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/TGFb_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_mouse_all, signaling = 'TGFb', targets.use= c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), title.space	= 1, layout = "chord")
# dev.off()


#########################################################
## signaling from PRN clusters to epithelium

# # pathway level
# png('./figures/cellchat_allCompartments/PRN_epith_pathways.png', width = 3500, height = 3000, res = 350)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      slot = 'netP', 
#                      sources.use= c('c5', 'c6', 'c7'), 
#                      targets.use = 'epithelium', 
#                      lab.cex = 1, 
#                      legend.pos.y = 180, 
#                      legend.pos.x = 50, 
#                      thresh = 0.01, small.gap = 5,
#                      #title.name = 'Signaling pathways from the PRN clusters to the epithelium',	
#                      #reduce = -1,
# )
# dev.off()
# 
# # LR interactions
# png('./figures/cellchat_allCompartments/PRN_epith_collagen.png', width = 2000, height = 1500, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      signaling = 'COLLAGEN', 
#                      sources.use= c('c5', 'c6', 'c7'), 
#                      targets.use = 'epithelium', 
#                      lab.cex = 0.5, 
#                      legend.pos.y = 50, 
#                      legend.pos.x = 20,
#                      reduce = -1, small.gap = 3,
#                      title.name = 'Collagen signaling network from the PRN clusters to the tumor epithelium'	
# )
# dev.off()

#########################################################
## signaling from WNT clusters to epithelium

# pathway level

# tiff('./figures/cellchat_allCompartments/epith_C3C4_pathways.tiff', width = 2500, height = 2000, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      slot = 'netP', 
#                      sources.use = 'epithelium', 
#                      targets.use= c('c3', 'c4'), 
#                      lab.cex = 0.5, 
#                      legend.pos.y = 165, 
#                      legend.pos.x = 35, 
#                      thresh = 0.05, small.gap = 3,
#                      title.name = 'Signaling pathways from the epithelium to c3 and c4',	
#                      reduce = -1,
# )
# dev.off()

#############################
# all pathways epithelium to stroma
# tiff('./figures/cellchat_allCompartments/epith_stroma_pathways.tiff', width = 2500, height = 2000, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      slot = 'netP', 
#                      sources.use = 'epithelium',
#                      targets.use= c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
#                      lab.cex = 0.5, 
#                      legend.pos.y = 50, 
#                      legend.pos.x = 20, 
#                      thresh = 0.05, small.gap = 3,
#                      title.name = 'Signaling pathways from the epithelium to the stroma',	
#                      reduce = -1,
# )
# dev.off()

#####################
# all pathways stroma to epithelium
# tiff('./figures/cellchat_allCompartments/stroma_epith_pathways.tiff', width = 2500, height = 2000, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      slot = 'netP', 
#                      targets.use = 'epithelium',
#                      sources.use= c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
#                      lab.cex = 0.3, 
#                      legend.pos.y = 50, 
#                      legend.pos.x = 20, 
#                      thresh = 0.05, small.gap = 2,
#                      title.name = 'Signaling pathways from the stroma to the epithelium',	
#                      reduce = -1,
# )
# dev.off()

########################
## LR interactions epthelium to stroma

# WNT
# tiff('./figures/cellchat_allCompartments/epith_stroma_wnt.tiff', width = 2000, height = 1500, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      signaling = c('WNT'), 
#                      sources.use = 'epithelium', 
#                      targets.use= c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
#                      lab.cex = 0.6, 
#                      legend.pos.y = 50, 
#                      legend.pos.x = 20,
#                      reduce = -1,
#                      thresh = 0.05,
#                      title.name = 'WNT signaling network from the epithelium to the stroma'	
# )
# dev.off()

## CCL Treg
# tiff('./figures/cellchat_allCompartments/immun_stroma_CCL.tiff', width = 2000, height = 1500, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      signaling = c('CCL'), 
#                      targets.use = 'Treg', 
#                      sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
#                      lab.cex = 0.6, 
#                      legend.pos.y = 50, 
#                      legend.pos.x = 15,
#                      reduce = -1,
#                      thresh = 0.05,
#                      title.name = 'CCL signaling network from the stroma to Treg'	
# )
# dev.off()

## CCL mono/macrophages
# tiff('./figures/cellchat_allCompartments/immun_macrophages_stroma_CCL.tiff', width = 2000, height = 1500, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      signaling = c('CCL'), 
#                      targets.use = 'monocytes/macrophages', 
#                      sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
#                      lab.cex = 0.6, 
#                      legend.pos.y = 50, 
#                      legend.pos.x = 1,
#                      reduce = -1,
#                      thresh = 0.05,
#                      title.name = 'CCL signaling network from the stroma to Macrophages'	
# )
# dev.off()

##############################################################
##############################################################
## for the figures
# color map
levels(cellchat_mouse_all@idents)
colMap <- c('powderblue', 'plum2', 'violetred4', 'slategray4', 'blue', 'darkorange1', 'forestgreen', 'red2', 'darkorchid', 'saddlebrown', 'hotpink1', 'olivedrab3', 'lemonchiffon4', 'bisque2', 'salmon3')



# Figure 3: common clusters

# to epithelium
png('./figures/cellchat_allCompartments/c0c1c2_epithelium.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                    slot.name = "netP",
                    targets.use = 'epithelium', 
                    sources.use = c('c0', 'c1', 'c2'), 
                    lab.cex = 1, 
                    legend.pos.y = 165, 
                    legend.pos.x = 25,
                    #reduce = -1,
                    thresh = 0.01, 
                    small.gap = 5,
                    color.use = colMap
                    #title.name = 'Signling networks from common clusters to the epithelium'	
)
dev.off()

# opposite
# tiff('./figures/cellchat_allCompartments/epithelium_c0c1c2.tiff', width = 2500, height = 3000, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      slot.name = "netP",
#                      sources.use = 'epithelium', 
#                      targets.use = c('c0', 'c1', 'c2'), 
#                      lab.cex = 0.8, 
#                      legend.pos.y = 100, 
#                      legend.pos.x = 8,
#                      reduce = -1,
#                      thresh = 0.01, 
#                      small.gap = 3,
#                      title.name = 'Signling networks from epithelium to common clusters'	
# )
# dev.off()

# to immune
png('./figures/cellchat_allCompartments/c0c1c2_immune.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "netP",
                     targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     sources.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.01, 
                     small.gap = 4,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()




##############################################################
##########
# Figure 4 : E
png('./figures/cellchat_allCompartments/new/c3c4_immune.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "netP",
                     #signaling = c('COMPLEMENT', 'CCL'),
                     targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     sources.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signaling Networks from c3 and c4 to immune cells'	
)
dev.off()


png('./figures/cellchat_allCompartments/new/C3C4_epith_pathways.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot = 'netP', 
                     sources.use= c('c3', 'c4'), 
                     targets.use = 'epithelium', 
                     lab.cex = 1, 
                     legend.pos.y = 165, 
                     legend.pos.x = 25, 
                     thresh = 0.05, small.gap = 5,
                     color.use = colMap,
                     #title.name = 'Signaling pathways from c3 and c4 to the epithelium',	
                     #reduce = -1,
)
dev.off()


##############################################################
# figure 4: F
png('./figures/cellchat_allCompartments/new/c5c6c7_immune.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "netP",
                     #signaling = c('COMPLEMENT', 'CCL'),
                     targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     sources.use = c('c5', 'c6', 'c7'), 
                     lab.cex = 0.9, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.01, 
                     small.gap = 3.5,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signaling Networks from PRN clusters to immune cells'	
)
dev.off()

png('./figures/cellchat_allCompartments/new/c5c6c7_epithelium.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "netP",
                     #signaling = c('COMPLEMENT', 'CCL'),
                     targets.use = c("epithelium"), 
                     sources.use = c('c5', 'c6', 'c7'), 
                     lab.cex = 1, 
                     legend.pos.y = 165, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signaling Networks from PRN clusters to epithelium'	
)
dev.off()


##############################################################
##############################################################
## Figure S3: 

# CCL signaling

png('./figures/cellchat_allCompartments/CCL_Chord.png', width = 3500, height = 3000, res = 350)
netVisual_aggregate(cellchat_mouse_all, 
                    signaling = 'CCL',
                    sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"),
                    targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"),
                    vertex.label.cex = 1, 
                    title.space	= 0,
                    #legend.pos.y = 170, 
                    #legend.pos.x = 20,
                    small.gap = 5,
                    color.use = colMap,
                    layout = "chord"
                    )
dev.off()


# Tregs to stroma
png('./figures/cellchat_allCompartments/Treg_stroma.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     sources.use = c('Treg'), 
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     lab.cex = 1, 
                     legend.pos.y = 166, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     color.use = colMap,
                     title.name = ''	
)
dev.off()

# macrophages to stroma
png('./figures/cellchat_allCompartments/macrophages_stroma.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c('monocytes/macrophages'), 
                     lab.cex = 0.9, 
                     legend.pos.y = 166, 
                     legend.pos.x = 4,
                     #reduce = -1,
                     thresh = 0.01, 
                     small.gap = 3.4,
                     #directional = 2,
                     color.use = colMap,
                     title.name = ''	
)
dev.off()

###############################################################
###############################################################
## Figure S4

# stroma-stroma interactions
png('./figures/cellchat_allCompartments/stroma_stroma.png', width = 4000, height = 4500, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "netP",
                     #signaling = c('COMPLEMENT', 'CCL'),
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     lab.cex = 0.5, 
                     legend.pos.y = 200, 
                     legend.pos.x = 4,
                     reduce = -1,
                     transparency = 0,
                     thresh = 0.01, 
                     small.gap = 1.2,
                     #directional = 2, 
                     color.use = colMap,
                     title.name = ''	
)
dev.off()

# WNT
png('./figures/cellchat_allCompartments/stroma_stroma_wnt.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     signaling = c('WNT', 'ncWNT'),
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     lab.cex = 0.8, 
                     legend.pos.y = 170, 
                     legend.pos.x = 20,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 4,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'WNT stroma stroma interactions'	
)
dev.off()

# POSTN
png('./figures/cellchat_allCompartments/Postn_stromaOnly.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     signaling = c('PERIOSTIN'),
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     lab.cex = 0.8, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     color.use = colMap,
                     title.name = ''	
)
dev.off()


###############################################################
# for S4 legend: immune cell types
levels(cellchat_mouse_all@idents)
colMap2 <- c('powderblue', 'plum2', 'violetred4', 'slategray4', 'blue', 'darkorange1', 'forestgreen', 'red2', 'darkorchid', 'saddlebrown', 'hotpink1', 'olivedrab3', 'lemonchiffon4', 'bisque2', 'salmon3')

png('./figures/cellchat_allCompartments/immuneLegend.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = 'netP',
                    signaling = 'CD96',
                    sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"),
                    targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"),
                    lab.cex = 1, 
                    legend.pos.y = 170, 
                    legend.pos.x = 10,
                    small.gap = 5,
                    color.use = colMap2,
                    title.name = ""
)
dev.off()
###############################################################
# figure S4
# complement and CCL signaling common clusters to immune
# png('./figures/cellchat_allCompartments/c0c1c2_immuneSpecific.png', width = 2500, height = 3000, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      slot.name = "net",
#                      signaling = c('COMPLEMENT', 'CCL'),
#                      targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
#                      sources.use = c('c0', 'c1', 'c2'), 
#                      lab.cex = 0.8, 
#                      legend.pos.y = 175, 
#                      legend.pos.x = 5,
#                      reduce = -1,
#                      thresh = 0.01, 
#                      small.gap = 5,
#                      #directional = 2,
#                      title.name = 'Complement and CCL Signling networks from the common clusters to immune cells'	
# )
# dev.off()
# 
# png('./figures/cellchat_allCompartments/immuneSpecific.png', width = 2500, height = 3000, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, 
#                      slot.name = "net",
#                      signaling = c('COMPLEMENT', 'CCL'),
#                      targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
#                      sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
#                      lab.cex = 0.8, 
#                      legend.pos.y = 175, 
#                      legend.pos.x = 5,
#                      reduce = -1,
#                      thresh = 0.01, 
#                      small.gap = 5,
#                      #directional = 2,
#                      title.name = 'Complement and CCL Signling networks from the stroma clusters to immune cells'	
# )
# dev.off()


##############################################################################################################################
##############################################################################################################################

netAnalysis_signalingRole_heatmap(cellchat_mouse_all, signaling = c("COMPLEMENT", "CCL"), slot.name = 'net')

#############################
# wnt signaling only stroma
png('./figures/cellchat_allCompartments/wnt_stromaOnly.png', width = 2500, height = 3000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     signaling = c('WNT', 'ncWNT'),
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     lab.cex = 0.8, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     title.name = ''	
)
dev.off()

#############################
# POSTN signaling only stroma


#############################
# THY1 signaling only stroma
png('./figures/cellchat_allCompartments/THY1_stromaOnly.png', width = 2500, height = 3000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     signaling = c('THY1'),
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     lab.cex = 0.8, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     title.name = ''	
)
dev.off()

#############################
##############################
# wnt signaling stroma-epithelium
png('./figures/cellchat_allCompartments/wnt_stroma_epith.png', width = 2500, height = 3000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     signaling = c('WNT', 'ncWNT'),
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c('epithelium'), 
                     lab.cex = 0.8, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     title.name = ''	
)
dev.off()

##############################
# wnt signaling stroma-immune
png('./figures/cellchat_allCompartments/wnt_stroma_immune.png', width = 2500, height = 3000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     signaling = c('WNT', 'ncWNT'),
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     lab.cex = 0.8, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     title.name = ''	
)
dev.off()


##############################
tiff('./figures/cellchat_allCompartments/c0c1c2_immuneSpecific.tiff', width = 3000, height = 3000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     signaling = c('COMPLEMENT', 'CCL'),
                     targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     sources.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 0.6, 
                     legend.pos.y = 12, 
                     legend.pos.x = 8,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 5,
                     #directional = 2,
                     title.name = 'Complement and CCL Signling networks from common clusters to immune cells'	
)
dev.off()

#################################


# CD8 to stroma
png('./figures/cellchat_allCompartments/CD8_stroma.png', width = 3500, height = 4000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c("NK/cytoxic T lymphocytes"), 
                     lab.cex = 0.8, 
                     legend.pos.y = 250, 
                     legend.pos.x = 5,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 2,
                     #directional = 2,
                     title.name = ''	
)
dev.off()

# Treg to stroma
png('./figures/cellchat_allCompartments/Treg_stroma.png', width = 3500, height = 4000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     sources.use = c("Treg"), 
                     lab.cex = 0.8, 
                     legend.pos.y = 250, 
                     legend.pos.x = 5,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 2,
                     #directional = 2,
                     title.name = ''	
)
dev.off()

# stroma to CD8
png('./figures/cellchat_allCompartments/stroma_CD8.png', width = 3500, height = 4000, res = 300)
netVisual_chord_gene(cellchat_mouse_all, 
                     slot.name = "net",
                     sources.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), 
                     targets.use = c("NK/cytoxic T lymphocytes"), 
                     lab.cex = 0.6, 
                     legend.pos.y = 170, 
                     legend.pos.x = 3,
                     reduce = -1,
                     thresh = 0.01, 
                     small.gap = 1,
                     #directional = 2,
                     title.name = ''	
)
dev.off()
##################################################
# Plot the signaling gene expression distribution using violin/dot plot

png('./figures/cellchat_allCompartments/PERIOSTIN_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_all, signaling = "PERIOSTIN")
dev.off()

png('./figures/cellchat_allCompartments/PDGF_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_all, signaling = "PDGF")
dev.off()

#################################
# Aggregated cell communication network
mat <- cellchat_mouse_all@net$weight
par(mfrow = c(3,5), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  png(filename = paste0('./figures/cellchat_allCompartments/', gsub('/', '', rownames(mat)[i]), '.png'), width = 2500, height = 3000, res = 300)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}


##############################################################

# Heatmap
png('./figures/cellchat_allCompartments/PERIOSTIN_heatmap.png', width = 2000, height = 2000, res = 300)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_mouse_all, signaling = 'PERIOSTIN', color.heatmap = "Reds")
dev.off()


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# png('./figures/cellchat/c3.png', width = 2000, height = 2000, res = 300)
# netVisual_chord_gene(cellchat_mouse_all, sources.use = 4, lab.cex = 0.2,legend.pos.y = 20)
# dev.off()

# Plot the signaling gene expression distribution using violin/dot plot
png('./figures/cellchat_allCompartments/Wnt_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_all, signaling = "WNT")
dev.off()

png('./figures/cellchat_allCompartments/PERIOSTIN_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_all, signaling = "PERIOSTIN")
dev.off()

png('./figures/cellchat_allCompartments/PDGF_signaling_expression.png', width = 2000, height = 2000, res = 300)
plotGeneExpression(cellchat_mouse_all, signaling = " v")
dev.off()

plotGeneExpression(cellchat_mouse_all, signaling = "THY1")

##############
# Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat_mouse_all <- netAnalysis_computeCentrality(cellchat_mouse_all, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

png('./figures/cellchat_allCompartments/complement_ccl_signaling.png', width = 3300, height = 1700, res = 300)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_mouse_all, pattern = "outgoing", signaling = c("COMPLEMENT", "CCL"))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_mouse_all, pattern = "incoming", signaling = c("COMPLEMENT", "CCL"))
ht1 + ht2
dev.off()

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
png('./figures/cellchat/WNT_signaling_role.png', width = 2000, height = 2000, res = 300)
netAnalysis_signalingRole_network(cellchat_mouse_all, signaling = 'WNT', width = 8, height = 2.5, font.size = 10)
dev.off()

png('./figures/cellchat/PERIOSTIN_signaling_role.png', width = 2000, height = 2000, res = 300)
netAnalysis_signalingRole_network(cellchat_mouse_all, signaling = 'PERIOSTIN', width = 8, height = 2.5, font.size = 10)
dev.off()

png('./figures/cellchat/PDGF_signaling_role.png', width = 2000, height = 2000, res = 300)
netAnalysis_signalingRole_network(cellchat_mouse_all, signaling = 'PDGF', width = 8, height = 2.5, font.size = 10)
dev.off()

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_mouse_all)
gg2 <- netAnalysis_signalingRole_scatter(cellchat_mouse_all, signaling = c("WNT", "ncWNT"))
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_mouse_all, pattern = "outgoing", font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_mouse_all, pattern = "incoming", font.size = 5)

png('./figures/cellchat/outgoing_signaling_role.png', width = 4000, height = 2000, res = 300)
ht1
dev.off()

#####################################
## Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
selectK(cellchat_mouse_all, pattern = "outgoing")

nPatterns = 4
cellchat_mouse_all <- identifyCommunicationPatterns(cellchat_mouse_all, pattern = "outgoing", k = nPatterns)

# river plot
png('./figures/cellchat/outgoing_signaling_pattern.png', width = 4000, height = 2000, res = 300)
netAnalysis_river(cellchat_mouse_all, pattern = "outgoing")
dev.off()

## Identify and visualize incoming communication pattern of secreting cells
selectK(cellchat_mouse_all, pattern = "incoming")

nPatterns = 4
cellchat_mouse_all <- identifyCommunicationPatterns(cellchat_mouse_all, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat_mouse_all, pattern = "incoming")

###############################################
## Manifold and classification learning analysis of signaling networks
# Identify signaling groups based on their functional similarity
cellchat_mouse_all <- computeNetSimilarity(cellchat_mouse_all, type = "functional")
cellchat_mouse_all <- netEmbedding(cellchat_mouse_all, type = "functional")
cellchat_mouse_all <- netClustering(cellchat_mouse_all, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat_mouse_all, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
cellchat_mouse_all <- computeNetSimilarity(cellchat_mouse_all, type = "structural")
cellchat_mouse_all <- netEmbedding(cellchat_mouse_all, type = "structural")
cellchat_mouse_all <- netClustering(cellchat_mouse_all, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat_mouse_all, type = "structural", label.size = 3.5)


netVisual_embeddingZoomIn(cellchat_mouse_all, type = "structural", nCol = 2)



############
# save
saveRDS(cellchat_mouse_all, file = "./forCellChat/cellchat_mouse_all_LS.rds")
cellchat_mouse_all <- readRDS('./forCellChat/cellchat_mouse_all_LS.rds')
