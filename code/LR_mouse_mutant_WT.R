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

netAnalysis_signalingRole_heatmap <- function (object, signaling = NULL, pattern = c("outgoing", "incoming", 
                                                "all"), slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
          title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
          cluster.rows = FALSE, cluster.cols = FALSE) 
{
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- title
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                             name = color.heatmap))))(100)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), ylim = c(0, 0.2),   
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = "Relative strength", bottom_annotation = col_annotation, 
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = 10),   
                column_names_gp = gpar(fontsize = 15), width = unit(width, 
                                                                           "cm"), height = unit(height, "cm"), column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 45, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10, 
                                                            fontface = "plain"), title_position = "leftcenter-rot", 
                                            border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                 "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(4, 
                                                                                                                                                          "mm")))
  return(ht1)
}


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
                                  title = "TME compartments", legend_height = unit(12, "cm"), labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize = 14))
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

##################################################################


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

# Add cell information into meta slot of the object
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
cellchat_mouse_all_wt <- computeCommunProb(cellchat_mouse_all_wt, population.size = TRUE)
cellchat_mouse_all_mutant <- computeCommunProb(cellchat_mouse_all_mutant, population.size = TRUE)

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

#####################################
## Merge the 2 objects

object.list <- list(WT = cellchat_mouse_all_wt, mutant = cellchat_mouse_all_mutant)

# change the names of different compartments


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

####################
# pathway level
png('./figures/cellchat_mutant_wt/new/up_mutants_pathways.png', width = 3500, height = 3000, res = 370)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', net = net.up, 
                     color.use = c('bisque2', 'darkblue', 'darkgreen'),
                     thresh = 0.05, lab.cex = 1, small.gap = 4,
                     legend.pos.y = 180, 
                     legend.pos.x = 20, 
                     link.visible	= TRUE,
                     title.name = "")
dev.off()

png('./figures/cellchat_mutant_wt/new/up_WT_pathways.png', width = 3500, height = 3000, res = 370)
netVisual_chord_gene(object.list[[1]], slot.name = 'netP', net = net.down, 
                     thresh = 0.05, lab.cex = 1, small.gap = 4, 
                     color.use = c('bisque2', 'darkblue', 'darkgreen'),
                     legend.pos.y = 180, 
                     legend.pos.x = 20, 
                     link.visible = T,
                     title.name = "")
dev.off()

######################################################
# change the names of different compartments
#levels(object.list$WT@idents) <- c('E', 'I', 'S')
#levels(object.list$mutant@idents) <- c('E', 'I', 'S')

# Compare outgoing (or incoming) signaling associated with each cell population
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = 'Wild types', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14, color.use = c('bisque2', 'darkblue', 'darkgreen'))
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = 'Mutants', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14, color.use = c('bisque2', 'darkblue', 'darkgreen'))


tiff('./figures/cellchat_mutant_wt/new/Diff_heatmap2.tiff', width = 4000, height = 3000, res = 370)
draw(ht2 + ht1, ht_gap = unit(2, "cm"))
dev.off()

############
# save
saveRDS(cellchat, file = "./forCellChat/cellchat_mouse_WT_vs_mutants.rds")

# load
cellchat = readRDS("./forCellChat/cellchat_mouse_WT_vs_mutants.rds")




