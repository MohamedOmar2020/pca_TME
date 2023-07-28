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
                row_names_rot = 0, row_names_gp = gpar(fontsize = 5),   
                column_names_gp = gpar(fontsize = 8), width = unit(width, 
                                                                   "cm"), height = unit(height, "cm"), column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 45, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10, 
                                                            fontface = "plain"), title_position = "leftcenter-rot", 
                                            border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                 "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(4, 
                                                                                                                                                          "mm")))
  return(ht1)
}


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
                         self.link = FALSE,
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
                                  title = "cell groups", legend_height = unit(12, "cm"), labels_gp = grid::gpar(fontsize=12), title_gp = grid::gpar(fontsize = 14))
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

# modify compareInteractions function
compareInteractions <- function (object, measure = c("count", "weight"), color.use = NULL, 
                                 group = NULL, group.levels = NULL, group.facet = NULL, group.facet.levels = NULL, 
                                 n.row = 1, color.alpha = 1, legend.title = NULL, width = 0.6, 
                                 title.name = NULL, digits = 3, xlabel = NULL, ylabel = NULL, 
                                 remove.xtick = FALSE, show.legend = TRUE, x.lab.rot = FALSE, 
                                 angle.x = 45, vjust.x = NULL, hjust.x = 1, size.text = 10, 
                                 sources.use = NULL, targets.use = NULL) 
{
  # Add validation for sources and targets
  if (!is.null(sources.use) && !all(sources.use %in% object@idents$joint))
    stop("Provided sources are not in the object@idents")
  if (!is.null(targets.use) && !all(targets.use %in% object@idents$joint))
    stop("Provided targets are not in the object@idents")
  
  measure <- match.arg(measure)
  if (measure == "count") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$count[sources.use, targets.use])))
    if (is.null(ylabel)) {
      ylabel = "Number of inferred interactions"
    }
  } else if (measure == "weight") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$weight[sources.use, targets.use])))
    df[,1] <- round(df[,1],digits)
    if (is.null(ylabel)) {
      ylabel = "Interaction strength"
    }
  }
  colnames(df) <- "count"
  
  df$dataset <- names(object@net)
  if (is.null(group)) {
    group <- 1
  }
  df$group <- group
  df$dataset <- factor(df$dataset, levels = names(object@net))
  if (is.null(group.levels)) {
    df$group <- factor(df$group)
  }
  else {
    df$group <- factor(df$group, levels = group.levels)
  }
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(group)))
  }
  if (!is.null(group.facet)) {
    if (all(group.facet %in% colnames(df))) {
      gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
        geom_bar(stat = "identity", width = width, position = position_dodge())
      gg <- gg + facet_wrap(group.facet, nrow = n.row)
    }
    else {
      df$group.facet <- group.facet
      if (is.null(group.facet.levels)) {
        df$group.facet <- factor(df$group.facet)
      }
      else {
        df$group.facet <- factor(df$group.facet, levels = group.facet.levels)
      }
      gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
        geom_bar(stat = "identity", width = width, position = position_dodge())
      gg <- gg + facet_wrap(~group.facet, nrow = n.row)
    }
  }
  else {
    gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
      geom_bar(stat = "identity", width = width, position = position_dodge())
  }
  gg <- gg + geom_text(aes(label = count), vjust = -0.3, size = 3, 
                       position = position_dodge(0.9))
  gg <- gg + ylab(ylabel) + xlab(xlabel) + theme_classic() + 
    labs(title = title.name) + theme(plot.title = element_text(size = 10, 
                                                               face = "bold", hjust = 0.5)) + theme(text = element_text(size = size.text), 
                                                                                                    axis.text = element_text(colour = "black"))
  gg <- gg + scale_fill_manual(values = alpha(color.use, alpha = color.alpha), 
                               drop = FALSE)
  if (remove.xtick) {
    gg <- gg + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  }
  else {
    gg <- gg + guides(fill = guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = angle.x, 
                                                hjust = hjust.x, vjust = vjust.x, size = size.text))
  }
  gg
  return(gg)
}


###################################
# set the active conda environment
use_condaenv(condaenv = "scutils", required = TRUE)

# load the anndata module
ad <- import("anndata", convert = FALSE)

# load the human h5ad object
adata_human_all <- ad$read_h5ad('forCellChat/human_all_annot_raw.h5ad')
# access normalized data matrix
human_all_dataInput <- t(py_to_r(adata_human_all$X))
range(human_all_dataInput)
rownames(human_all_dataInput) <- rownames(py_to_r(adata_human_all$var))
colnames(human_all_dataInput) <- rownames(py_to_r(adata_human_all$obs))

## access meta data
human_all_metaData <- py_to_r(adata_human_all$obs)
human_all_meta <- human_all_metaData
# some cleaning
table(human_all_meta$cluster)
human_all_meta <- human_all_meta[!is.na(human_all_meta$cluster), ]
#levels(human_all_meta$cluster) <- c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages")
human_all_dataInput <- human_all_dataInput[, rownames(human_all_meta)]

############################################
## separate into ERG+ and ERG-
table(human_all_meta$erg)

# ERG pos
human_all_meta_ERGpos <- human_all_meta[human_all_meta$erg == 'positive', ]
human_all_dataInput_ERGpos <- human_all_dataInput[, rownames(human_all_meta_ERGpos)]
human_all_meta_ERGpos$cluster <- droplevels(human_all_meta_ERGpos$cluster)

# ERG neg
human_all_meta_ERGneg <- human_all_meta[human_all_meta$erg == 'negative', ]
human_all_dataInput_ERGneg <- human_all_dataInput[, rownames(human_all_meta_ERGneg)]
human_all_meta_ERGneg$cluster <- droplevels(human_all_meta_ERGneg$cluster)

#############################################
# creat a cellchat object
cellchat_human_all <- createCellChat(object = human_all_dataInput, meta = human_all_meta, group.by = "cluster")
cellchat_human_all_ERGpos <- createCellChat(object = human_all_dataInput_ERGpos, meta = human_all_meta_ERGpos, group.by = "cluster")
cellchat_human_all_ERGneg <- createCellChat(object = human_all_dataInput_ERGneg, meta = human_all_meta_ERGneg, group.by = "cluster")

# Add cell information into meta slot of the object (Optional)
cellchat_human_all <- addMeta(cellchat_human_all, meta = human_all_meta)
cellchat_human_all <- setIdent(cellchat_human_all, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_human_all@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_human_all@idents)) # number of cells in each cell group

cellchat_human_all_ERGpos <- addMeta(cellchat_human_all_ERGpos, meta = human_all_meta_ERGpos)
cellchat_human_all_ERGpos <- setIdent(cellchat_human_all_ERGpos, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_human_all_ERGpos@idents) # show factor levels of the cell labels
groupSize_ERGpos <- as.numeric(table(cellchat_human_all_ERGpos@idents)) # number of cells in each cell group

cellchat_human_all_ERGneg <- addMeta(cellchat_human_all_ERGneg, meta = human_all_meta_ERGneg)
cellchat_human_all_ERGneg <- setIdent(cellchat_human_all_ERGneg, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_human_all_ERGneg@idents) # show factor levels of the cell labels
groupSize_ERGneg <- as.numeric(table(cellchat_human_all_ERGneg@idents)) # number of cells in each cell group


# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB

cellchat_human_all@DB <- CellChatDB.use
cellchat_human_all_ERGpos@DB <- CellChatDB.use
cellchat_human_all_ERGneg@DB <- CellChatDB.use

######################################################
# clean
rm(human_all_dataInput, human_all_dataInput_ERGpos, human_all_dataInput_ERGneg)

#################################################
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat_human_all <- subsetData(cellchat_human_all) # This step is necessary even if using the whole database
cellchat_human_all_ERGpos <- subsetData(cellchat_human_all_ERGpos) # This step is necessary even if using the whole database
cellchat_human_all_ERGneg <- subsetData(cellchat_human_all_ERGneg) # This step is necessary even if using the whole database

future::plan("multicore", workers = 8) # do parallel
cellchat_human_all <- identifyOverExpressedGenes(cellchat_human_all)
cellchat_human_all_ERGpos <- identifyOverExpressedGenes(cellchat_human_all_ERGpos)
cellchat_human_all_ERGneg <- identifyOverExpressedGenes(cellchat_human_all_ERGneg)

cellchat_human_all <- identifyOverExpressedInteractions(cellchat_human_all)
cellchat_human_all_ERGpos <- identifyOverExpressedInteractions(cellchat_human_all_ERGpos)
cellchat_human_all_ERGneg <- identifyOverExpressedInteractions(cellchat_human_all_ERGneg)

# project gene expression data onto PPI network (optional)
cellchat_human_all <- projectData(cellchat_human_all, PPI.human)
cellchat_human_all_ERGpos <- projectData(cellchat_human_all_ERGpos, PPI.human)
cellchat_human_all_ERGneg <- projectData(cellchat_human_all_ERGneg, PPI.human)

#####################################################
## Inference of cell-cell communication network
options(future.globals.maxSize=20*1024^3)


# Compute the communication probability and infer cellular communication network
cellchat_human_all <- computeCommunProb(cellchat_human_all, population.size = TRUE)
cellchat_human_all_ERGpos <- computeCommunProb(cellchat_human_all_ERGpos, population.size = TRUE)
cellchat_human_all_ERGneg <- computeCommunProb(cellchat_human_all_ERGneg, population.size = TRUE)


# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_human_all <- filterCommunication(cellchat_human_all, min.cells = 10)
cellchat_human_all_ERGpos <- filterCommunication(cellchat_human_all_ERGpos, min.cells = 10)
cellchat_human_all_ERGneg <- filterCommunication(cellchat_human_all_ERGneg, min.cells = 10)

######################
# Infer the cell-cell communication at a signaling pathway level
cellchat_human_all <- computeCommunProbPathway(cellchat_human_all)
cellchat_human_all_ERGpos <- computeCommunProbPathway(cellchat_human_all_ERGpos)
cellchat_human_all_ERGneg <- computeCommunProbPathway(cellchat_human_all_ERGneg)

# Calculate the aggregated cell-cell communication network
cellchat_human_all <- aggregateNet(cellchat_human_all)
cellchat_human_all_ERGpos <- aggregateNet(cellchat_human_all_ERGpos)
cellchat_human_all_ERGneg <- aggregateNet(cellchat_human_all_ERGneg)

## Get the dataframe of communications at the L/R level
df.net_human <- subsetCommunication(cellchat_human_all)
df.net_human_ERGpos <- subsetCommunication(cellchat_human_all_ERGpos)
df.net_human_ERGneg <- subsetCommunication(cellchat_human_all_ERGneg)

## Get the dataframe of communications at the pathway level
df.net_human_pathway <- subsetCommunication(cellchat_human_all, slot.name = "netP")
df.net_human_pathway_ERGpos <- subsetCommunication(cellchat_human_all_ERGpos, slot.name = "netP")
df.net_human_pathway_ERGneg <- subsetCommunication(cellchat_human_all_ERGneg, slot.name = "netP")

# all the pathways
table(df.net_human_pathway$pathway_name)

save(cellchat_human_all, df.net_human, df.net_human_pathway, 
     cellchat_human_all_ERGpos, df.net_human_ERGpos, df.net_human_pathway_ERGpos, 
     cellchat_human_all_ERGneg, df.net_human_ERGneg, df.net_human_pathway_ERGneg, 
     file = 'forCellChat/cellchat_human_all.rda')

write.csv(df.net_human, file = 'forCellChat/df.net_human.csv')

#####################################################
#####################################################
#####################################################

load('forCellChat/cellchat_human_all.rda')

#####################################################
#  visualize the aggregated cell-cell communication network:  the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat_human_all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_human_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_human_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

###########
# examine the signaling sent from each cell group.
mat_human <- cellchat_human_all@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat_human)) {
  mat2 <- matrix(0, nrow = nrow(mat_human), ncol = ncol(mat_human), dimnames = dimnames(mat_human))
  mat2[i, ] <- mat_human[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_human), title.name = rownames(mat_human)[i])
}

##########################################################
levels(cellchat_human_all@idents)


################################################################
# Chord diagram
# png('./figures/cellchat_allCompartments/WNT_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'WNT', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/PERIOSTIN_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'PERIOSTIN', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/TGFb_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'TGFb', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/ncWNT_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'ncWNT', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/PDGF_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'PDGF', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/PTN_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'PTN', layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/COLLAGEN_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'COLLAGEN', sources.use= c('c5', 'c6', 'c7'), targets.use = 'epithelium', title.space	= 1, layout = "chord")
# dev.off()
# 
# png('./figures/cellchat_allCompartments/TGFb_chrord_heatmap.png', width = 2000, height = 2000, res = 300)
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_human_all, signaling = 'TGFb', targets.use= c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), title.space	= 1, layout = "chord")
# dev.off()


#########################################################
## signaling from PRN clusters to epithelium

# # pathway level
# png('./figures/cellchat_allCompartments/PRN_epith_pathways.png', width = 3500, height = 3000, res = 350)
# netVisual_chord_gene(cellchat_human_all, 
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
# netVisual_chord_gene(cellchat_human_all, 
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
# netVisual_chord_gene(cellchat_human_all, 
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
# netVisual_chord_gene(cellchat_human_all, 
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
# netVisual_chord_gene(cellchat_human_all, 
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
# netVisual_chord_gene(cellchat_human_all, 
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
# netVisual_chord_gene(cellchat_human_all, 
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
# netVisual_chord_gene(cellchat_human_all, 
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
levels(cellchat_human_all@idents)
colMap <- c('powderblue', 'plum2', 'violetred4', 'slategray4', 'blue', 'darkorange1', 'forestgreen', 'red2', 'darkorchid', 'saddlebrown', 'hotpink1', 'olivedrab3', 'lemonchiffon4', 'bisque2', 'salmon3', 'gold3', 'violetred')


##############################################################
# common clusters
##############################################################

## from epithelium

# all 

# Now plot the chord diagram

png('./figures/cellchat_allCompartments/human/epithelium_c0c1c2.png', width = 3500, height = 3500, res = 350)
CellChat::netVisual_chord_gene(cellchat_human_all, 
                     slot.name = "net",
                     sources.use = c('luminal', 'basal'), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 0.7, 
                     legend.pos.y = 190, 
                     legend.pos.x = 25,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 5,
                     color.use = colMap,
                     #title.name = 'Signling networks from common clusters to the epithelium'	
)
dev.off()

# ERG pos
png('./figures/cellchat_allCompartments/human/epithelium_c0c1c2_ERGpos.png', width = 3500, height = 3500, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGpos, 
                     slot.name = "net",
                     sources.use = c('luminal', 'basal'), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 0.7, 
                     legend.pos.y = 190, 
                     legend.pos.x = 25,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 5,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to the epithelium'	
)
dev.off()

# ERG neg
png('./figures/cellchat_allCompartments/human/epithelium_c0c1c2_ERGneg.png', width = 3500, height = 3500, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGneg, 
                     slot.name = "net",
                     sources.use = c('luminal', 'basal'), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 0.7, 
                     legend.pos.y = 190, 
                     legend.pos.x = 25,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 5,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to the epithelium'	
)
dev.off()

############################
## from immune

# all
png('./figures/cellchat_allCompartments/human/immune_c0c1c2.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all, 
                     slot.name = "net",
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 0.7, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 1,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

# ERG pos
png('./figures/cellchat_allCompartments/human/immune_c0c1c2_ERGpos.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGpos, 
                     slot.name = "net",
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 0.7, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 1,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

# ERG neg
png('./figures/cellchat_allCompartments/human/immune_c0c1c2_ERGneg.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGneg, 
                     slot.name = "net",
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 0.7, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 1,
                     self.link = FALSE,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

##############################################################
##########
# to c3 and c4
##############################################################

############################
## from immune

# all
png('./figures/cellchat_allCompartments/human/immune_c3c4.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all, 
                     slot.name = "net",
                     #signaling = c('COMPLEMENT', 'CCL'),
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signaling Networks from c3 and c4 to immune cells'	
)
dev.off()

# ERG pos
png('./figures/cellchat_allCompartments/human/immune_c3c4_ERGpos.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGpos, 
                     slot.name = "net",
                     #signaling = c('COMPLEMENT', 'CCL'),
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signaling Networks from c3 and c4 to immune cells'	
)
dev.off()

# ERG neg
png('./figures/cellchat_allCompartments/human/immune_c3c4_ERGneg.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGneg, 
                     slot.name = "net",
                     #signaling = c('COMPLEMENT', 'CCL'),
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 175, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signaling Networks from c3 and c4 to immune cells'	
)
dev.off()

############################
## from epithelium

# all
png('./figures/cellchat_allCompartments/human/epith_c3c4.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all, 
                     slot = 'net', 
                     targets.use= c('c3', 'c4'), 
                     sources.use = 'epithelium', 
                     lab.cex = 1, 
                     legend.pos.y = 165, 
                     legend.pos.x = 25, 
                     thresh = 0.05, small.gap = 3,
                     color.use = colMap,
                     #title.name = 'Signaling pathways from c3 and c4 to the epithelium',	
                     #reduce = -1,
)
dev.off()

# ERG pos
png('./figures/cellchat_allCompartments/human/epith_c3c4_ERGpos.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGpos, 
                     slot = 'net', 
                     targets.use= c('c3', 'c4'), 
                     sources.use = 'epithelium', 
                     lab.cex = 1, 
                     legend.pos.y = 165, 
                     legend.pos.x = 25, 
                     thresh = 0.05, small.gap = 3,
                     color.use = colMap,
                     #title.name = 'Signaling pathways from c3 and c4 to the epithelium',	
                     #reduce = -1,
)
dev.off()

# ERG neg
png('./figures/cellchat_allCompartments/human/epith_c3c4_ERGneg.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_human_all_ERGneg, 
                     slot = 'net', 
                     targets.use= c('c3', 'c4'), 
                     sources.use = 'epithelium', 
                     lab.cex = 1, 
                     legend.pos.y = 165, 
                     legend.pos.x = 25, 
                     thresh = 0.05, small.gap = 3,
                     color.use = colMap,
                     #title.name = 'Signaling pathways from c3 and c4 to the epithelium',	
                     #reduce = -1,
)
dev.off()
##########################################################################################################################
##########################################################################################################################

#############################################################
# comparison
##############################################################

##########################################
# compute centrality
##########################################
cellchat_human_all_ERGpos <- netAnalysis_computeCentrality(cellchat_human_all_ERGpos, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_human_all_ERGneg <- netAnalysis_computeCentrality(cellchat_human_all_ERGneg, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#####################################
#Lift up CellChat object and merge together
#####################################


# Define the cell labels to lift up
group.new = levels(cellchat_human_all@idents)

cellchat_human_all_ERGpos <- liftCellChat(cellchat_human_all_ERGpos, group.new)
cellchat_human_all_ERGneg <- liftCellChat(cellchat_human_all_ERGneg, group.new)


object.list <- list(ERGpos = cellchat_human_all_ERGpos, 
                               ERGneg = cellchat_human_all_ERGneg
)
cellchat_ERGstatus <- mergeCellChat(object.list, add.names = names(object.list))


##########################################################################################################
# Compare the total number of interactions and interaction strength
##########################################################################################################

# from epithelium
gg_number_fromEpithelium <- compareInteractions(cellchat_ERGstatus, show.legend = F, group = c(1,2), sources.use = 'epithelium', targets.use = c('c3', 'c4'), title.name = 'N of interactions')
gg_weight_fromEpithelium <- compareInteractions(cellchat_ERGstatus, show.legend = F, group = c(1,2), measure = "weight", sources.use = 'epithelium', targets.use = c('c3', 'c4'), digits = 4, title.name = 'Weight of interactions')

png('./figures/cellchat_allCompartments/human/Diff_fromEpith.png', width = 3500, height = 3000, res = 350)
gg_number_fromEpithelium + gg_weight_fromEpithelium
dev.off()

# from immune 
gg_number_fromImmune <- compareInteractions(cellchat_ERGstatus, show.legend = F, group = c(1,2), sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c("c3", 'c4'))
gg_weight_fromImmune <- compareInteractions(cellchat_ERGstatus, show.legend = F, group = c(1,2), measure = "weight", sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c("c3", 'c4'), digits = 4)
png('./figures/cellchat_allCompartments/human/Diff_fromImmune.png', width = 3500, height = 3000, res = 350)
gg_number_fromImmune + gg_weight_fromImmune
dev.off()

#############################################
# heatmap of differential number of interactions and strength
#############################################
gg_heatmap <- netVisual_heatmap(cellchat_ERGstatus)
gg_heatmap_weight <- netVisual_heatmap(cellchat_ERGstatus, measure = "weight")
gg_heatmap + gg_heatmap_weight

#############################################
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

## Number of interactions between ERG pos and neg
png('./figures/cellchat_allCompartments/human/Diff_N_interactions.png', width = 3500, height = 3000, res = 350)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, 
                   sources.use = c('epithelium', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"),
                   targets.use = c('c3', 'c4'),
                   remove.isolate = T,
                   weight.scale = T, 
                   label.edge= F, 
                   edge.weight.max = weight.max[2], 
                   #edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", 
                                       names(object.list)[i])
  )
}
dev.off()


# Compare the overall information flow of each signaling pathway
gg_overall <- rankNet(cellchat_ERGstatus, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('cyan3', 'indianred1'), font.size = 6)

png('./figures/cellchat_allCompartments/human/Diff_informationFlow.png', width = 1500, height = 2500, res = 300)
gg_overall 
dev.off()


#####################################################
# Identify dysfunctional signaling by using differential expression analysis

# perform differential expression analysis
cellchat_ERGstatus <- identifyOverExpressedGenes(cellchat_ERGstatus, group.dataset = "datasets", pos.dataset = 'ERGpos', features.name = 'ERGpos', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_ERGstatus, features.name = 'ERGpos')

# extract the ligand-receptor pairs with upregulated ligands in mutant
net.up <- subsetCommunication(cellchat_ERGstatus, net = net, datasets = "ERGpos", ligand.logFC = 0.1, receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in wildtypes, i.e.,downregulated in mutants
net.down <- subsetCommunication(cellchat_ERGstatus, net = net, datasets = "ERGpos", ligand.logFC = -0.1, receptor.logFC = -0.1)

# do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_ERGstatus)

gene.down <- extractGeneSubsetFromPair(net.down, cellchat_ERGstatus)

# visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up = net.up[, "interaction_name", drop = F]

pairLR.use.down = net.down[, "interaction_name", drop = F]

######################################################################
# from epith to common clusters
######################################################################
colMap <- c('powderblue', 'plum2', 'violetred4', 'slategray4', 'blue', 'darkorange1', 'forestgreen', 'red2', 'darkorchid', 'saddlebrown', 'hotpink1', 'olivedrab3', 'lemonchiffon4', 'CadetBlue3', 'salmon3', 'gold3', 'violetred')

png('./figures/cellchat_allCompartments/human/up_epith_commonClusters.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list[[1]], 
                     slot.name = 'net', 
                     net = net.up, 
                     sources.use = c('basal', 'luminal'), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     thresh = 0.05, 
                     lab.cex = 0.7,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 25, 
                     color.use = colMap,
                     title.name = "")
dev.off()

######################################################################
# from epith to c3,c4 
######################################################################
png('./figures/cellchat_allCompartments/human/up_epith_c3c4.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list[[1]], 
                     slot.name = 'net', 
                     net = net.up, 
                     sources.use = c('basal', 'luminal'), 
                     targets.use = c('c3', 'c4'), 
                     thresh = 0.05, 
                     lab.cex = 0.7,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 25, 
                     color.use = colMap,
                     title.name = "")
dev.off()



###########################################################################################

######################################################################
# from immune to c0,c1,c2
######################################################################
png('./figures/cellchat_allCompartments/human/up_immune_CommonClusters.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list[[1]], 
                     slot.name = 'net', 
                     net = net.up, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     thresh = 0.05, 
                     lab.cex = 0.7,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 175, 
                     legend.pos.x = 5, 
                     color.use = colMap,
                     title.name = "")
dev.off()

######################################################################
# from immune to c3,c4 in PN, HiMYC, Terg vs wt
######################################################################
png('./figures/cellchat_allCompartments/human/up_immune_c3c4.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list[[1]], 
                     slot.name = 'net', 
                     net = net.up, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     thresh = 0.05, 
                     lab.cex = 0.7,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 175, 
                     legend.pos.x = 5, 
                     color.use = colMap,
                     title.name = "")
dev.off()

