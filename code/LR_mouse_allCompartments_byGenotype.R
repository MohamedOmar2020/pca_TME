# clean workspace
rm(list = ls())

# load libraries
library(data.table)
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(reticulate)
library(xlsx)

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

##############################

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

# load the mouse h5ad object
adata_mouse_all <- ad$read_h5ad("./forCellChat/mouse_all_raw2.h5ad")
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


############################################
## separate into genotypes
table(mouse_all_meta$key_new)

# PN
mouse_all_meta_PN <- mouse_all_meta[mouse_all_meta$key_new == 'NP', ]
mouse_all_dataInput_PN <- mouse_all_dataInput[, rownames(mouse_all_meta_PN)]

# HiMYC
mouse_all_meta_HiMYC <- mouse_all_meta[mouse_all_meta$key_new == 'Hi-MYC', ]
mouse_all_dataInput_HiMYC <- mouse_all_dataInput[, rownames(mouse_all_meta_HiMYC)]

# T-ERG
mouse_all_meta_Terg <- mouse_all_meta[mouse_all_meta$key_new == 'T-ERG', ]
mouse_all_dataInput_Terg <- mouse_all_dataInput[, rownames(mouse_all_meta_Terg)]

# PRN
mouse_all_meta_PRN <- mouse_all_meta[mouse_all_meta$key_new == 'PRN', ]
mouse_all_dataInput_PRN <- mouse_all_dataInput[, rownames(mouse_all_meta_PRN)]

# FVBN
mouse_all_meta_FVBN <- mouse_all_meta[mouse_all_meta$key_new == 'FVBN', ]
mouse_all_dataInput_FVBN <- mouse_all_dataInput[, rownames(mouse_all_meta_FVBN)]

# B6
mouse_all_meta_B6 <- mouse_all_meta[mouse_all_meta$key_new == 'B6', ]
mouse_all_dataInput_B6 <- mouse_all_dataInput[, rownames(mouse_all_meta_B6)]

# B6.129
mouse_all_meta_B6_129 <- mouse_all_meta[mouse_all_meta$key_new == 'B6.129', ]
mouse_all_dataInput_B6_129 <- mouse_all_dataInput[, rownames(mouse_all_meta_B6_129)]

# PRNwt
mouse_all_meta_PRN_wt <- mouse_all_meta[mouse_all_meta$key_new == 'WT for PRN', ]
mouse_all_dataInput_PRN_wt <- mouse_all_dataInput[, rownames(mouse_all_meta_PRN_wt)]

# WT for PN
mouse_all_meta_PN_wt <- mouse_all_meta[mouse_all_meta$key_new == 'WT for PN', ]
mouse_all_dataInput_PN_wt <- mouse_all_dataInput[, rownames(mouse_all_meta_PN_wt)]

# all wt for c3,c4
mouse_all_meta_c3c4_wt <- mouse_all_meta[mouse_all_meta$key_new %in% c('B6', 'B6.129', 'WT for PN', 'FVBN'), ]
mouse_all_dataInput_c3c4_wt <- mouse_all_dataInput[, rownames(mouse_all_meta_c3c4_wt)]

# all c3,c4
mouse_all_meta_c3c4 <- mouse_all_meta[mouse_all_meta$key_new %in% c('T-ERG', 'Hi-MYC', 'NP'), ]
mouse_all_dataInput_c3c4 <- mouse_all_dataInput[, rownames(mouse_all_meta_c3c4)]

# all mutants
mouse_all_meta_mutants <- mouse_all_meta[mouse_all_meta$condition =='mutant', ]
mouse_all_dataInput_mutants <- mouse_all_dataInput[, rownames(mouse_all_meta_mutants)]

# all wildtype
mouse_all_meta_WT <- mouse_all_meta[mouse_all_meta$condition =='wildtype', ]
mouse_all_dataInput_WT <- mouse_all_dataInput[, rownames(mouse_all_meta_WT)]

#############################################
# creat cellchat objects
#############################################
cellchat_mouse_all_PN <- createCellChat(object = mouse_all_dataInput_PN, meta = mouse_all_meta_PN, group.by = "cluster")
cellchat_mouse_all_HiMYC <- createCellChat(object = mouse_all_dataInput_HiMYC, meta = mouse_all_meta_HiMYC, group.by = "cluster")
cellchat_mouse_all_Terg <- createCellChat(object = mouse_all_dataInput_Terg, meta = mouse_all_meta_Terg, group.by = "cluster")
cellchat_mouse_all_PRN <- createCellChat(object = mouse_all_dataInput_PRN, meta = mouse_all_meta_PRN, group.by = "cluster")
cellchat_mouse_all_FVBN <- createCellChat(object = mouse_all_dataInput_FVBN, meta = mouse_all_meta_FVBN, group.by = "cluster")
cellchat_mouse_all_B6 <- createCellChat(object = mouse_all_dataInput_B6, meta = mouse_all_meta_B6, group.by = "cluster")
cellchat_mouse_all_B6_129 <- createCellChat(object = mouse_all_dataInput_B6_129, meta = mouse_all_meta_B6_129, group.by = "cluster")
cellchat_mouse_all_PRN_wt <- createCellChat(object = mouse_all_dataInput_PRN_wt, meta = mouse_all_meta_PRN_wt, group.by = "cluster")
cellchat_mouse_all_PN_wt <- createCellChat(object = mouse_all_dataInput_PN_wt, meta = mouse_all_meta_PN_wt, group.by = "cluster")
cellchat_mouse_all_c3c4 <- createCellChat(object = mouse_all_dataInput_c3c4, meta = mouse_all_meta_c3c4, group.by = "cluster")
cellchat_mouse_all_c3c4_wt <- createCellChat(object = mouse_all_dataInput_c3c4_wt, meta = mouse_all_meta_c3c4_wt, group.by = "cluster")
cellchat_mouse_all_mutants <- createCellChat(object = mouse_all_dataInput_mutants, meta = mouse_all_meta_mutants, group.by = "cluster")
cellchat_mouse_all_WT <- createCellChat(object = mouse_all_dataInput_WT, meta = mouse_all_meta_WT, group.by = "cluster")

#############################################
# Add cell information into meta slot of the object
#############################################
mouse_all_meta_PN$cluster <- droplevels(mouse_all_meta_PN$cluster)
cellchat_mouse_all_PN <- addMeta(cellchat_mouse_all_PN, meta = mouse_all_meta_PN)
cellchat_mouse_all_PN <- setIdent(cellchat_mouse_all_PN, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_PN@idents) # show factor levels of the cell labels
groupSize_PN <- as.numeric(table(cellchat_mouse_all_PN@idents)) # number of cells in each cell group

####
mouse_all_meta_HiMYC$cluster <- droplevels(mouse_all_meta_HiMYC$cluster)
cellchat_mouse_all_HiMYC <- addMeta(cellchat_mouse_all_HiMYC, meta = mouse_all_meta_HiMYC)
cellchat_mouse_all_HiMYC <- setIdent(cellchat_mouse_all_HiMYC, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_HiMYC@idents) # show factor levels of the cell labels
groupSize_HiMYC <- as.numeric(table(cellchat_mouse_all_HiMYC@idents)) # number of cells in each cell group

####
mouse_all_meta_Terg$cluster <- droplevels(mouse_all_meta_Terg$cluster)
cellchat_mouse_all_Terg <- addMeta(cellchat_mouse_all_Terg, meta = mouse_all_meta_Terg)
cellchat_mouse_all_Terg <- setIdent(cellchat_mouse_all_Terg, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_Terg@idents) # show factor levels of the cell labels
groupSize_Terg <- as.numeric(table(cellchat_mouse_all_Terg@idents)) # number of cells in each cell group

####
mouse_all_meta_PRN$cluster <- droplevels(mouse_all_meta_PRN$cluster)
cellchat_mouse_all_PRN <- addMeta(cellchat_mouse_all_PRN, meta = mouse_all_meta_PRN)
cellchat_mouse_all_PRN <- setIdent(cellchat_mouse_all_PRN, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_PRN@idents) # show factor levels of the cell labels
groupSize_PRN <- as.numeric(table(cellchat_mouse_all_PRN@idents)) # number of cells in each cell group

####
mouse_all_meta_FVBN$cluster <- droplevels(mouse_all_meta_FVBN$cluster)
cellchat_mouse_all_FVBN <- addMeta(cellchat_mouse_all_FVBN, meta = mouse_all_meta_FVBN)
cellchat_mouse_all_FVBN <- setIdent(cellchat_mouse_all_FVBN, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_FVBN@idents) # show factor levels of the cell labels
groupSize_FVBN <- as.numeric(table(cellchat_mouse_all_FVBN@idents)) # number of cells in each cell group

####
mouse_all_meta_B6$cluster <- droplevels(mouse_all_meta_B6$cluster)
cellchat_mouse_all_B6 <- addMeta(cellchat_mouse_all_B6, meta = mouse_all_meta_B6)
cellchat_mouse_all_B6 <- setIdent(cellchat_mouse_all_B6, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_B6@idents) # show factor levels of the cell labels
groupSize_B6 <- as.numeric(table(cellchat_mouse_all_B6@idents)) # number of cells in each cell group

####
mouse_all_meta_B6_129$cluster <- droplevels(mouse_all_meta_B6_129$cluster)
cellchat_mouse_all_B6_129 <- addMeta(cellchat_mouse_all_B6_129, meta = mouse_all_meta_B6_129)
cellchat_mouse_all_B6_129 <- setIdent(cellchat_mouse_all_B6_129, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_B6_129@idents) # show factor levels of the cell labels
groupSize_B6_129 <- as.numeric(table(cellchat_mouse_all_B6_129@idents)) # number of cells in each cell group

####
mouse_all_meta_PRN_wt$cluster <- droplevels(mouse_all_meta_PRN_wt$cluster)
cellchat_mouse_all_PRN_wt <- addMeta(cellchat_mouse_all_PRN_wt, meta = mouse_all_meta_PRN_wt)
cellchat_mouse_all_PRN_wt <- setIdent(cellchat_mouse_all_PRN_wt, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_PRN_wt@idents) # show factor levels of the cell labels
groupSize_PRN_WT <- as.numeric(table(cellchat_mouse_all_PRN_wt@idents)) # number of cells in each cell group

####
mouse_all_meta_PN_wt$cluster <- droplevels(mouse_all_meta_PN_wt$cluster)
cellchat_mouse_all_PN_wt <- addMeta(cellchat_mouse_all_PN_wt, meta = mouse_all_meta_PN_wt)
cellchat_mouse_all_PN_wt <- setIdent(cellchat_mouse_all_PN_wt, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_PN_wt@idents) # show factor levels of the cell labels
groupSize_PN_wt <- as.numeric(table(cellchat_mouse_all_PN_wt@idents)) # number of cells in each cell group

####
mouse_all_meta_c3c4$cluster <- droplevels(mouse_all_meta_c3c4$cluster)
cellchat_mouse_all_c3c4 <- addMeta(cellchat_mouse_all_c3c4, meta = mouse_all_meta_c3c4)
cellchat_mouse_all_c3c4 <- setIdent(cellchat_mouse_all_c3c4, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_c3c4@idents) # show factor levels of the cell labels
groupSize_c3c4 <- as.numeric(table(cellchat_mouse_all_c3c4@idents)) # number of cells in each cell group

####
mouse_all_meta_c3c4_wt$cluster <- droplevels(mouse_all_meta_c3c4_wt$cluster)
cellchat_mouse_all_c3c4_wt <- addMeta(cellchat_mouse_all_c3c4_wt, meta = mouse_all_meta_c3c4_wt)
cellchat_mouse_all_c3c4_wt <- setIdent(cellchat_mouse_all_c3c4_wt, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_c3c4_wt@idents) # show factor levels of the cell labels
groupSize_c3c4_wt <- as.numeric(table(cellchat_mouse_all_c3c4_wt@idents)) # number of cells in each cell group

####
mouse_all_meta_mutants$cluster <- droplevels(mouse_all_meta_mutants$cluster)
cellchat_mouse_all_mutants <- addMeta(cellchat_mouse_all_mutants, meta = mouse_all_meta_mutants)
cellchat_mouse_all_mutants <- setIdent(cellchat_mouse_all_mutants, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_mutants@idents) # show factor levels of the cell labels
groupSize_mutants <- as.numeric(table(cellchat_mouse_all_mutants@idents)) # number of cells in each cell group

####
mouse_all_meta_WT$cluster <- droplevels(mouse_all_meta_WT$cluster)
cellchat_mouse_all_WT <- addMeta(cellchat_mouse_all_WT, meta = mouse_all_meta_WT)
cellchat_mouse_all_WT <- setIdent(cellchat_mouse_all_WT, ident.use = "cluster") # set "labels" as default cell identity
levels(cellchat_mouse_all_WT@idents) # show factor levels of the cell labels
groupSize_WT <- as.numeric(table(cellchat_mouse_all_WT@idents)) # number of cells in each cell group

######################################################
# clean
######################################################
rm(adata_mouse_all, 
   mouse_all_dataInput, 
   mouse_all_dataInput_B6, 
   mouse_all_dataInput_B6_129, 
   mouse_all_dataInput_FVBN,
   mouse_all_dataInput_HiMYC, 
   mouse_all_dataInput_PN, 
   mouse_all_dataInput_PN_wt, 
   mouse_all_dataInput_PRN, 
   mouse_all_dataInput_PRN_wt,
   mouse_all_dataInput_c3c4,
   mouse_all_dataInput_c3c4_wt,
   mouse_all_dataInput_mutants,
   mouse_all_dataInput_WT
   )



#######################################################################################################################################
# Set the ligand-receptor interaction database
#######################################################################################################################################
CellChatDB <- CellChatDB.mouse
#showDatabaseCategory(CellChatDB)

# Show the structure of the database
#dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB
cellchat_mouse_all_PN@DB <- CellChatDB.use
cellchat_mouse_all_HiMYC@DB <- CellChatDB.use
cellchat_mouse_all_Terg@DB <- CellChatDB.use
cellchat_mouse_all_PRN@DB <- CellChatDB.use
cellchat_mouse_all_FVBN@DB <- CellChatDB.use
cellchat_mouse_all_B6@DB <- CellChatDB.use
cellchat_mouse_all_B6_129@DB <- CellChatDB.use
cellchat_mouse_all_PRN_wt@DB <- CellChatDB.use
cellchat_mouse_all_PN_wt@DB <- CellChatDB.use
cellchat_mouse_all_c3c4@DB <- CellChatDB.use
cellchat_mouse_all_c3c4_wt@DB <- CellChatDB.use
cellchat_mouse_all_mutants@DB <- CellChatDB.use
cellchat_mouse_all_WT@DB <- CellChatDB.use

#################################################
# Preprocessing the expression data for cell-cell communication analysis
#################################################

# subset the expression data of signaling genes for saving computation cost
cellchat_mouse_all_PN <- subsetData(cellchat_mouse_all_PN) # This step is necessary even if using the whole database
future::plan("multicore", workers = 8) # do parallel
cellchat_mouse_all_PN <- identifyOverExpressedGenes(cellchat_mouse_all_PN)
cellchat_mouse_all_PN <- identifyOverExpressedInteractions(cellchat_mouse_all_PN)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_PN <- projectData(cellchat_mouse_all_PN, PPI.mouse)

#####
cellchat_mouse_all_HiMYC <- subsetData(cellchat_mouse_all_HiMYC) # This step is necessary even if using the whole database
cellchat_mouse_all_HiMYC <- identifyOverExpressedGenes(cellchat_mouse_all_HiMYC)
cellchat_mouse_all_HiMYC <- identifyOverExpressedInteractions(cellchat_mouse_all_HiMYC)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_HiMYC <- projectData(cellchat_mouse_all_HiMYC, PPI.mouse)

#####
cellchat_mouse_all_Terg <- subsetData(cellchat_mouse_all_Terg) # This step is necessary even if using the whole database
cellchat_mouse_all_Terg <- identifyOverExpressedGenes(cellchat_mouse_all_Terg)
cellchat_mouse_all_Terg <- identifyOverExpressedInteractions(cellchat_mouse_all_Terg)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_Terg <- projectData(cellchat_mouse_all_Terg, PPI.mouse)

#####
cellchat_mouse_all_PRN <- subsetData(cellchat_mouse_all_PRN) # This step is necessary even if using the whole database
cellchat_mouse_all_PRN <- identifyOverExpressedGenes(cellchat_mouse_all_PRN)
cellchat_mouse_all_PRN <- identifyOverExpressedInteractions(cellchat_mouse_all_PRN)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_PRN <- projectData(cellchat_mouse_all_PRN, PPI.mouse)

#####
cellchat_mouse_all_FVBN <- subsetData(cellchat_mouse_all_FVBN) # This step is necessary even if using the whole database
cellchat_mouse_all_FVBN <- identifyOverExpressedGenes(cellchat_mouse_all_FVBN)
cellchat_mouse_all_FVBN <- identifyOverExpressedInteractions(cellchat_mouse_all_FVBN)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_FVBN <- projectData(cellchat_mouse_all_FVBN, PPI.mouse)

#####
cellchat_mouse_all_B6 <- subsetData(cellchat_mouse_all_B6) # This step is necessary even if using the whole database
cellchat_mouse_all_B6 <- identifyOverExpressedGenes(cellchat_mouse_all_B6)
cellchat_mouse_all_B6 <- identifyOverExpressedInteractions(cellchat_mouse_all_B6)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_B6 <- projectData(cellchat_mouse_all_B6, PPI.mouse)

#####
cellchat_mouse_all_B6_129 <- subsetData(cellchat_mouse_all_B6_129) # This step is necessary even if using the whole database
cellchat_mouse_all_B6_129 <- identifyOverExpressedGenes(cellchat_mouse_all_B6_129)
cellchat_mouse_all_B6_129 <- identifyOverExpressedInteractions(cellchat_mouse_all_B6_129)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_B6_129 <- projectData(cellchat_mouse_all_B6_129, PPI.mouse)

#####
cellchat_mouse_all_PRN_wt <- subsetData(cellchat_mouse_all_PRN_wt) # This step is necessary even if using the whole database
cellchat_mouse_all_PRN_wt <- identifyOverExpressedGenes(cellchat_mouse_all_PRN_wt)
cellchat_mouse_all_PRN_wt <- identifyOverExpressedInteractions(cellchat_mouse_all_PRN_wt)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_PRN_wt <- projectData(cellchat_mouse_all_PRN_wt, PPI.mouse)

#####
cellchat_mouse_all_PN_wt <- subsetData(cellchat_mouse_all_PN_wt) # This step is necessary even if using the whole database
cellchat_mouse_all_PN_wt <- identifyOverExpressedGenes(cellchat_mouse_all_PN_wt)
cellchat_mouse_all_PN_wt <- identifyOverExpressedInteractions(cellchat_mouse_all_PN_wt)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_PN_wt <- projectData(cellchat_mouse_all_PN_wt, PPI.mouse)

#####
cellchat_mouse_all_c3c4 <- subsetData(cellchat_mouse_all_c3c4) # This step is necessary even if using the whole database
cellchat_mouse_all_c3c4 <- identifyOverExpressedGenes(cellchat_mouse_all_c3c4)
cellchat_mouse_all_c3c4 <- identifyOverExpressedInteractions(cellchat_mouse_all_c3c4)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_c3c4 <- projectData(cellchat_mouse_all_c3c4, PPI.mouse)

#####
cellchat_mouse_all_c3c4_wt <- subsetData(cellchat_mouse_all_c3c4_wt) # This step is necessary even if using the whole database
cellchat_mouse_all_c3c4_wt <- identifyOverExpressedGenes(cellchat_mouse_all_c3c4_wt)
cellchat_mouse_all_c3c4_wt <- identifyOverExpressedInteractions(cellchat_mouse_all_c3c4_wt)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_c3c4_wt <- projectData(cellchat_mouse_all_c3c4_wt, PPI.mouse)

#####
cellchat_mouse_all_mutants <- subsetData(cellchat_mouse_all_mutants) # This step is necessary even if using the whole database
cellchat_mouse_all_mutants <- identifyOverExpressedGenes(cellchat_mouse_all_mutants)
cellchat_mouse_all_mutants <- identifyOverExpressedInteractions(cellchat_mouse_all_mutants)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_mutants <- projectData(cellchat_mouse_all_mutants, PPI.mouse)

#####
cellchat_mouse_all_WT <- subsetData(cellchat_mouse_all_WT) # This step is necessary even if using the whole database
cellchat_mouse_all_WT <- identifyOverExpressedGenes(cellchat_mouse_all_WT)
cellchat_mouse_all_WT <- identifyOverExpressedInteractions(cellchat_mouse_all_WT)
# project gene expression data onto PPI network (optional)
cellchat_mouse_all_WT <- projectData(cellchat_mouse_all_WT, PPI.mouse)

#####################################################
## Inference of cell-cell communication network
#####################################################
options(future.globals.maxSize=20*1024^3)

# Compute the communication probability and infer cellular communication network
cellchat_mouse_all_PN <- computeCommunProb(cellchat_mouse_all_PN, population.size = TRUE)
cellchat_mouse_all_HiMYC <- computeCommunProb(cellchat_mouse_all_HiMYC, population.size = TRUE)
cellchat_mouse_all_Terg <- computeCommunProb(cellchat_mouse_all_Terg, population.size = TRUE)
cellchat_mouse_all_PRN <- computeCommunProb(cellchat_mouse_all_PRN, population.size = TRUE)
cellchat_mouse_all_FVBN <- computeCommunProb(cellchat_mouse_all_FVBN, population.size = TRUE)
cellchat_mouse_all_B6 <- computeCommunProb(cellchat_mouse_all_B6, population.size = TRUE)
cellchat_mouse_all_B6_129 <- computeCommunProb(cellchat_mouse_all_B6_129, population.size = TRUE)
cellchat_mouse_all_PRN_wt <- computeCommunProb(cellchat_mouse_all_PRN_wt, population.size = TRUE)
cellchat_mouse_all_PN_wt <- computeCommunProb(cellchat_mouse_all_PN_wt, population.size = TRUE)
cellchat_mouse_all_c3c4 <- computeCommunProb(cellchat_mouse_all_c3c4, population.size = TRUE)
cellchat_mouse_all_c3c4_wt <- computeCommunProb(cellchat_mouse_all_c3c4_wt, population.size = TRUE)
cellchat_mouse_all_mutants <- computeCommunProb(cellchat_mouse_all_mutants, population.size = TRUE)
cellchat_mouse_all_WT <- computeCommunProb(cellchat_mouse_all_WT, population.size = TRUE)



# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_mouse_all_PN <- filterCommunication(cellchat_mouse_all_PN, min.cells = 10)
cellchat_mouse_all_HiMYC <- filterCommunication(cellchat_mouse_all_HiMYC, min.cells = 10)
cellchat_mouse_all_Terg <- filterCommunication(cellchat_mouse_all_Terg, min.cells = 10)
cellchat_mouse_all_PRN <- filterCommunication(cellchat_mouse_all_PRN, min.cells = 10)
cellchat_mouse_all_FVBN <- filterCommunication(cellchat_mouse_all_FVBN, min.cells = 10)
cellchat_mouse_all_B6 <- filterCommunication(cellchat_mouse_all_B6, min.cells = 10)
cellchat_mouse_all_B6_129 <- filterCommunication(cellchat_mouse_all_B6_129, min.cells = 10)
cellchat_mouse_all_PRN_wt <- filterCommunication(cellchat_mouse_all_PRN_wt, min.cells = 10)
cellchat_mouse_all_PN_wt <- filterCommunication(cellchat_mouse_all_PN_wt, min.cells = 10)
cellchat_mouse_all_c3c4 <- filterCommunication(cellchat_mouse_all_c3c4, min.cells = 10)
cellchat_mouse_all_c3c4_wt <- filterCommunication(cellchat_mouse_all_c3c4_wt, min.cells = 10)
cellchat_mouse_all_mutants <- filterCommunication(cellchat_mouse_all_mutants, min.cells = 10)
cellchat_mouse_all_WT <- filterCommunication(cellchat_mouse_all_WT, min.cells = 10)


############################################
# Infer the cell-cell communication at a signaling pathway level
############################################
cellchat_mouse_all_PN <- computeCommunProbPathway(cellchat_mouse_all_PN)
cellchat_mouse_all_HiMYC <- computeCommunProbPathway(cellchat_mouse_all_HiMYC)
cellchat_mouse_all_Terg <- computeCommunProbPathway(cellchat_mouse_all_Terg)
cellchat_mouse_all_PRN <- computeCommunProbPathway(cellchat_mouse_all_PRN)
cellchat_mouse_all_FVBN <- computeCommunProbPathway(cellchat_mouse_all_FVBN)
cellchat_mouse_all_B6 <- computeCommunProbPathway(cellchat_mouse_all_B6)
cellchat_mouse_all_B6_129 <- computeCommunProbPathway(cellchat_mouse_all_B6_129)
cellchat_mouse_all_PRN_wt <- computeCommunProbPathway(cellchat_mouse_all_PRN_wt)
cellchat_mouse_all_PN_wt <- computeCommunProbPathway(cellchat_mouse_all_PN_wt)
cellchat_mouse_all_c3c4 <- computeCommunProbPathway(cellchat_mouse_all_c3c4)
cellchat_mouse_all_c3c4_wt <- computeCommunProbPathway(cellchat_mouse_all_c3c4_wt)
cellchat_mouse_all_mutants <- computeCommunProbPathway(cellchat_mouse_all_mutants)
cellchat_mouse_all_WT <- computeCommunProbPathway(cellchat_mouse_all_WT)

############################################
# Calculate the aggregated cell-cell communication network
############################################
cellchat_mouse_all_PN <- aggregateNet(cellchat_mouse_all_PN)
cellchat_mouse_all_HiMYC <- aggregateNet(cellchat_mouse_all_HiMYC)
cellchat_mouse_all_Terg <- aggregateNet(cellchat_mouse_all_Terg)
cellchat_mouse_all_PRN <- aggregateNet(cellchat_mouse_all_PRN)
cellchat_mouse_all_FVBN <- aggregateNet(cellchat_mouse_all_FVBN)
cellchat_mouse_all_B6 <- aggregateNet(cellchat_mouse_all_B6)
cellchat_mouse_all_B6_129 <- aggregateNet(cellchat_mouse_all_B6_129)
cellchat_mouse_all_PRN_wt <- aggregateNet(cellchat_mouse_all_PRN_wt)
cellchat_mouse_all_PN_wt <- aggregateNet(cellchat_mouse_all_PN_wt)
cellchat_mouse_all_c3c4 <- aggregateNet(cellchat_mouse_all_c3c4)
cellchat_mouse_all_c3c4_wt <- aggregateNet(cellchat_mouse_all_c3c4_wt)
cellchat_mouse_all_mutants <- aggregateNet(cellchat_mouse_all_mutants)
cellchat_mouse_all_WT <- aggregateNet(cellchat_mouse_all_WT)

############################################
# save
############################################
save(cellchat_mouse_all_PN, 
     cellchat_mouse_all_HiMYC, 
     cellchat_mouse_all_Terg, 
     cellchat_mouse_all_PRN,
     cellchat_mouse_all_FVBN,
     cellchat_mouse_all_B6,
     cellchat_mouse_all_B6_129,
     cellchat_mouse_all_PRN_wt,
     cellchat_mouse_all_PN_wt,
     cellchat_mouse_all_c3c4,
     cellchat_mouse_all_c3c4_wt,
     cellchat_mouse_all_mutants,
     cellchat_mouse_all_WT,
     file = './objs/LR_genotypes.rda'
)


load('./objs/LR_genotypes.rda')

############################################
## Get the dataframe of communications at the L/R level
############################################
df.net_mouse_PN <- subsetCommunication(cellchat_mouse_all_PN)
df.net_mouse_HiMyc <- subsetCommunication(cellchat_mouse_all_HiMYC)
df.net_mouse_Terg <- subsetCommunication(cellchat_mouse_all_Terg)
df.net_mouse_PRN <- subsetCommunication(cellchat_mouse_all_PRN)
df.net_mouse_FVBN <- subsetCommunication(cellchat_mouse_all_FVBN)
df.net_mouse_B6 <- subsetCommunication(cellchat_mouse_all_B6)
df.net_mouse_B6_129 <- subsetCommunication(cellchat_mouse_all_B6_129)
df.net_mouse_PRN_wt <- subsetCommunication(cellchat_mouse_all_PRN_wt)
df.net_mouse_PN_wt <- subsetCommunication(cellchat_mouse_all_PN_wt)
df.net_mouse_c3c4 <- subsetCommunication(cellchat_mouse_all_c3c4)
df.net_mouse_c3c4_wt <- subsetCommunication(cellchat_mouse_all_c3c4_wt)
df.net_mouse_mutants <- subsetCommunication(cellchat_mouse_all_mutants)
df.net_mouse_WT <- subsetCommunication(cellchat_mouse_all_WT)

############################################
## Get the dataframe of communications at the pathway level
############################################
df.net_mouse_pathway_PN <- subsetCommunication(cellchat_mouse_all_PN, slot.name = "netP")
df.net_mouse_pathway_HiMYC <- subsetCommunication(cellchat_mouse_all_HiMYC, slot.name = "netP")
df.net_mouse_pathway_Terg <- subsetCommunication(cellchat_mouse_all_Terg, slot.name = "netP")
df.net_mouse_pathway_PRN <- subsetCommunication(cellchat_mouse_all_PRN, slot.name = "netP")
df.net_mouse_pathway_FVBN <- subsetCommunication(cellchat_mouse_all_FVBN, slot.name = "netP")
df.net_mouse_pathway_B6 <- subsetCommunication(cellchat_mouse_all_B6, slot.name = "netP")
df.net_mouse_pathway_B6_129 <- subsetCommunication(cellchat_mouse_all_B6_129, slot.name = "netP")
df.net_mouse_pathway_PRN_wt <- subsetCommunication(cellchat_mouse_all_PRN_wt, slot.name = "netP")
df.net_mouse_pathway_PN_wt <- subsetCommunication(cellchat_mouse_all_PN_wt, slot.name = "netP")
df.net_mouse_pathway_c3c4 <- subsetCommunication(cellchat_mouse_all_c3c4, slot.name = "netP")
df.net_mouse_pathway_c3c4_wt <- subsetCommunication(cellchat_mouse_all_c3c4_wt, slot.name = "netP")
df.net_mouse_pathway_mutants <- subsetCommunication(cellchat_mouse_all_mutants, slot.name = "netP")
df.net_mouse_pathway_WT <- subsetCommunication(cellchat_mouse_all_WT, slot.name = "netP")


##########################################
# compute centrality
##########################################
cellchat_mouse_all_PN <- netAnalysis_computeCentrality(cellchat_mouse_all_PN, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_HiMYC <- netAnalysis_computeCentrality(cellchat_mouse_all_HiMYC, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_Terg <- netAnalysis_computeCentrality(cellchat_mouse_all_Terg, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_PRN <- netAnalysis_computeCentrality(cellchat_mouse_all_PRN, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_FVBN <- netAnalysis_computeCentrality(cellchat_mouse_all_FVBN, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_B6 <- netAnalysis_computeCentrality(cellchat_mouse_all_B6, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_B6_129 <- netAnalysis_computeCentrality(cellchat_mouse_all_B6_129, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_PRN_wt <- netAnalysis_computeCentrality(cellchat_mouse_all_PRN_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_PN_wt <- netAnalysis_computeCentrality(cellchat_mouse_all_PN_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_c3c4 <- netAnalysis_computeCentrality(cellchat_mouse_all_c3c4, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_c3c4_wt <- netAnalysis_computeCentrality(cellchat_mouse_all_c3c4_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_mutants <- netAnalysis_computeCentrality(cellchat_mouse_all_mutants, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_mouse_all_WT <- netAnalysis_computeCentrality(cellchat_mouse_all_WT, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#####################################
#Lift up CellChat object and merge together
#####################################


# Define the cell labels to lift up
group.new = levels(cellchat_mouse_all_mutants@idents)

cellchat_mouse_all_FVBN <- liftCellChat(cellchat_mouse_all_FVBN, group.new)
cellchat_mouse_all_B6 <- liftCellChat(cellchat_mouse_all_B6, group.new)
cellchat_mouse_all_B6_129 <- liftCellChat(cellchat_mouse_all_B6_129, group.new)
cellchat_mouse_all_PN <- liftCellChat(cellchat_mouse_all_PN, group.new)
cellchat_mouse_all_PRN_wt <- liftCellChat(cellchat_mouse_all_PRN_wt, group.new)
cellchat_mouse_all_PRN <- liftCellChat(cellchat_mouse_all_PRN, group.new)
cellchat_mouse_all_HiMYC <- liftCellChat(cellchat_mouse_all_HiMYC, group.new)
cellchat_mouse_all_Terg <- liftCellChat(cellchat_mouse_all_Terg, group.new)
cellchat_mouse_all_WT <- liftCellChat(cellchat_mouse_all_WT, group.new)


#####################################
## Merge the objects

# object.list_all <- list(PN = cellchat_mouse_all_PN, 
#                     HiMYC = cellchat_mouse_all_HiMYC,
#                     Terg = cellchat_mouse_all_Terg,
#                     PRN = cellchat_mouse_all_PRN,
#                     FVBN = cellchat_mouse_all_FVBN,
#                     B6 = cellchat_mouse_all_B6,
#                     B6_129 = cellchat_mouse_all_B6_129,
#                     PRN_wt = cellchat_mouse_all_PRN_wt,
#                     PN_wt = cellchat_mouse_all_PN_wt
#                     )
# 
# 
# cellchat <- mergeCellChat(object.list_all, add.names = names(object.list))

#####################################
## For PN

object.list_NP <- list(NP = cellchat_mouse_all_PN, 
                    NP_wt = cellchat_mouse_all_PN_wt
)


cellchat_NP_with_wt <- mergeCellChat(object.list_NP, add.names = names(object.list_NP))

#####################################
## For Terg

object.list_Terg <- list(Terg = cellchat_mouse_all_Terg, 
                       WT = cellchat_mouse_all_FVBN
)


cellchat_Terg_with_wt <- mergeCellChat(object.list_Terg, add.names = names(object.list_Terg))


#####################################
## For HiMYC

object.list_HiMyc <- list(HiMYC = cellchat_mouse_all_HiMYC, 
                         WT = cellchat_mouse_all_B6
)


cellchat_HiMyc_with_wt <- mergeCellChat(object.list_HiMyc, add.names = names(object.list_HiMyc))

#####################################
## For c3-c4

object.list_c3c4 <- list(c3c4 = cellchat_mouse_all_c3c4, 
                          WT = cellchat_mouse_all_c3c4_wt
)


cellchat_c3c4_with_wt <- mergeCellChat(object.list_c3c4, add.names = c('PN / HiMYC / T-ERG', 'WT'))

#####################################
## For PRN

object.list_PRN <- list(PRN = cellchat_mouse_all_PRN, 
                        WT = cellchat_mouse_all_PRN_wt
)


cellchat_PRN_with_wt <- mergeCellChat(object.list_PRN, add.names = names(object.list_PRN))


#####################################
## PRN and HiMYC

# object.list_PRN_HiMYC <- list(PRN = cellchat_mouse_all_PRN, 
#                         HiMYC = cellchat_mouse_all_HiMYC
# )
# 
# 
# cellchat_PRN_with_HiMYC <- mergeCellChat(object.list_PRN_HiMYC, add.names = names(object.list_PRN_HiMYC))

#####################################
## PRN and PN

# object.list_PRN_PN <- list(PRN = cellchat_mouse_all_PRN, 
#                               HiMYC = cellchat_mouse_all_PN
# )
# 
# 
# cellchat_PRN_with_PN <- mergeCellChat(object.list_PRN_PN, add.names = names(object.list_PRN_PN))

#####################################
## HiMYC and Terg

object.list_HiMYC_Terg <- list(HiMYC = cellchat_mouse_all_HiMYC, 
                             Terg = cellchat_mouse_all_Terg
  )
cellchat_HiMYC_with_Terg <- mergeCellChat(object.list_HiMYC_Terg, add.names = names(object.list_HiMYC_Terg))


object.list_Terg_HiMYC <- list( Terg = cellchat_mouse_all_Terg,
                                HiMYC = cellchat_mouse_all_HiMYC
)
cellchat_Terg_with_HiMyc <- mergeCellChat(object.list_Terg_HiMYC, add.names = names(object.list_Terg_HiMYC))

#####################################
## NP and Terg
object.list_Terg_PN <- list( Terg = cellchat_mouse_all_Terg,
                                PN = cellchat_mouse_all_PN
)
cellchat_Terg_with_PN <- mergeCellChat(object.list_Terg_PN, add.names = names(object.list_Terg_PN))

#####################################
## NP and HiMYC

object.list_NP_HiMYC <- list(NP = cellchat_mouse_all_PN,
  HiMYC = cellchat_mouse_all_HiMYC
)
cellchat_NP_with_HiMYC <- mergeCellChat(object.list_NP_HiMYC, add.names = names(object.list_NP_HiMYC))


#####################################
## PRN and c3c4

object.list_PRN_c3c4 <- list(PRN = cellchat_mouse_all_PRN, 
                             c3c4 = cellchat_mouse_all_c3c4
)


cellchat_PRN_with_c3c4 <- mergeCellChat(object.list_PRN_c3c4, add.names = c('PRN', 'PN / HiMYC / T-ERG'))

#####################################
## all mutants and all WT

object.list_mutants_WT <- list(mutants = cellchat_mouse_all_mutants, 
                             WT = cellchat_mouse_all_WT
)


cellchat_mutants_with_WT <- mergeCellChat(object.list_mutants_WT, add.names = c('mutants', 'WT'))

##########################################################################################################
# Compare the total number of interactions and interaction strength
##########################################################################################################
# NP vs wt

# # from epithelium
# gg_NP_number_fromEpithelium <- compareInteractions(cellchat_NP_with_wt, show.legend = F, group = c(1,2), sources.use = 'epithelium', targets.use = c('c3', 'c4'), title.name = 'N of interactions')
# gg_NP_weight_fromEpithelium <- compareInteractions(cellchat_NP_with_wt, show.legend = F, group = c(1,2), measure = "weight", sources.use = 'epithelium', targets.use = c('c3', 'c4'), digits = 4, title.name = 'Weight of interactions')
# png('./figures/LR_byGenotype/Diff_NP_vs_wt_fromEpith.png', width = 2000, height = 1500, res = 300)
# gg_NP_number_fromEpithelium + gg_NP_weight_fromEpithelium
# dev.off()
# 
# # from immune 
# gg_NP_number_fromImmune <- compareInteractions(cellchat_NP_with_wt, show.legend = F, group = c(1,2), sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"))
# gg_NP_weight_fromImmune <- compareInteractions(cellchat_NP_with_wt, show.legend = F, group = c(1,2), measure = "weight", sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), digits = 4)
# png('./figures/LR_byGenotype/Diff_NP_vs_wt_fromImmune.png', width = 2000, height = 1500, res = 300)
# gg_NP_number_fromImmune + gg_NP_weight_fromImmune
# dev.off()


# all
# gg_NP_number_all <- compareInteractions(cellchat_NP_with_wt, show.legend = F, group = c(1,2), sources.use = 'monocytes/macrophages', targets.use = levels(cellchat_mouse_all_PN@idents))
# gg_NP_weight_all <- compareInteractions(cellchat_NP_with_wt, show.legend = F, group = c(1,2), measure = "weight", digits = 4, sources.use = 'monocytes/macrophages', targets.use = levels(cellchat_mouse_all_PN@idents))
# png('./figures/LR_byGenotype/Diff_NP_vs_wt_all.png', width = 2000, height = 1500, res = 300)
# gg_NP_number_all + gg_NP_weight_all
# dev.off()
# 
# ######################
# # c3c4 vs wt
# 
# # from epithelium
# gg_c3c4_number_fromEpithelium <- compareInteractions(cellchat_c3c4_with_wt, show.legend = F, group = c(1,2), sources.use = 'epithelium', targets.use = c('c3', 'c4'), title.name = 'N of interactions')
# gg_c3c4_weight_fromEpithelium <- compareInteractions(cellchat_c3c4_with_wt, show.legend = F, group = c(1,2), measure = "weight", sources.use = 'epithelium', targets.use = c('c3', 'c4'), digits = 4, title.name = 'Weight of interactions')
# png('./figures/LR_byGenotype/Diff_c3c4_vs_wt_fromEpith.png', width = 2000, height = 1500, res = 300)
# gg_c3c4_number_fromEpithelium + gg_c3c4_weight_fromEpithelium
# dev.off()
# 
# # from immune
# gg_c3c4_number_fromImmune <- compareInteractions(cellchat_c3c4_with_wt, show.legend = F, group = c(1,2), sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c3', 'c4'))
# gg_c3c4_weight_fromImmune <- compareInteractions(cellchat_c3c4_with_wt, show.legend = F, group = c(1,2), measure = "weight", sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c3', 'c4'), digits = 4)
# png('./figures/LR_byGenotype/Diff_c3c4_vs_wt_fromImmune.png', width = 2000, height = 1500, res = 300)
# gg_c3c4_number_fromImmune + gg_c3c4_weight_fromImmune
# dev.off()


################
# PRN vs wt

# from epithelium
gg_PRN_number_fromEpithelium <- compareInteractions(cellchat_PRN_with_wt, show.legend = F, group = c(1,2), sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c5', 'c6', 'c7'))
gg_PRN_weight_fromEpithelium <- compareInteractions(cellchat_PRN_with_wt, show.legend = F, group = c(1,2), measure = "weight", sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c5', 'c6', 'c7'), digits = 4)
png('./figures/LR_byGenotype/Diff_PRN_vs_wt_fromEpithelium2.png', width = 2000, height = 1500, res = 300)
gg_PRN_number_fromEpithelium + gg_PRN_weight_fromEpithelium
dev.off()

# from immune
gg_PRN_number_fromImmune <- compareInteractions(cellchat_PRN_with_wt, show.legend = F, group = c(1,2), sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c5', 'c6', 'c7'))
gg_PRN_weight_fromImmune <- compareInteractions(cellchat_PRN_with_wt, show.legend = F, group = c(1,2), measure = "weight", sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c5', 'c6', 'c7'), digits = 4)
png('./figures/LR_byGenotype/Diff_PRN_vs_wt_fromImmune.png', width = 2000, height = 1500, res = 300)
gg_PRN_number_fromImmune + gg_PRN_weight_fromImmune
dev.off()

################
# NP vs HiMYC

# from epithelium
gg_NP_vs_HiMYC_number_fromEpithelium <- compareInteractions(cellchat_NP_with_HiMYC, show.legend = F, group = c(1,2), sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c3', 'c4'))
gg_NP_vs_HiMYC_weight_fromEpithelium <- compareInteractions(cellchat_NP_with_HiMYC, show.legend = F, group = c(1,2), measure = "weight", sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c3', 'c4'), digits = 4)
png('./figures/LR_byGenotype/Diff_NP_vs_HiMYC_fromEpithelium2.png', width = 2000, height = 1500, res = 300)
gg_NP_vs_HiMYC_number_fromEpithelium + gg_NP_vs_HiMYC_weight_fromEpithelium
dev.off()

# from immune
gg_NP_vs_HiMYC_number_fromImmune <- compareInteractions(cellchat_NP_with_HiMYC, show.legend = F, group = c(1,2), sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c3', 'c4'))
gg_NP_vs_HiMYC_weight_fromImmune <- compareInteractions(cellchat_NP_with_HiMYC, show.legend = F, group = c(1,2), measure = "weight", sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c3', 'c4'), digits = 4)
png('./figures/LR_byGenotype/Diff_NP_vs_HiMYC_fromImmune.png', width = 2000, height = 1500, res = 300)
gg_NP_vs_HiMYC_number_fromImmune + gg_NP_vs_HiMYC_weight_fromImmune
dev.off()


################
# PRN vs c3c4

# from epithelium
gg_PRN_vs_c3c4_number_fromEpitheliumToAllstroma <- compareInteractions(cellchat_PRN_with_c3c4, show.legend = F, group = c(1,2), sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c5', 'c6', 'c7'))
gg_PRN_vs_c3c4_weight_fromEpithelium_toAllstroma <- compareInteractions(cellchat_PRN_with_c3c4, show.legend = F, group = c(1,2), measure = "weight", sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c5', 'c6', 'c7'), digits = 4)
png('./figures/LR_byGenotype/Diff_PRN_vs_c3c4_fromEpitheliumTo_c5c6c7.png', width = 2000, height = 1500, res = 300)
gg_PRN_vs_c3c4_number_fromEpitheliumToAllstroma + gg_PRN_vs_c3c4_weight_fromEpithelium_toAllstroma
dev.off()

# from immune
gg_PRN_vs_c3c4_number_fromImmuneToAllstroma <- compareInteractions(cellchat_PRN_with_c3c4, show.legend = F, group = c(1,2), sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c5', 'c6', 'c7'))
gg_PRN_vs_c3c4_weight_fromImmuneToAllstroma <- compareInteractions(cellchat_PRN_with_c3c4, show.legend = F, group = c(1,2), measure = "weight", sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c5', 'c6', 'c7'), digits = 4)
png('./figures/LR_byGenotype/Diff_PRN_vs_c3c4_fromImmuneTo_c5c6c7.png', width = 2000, height = 1500, res = 300)
gg_PRN_vs_c3c4_number_fromImmuneToAllstroma + gg_PRN_vs_c3c4_weight_fromImmuneToAllstroma
dev.off()

################
# all mutants vs all WT

# from epithelium
gg_mutants_vs_WT_number_fromEpitheliumToAllstroma <- compareInteractions(cellchat_mutants_with_WT, show.legend = F, group = c(1,2), sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'))
gg_mutants_vs_WT_weight_fromEpithelium_toAllstroma <- compareInteractions(cellchat_mutants_with_WT, show.legend = F, group = c(1,2), measure = "weight", sources.use = c('luminal', 'basal', 'neuroendocrine'), targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), digits = 4)
png('./figures/LR_byGenotype/Diff_mutants_vs_WT_fromEpitheliumToAllStroma2.png', width = 2000, height = 1500, res = 300)
gg_mutants_vs_WT_number_fromEpitheliumToAllstroma + gg_mutants_vs_WT_weight_fromEpithelium_toAllstroma
dev.off()

# from immune
gg_mutants_vs_WT_number_fromImmuneToAllstroma <- compareInteractions(cellchat_mutants_with_WT, show.legend = F, group = c(1,2), sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'))
gg_mutants_vs_WT_weight_fromImmuneToAllstroma <- compareInteractions(cellchat_mutants_with_WT, show.legend = F, group = c(1,2), measure = "weight", sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), targets.use = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'), digits = 4)
png('./figures/LR_byGenotype/Diff_mutants_vs_WT_fromImmuneToAllStroma.png', width = 2000, height = 1500, res = 300)
gg_mutants_vs_WT_number_fromImmuneToAllstroma + gg_mutants_vs_WT_weight_fromImmuneToAllstroma
dev.off()



#############################################
# heatmap of differential number of interactions and strength
#############################################
gg_c3c4_vs_wt_heatmap <- netVisual_heatmap(cellchat_c3c4_with_wt)
gg_c3c4_vs_wt_heatmap_weight <- netVisual_heatmap(cellchat_c3c4_with_wt, measure = "weight")
gg_c3c4_vs_wt_heatmap + gg_c3c4_vs_wt_heatmap_weight



#############################################
weight.max_c3c4 <- getMaxWeight(object.list_c3c4, attribute = c("idents","count"))
weight.max_PRN <- getMaxWeight(object.list_PRN, attribute = c("idents","count"))

## Number of interactions between wt and mutant
# C3c4 vs wt
png('./figures/LR_byGenotype/c3c4_vs_wt_Diff_N_interactions.png', width = 2000, height = 1500, res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_c3c4)) {
  netVisual_circle(object.list_c3c4[[i]]@net$count, 
                   sources.use = c('epithelium', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"),
                   targets.use = c('c3', 'c4'),
                   remove.isolate = T,
                   weight.scale = T, 
                   label.edge= F, 
                   edge.weight.max = weight.max_c3c4[2], 
                   #edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", 
                                       names(object.list_c3c4)[i])
  )
}
dev.off()

# PRN vs wt
png('./figures/LR_byGenotype/PRN_vs_wt_Diff_N_interactions.png', width = 2000, height = 1500, res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_PRN)) {
  netVisual_circle(object.list_PRN[[i]]@net$count, 
                   sources.use = c('epithelium', "B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"),
                   targets.use = c('c5', 'c6', 'c7'),
                   remove.isolate = T,
                   weight.scale = T, 
                   label.edge= F, 
                   edge.weight.max = weight.max_PRN[2], 
                   #edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", 
                                       names(object.list_PRN)[i])
                   )
}
dev.off()






# Compute and visualize the pathway distance in the learned joint manifold
#rankSimilarity(cellchat, type = "functional")

# Compare the overall information flow of each signaling pathway
gg_overall_c3c4_vs_wt <- rankNet(cellchat_c3c4_with_wt, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('cyan3', 'indianred1'), font.size = 6)
gg_overall_PRN_vs_wt <- rankNet(cellchat_PRN_with_wt, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('cyan3', 'indianred1'), font.size = 6)
gg_overall_mutants_vs_WT <- rankNet(cellchat_mutants_with_WT, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('cyan3', 'indianred1'), font.size = 6)

png('./figures/LR_byGenotype/c3c4_vs_WT_Diff_informationFlow.png', width = 1500, height = 2500, res = 300)
gg_overall_c3c4_vs_wt 
dev.off()

png('./figures/LR_byGenotype/PRN_vs_WT_Diff_informationFlow.png', width = 1500, height = 2500, res = 300)
gg_overall_PRN_vs_wt 
dev.off()

png('./figures/LR_byGenotype/mutants_vs_WT_Diff_informationFlow.png', width = 1500, height = 2500, res = 300)
gg_overall_mutants_vs_WT 
dev.off()

#####################################################
# Identify the upgulated and down-regulated signaling ligand-receptor pairs
#gg1 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in mutants", angle.x = 45, remove.isolate = T, n.colors = 5, min.quantile = 0.25, line.on= T)
#gg2 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in mutants", angle.x = 45, remove.isolate = T, n.colors = 5)

#gg1 + gg2

#####################################################
# Identify dysfunctional signaling by using differential expression analysis

# perform differential expression analysis
cellchat_PRN_with_wt <- identifyOverExpressedGenes(cellchat_PRN_with_wt, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.2, thresh.fc = 0.1, thresh.p = 1)
cellchat_NP_with_wt <- identifyOverExpressedGenes(cellchat_NP_with_wt, group.dataset = "datasets", pos.dataset = 'NP', features.name = 'NP', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_NP_with_HiMYC <- identifyOverExpressedGenes(cellchat_NP_with_HiMYC, group.dataset = "datasets", pos.dataset = 'NP', features.name = 'NP', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_HiMYC_with_Terg <- identifyOverExpressedGenes(cellchat_HiMYC_with_Terg, group.dataset = "datasets", pos.dataset = 'HiMYC', features.name = 'HiMYC', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_Terg_with_HiMyc <- identifyOverExpressedGenes(cellchat_Terg_with_HiMyc, group.dataset = "datasets", pos.dataset = 'Terg', features.name = 'Terg', only.pos = FALSE, thresh.pc = 0.2, thresh.fc = 0.1, thresh.p = 1)
cellchat_Terg_with_PN <- identifyOverExpressedGenes(cellchat_Terg_with_PN, group.dataset = "datasets", pos.dataset = 'Terg', features.name = 'Terg', only.pos = FALSE, thresh.pc = 0.2, thresh.fc = 0.1, thresh.p = 1)
#cellchat_PRN_with_Terg <- identifyOverExpressedGenes(cellchat_PRN_with_Terg, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#cellchat_PRN_with_PN <- identifyOverExpressedGenes(cellchat_PRN_with_PN, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_c3c4_with_wt <- identifyOverExpressedGenes(cellchat_c3c4_with_wt, group.dataset = "datasets", pos.dataset = 'PN / HiMYC / T-ERG', features.name = 'PN / HiMYC / T-ERG', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_PRN_with_c3c4 <- identifyOverExpressedGenes(cellchat_PRN_with_c3c4, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_mutants_with_WT <- identifyOverExpressedGenes(cellchat_mutants_with_WT, group.dataset = "datasets", pos.dataset = 'mutants', features.name = 'mutants', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net_PRN_vs_wt <- netMappingDEG(cellchat_PRN_with_wt, features.name = 'PRN')
net_NP_vs_wt <- netMappingDEG(cellchat_NP_with_wt, features.name = 'NP')
net_NP_vs_HiMYC <- netMappingDEG(cellchat_NP_with_HiMYC, features.name = 'NP')
net_HiMYC_vs_Terg <- netMappingDEG(cellchat_HiMYC_with_Terg, features.name = 'HiMYC')
net_Terg_vs_HiMYC <- netMappingDEG(cellchat_Terg_with_HiMyc, features.name = 'Terg')
net_Terg_vs_PN <- netMappingDEG(cellchat_Terg_with_PN, features.name = 'Terg')
#net_PRNvsHiMYC <- netMappingDEG(cellchat_PRN_with_HiMYC, features.name = 'PRN')
#net_PRNvsTerg <- netMappingDEG(cellchat_PRN_with_Terg, features.name = 'PRN')
net_c3c4_vs_wt <- netMappingDEG(cellchat_c3c4_with_wt, features.name = 'PN / HiMYC / T-ERG')
net_PRN_vs_c3c4 <- netMappingDEG(cellchat_PRN_with_c3c4, features.name = 'PRN')
net_mutants_vs_WT <- netMappingDEG(cellchat_mutants_with_WT, features.name = 'mutants')

# extract the ligand-receptor pairs with upregulated ligands in mutant
net.up_PRN_vs_wt <- subsetCommunication(cellchat_PRN_with_wt, net = net_PRN_vs_wt, datasets = "PRN", ligand.logFC = 0.1, receptor.logFC = NULL)
net.up_NP_vs_wt <- subsetCommunication(cellchat_NP_with_wt, net = net_NP_vs_wt, datasets = "NP", ligand.logFC = 0.1, receptor.logFC = NULL)
net.up_NP_vs_HiMYC <- subsetCommunication(cellchat_NP_with_HiMYC, net = net_NP_vs_HiMYC, datasets = "NP", ligand.logFC = 0.1, receptor.logFC = NULL)
net.up_HiMYC_vs_Terg <- subsetCommunication(cellchat_HiMYC_with_Terg, net = net_HiMYC_vs_Terg, datasets = "HiMYC", ligand.logFC = 0.1, receptor.logFC = NULL)
net.up_Terg_vs_HiMYC <- subsetCommunication(cellchat_Terg_with_HiMyc, net = net_Terg_vs_HiMYC, datasets = "Terg", ligand.logFC = 0.1, receptor.logFC = NULL)
net.up_Terg_vs_PN <- subsetCommunication(cellchat_Terg_with_PN, net = net_Terg_vs_PN, datasets = "Terg", ligand.logFC = 0.1, receptor.logFC = NULL)
#net.up_PRNvshiMYC <- subsetCommunication(cellchat_PRN_with_HiMYC, net = net_PRNvsHiMYC, datasets = "PRN", ligand.logFC = 0.2, receptor.logFC = NULL)
#net.up_PRNvsTerg <- subsetCommunication(cellchat_PRN_with_Terg, net = net_PRNvsTerg, datasets = "PRN", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_c3c4_vs_wt <- subsetCommunication(cellchat_c3c4_with_wt, net = net_c3c4_vs_wt, datasets = "PN / HiMYC / T-ERG", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_PRN_vs_c3c4 <- subsetCommunication(cellchat_PRN_with_c3c4, net = net_PRN_vs_c3c4, datasets = "PRN", ligand.logFC = 0.1, receptor.logFC = NULL)
net.up_mutants_vs_WT <- subsetCommunication(cellchat_mutants_with_WT, net = net_mutants_vs_WT, datasets = "mutants", ligand.logFC = 0.1, receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in wildtypes, i.e.,downregulated in mutants
net.down_PRN_vs_wt <- subsetCommunication(cellchat_PRN_with_wt, net = net_PRN_vs_wt, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_NP_vs_wt <- subsetCommunication(cellchat_NP_with_wt, net = net_NP_vs_wt, datasets = "NP", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_NP_vs_HiMYC <- subsetCommunication(cellchat_NP_with_HiMYC, net = net_NP_vs_HiMYC, datasets = "NP", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_HiMYC_vs_Terg <- subsetCommunication(cellchat_HiMYC_with_Terg, net = net_HiMYC_vs_Terg, datasets = "HiMYC", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_Terg_vs_HiMYC <- subsetCommunication(cellchat_Terg_with_HiMyc, net = net_Terg_vs_HiMYC, datasets = "HiMYC", ligand.logFC = -0.1, receptor.logFC = -0.1)
#net.down_PRNvsHiMYC <- subsetCommunication(cellchat_PRN_with_HiMYC, net = net_PRNvsHiMYC, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)
#net.down_PRNvsTerg <- subsetCommunication(cellchat_PRN_with_Terg, net = net_PRNvsTerg, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_c3c4_vs_wt <- subsetCommunication(cellchat_c3c4_with_wt, net = net_c3c4_vs_wt, datasets = "PN / HiMYC / T-ERG", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_PRN_vs_c3c4 <- subsetCommunication(cellchat_PRN_with_c3c4, net = net_PRN_vs_c3c4, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_mutants_vs_WT <- subsetCommunication(cellchat_mutants_with_WT, net = net_mutants_vs_WT, datasets = "mutants", ligand.logFC = -0.1, receptor.logFC = -0.1)

# do further deconvolution to obtain the individual signaling genes.
gene.up_PRN_vs_wt <- extractGeneSubsetFromPair(net.up_PRN_vs_wt, cellchat_PRN_with_wt)
gene.up_NP_vs_wt <- extractGeneSubsetFromPair(net.up_NP_vs_wt, cellchat_NP_with_wt)
#gene.up_Terg <- extractGeneSubsetFromPair(net.up_Terg, cellchat_Terg_with_wt)
gene.up_HiMYC_vs_Terg <- extractGeneSubsetFromPair(net.up_HiMYC_vs_Terg, cellchat_HiMYC_with_Terg)
gene.up_Terg_vs_HiMYC <- extractGeneSubsetFromPair(net.up_Terg_vs_HiMYC, cellchat_Terg_with_HiMyc)
#gene.up_PRNvsHiMYC <- extractGeneSubsetFromPair(net.up_PRNvshiMYC, cellchat_PRN_with_HiMYC)
#gene.up_PRNvsTerg <- extractGeneSubsetFromPair(net.up_PRNvsTerg, cellchat_PRN_with_Terg)
gene.up_c3c4_vs_wt <- extractGeneSubsetFromPair(net.up_c3c4_vs_wt, cellchat_c3c4_with_wt)
gene.up_PRN_vs_c3c4 <- extractGeneSubsetFromPair(net.up_PRN_vs_c3c4, cellchat_PRN_with_c3c4)
gene.up_mutants_vs_WT <- extractGeneSubsetFromPair(net.up_mutants_vs_WT, cellchat_mutants_with_WT)

gene.down_PRN_vs_wt <- extractGeneSubsetFromPair(net.down_PRN_vs_wt, cellchat_PRN_with_wt)
gene.down_NP_vs_wt <- extractGeneSubsetFromPair(net.down_NP_vs_wt, cellchat_NP_with_wt)
#gene.down_Terg <- extractGeneSubsetFromPair(net.down_Terg, cellchat_Terg_with_wt)
gene.down_HiMYC_vs_Terg <- extractGeneSubsetFromPair(net.down_HiMYC_vs_Terg, cellchat_HiMYC_with_Terg)
gene.down_Terg_vs_HiMYC <- extractGeneSubsetFromPair(net.down_Terg_vs_HiMYC, cellchat_Terg_with_HiMyc)
#gene.down_PRNvsHiMYC <- extractGeneSubsetFromPair(net.down_PRNvsHiMYC, cellchat_PRN_with_HiMYC)
#gene.down_PRNvsTerg <- extractGeneSubsetFromPair(net.down_PRNvsTerg, cellchat_PRN_with_Terg)
gene.down_c3c4_vs_wt <- extractGeneSubsetFromPair(net.down_c3c4_vs_wt, cellchat_c3c4_with_wt)
gene.down_PRN_vs_c3c4 <- extractGeneSubsetFromPair(net.down_PRN_vs_c3c4, cellchat_PRN_with_c3c4)
gene.down_mutants_vs_WT <- extractGeneSubsetFromPair(net.down_mutants_vs_WT, cellchat_mutants_with_WT)

# visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up_PRN_vs_wt = net.up_PRN_vs_wt[, "interaction_name", drop = F]
pairLR.use.up_NP_vs_wt = net.up_NP_vs_wt[, "interaction_name", drop = F]
#pairLR.use.up_Terg = net.up_Terg[, "interaction_name", drop = F]
pairLR.use.up_HiMYC_vs_Terg = net.up_HiMYC_vs_Terg[, "interaction_name", drop = F]
pairLR.use.up_Terg_vs_HiMYC = net.up_Terg_vs_HiMYC[, "interaction_name", drop = F]
#pairLR.use.up_PRNvsHiMYC = net.up_PRNvshiMYC[, "interaction_name", drop = F]
#pairLR.use.up_PRNvsTerg = net.up_PRNvsTerg[, "interaction_name", drop = F]
pairLR.use.up_c3c4_vs_wt = net.up_c3c4_vs_wt[, "interaction_name", drop = F]
pairLR.use.up_PRN_vs_c3c4 = net.up_PRN_vs_c3c4[, "interaction_name", drop = F]
pairLR.use.up_mutants_vs_WT = net.up_mutants_vs_WT[, "interaction_name", drop = F]

pairLR.use.down_PRN_vs_wt = net.down_PRN_vs_wt[, "interaction_name", drop = F]
pairLR.use.down_NP_vs_wt = net.down_NP_vs_wt[, "interaction_name", drop = F]
#pairLR.use.down_Terg = net.down_Terg[, "interaction_name", drop = F]
pairLR.use.down_HiMYC_vs_Terg = net.down_HiMYC_vs_Terg[, "interaction_name", drop = F]
pairLR.use.down_Terg_vs_HiMYC = net.down_Terg_vs_HiMYC[, "interaction_name", drop = F]
#pairLR.use.down_PRNvsHiMYC = net.down_PRNvsHiMYC[, "interaction_name", drop = F]
#pairLR.use.down_PRNvsTerg = net.down_PRNvsTerg[, "interaction_name", drop = F]
pairLR.use.down_c3c4_vs_wt = net.down_c3c4_vs_wt[, "interaction_name", drop = F]
pairLR.use.down_PRN_vs_c3c4 = net.down_PRN_vs_c3c4[, "interaction_name", drop = F]
pairLR.use.down_mutants_vs_WT = net.down_mutants_vs_WT[, "interaction_name", drop = F]

######################################################################
# from epith to common clusters in all mutants vs WT
######################################################################
colMap <- c('powderblue', 'plum2', 'violetred4', 'slategray4', 'blue', 'darkorange1', 'forestgreen', 'red2', 'darkorchid', 'saddlebrown', 'hotpink1', 'olivedrab3', 'lemonchiffon4', 'CadetBlue3', 'salmon3', 'mediumvioletred', 'cyan4')

######################################################################
# from epith to c3,c4 in Terg vs HiMYC
######################################################################
png('./figures/LR_byGenotype/up_epith_Terg_vs_HiMYC2.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_Terg_HiMYC[[1]], 
                     slot.name = 'net',
                     net = net.up_Terg_vs_HiMYC,
                     sources.use = c('luminal', 'basal', 'neuroendocrine'),
                     targets.use = c('c3', 'c4'),
                     thresh = 0.05,
                     lab.cex = 1,
                     small.gap = 5,
                     #big.gap = 20,
                     legend.pos.y = 165,
                     legend.pos.x = 25,
                     title.name = "",
                     color.use = colMap
                     )
dev.off()

######################################################################
# from epith to c3,c4 in Terg vs PN
######################################################################
png('./figures/LR_byGenotype/Terg_vs_PN/up_epith_Terg_vs_PN2.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_Terg_PN[[1]], 
                     slot.name = 'net',
                     net = net.up_Terg_vs_PN,
                     sources.use = c('luminal', 'basal', 'neuroendocrine'),
                     targets.use = c('c3', 'c4'),
                     thresh = 0.05,
                     lab.cex = 1,
                     small.gap = 5,
                     #big.gap = 20,
                     legend.pos.y = 165,
                     legend.pos.x = 25,
                     title.name = "",
                     color.use = colMap
                     
                     )
dev.off()
######################################################################
# from epith to c3,c4 in NP vs HiMYC
######################################################################
png('./figures/LR_byGenotype/up_epith_NP_vs_HiMYC2.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_NP_HiMYC[[1]], 
                     slot.name = 'net',
                     net = net.up_NP_vs_HiMYC,
                     sources.use = c('luminal', 'basal', 'neuroendocrine'),
                     targets.use = c('c3', 'c4'),
                     thresh = 0.05,
                     lab.cex = 1,
                     small.gap = 5,
                     #big.gap = 20,
                     legend.pos.y = 165,
                     legend.pos.x = 25,
                     title.name = "",
                     color.use = colMap
                     )
dev.off()

######################################################################
# from epith to c5,c6,c7 in PRN vs wt
######################################################################
png('./figures/LR_byGenotype/up_epith_PRN_vs_wt2.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRN_vs_wt, 
                     sources.use = c('luminal', 'basal', 'neuroendocrine'), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.01, 
                     lab.cex = 0.8, 
                     small.gap = 4, 
                     #big.gap = 20,
                     legend.pos.y = 180, 
                     legend.pos.x = 10, 
                     color.use = colMap,
                     title.name = "")
dev.off()

png('./figures/LR_byGenotype/up_epith_PRN.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(cellchat_mouse_all_PRN, 
                     slot.name = 'net', 
                     sources.use = c('luminal', 'basal', 'neuroendocrine'), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.01, 
                     lab.cex = 0.5, 
                     small.gap = 4, 
                     #big.gap = 20,
                     legend.pos.y = 180, 
                     legend.pos.x = 10, 
                     color.use = colMap,
                     title.name = "")
dev.off()

######################################################################
# from epith to c5,c6,c7 in PRN vs c3c4  
######################################################################
png('./figures/LR_byGenotype/up_epith_PRN_vs_c3c4.png', width = 3500, height = 3500, res = 300)
netVisual_chord_gene(object.list_PRN_c3c4[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRN_vs_c3c4, 
                     sources.use = c('luminal', 'basal', 'neuroendocrine'), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.01, 
                     lab.cex = 0.55, 
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 210, 
                     legend.pos.x = 10, 
                     color.use = colMap,
                     title.name = "")
dev.off()



###########################################################################################
######################################################################
# from immune to c0,c1,c2 in all mutants vs all WT
######################################################################
png('./figures/LR_byGenotype/up_epithelium_commonClusters_mutants_vs_WT2.png',  width = 3500, height = 3500, res = 350)
netVisual_chord_gene(object.list_mutants_WT[[1]], 
                     slot.name = 'net', 
                     net = net.up_mutants_vs_WT, 
                     sources.use = c('luminal', 'basal', 'neuroendocrine'), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     thresh = 0.01, 
                     lab.cex = 0.6,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 175, 
                     legend.pos.x = 5, 
                     color.use = colMap,
                     title.name = "")
dev.off()

######################################################################
# from immune to c0,c1,c2 in all mutants vs all WT
######################################################################
png('./figures/LR_byGenotype/up_immune_commonClusters_mutants_vs_WT2.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_mutants_WT[[1]], 
                     slot.name = 'net', 
                     net = net.up_mutants_vs_WT, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c0', 'c1', 'c2'), 
                     thresh = 0.05, 
                     lab.cex = 1,
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
png('./figures/LR_byGenotype/up_immune_c3c4_vs_wt.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_c3c4[[1]], 
                     slot.name = 'net', 
                     net = net.up_c3c4_vs_wt, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     thresh = 0.05, 
                     lab.cex = 1,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 175, 
                     legend.pos.x = 5, 
                     color.use = colMap,
                     title.name = "")
dev.off()

######################################################################
# from immune to c3,c4 in HiMYC vs Terg
######################################################################
png('./figures/LR_byGenotype/up_immune_HiMYC_vs_Terg.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_HiMYC_Terg[[1]], 
                     slot.name = 'net',
                     net = net.up_HiMYC_vs_Terg,
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'),
                     thresh = 0.05,
                     lab.cex = 1,
                     small.gap = 5,
                     #big.gap = 20,
                     legend.pos.y = 165,
                     legend.pos.x = 25,
                     title.name = "")
dev.off()

######################################################################
# from immune to c3,c4 in Terg vs HiMYC
######################################################################
png('./figures/LR_byGenotype/up_immune_Terg_vs_HiMYC.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_Terg_HiMYC[[1]], 
                     slot.name = 'net',
                     net = net.up_Terg_vs_HiMYC,
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'),
                     thresh = 0.05,
                     lab.cex = 1,
                     small.gap = 5,
                     #big.gap = 20,
                     legend.pos.y = 165,
                     legend.pos.x = 10,
                     title.name = "",
                     color.use = colMap
                     )
dev.off()

######################################################################
# from immune to c3,c4 in Terg vs HiMYC
######################################################################
png('./figures/LR_byGenotype/up_immune_Terg_vs_PN.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_Terg_PN[[1]], 
                     slot.name = 'net',
                     net = net.up_Terg_vs_PN,
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'),
                     thresh = 0.05,
                     lab.cex = 1,
                     small.gap = 5,
                     #big.gap = 20,
                     legend.pos.y = 165,
                     legend.pos.x = 10,
                     title.name = "")
dev.off()

######################################################################
# from immune to c3,c4 in NP vs HiMYC
######################################################################
png('./figures/LR_byGenotype/up_immune_NP_vs_HiMYC.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_NP_HiMYC[[1]], 
                     slot.name = 'net',
                     net = net.up_NP_vs_HiMYC,
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'),
                     thresh = 0.05,
                     lab.cex = 1,
                     small.gap = 5,
                     #big.gap = 20,
                     legend.pos.y = 165,
                     legend.pos.x = 25,
                     title.name = "")
dev.off()

######################################################################
# from immune to c5,c6,c7 in PRN vs wt
######################################################################
png('./figures/LR_byGenotype/up_immune_PRN_vs_wt.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRN_vs_wt, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 1, 
                     small.gap = 3, 
                     #big.gap = 20,
                     legend.pos.y = 200, 
                     legend.pos.x = 5, 
                     color.use = colMap,
                     title.name = "")
dev.off()

######################################################################
# immune-immune interactions in PRN
######################################################################
png('./figures/LR_byGenotype/up_immuneWithImmune_PRN_vs_wt.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRN_vs_wt, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     thresh = 0.05, 
                     lab.cex = 0.6, 
                     small.gap = 1.4, 
                     #big.gap = 20,
                     legend.pos.y = 200, 
                     legend.pos.x = 5, 
                     color.use = colMap,
                     title.name = "")
dev.off()

######################################################################
# from immune to c5,c6,c7 in PRN vs c3c4
######################################################################
png('./figures/LR_byGenotype/up_immune_PRN_vs_c3c4.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN_c3c4[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRN_vs_c3c4, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 1, 
                     small.gap = 3, 
                     #big.gap = 20,
                     legend.pos.y = 200, 
                     legend.pos.x = 5, 
                     color.use = colMap,
                     title.name = "")
dev.off()


######################################################
# Compare outgoing signaling associated with each cell population: c3c4 vs wt
######################################################
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union_c3c4_vs_wt <- union(object.list_c3c4[[i]]@netP$pathways, object.list_c3c4[[i+1]]@netP$pathways)

ht1_c3c4_vs_wt = netAnalysis_signalingRole_heatmap(object.list_c3c4[[i]], pattern = "outgoing", signaling = pathway.union_c3c4_vs_wt, title = 'PN / HiMYC / T-ERG', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14, color.use = colMap)
ht2_c3c4_vs_wt = netAnalysis_signalingRole_heatmap(object.list_c3c4[[i+1]], pattern = "outgoing", signaling = pathway.union_c3c4_vs_wt, title = 'WT', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14, color.use = colMap)


tiff('./figures/LR_byGenotype/Diff_heatmap_c3c4_vs_wt.tiff', width = 4000, height = 3000, res = 350)
draw(ht1_c3c4_vs_wt + ht2_c3c4_vs_wt, ht_gap = unit(2, "cm"))
dev.off()


######################################################
# Compare outgoing signaling associated with each cell population: PRN vs wt
######################################################
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union_PRN_vs_wt <- union(object.list_PRN[[i]]@netP$pathways, object.list_PRN[[i+1]]@netP$pathways)

ht1_PRN_vs_wt = netAnalysis_signalingRole_heatmap(object.list_PRN[[i]], pattern = "outgoing", signaling = pathway.union_PRN_vs_wt, title = 'PRN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)
ht2_PRN_vs_wt = netAnalysis_signalingRole_heatmap(object.list_PRN[[i+1]], pattern = "outgoing", signaling = pathway.union_PRN_vs_wt, title = 'WT', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)


tiff('./figures/LR_byGenotype/Diff_heatmap_PRN_vs_wt.tiff', width = 4000, height = 3000, res = 350)
draw(ht1_PRN_vs_wt + ht2_PRN_vs_wt, ht_gap = unit(2, "cm"))
dev.off()


######################################################
# Compare outgoing signaling associated with each cell population: PRN vs c3c4
######################################################
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union_PRN_vs_c3c4 <- union(object.list_PRN_c3c4[[i]]@netP$pathways, object.list_PRN_c3c4[[i+1]]@netP$pathways)

ht1_PRN_vs_c3c4 = netAnalysis_signalingRole_heatmap(object.list_PRN_c3c4[[i]], pattern = "outgoing", signaling = pathway.union_PRN_vs_c3c4, title = 'PRN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)
ht2_PRN_vs_c3c4 = netAnalysis_signalingRole_heatmap(object.list_PRN_c3c4[[i+1]], pattern = "outgoing", signaling = pathway.union_PRN_vs_c3c4, title = 'PN / HiMYC / T-ERG', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)


tiff('./figures/LR_byGenotype/Diff_heatmap_PRN_vs_c3c4.tiff', width = 4000, height = 3000, res = 350)
draw(ht1_PRN_vs_c3c4 + ht2_PRN_vs_c3c4, ht_gap = unit(2, "cm"))
dev.off()



######################################################
# HiMYC interactions
######################################################
# from epithelium
png('./figures/LR_byGenotype/epithelium_c3c4_hiMYC2.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all_HiMYC, 
                     slot.name = "net",
                     sources.use = c('luminal', 'basal', 'neuroendocrine'), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     #color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

# from immune
png('./figures/LR_byGenotype/immune_c3c4_hiMYC.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all_HiMYC, 
                     slot.name = "net",
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

######################################################
# Terg interactions
######################################################
# from epithelium
png('./figures/LR_byGenotype/epithelium_c3c4_Terg2.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all_Terg, 
                     slot.name = "net",
                     sources.use = c('luminal', 'basal', 'neuroendocrine'), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     #color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

# from immune
png('./figures/LR_byGenotype/immune_c3c4_Terg.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all_Terg, 
                     slot.name = "net",
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()


######################################################
# NP interactions
######################################################
# from epithelium
png('./figures/LR_byGenotype/PN/epithelium_c3c4_PN2.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all_PN, 
                     slot.name = "net",
                     sources.use = c('luminal', 'basal', 'neuroendocrine'), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     #color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

# from immune
png('./figures/LR_byGenotype/PN/immune_c3c4_PN.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all_PN, 
                     slot.name = "net",
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.05, 
                     small.gap = 3,
                     #directional = 2,
                     #color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()


###########################################################
# Stroma to epithelium interactions in WTs

png('./figures/LR_byGenotype/PN/stroma_to_epithelium_WTs.png', width = 3500, height = 3000, res = 350)
netVisual_chord_gene(cellchat_mouse_all_WT, 
                     slot.name = "net",
                     targets.use = c('luminal', 'basal'), 
                     sources.use = c('c0', 'c1', 'c2'), 
                     lab.cex = 1, 
                     legend.pos.y = 170, 
                     legend.pos.x = 5,
                     #reduce = -1,
                     thresh = 0.01, 
                     small.gap = 3,
                     #directional = 2,
                     #color.use = colMap
                     #title.name = 'Signling networks from common clusters to immune cells'	
)
dev.off()

#############
# save tables
write.xlsx(net.up_Terg_vs_HiMYC, 
           file = 'tables/Table_S3.xlsx', 
           sheetName = 'LR_TRG_vs_HiMYC',
           row.names = F, append = TRUE
           )

write.xlsx(net.up_HiMYC_vs_Terg, 
           file = 'tables/Table_S3.xlsx', 
           sheetName = 'LR_HiMYC_vs_TRG',
           row.names = F,  append = TRUE
           )

write.xlsx(net.up_NP_vs_HiMYC, 
           file = 'tables/Table_S3.xlsx', 
           sheetName = 'LR_PN_vs_HiMYC',
           row.names = F,  append = TRUE
)

write.xlsx(net.up_Terg_vs_PN, 
           file = 'tables/Table_S3.xlsx', 
           sheetName = 'LR_TRG_vs_PN',
           row.names = F,  append = TRUE
)

write.xlsx(net.up_PRN_vs_wt, 
           file = 'tables/Table_S3.xlsx', 
           sheetName = 'LR_PRN_vs_WT',
           row.names = F,  append = TRUE
)



