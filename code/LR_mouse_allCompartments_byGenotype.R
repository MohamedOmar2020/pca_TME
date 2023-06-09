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

######################################################
# clean
######################################################
rm(adata_mouse_all, mouse_all_dataInput, mouse_all_dataInput_B6, mouse_all_dataInput_B6_129, mouse_all_dataInput_FVBN, mouse_all_dataInput_HiMYC, mouse_all_dataInput_PN, mouse_all_dataInput_PN_wt, mouse_all_dataInput_PRN, mouse_all_dataInput_PRN_wt)



#######################################################################################################################################
# Set the ligand-receptor interaction database
#######################################################################################################################################
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

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
     file = './objs/LR_genotypes.rda'
)

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

#####################################
#Lift up CellChat object and merge together
#####################################


# Define the cell labels to lift up
group.new = levels(cellchat_mouse_all_HiMYC@idents)

cellchat_mouse_all_B6_129 <- liftCellChat(cellchat_mouse_all_B6_129, group.new)
cellchat_mouse_all_PN <- liftCellChat(cellchat_mouse_all_PN, group.new)
cellchat_mouse_all_PRN_wt <- liftCellChat(cellchat_mouse_all_PRN_wt, group.new)
cellchat_mouse_all_Terg <- liftCellChat(cellchat_mouse_all_Terg, group.new)


#####################################
## Merge the objects

object.list_all <- list(PN = cellchat_mouse_all_PN, 
                    HiMYC = cellchat_mouse_all_HiMYC,
                    Terg = cellchat_mouse_all_Terg,
                    PRN = cellchat_mouse_all_PRN,
                    FVBN = cellchat_mouse_all_FVBN,
                    B6 = cellchat_mouse_all_B6,
                    B6_129 = cellchat_mouse_all_B6_129,
                    PRN_wt = cellchat_mouse_all_PRN_wt,
                    PN_wt = cellchat_mouse_all_PN_wt
                    )


cellchat <- mergeCellChat(object.list_all, add.names = names(object.list))

#####################################
## For PN

object.list_PN <- list(PN = cellchat_mouse_all_PN, 
                    PN_wt = cellchat_mouse_all_PN_wt
)


cellchat_PN_with_wt <- mergeCellChat(object.list_PN, add.names = names(object.list_PN))

#####################################
## For Terg

object.list_Terg <- list(Terg = cellchat_mouse_all_Terg, 
                       WT = cellchat_mouse_all_FVBN
)


cellchat_Terg_with_wt <- mergeCellChat(object.list_Terg, add.names = names(object.list_Terg))

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
## For PRN

object.list_PRN <- list(PRN = cellchat_mouse_all_PRN, 
                       PRN_wt = cellchat_mouse_all_PRN_wt
)


cellchat_PRN_with_wt <- mergeCellChat(object.list_PRN, add.names = names(object.list_PRN))


#####################################
## PRN and HiMYC

object.list_PRN_HiMYC <- list(PRN = cellchat_mouse_all_PRN, 
                        HiMYC = cellchat_mouse_all_HiMYC
)


cellchat_PRN_with_HiMYC <- mergeCellChat(object.list_PRN_HiMYC, add.names = names(object.list_PRN_HiMYC))

#####################################
## PRN and PN

object.list_PRN_PN <- list(PRN = cellchat_mouse_all_PRN, 
                              HiMYC = cellchat_mouse_all_PN
)


cellchat_PRN_with_PN <- mergeCellChat(object.list_PRN_PN, add.names = names(object.list_PRN_PN))

#####################################
## PRN and Terg

object.list_PRN_Terg <- list(PRN = cellchat_mouse_all_PRN, 
                           HiMYC = cellchat_mouse_all_Terg
)


cellchat_PRN_with_Terg <- mergeCellChat(object.list_PRN_Terg, add.names = names(object.list_PRN_Terg))

##########################################################################################################
# Compare the total number of interactions and interaction strength
##########################################################################################################
# For PN
gg_PN_number <- compareInteractions(cellchat_PN_with_wt, show.legend = F, group = c(1,2))
gg_PN_weight <- compareInteractions(cellchat_PN_with_wt, show.legend = F, group = c(1,2), measure = "weight")
gg_PN_number + gg_PN_weight

# For PRN
gg_PRN_number <- compareInteractions(cellchat_PRN_with_wt, show.legend = F, group = c(1,2))
gg_PRN_weight <- compareInteractions(cellchat_PRN_with_wt, show.legend = F, group = c(1,2), measure = "weight")
gg_PRN_number + gg_PRN_weight


gg_PN_heatmap <- netVisual_heatmap(cellchat_PN_with_wt)
#> Do heatmap based on a merged object
gg_PN_heatmap_weight <- netVisual_heatmap(cellchat_PN_with_wt, measure = "weight")
#> Do heatmap based on a merged object
gg_PN_heatmap + gg_PN_heatmap_weight


#############################################
weight.max_PN <- getMaxWeight(object.list_PN, attribute = c("idents","count"))
weight.max_PRN <- getMaxWeight(object.list_PRN, attribute = c("idents","count"))

# Number of interactions between wt and mutant
png('./figures/LR_byGenotype/PN_Diff_N_interactions.png', width = 2000, height = 1500, res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_PN)) {
  netVisual_circle(object.list_PN[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max_PN[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list_PN)[i]))
}
dev.off()

png('./figures/LR_byGenotype/PRN_Diff_N_interactions.png', width = 2000, height = 1500, res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_PRN)) {
  netVisual_circle(object.list_PRN[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max_PRN[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list_PRN)[i]))
}
dev.off()


# Compute and visualize the pathway distance in the learned joint manifold
#rankSimilarity(cellchat, type = "functional")

# Compare the overall information flow of each signaling pathway
gg_overall_PN <- rankNet(cellchat_PN_with_wt, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('cyan3', 'indianred1'))
gg_overall_PRN <- rankNet(cellchat_PRN_with_wt, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c('cyan3', 'indianred1'))

png('./figures/LR_byGenotype/PN_Diff_informationFlow.png', width = 1500, height = 2000, res = 300)
gg_overall_PN 
dev.off()

png('./figures/LR_byGenotype/PRN_Diff_informationFlow.png', width = 1500, height = 2000, res = 300)
gg_overall_PRN 
dev.off()


#####################################################
# Identify the upgulated and down-regulated signaling ligand-receptor pairs
#gg1 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in mutants", angle.x = 45, remove.isolate = T, n.colors = 5, min.quantile = 0.25, line.on= T)
#gg2 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in mutants", angle.x = 45, remove.isolate = T, n.colors = 5)

#gg1 + gg2

#####################################################
# Identify dysfunctional signaling by using differential expression analysis

# perform differential expression analysis
cellchat_PN_with_wt <- identifyOverExpressedGenes(cellchat_PN_with_wt, group.dataset = "datasets", pos.dataset = 'PN', features.name = 'PN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_PRN_with_wt <- identifyOverExpressedGenes(cellchat_PRN_with_wt, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_Terg_with_wt <- identifyOverExpressedGenes(cellchat_Terg_with_wt, group.dataset = "datasets", pos.dataset = 'Terg', features.name = 'Terg', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_HiMYC_with_wt <- identifyOverExpressedGenes(cellchat_HiMyc_with_wt, group.dataset = "datasets", pos.dataset = 'HiMYC', features.name = 'HiMYC', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_PRN_with_HiMYC <- identifyOverExpressedGenes(cellchat_PRN_with_HiMYC, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_PRN_with_Terg <- identifyOverExpressedGenes(cellchat_PRN_with_Terg, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
cellchat_PRN_with_PN <- identifyOverExpressedGenes(cellchat_PRN_with_PN, group.dataset = "datasets", pos.dataset = 'PRN', features.name = 'PRN', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net_PN <- netMappingDEG(cellchat_PN_with_wt, features.name = 'PN')
net_PRN <- netMappingDEG(cellchat_PRN_with_wt, features.name = 'PRN')
net_Terg <- netMappingDEG(cellchat_Terg_with_wt, features.name = 'Terg')
net_HiMYC <- netMappingDEG(cellchat_HiMYC_with_wt, features.name = 'HiMYC')
net_PRNvsHiMYC <- netMappingDEG(cellchat_PRN_with_HiMYC, features.name = 'PRN')
net_PRNvsTerg <- netMappingDEG(cellchat_PRN_with_Terg, features.name = 'PRN')
net_PRNvsPN <- netMappingDEG(cellchat_PRN_with_PN, features.name = 'PRN')

# extract the ligand-receptor pairs with upregulated ligands in mutant
net.up_PN <- subsetCommunication(cellchat_PN_with_wt, net = net_PN, datasets = "PN", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_PRN <- subsetCommunication(cellchat_PRN_with_wt, net = net_PRN, datasets = "PRN", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_Terg <- subsetCommunication(cellchat_Terg_with_wt, net = net_Terg, datasets = "Terg", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_HiMYC <- subsetCommunication(cellchat_HiMYC_with_wt, net = net_HiMYC, datasets = "HiMYC", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_PRNvshiMYC <- subsetCommunication(cellchat_PRN_with_HiMYC, net = net_PRNvsHiMYC, datasets = "PRN", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_PRNvsTerg <- subsetCommunication(cellchat_PRN_with_Terg, net = net_PRNvsTerg, datasets = "PRN", ligand.logFC = 0.2, receptor.logFC = NULL)
net.up_PRNvsPN <- subsetCommunication(cellchat_PRN_with_PN, net = net_PRNvsPN, datasets = "PRN", ligand.logFC = 0.2, receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in wildtypes, i.e.,downregulated in mutants
net.down_PN <- subsetCommunication(cellchat_PN_with_wt, net = net_PN, datasets = "PN", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_PRN <- subsetCommunication(cellchat_PRN_with_wt, net = net_PRN, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_Terg <- subsetCommunication(cellchat_Terg_with_wt, net = net_Terg, datasets = "Terg", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_HiMYC <- subsetCommunication(cellchat_HiMYC_with_wt, net = net_HiMYC, datasets = "HiMYC", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_PRNvsHiMYC <- subsetCommunication(cellchat_PRN_with_HiMYC, net = net_PRNvsHiMYC, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_PRNvsTerg <- subsetCommunication(cellchat_PRN_with_Terg, net = net_PRNvsTerg, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down_PRNvsPN <- subsetCommunication(cellchat_PRN_with_PN, net = net_PRNvsPN, datasets = "PRN", ligand.logFC = -0.1, receptor.logFC = -0.1)

# do further deconvolution to obtain the individual signaling genes.
gene.up_PN <- extractGeneSubsetFromPair(net.up_PN, cellchat_PN_with_wt)
gene.up_PRN <- extractGeneSubsetFromPair(net.up_PRN, cellchat_PRN_with_wt)
gene.up_Terg <- extractGeneSubsetFromPair(net.up_Terg, cellchat_Terg_with_wt)
gene.up_HiMYC <- extractGeneSubsetFromPair(net.up_HiMYC, cellchat_HiMYC_with_wt)
gene.up_PRNvsHiMYC <- extractGeneSubsetFromPair(net.up_PRNvshiMYC, cellchat_PRN_with_HiMYC)
gene.up_PRNvsTerg <- extractGeneSubsetFromPair(net.up_PRNvsTerg, cellchat_PRN_with_Terg)
gene.up_PRNvsPN <- extractGeneSubsetFromPair(net.up_PRNvsPN, cellchat_PRN_with_PN)

gene.down_PN <- extractGeneSubsetFromPair(net.down_PN, cellchat_PN_with_wt)
gene.down_PRN <- extractGeneSubsetFromPair(net.down_PRN, cellchat_PRN_with_wt)
gene.down_Terg <- extractGeneSubsetFromPair(net.down_Terg, cellchat_Terg_with_wt)
gene.down_HiMYC <- extractGeneSubsetFromPair(net.down_HiMYC, cellchat_HiMYC_with_wt)
gene.down_PRNvsHiMYC <- extractGeneSubsetFromPair(net.down_PRNvsHiMYC, cellchat_PRN_with_HiMYC)
gene.down_PRNvsTerg <- extractGeneSubsetFromPair(net.down_PRNvsTerg, cellchat_PRN_with_Terg)
gene.down_PRNvsPN <- extractGeneSubsetFromPair(net.down_PRNvsPN, cellchat_PRN_with_PN)

# visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up_PN = net.up_PN[, "interaction_name", drop = F]
pairLR.use.up_PRN = net.up_PRN[, "interaction_name", drop = F]
pairLR.use.up_Terg = net.up_Terg[, "interaction_name", drop = F]
pairLR.use.up_HiMYC = net.up_HiMYC[, "interaction_name", drop = F]
pairLR.use.up_PRNvsHiMYC = net.up_PRNvshiMYC[, "interaction_name", drop = F]
pairLR.use.up_PRNvsTerg = net.up_PRNvsTerg[, "interaction_name", drop = F]
pairLR.use.up_PRNvsPN = net.up_PRNvsPN[, "interaction_name", drop = F]

pairLR.use.down_PN = net.down_PN[, "interaction_name", drop = F]
pairLR.use.down_PRN = net.down_PRN[, "interaction_name", drop = F]
pairLR.use.down_Terg = net.down_Terg[, "interaction_name", drop = F]
pairLR.use.down_HiMYC = net.down_HiMYC[, "interaction_name", drop = F]
pairLR.use.down_PRNvsHiMYC = net.down_PRNvsHiMYC[, "interaction_name", drop = F]
pairLR.use.down_PRNvsTerg = net.down_PRNvsTerg[, "interaction_name", drop = F]
pairLR.use.down_PRNvsPN = net.down_PRNvsPN[, "interaction_name", drop = F]

######################################################################
# from epith to c3,c4 in PN
######################################################################
png('./figures/LR_byGenotype/up_epith_PN.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_PN[[2]], 
                     slot.name = 'net', 
                     net = net.up_PN, 
                     sources.use = 'epithelium', 
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
# from epith to c3,c4 in Terg
######################################################################
png('./figures/LR_byGenotype/up_epith_Terg.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_Terg[[1]], 
                     slot.name = 'net', 
                     net = net.up_Terg, 
                     sources.use = 'epithelium', 
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
# from epith to c3,c4 in HiMYC
######################################################################
png('./figures/LR_byGenotype/up_epith_HiMYC.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_HiMyc[[1]], 
                     slot.name = 'net', 
                     net = net.up_HiMYC, 
                     sources.use = 'epithelium', 
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
# from epith to c5,c6,c7 in PRN
######################################################################
png('./figures/LR_byGenotype/up_epith_PRN.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRN, 
                     sources.use = 'epithelium', 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 0.8, 
                     small.gap = 4, 
                     #big.gap = 20,
                     legend.pos.y = 180, 
                     legend.pos.x = 10, 
                     title.name = "")
dev.off()

######################################################################
# from immune to c3,c4 in PN
######################################################################
png('./figures/LR_byGenotype/up_immune_PN.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_PN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PN, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     thresh = 0.05, 
                     lab.cex = 1,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 5, 
                     title.name = "")
dev.off()

######################################################################
# from immune to c3,c4 in Terg
######################################################################
png('./figures/LR_byGenotype/up_immune_Terg.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_Terg[[1]], 
                     slot.name = 'net', 
                     net = net.up_Terg, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     thresh = 0.05, 
                     lab.cex = 1,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 5, 
                     title.name = "")
dev.off()

######################################################################
# from immune to c3,c4 in HiMYC
######################################################################
png('./figures/LR_byGenotype/up_immune_HiMYC.png',  width = 3500, height = 3000, res = 350)
netVisual_chord_gene(object.list_HiMyc[[1]], 
                     slot.name = 'net', 
                     net = net.up_HiMYC, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c3', 'c4'), 
                     thresh = 0.05, 
                     lab.cex = 1,
                     small.gap = 5, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 5, 
                     title.name = "")
dev.off()

######################################################################
# from immune to c5,c6,c7 in PRN
######################################################################
png('./figures/LR_byGenotype/up_immune_PRN.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRN, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 1, 
                     small.gap = 3, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 5, 
                     title.name = "")
dev.off()



######################################################
# Signaling up in PRN compared to other models
######################################################

# from epithelium

# VS HiMYC
png('./figures/LR_byGenotype/up_epith_PRNvsHiMYC.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN_HiMYC[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRNvshiMYC, 
                     sources.use = 'epithelium', 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 0.8, 
                     small.gap = 4, 
                     #big.gap = 20,
                     legend.pos.y = 180, 
                     legend.pos.x = 10, 
                     title.name = "")
dev.off()

# VS PN
png('./figures/LR_byGenotype/up_epith_PRNvsPN.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN_PN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRNvsPN, 
                     sources.use = 'epithelium', 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 0.8, 
                     small.gap = 4, 
                     #big.gap = 20,
                     legend.pos.y = 180, 
                     legend.pos.x = 10, 
                     title.name = "")
dev.off()


# VS Terg
png('./figures/LR_byGenotype/up_epith_PRNvsTerg.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN_Terg[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRNvsTerg, 
                     sources.use = 'epithelium', 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 0.8, 
                     small.gap = 4, 
                     #big.gap = 20,
                     legend.pos.y = 180, 
                     legend.pos.x = 10, 
                     title.name = "")
dev.off()

################
# from immune

# VS HiMYC
png('./figures/LR_byGenotype/up_immune_PRNvsHiMYC.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN_HiMYC[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRNvshiMYC, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 0.8, 
                     small.gap = 3, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 5, 
                     title.name = "")
dev.off()

# VS PN
png('./figures/LR_byGenotype/up_immune_PRNvsPN.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN_PN[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRNvsPN, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 0.8, 
                     small.gap = 3, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 5, 
                     title.name = "")
dev.off()

# VS Terg
png('./figures/LR_byGenotype/up_immune_PRNvsTerg.png', width = 3500, height = 3000, res = 300)
netVisual_chord_gene(object.list_PRN_Terg[[1]], 
                     slot.name = 'net', 
                     net = net.up_PRNvsTerg, 
                     sources.use = c("B cells", "CD4+ T lymphocytes", "NK/cytoxic T lymphocytes", "Treg", "dendritic cells", "monocytes/macrophages"), 
                     targets.use = c('c5', 'c6', 'c7'), 
                     thresh = 0.05, 
                     lab.cex = 0.8, 
                     small.gap = 3, 
                     #big.gap = 20,
                     legend.pos.y = 165, 
                     legend.pos.x = 5, 
                     title.name = "")
dev.off()


######################################################
# Compare outgoing (or incoming) signaling associated with each cell population
######################################################
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union_PN <- union(object.list_PN[[i]]@netP$pathways, object.list_PN[[i+1]]@netP$pathways)

ht1_PN = netAnalysis_signalingRole_heatmap(object.list_PN[[i]], pattern = "outgoing", signaling = pathway.union_PN, title = 'PN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)
ht2_PN = netAnalysis_signalingRole_heatmap(object.list_PN[[i+1]], pattern = "outgoing", signaling = pathway.union_PN, title = 'WT for PN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)


tiff('./figures/LR_byGenotype/Diff_heatmap_PN.tiff', width = 4000, height = 3000, res = 350)
draw(ht2_PN + ht1_PN, ht_gap = unit(2, "cm"))
dev.off()


######################################################
# Compare outgoing (or incoming) signaling associated with each cell population
######################################################
i = 1

# PRN vs PN
pathway.union_PRN_vs_PN <- union(object.list_PRN_PN[[i]]@netP$pathways, object.list_PRN_PN[[i+1]]@netP$pathways)

ht1_PRN_vs_PN = netAnalysis_signalingRole_heatmap(object.list_PRN_PN[[i]], pattern = "outgoing", signaling = pathway.union_PRN_vs_PN, title = 'PRN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)
ht2_PRN_vs_PN = netAnalysis_signalingRole_heatmap(object.list_PRN_PN[[i+1]], pattern = "outgoing", signaling = pathway.union_PRN_vs_PN, title = 'PN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)


tiff('./figures/LR_byGenotype/Diff_heatmap_PRN_vs_PN.tiff', width = 4000, height = 3000, res = 350)
draw(ht1_PRN_vs_PN + ht2_PRN_vs_PN, ht_gap = unit(2, "cm"))
dev.off()

########
# PRN vs Terg
pathway.union_PRN_vs_Terg <- union(object.list_PRN_Terg[[i]]@netP$pathways, object.list_PRN_Terg[[i+1]]@netP$pathways)

ht1_PRN_vs_Terg = netAnalysis_signalingRole_heatmap(object.list_PRN_Terg[[i]], pattern = "outgoing", signaling = pathway.union_PRN_vs_Terg, title = 'PRN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)
ht2_PRN_vs_Terg = netAnalysis_signalingRole_heatmap(object.list_PRN_Terg[[i+1]], pattern = "outgoing", signaling = pathway.union_PRN_vs_Terg, title = 'Terg', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)


tiff('./figures/LR_byGenotype/Diff_heatmap_PRN_vs_Terg.tiff', width = 4000, height = 3000, res = 350)
draw(ht1_PRN_vs_Terg + ht2_PRN_vs_Terg, ht_gap = unit(2, "cm"))
dev.off()

########
# PRN vs HiMYC
pathway.union_PRN_vs_HiMYC <- union(object.list_PRN_Terg[[i]]@netP$pathways, object.list_PRN_Terg[[i+1]]@netP$pathways)

ht1_PRN_vs_HiMYC = netAnalysis_signalingRole_heatmap(object.list_PRN_HiMYC[[i]], pattern = "outgoing", signaling = pathway.union_PRN_vs_HiMYC, title = 'PRN', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)
ht2_PRN_vs_HiMYC = netAnalysis_signalingRole_heatmap(object.list_PRN_HiMYC[[i+1]], pattern = "outgoing", signaling = pathway.union_PRN_vs_HiMYC, title = 'HiMYC', width = 10, height = 15, color.heatmap = "OrRd", font.size.title = 14)


tiff('./figures/LR_byGenotype/Diff_heatmap_PRN_vs_HiMYC.tiff', width = 4000, height = 3000, res = 350)
draw(ht1_PRN_vs_HiMYC + ht2_PRN_vs_HiMYC, ht_gap = unit(2, "cm"))
dev.off()



############
# save
saveRDS(cellchat, file = "./forCellChat/cellchat_mouse_WT_vs_mutants.rds")

# load
cellchat = readRDS("./forCellChat/cellchat_mouse_WT_vs_mutants.rds")




