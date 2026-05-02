library(Seurat)
library(cowplot)
library(magrittr)
library(dplyr)
library(ggplot2)

library(SeuratWrappers)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

cellchat_human_young_ss <- readRDS("C:/Users/annub/Downloads/Aging/human_cellchat/cellchat_human_young.rds")
cellchat_human_old_ss <- readRDS("C:/Users/annub/Downloads/Aging/human_cellchat/cellchat_human_old.rds")

object.list <- list(HY = cellchat_human_young, HO =  cellchat_human_old )
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Define the correct order of cell types
cellchat@idents$joint <- factor(
  cellchat@idents$joint, 
  levels = c("NE", "OL", "OPC", "MI", "AS", "EC")  # Your desired order
)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, label.edge = FALSE, vertex.label.cex = 0.02)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = FALSE, vertex.label.cex = 0.02)

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)



#identify signaling groups based on fuctional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2


#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)


rankSimilarity(cellchat, type = "functional")

#compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2


library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
custom_order <- c("Neuron", "OPC", "Oligodendrocyte", "Microglia", "Astrocyte", "Endothelial cells")  # Replace with your desired order
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 14, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 14, height = 16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 14, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 14, height = 16, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 14, height = 16, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 14, height = 16, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)

save(object.list, file = "cellchat_object.list_pos_neg_ss.RData")
save(cellchat, file = "cellchat_merged_pos_neg_ss.RData")
