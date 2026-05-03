library(CellChat)
library(patchwork)
library(dplyr)
options(stringsAsFactors = FALSE)

#load data
mouse <- readRDS("C:/Users/annub/Downloads/Aging_YT/mouse/mouse.rds")
Idents(mouse) <- as.factor(mouse$CellType)

Young<- subset(x=mouse, subset = Type == "Young")

data.input <- Young[["RNA"]]$data # normalized data matrix`
labels <- Idents(Young)
meta <- data.frame(labels = labels, row.names = names(labels))

cellchat <- createCellChat(object = Young, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

#use secreted signaling and cell-cell contact 
CellChatDB.use <- subsetDB(CellChatDB,  search = c("Secreted Signaling", "Cell-Cell Contact"), key = "annotation") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
options(future.globals.maxSize = 5 * 1024^3)  # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)

# different diagrams for pathways
pathways.show <- c("APP") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

netAnalysis_contribution(cellchat, signaling = pathways.show)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Dominant sender and receiver
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 20, height = 12)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 20, height = 12)
ht1 + ht2

library(NMF)
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster

library(ggalluvial)
#outgoing
selectK(cellchat, pattern = "outgoing")

nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 15, height = 11)

#river plot
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

#incoming
selectK(cellchat, pattern = "incoming")

nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, width = 15, height = 11)

netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchat_mouse_young.rds")

