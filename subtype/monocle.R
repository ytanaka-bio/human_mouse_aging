library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
source("~/R_commands/box_strip_v2_large.R")

#human
load("OL_human.dat")
DimPlot(OL_human,group.by="Type") & scale_color_manual(values=hue_pal()(3)[c(2,1)])

cds_OL_human <- as.cell_data_set(OL_human)
cds_OL_human <- cluster_cells(cds_OL_human,reduction_method="UMAP")
cds_OL_human <- learn_graph(cds_OL_human)

#cells <- colnames(OL_human)[which(OL_human@active.ident == 9)]
#cds_OL_human <- order_cells(cds_OL_human,reduction_method="UMAP",root_cells=cells)
cells <- colnames(OL_human)[which(OL_human@active.ident == 5)]
cds_OL_human <- order_cells(cds_OL_human,reduction_method="UMAP",root_cells=cells)
samp <- pseudotime(cds_OL_human)
samp[is.infinite(samp)] <- NA
OL_human@meta.data$Pseudotime <- samp
FeaturePlot(OL_human,"Pseudotime",min.cutoff=0,max.cutoff=8) + scale_color_viridis()

data <- read.table("/Volumes/Backup_atlas/snorkel_v16/weasel/github/human_aging_train_wsl2_result.csv",sep=",",header=T,row.names=1)
OL_human@meta.data$scIDST <- data[colnames(OL_human),1]
FeaturePlot(OL_human,"scIDST",min.cutoff=0,max.cutoff=1) + scale_color_viridis(option="magma")

young <- which(OL_human$Type == "Young")
old <- which(OL_human$Type == "Old")
vec <- vector("list",2)
vec[[1]] <- OL_human$Pseudotime[young]
vec[[2]] <- OL_human$Pseudotime[old]
wilcox.test(vec[[1]],vec[[2]])
box_strip_v2_large(vec,c(0,8),hue_pal()(3)[c(1,2)])

young <- which(OL_human$Type == "Young")
old <- which(OL_human$Type == "Old")
vec <- vector("list",2)
vec[[1]] <- OL_human$scIDST[young]
vec[[2]] <- OL_human$scIDST[old]
wilcox.test(vec[[1]],vec[[2]])
box_strip_v2_large(vec,c(0,1),hue_pal()(3)[c(1,2)])

#mouse
load("OL_mouse.dat")
DimPlot(OL_mouse,group.by="Type") & scale_color_manual(values=hue_pal()(3)[c(3,2,1)])

cds_OL_mouse <- as.cell_data_set(OL_mouse)
cds_OL_mouse <- cluster_cells(cds_OL_mouse,reduction_method="UMAP")
cds_OL_mouse <- learn_graph(cds_OL_mouse)

cells <- colnames(OL_mouse)[which(OL_mouse@active.ident == 3)]
cds_OL_mouse <- order_cells(cds_OL_mouse,reduction_method="UMAP",root_cells=cells)
samp <- pseudotime(cds_OL_mouse)
samp[is.infinite(samp)] <- NA
OL_mouse@meta.data$Pseudotime <- samp
FeaturePlot(OL_mouse,"Pseudotime",min.cutoff=0,max.cutoff=8) + scale_color_viridis()
data <- read.table("/Volumes/Backup_atlas/snorkel_v16/weasel/github/mouse_aging_train_wsl3_result.csv",sep=",",header=T,row.names=1)
data2 <- read.table("/Volumes/Backup_atlas/snorkel_v16/weasel/github/mouse_aging_class_middle_predict_wsl3.csv",sep=",",header=T,row.names=1)
data <- rbind(data,data2)
data <- data[colnames(OL_mouse),]
data[,1] <- (data[,1] - min(data[,1]))/(max(data[,1]) - min(data[,1]))
OL_mouse@meta.data$scIDST <- data[colnames(OL_mouse),1]
FeaturePlot(OL_mouse,"scIDST",min.cutoff=0,max.cutoff=1) + scale_color_viridis(option="magma")

young <- which(OL_mouse$Type == "Young")
middle <- which(OL_mouse$Type == "Middle")
old <- which(OL_mouse$Type == "Old")
vec <- vector("list",2)
vec[[1]] <- OL_mouse$Pseudotime[young]
vec[[2]] <- OL_mouse$Pseudotime[middle]
vec[[3]] <- OL_mouse$Pseudotime[old]
wilcox.test(vec[[1]],vec[[2]])
box_strip_v2_large(vec,c(0,8),hue_pal()(3)[c(1,3,2)])

young <- which(OL_mouse$Type == "Young")
middle <- which(OL_mouse$Type == "Middle")
old <- which(OL_mouse$Type == "Old")
vec <- vector("list",2)
vec[[1]] <- OL_mouse$scIDST[young]
vec[[2]] <- OL_mouse$scIDST[middle]
vec[[3]] <- OL_mouse$scIDST[old]
wilcox.test(vec[[1]],vec[[2]])
box_strip_v2_large(vec,c(0,1),hue_pal()(3)[c(1,3,2)])


save(OL_human,file="OL_human_v2.dat")
save(OL_mouse,file="OL_mouse_v2.dat")

cor(t(as.matrix(OL_human@assays$RNA@data[,samp])),OL_human$Pseudotime[samp])[,1] -> cor_pseudotime_human
cor(t(as.matrix(OL_human@assays$RNA@data)),OL_human$scIDST)[,1] -> cor_scIDST_human
cor(t(as.matrix(OL_mouse@assays$RNA@data[,samp])),OL_mouse$Pseudotime[samp])[,1] -> cor_pseudotime_mouse
cor(t(as.matrix(OL_mouse@assays$RNA@data)),OL_mouse$scIDST)[,1] -> cor_scIDST_mouse

mart <- read.table("mart_export.txt",sep="\t",header=F)
mart <- mart[which(mart[,2]!=""),]
mart <- mart[which(mart[,3]!=""),]
mart <- mart[which(mart[,4]!=""),]

mart <- mart[duplicated(mart[,2])==F,]
mart <- mart[duplicated(mart[,4])==F,]

df <- data.frame(cor_scIDST_human[mart[,2]],cor_scIDST_mouse[mart[,4]])

samp <- cor_scIDST_human
samp[is.na(samp) ==F] -> samp
write.table(sort(samp,decreasing=T),"cor_scIDST_human.rnk",sep="\t",quote=F,col.names=F)

samp <- cor_scIDST_mouse
samp[is.na(samp) ==F] -> samp
write.table(sort(samp,decreasing=T),"cor_scIDST_mouse.rnk",sep="\t",quote=F,col.names=F)

samp <- sort(cor_scIDST_human)
plot(rank(samp),samp,pch=20,ylim=c(-0.5,0.5),xlim=c(0,length(samp)))

abline(h=0,lty=2)
select <- c("STON2","SYNJ2","RNF13","RASSF2","MOBP","MOG","MALAT1","NEAT1","CLU","FKBP5","HSP90AA1","APOE","CST3")
par(new=T)
plot(rank(samp)[select],samp[select],pch=19,ylim=c(-0.5,0.5),xlim=c(0,length(samp)),col="red",ann=F,axes=F)

samp <- sort(cor_scIDST_mouse)
plot(rank(samp),samp,pch=20,ylim=c(-0.3,0.3),xlim=c(0,length(samp)))
abline(h=0,lty=2)
select <- c("Ston2","Synj2","Rnf13","Rassf2","Mobp","Mog","Malat1","Neat1","Clu","Fkbp5","Hsp90aa1","Apoe","Cst3")
par(new=T)
plot(rank(samp)[select],samp[select],pch=19,ylim=c(-0.3,0.3),xlim=c(0,length(samp)),col="red",ann=F,axes=F)


source("~/R_commands/plot_circle_posneg.R")
samp <- read.table("GO_aging.txt",sep="\t",header=T,row.names=1)
plot_circle_posneg(samp,3,bluered(128))