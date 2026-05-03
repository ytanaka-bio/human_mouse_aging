library(Seurat)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(scales)

#OL
#human
exp <- OL_human@assays$RNA@data
#cor(as.matrix(t(exp)),exp["OPALIN",]) -> cor_gene
#cor_gene[,1] -> cor_gene
#save(cor_gene,file="cor_gene_human.dat")
#cor_gene -> cor_gene_human
OL_pol_human <- FindMarkers(OL_human,ident.1="12",ident.2="2")

OL_rep_human <- vector("list",2)
names(OL_rep_human) <- 1:2
OL_rep_human[[1]] <- rownames(OL_pol_human[OL_pol_human[,5]<=1e-200 & OL_pol_human[,2] >= 2,])
OL_rep_human[[2]] <- rownames(OL_pol_human[OL_pol_human[,5]<=1e-200 & OL_pol_human[,2] <= -2,])

OL_pol_score <- colMeans(exp[OL_rep_human[[1]],]) - colMeans(exp[OL_rep_human[[2]],])
OL_pol_score <- scale(OL_pol_score)
OL_human@meta.data$OL_pol_score <- OL_pol_score
OL_human$OL_pol_score <- OL_pol_score
cor(t(as.matrix(exp)),OL_human$OL_pol_score)[,1] -> cor_gene
write.table(data.frame(sort(cor_gene,decreasing=T)),"cor_gene_OL_pol_human.rnk",sep="\t",quote=F,col.names=F)
save(cor_gene,file="cor_gene_human.dat")
cor_gene -> cor_gene_human
FeaturePlot(OL_human,features=c("RBFOX1","AFF3","OPALIN","FCHSD2"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_human,features=c("MBP","SLC13A1","MOBP","CNTN3"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))

#mouse
exp <- OL_mouse@assays$RNA@data
#cor(as.matrix(t(exp)),exp["Opalin",]) -> cor_gene
#cor_gene[,1] -> cor_gene
#save(cor_gene,file="cor_gene_mouse.dat")
#cor_gene -> cor_gene_mouse
OL_pol_mouse <- FindMarkers(OL_mouse,ident.1="7",ident.2="10")

OL_rep_mouse <- vector("list",2)
names(OL_rep_mouse) <- 1:2
OL_rep_mouse[[1]] <- rownames(OL_pol_mouse[OL_pol_mouse[,5]<=1e-200 & OL_pol_mouse[,2] >= 2,])
OL_rep_mouse[[2]] <- rownames(OL_pol_mouse[OL_pol_mouse[,5]<=1e-200 & OL_pol_mouse[,2] <= -2,])
OL_pol_score <- colMeans(exp[OL_rep_mouse[[1]],]) - colMeans(exp[OL_rep_mouse[[2]],])
OL_pol_score <- scale(OL_pol_score)
OL_mouse$OL_pol_score <- OL_pol_score
cor(t(as.matrix(exp)),OL_mouse$OL_pol_score)[,1] -> cor_gene
write.table(data.frame(sort(cor_gene,decreasing=T)),"cor_gene_OL_pol_mouse.rnk",sep="\t",quote=F,col.names=F)
samp <- sort(cor_gene,decreasing=T)
write.table(data.frame(GENE=toupper(names(samp)),SCORE=samp),"cor_gene_OL_pol_mouse_humangene.rnk",sep="\t",quote=F,row.names=F,col.names=F)
save(cor_gene,file="cor_gene_mouse.dat")
cor_gene -> cor_gene_mouse
FeaturePlot(OL_mouse,features=c("Rbfox1","Aff3","Opalin","Fchsd2"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_mouse,features=c("Mbp","Slc13a1","Mobp","Cntn3"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))

#cor comp
names(cor_gene_mouse) <- toupper(names(cor_gene_mouse))
gene <- intersect(names(cor_gene_human),names(cor_gene_mouse))
data.frame(Human=cor_gene_human[gene],Mouse=cor_gene_mouse[gene]) -> df
write.table(df,"cor_gene_OL_comp.txt",sep="\t",quote=F)
my.col <- rep("black",nrow(df))
names(my.col) <- rownames(df)
my.col[c("OPALIN","FCHSD2","SLC13A1","RBFOX1","CNTN3","AFF3")] <- "red"
plot(df[,1:2],pch=20,xlim=c(-0.8,0.8),ylim=c(-0.8,0.8),col=my.col)
abline(h=0,v=0,lty=2)

#DEG outlier
DEG_clust14 <- FindMarkers(OL_human,ident.1="14")
DEG_clust14_1p25fold_pval005 <- vector("list",2)
names(DEG_clust14_1p25fold_pval005) <- c("UP","DOWN")
DEG_clust14_1p25fold_pval005[[1]] <- rownames(DEG_clust14)[which(DEG_clust14[,1] < 0.05 & DEG_clust14[,2] > log2(1.25))]
DEG_clust14_1p25fold_pval005[[2]] <- rownames(DEG_clust14)[which(DEG_clust14[,1] < 0.05 & DEG_clust14[,2] < -log2(1.25))]

DEG_clust14_1p25fold_pval005 -> DEG_clust14_1p25fold_pval005_GO
source("~/R_command/all_go_analysis_v2_symbol2.R")
source("~/R_command/all_GO_heat2.R")
for(i in 1:2){ 
all_go_analysis_v2_symbol2(DEG_clust14_1p25fold_pval005[[i]],"human") -> DEG_clust14_1p25fold_pval005_GO[[i]]
}
 
all_GO_heat2(DEG_clust14_1p25fold_pval005,1) -> samp
write.table(samp,"DEG_clust14_1p25fold_pval005_GO.txt",sep="\t",quote=F)

subtype <- rep(1,ncol(OL_human))
subtype[which(OL_human@active.ident == 14)] <- 2
OL_human@meta.data$SubType <- factor(subtype)
FeaturePlot(OL_human,features=c("HSPB1","HSP90AA1","HSPA6","HSPA1B"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))

df <- OL_human@meta.data
df2 <- OL_mouse@meta.data

ggplot(df,aes(OL_pol_score,fill=Type)) + geom_density(alpha=0.5) + theme_bw() + xlim(c(-3,3)) + scale_fill_manual( values = hue_pal()(3)[c(2,1)])
ggplot(df2,aes(OL_pol_score,fill=Type)) + geom_density(alpha=0.5) + theme_bw() + xlim(c(-3,3)) + scale_fill_manual( values = hue_pal()(3)[c(3,2,1)])

#MG
#human
MG_human_group <- rep(1,ncol(MG_human))
names(MG_human_group) <- colnames(MG_human)
MG_human_group[which(MG_human@active.ident == 5)] <- 2
MG_human_group[which(MG_human@active.ident == 7)] <- 2
MG_human_group[which(MG_human@active.ident == 9)] <- 2
MG_human_group[which(MG_human@active.ident == 1)] <- 3
MG_human_group[which(MG_human@active.ident == 12)] <- 4
MG_human@active.ident <- factor(MG_human_group)
MG_dif_human <- FindAllMarkers(MG_human)
MG_DEG_human <- vector("list",4)
names(MG_DEG_human) <- 1:4
for(i in 1:4){
 MG_dif_human[which(MG_dif_human$avg_log2FC > 1 & MG_dif_human$p_val_adj < 1e-20 & MG_dif_human$cluster == i),7] -> MG_DEG_human[[i]]
}

#mouse
MG_mouse_group <- rep(1,ncol(MG_mouse))
names(MG_mouse_group) <- colnames(MG_mouse)
MG_mouse_group[which(MG_mouse@active.ident == 1)] <- 2
MG_mouse_group[which(MG_mouse@active.ident == 4)] <- 2
MG_mouse_group[which(MG_mouse@active.ident == 5)] <- 2
MG_mouse_group[which(MG_mouse@active.ident == 6)] <- 3
MG_mouse_group[which(MG_mouse@active.ident == 7)] <- 4
MG_mouse@active.ident <- factor(MG_mouse_group)
MG_dif_mouse <- FindAllMarkers(MG_mouse)
MG_DEG_mouse <- vector("list",4)
names(MG_DEG_mouse) <- 1:4
for(i in 1:4){
 MG_dif_mouse[which(MG_dif_mouse$avg_log2FC > 1 & MG_dif_mouse$p_val_adj < 1e-20 & MG_dif_mouse$cluster == i),7] -> MG_DEG_mouse[[i]]
}
MG_DEG_mouse -> MG_DEG_mouse_UL
for(i in 1:4){
 MG_DEG_mouse_UL[[i]] <- toupper(MG_DEG_mouse[[i]])
}
for(i in 1:4){
for(j in 1:4){
 print(c(i,j))
 print(length(intersect(MG_DEG_human[[i]],MG_DEG_mouse_UL[[j]])))
}
}

MG_DEG_human_GO <- MG_DEG_human
MG_DEG_mouse_GO <- MG_DEG_mouse
for(i in 1:4){ 
   all_go_analysis_v2_symbol2(MG_DEG_human[[i]],"human") -> MG_DEG_human_GO[[i]]
   all_go_analysis_v2_symbol2(MG_DEG_mouse[[i]],"mouse") -> MG_DEG_mouse_GO[[i]]
}
 
all_GO_heat2(MG_DEG_human_GO,1) -> samp
write.table(samp,"MG_DEG_human_GO.txt",sep="\t",quote=F)
all_GO_heat2(MG_DEG_mouse_GO,1) -> samp
write.table(samp,"MG_DEG_mouse_GO.txt",sep="\t",quote=F)

#OPC
#human
OPC_human_group <- rep(1,ncol(OPC_human))
names(OPC_human_group) <- colnames(OPC_human)
OPC_human_group[which(OPC_human@active.ident == 9)] <- 1
OPC_human_group[which(OPC_human@active.ident == 1)] <- 1
OPC_human_group[which(OPC_human@active.ident == 5)] <- 1
OPC_human_group[which(OPC_human@active.ident == 8)] <- 1
OPC_human_group[which(OPC_human@active.ident == 10)] <- 1
OPC_human_group[which(OPC_human@active.ident == 4)] <- 2
OPC_human_group[which(OPC_human@active.ident == 12)] <- 3
OPC_human_group[which(OPC_human@active.ident == 11)] <- 4
OPC_human_group[which(OPC_human@active.ident == 13)] <- 5
OPC_human@active.ident <- factor(OPC_human_group)
OPC_dif_human <- FindAllMarkers(OPC_human)
OPC_DEG_human <- vector("list",5)
names(OPC_DEG_human) <- 1:5
for(i in 1:5){
 OPC_dif_human[which(OPC_dif_human$avg_log2FC > 1 & OPC_dif_human$p_val_adj < 0.05 & OPC_dif_human$cluster == i),7] -> OPC_DEG_human[[i]]
}

#mouse
OPC_mouse_group <- rep(1,ncol(OPC_mouse))
names(OPC_mouse_group) <- colnames(OPC_mouse)
OPC_mouse_group[which(OPC_mouse@active.ident == 0)] <- 1
OPC_mouse_group[which(OPC_mouse@active.ident == 1)] <- 1
OPC_mouse_group[which(OPC_mouse@active.ident == 2)] <- 2
OPC_mouse_group[which(OPC_mouse@active.ident == 3)] <- 1
OPC_mouse@active.ident <- factor(OPC_mouse_group)
OPC_dif_mouse <- FindAllMarkers(OPC_mouse)
OPC_DEG_mouse <- vector("list",2)
names(OPC_DEG_mouse) <- 1:2
for(i in 1:2){
 OPC_dif_mouse[which(OPC_dif_mouse$avg_log2FC > 1 & OPC_dif_mouse$p_val_adj < 0.05 & OPC_dif_mouse$cluster == i),7] -> OPC_DEG_mouse[[i]]
}
OPC_DEG_mouse -> OPC_DEG_mouse_UL
for(i in 1:2){
 OPC_DEG_mouse_UL[[i]] <- toupper(OPC_DEG_mouse[[i]])
}
for(i in 1:5){
for(j in 1:2){
 print(c(i,j))
 print(length(intersect(OPC_DEG_human[[i]],OPC_DEG_mouse_UL[[j]])))
}
}

OPC_DEG_human_GO <- OPC_DEG_human
OPC_DEG_mouse_GO <- OPC_DEG_mouse
for(i in 1:5){ 
   all_go_analysis_v2_symbol2(OPC_DEG_human[[i]],"human") -> OPC_DEG_human_GO[[i]]
}
for(i in 1:2){
   all_go_analysis_v2_symbol2(OPC_DEG_mouse[[i]],"mouse") -> OPC_DEG_mouse_GO[[i]]
}
 
all_GO_heat2(OPC_DEG_human_GO,1) -> samp
write.table(samp,"OPC_DEG_human_GO.txt",sep="\t",quote=F)
all_GO_heat2(OPC_DEG_mouse_GO,1) -> samp
write.table(samp,"OPC_DEG_mouse_GO.txt",sep="\t",quote=F)


#plot
FeaturePlot(OPC_human,c("PDGFRA","OLIG2","MAG","MOBP"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OPC_mouse,c("Pdgfra","Olig2","Mag","Mobp"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
DimPlot(OPC_human,cols=c(hue_pal()(2),"gray","gray","gray"))
DimPlot(OPC_mouse,cols=c(hue_pal()(2)))

#AS
#human
AS_human_group <- rep(1,ncol(AS_human))
names(AS_human_group) <- colnames(AS_human)
AS_human_group[which(AS_human@active.ident == 2)] <- 2
AS_human_group[which(AS_human@active.ident == 5)] <- 2
AS_human_group[which(AS_human@active.ident == 8)] <- 2
AS_human_group[which(AS_human@active.ident == 13)] <- 2
AS_human_group[which(AS_human@active.ident == 12)] <- 3
AS_human@active.ident <- factor(AS_human_group)
AS_dif_human <- FindAllMarkers(AS_human)
AS_DEG_human <- vector("list",3)
names(AS_DEG_human) <- 1:3
for(i in 1:3){
 AS_dif_human[which(AS_dif_human$avg_log2FC > 1 & AS_dif_human$p_val_adj < 0.05 & AS_dif_human$cluster == i),7] -> AS_DEG_human[[i]]
}

#mouse
AS_mouse_group <- rep(1,ncol(AS_mouse))
names(AS_mouse_group) <- colnames(AS_mouse)
AS_mouse_group[which(AS_mouse@active.ident == 4)] <- 2
AS_mouse_group[which(AS_mouse@active.ident == 6)] <- 3
AS_mouse@active.ident <- factor(AS_mouse_group)
AS_dif_mouse <- FindAllMarkers(AS_mouse)
AS_DEG_mouse <- vector("list",3)
names(AS_DEG_mouse) <- 1:3
for(i in 1:3){
 AS_dif_mouse[which(AS_dif_mouse$avg_log2FC > 1 & AS_dif_mouse$p_val_adj < 0.05 & AS_dif_mouse$cluster == i),7] -> AS_DEG_mouse[[i]]
}
AS_DEG_mouse -> AS_DEG_mouse_UL
for(i in 1:3){
 AS_DEG_mouse_UL[[i]] <- toupper(AS_DEG_mouse[[i]])
}
for(i in 1:3){
for(j in 1:3){
 print(c(i,j))
 print(length(intersect(AS_DEG_human[[i]],AS_DEG_mouse_UL[[j]])))
}
}

AS_DEG_human_GO <- AS_DEG_human
AS_DEG_mouse_GO <- AS_DEG_mouse
for(i in 1:3){ 
   all_go_analysis_v2_symbol2(AS_DEG_human[[i]],"human") -> AS_DEG_human_GO[[i]]
}
for(i in 1:3){
   all_go_analysis_v2_symbol2(AS_DEG_mouse[[i]],"mouse") -> AS_DEG_mouse_GO[[i]]
}
 
all_GO_heat2(AS_DEG_human_GO,1) -> samp
write.table(samp,"AS_DEG_human_GO.txt",sep="\t",quote=F)
all_GO_heat2(AS_DEG_mouse_GO,1) -> samp
write.table(samp,"AS_DEG_mouse_GO.txt",sep="\t",quote=F)

exp <- AS_human@assays$RNA$data
AS_human_ratio <- rowMeans(exp[,which(AS_human@active.ident == 1)]) - rowMeans(exp[,which(AS_human@active.ident == 2)])
exp <- AS_mouse@assays$RNA$data
AS_mouse_ratio <- rowMeans(exp[,which(AS_mouse@active.ident == 1)]) - rowMeans(exp[,which(AS_mouse@active.ident == 2)])

write.table(data.frame(sort(AS_human_ratio,decreasing=T)),"AS_human_ratio.rnk",sep="\t",quote=F,col.names=F)
#names(AS_mouse_ratio) <- toupper(names(AS_mouse_ratio))
write.table(data.frame(sort(AS_mouse_ratio,decreasing=T)),"AS_mouse_ratio.rnk",sep="\t",quote=F,col.names=F)

DimPlot(AS_human,cols=c(hue_pal()(2),"gray"))
DimPlot(AS_mouse,cols=c(hue_pal()(2),"gray"))