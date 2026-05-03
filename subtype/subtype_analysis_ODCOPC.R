library(Seurat)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(SeuratWrappers)
library(monocle3)

#OL_OPC
#human
load("OL_OPC_human.dat")
OL_OPC_human <- FindClusters(OL_OPC_human,resolution=1)

FeaturePlot(OL_OPC_human,c("HSPB1","HSPA1B","OLIG1","OLIG2"),min.cutoff=0,max.cutoff=1,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("PDGFRA","CSPG4","SOX10","OLIG3"),min.cutoff=0,max.cutoff=1,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("MBP","MOBP","OPALIN","FCHSD2"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("SLC13A1","RBFOX1","CNTN3","AFF3"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("FTH1","PTGDS","MALAT1","NKAIN2"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("FTH1","PTGDS","MALAT1","NKAIN2"),min.cutoff=2,max.cutoff=4,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("FTH1","PTGDS","MALAT1","NKAIN2"),min.cutoff=4.5,max.cutoff=6,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("PCDH9","LRP1B","RPL10","EEF1A1"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("PCDH9","LRP1B","RPL10","EEF1A1"),min.cutoff=1.5,max.cutoff=3.5,cols=c("gray","red"))
FeaturePlot(OL_OPC_human,c("PCDH9","LRP1B","RPL10","EEF1A1"),min.cutoff=3.5,max.cutoff=5.5,cols=c("gray","red"))

group <- rep(1,ncol(OL_OPC_human))
names(group) <- colnames(OL_OPC_human)
group[which(OL_OPC_human@active.ident == 20)] <- 3
group[which(OL_OPC_human@active.ident == 23)] <- 2
group[which(OL_OPC_human@active.ident == 21)] <- 2
group[which(OL_OPC_human@active.ident == 11)] <- 2
group[which(OL_OPC_human@active.ident == 15)] <- 2
DimPlot(OL_OPC_human,label=T) -> p1
plot(p1$data[,1:2],pch=".",col=hue_pal()(3)[group])

OL_OPC_human@meta.data$Group <- group

DEG_human_c14_c4 <- FindMarkers(OL_OPC_human,ident.1=14,ident.2=4)
DEG_human_c3_c9 <- FindMarkers(OL_OPC_human,ident.1=3,ident.2=9)

OL_rep_human_1 <- vector("list",2)
names(OL_rep_human_1) <- 1:2 
OL_rep_human_1[[1]] <- rownames(DEG_human_c14_c4[DEG_human_c14_c4[,5]<=1e-50 & DEG_human_c14_c4[,2] >= 1,])
OL_rep_human_1[[2]] <- rownames(DEG_human_c14_c4[DEG_human_c14_c4[,5]<=1e-50 & DEG_human_c14_c4[,2] <= -1,])

OL_rep_human_2 <- vector("list",2)
names(OL_rep_human_2) <- 1:2
OL_rep_human_2[[1]] <- rownames(DEG_human_c3_c9[DEG_human_c3_c9[,5]<=1e-50 & DEG_human_c3_c9[,2] >= 1,])
OL_rep_human_2[[2]] <- rownames(DEG_human_c3_c9[DEG_human_c3_c9[,5]<=1e-50 & DEG_human_c3_c9[,2] <= -1,])

exp <- OL_OPC_human@assays$RNA@data
OL_score_1 <- colMeans(exp[OL_rep_human_1[[1]],]) - colMeans(exp[OL_rep_human_1[[2]],])
OL_score_2 <- colMeans(exp[OL_rep_human_2[[1]],]) - colMeans(exp[OL_rep_human_2[[2]],])
OL_score_1 <- scale(OL_score_1)
OL_score_2 <- scale(OL_score_2)
OL_OPC_human@meta.data$OL_score_1 <- OL_score_1
OL_OPC_human@meta.data$OL_score_2 <- OL_score_2
FeaturePlot(OL_OPC_human,"OL_score_1",min.cutoff=-2,max.cutoff=2) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PRGn")))
FeaturePlot(OL_OPC_human,"OL_score_2",min.cutoff=-2,max.cutoff=2) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

cor(t(as.matrix(exp)),OL_OPC_human$OL_score_1)[,1] -> cor_gene_human_1
cor(t(as.matrix(exp)),OL_OPC_human$OL_score_2)[,1] -> cor_gene_human_2
save(cor_gene_human_1,file="cor_gene_human_1.dat")
save(cor_gene_human_2,file="cor_gene_human_2.dat")
write.table(data.frame(sort(cor_gene_human_1,decreasing=T)),"cor_gene_OL_OPC_human1.rnk",sep="\t",quote=F,col.names=F)
write.table(data.frame(sort(cor_gene_human_2,decreasing=T)),"cor_gene_OL_OPC_human2.rnk",sep="\t",quote=F,col.names=F)

#Heatshock group
HS_OPC <- FindMarkers(OL_OPC_human,ident.1=3,ident.2=1,group.by="Group")
dif_gene_2fold_pval005_HS_OPC <- vector("list",2)
names(dif_gene_2fold_pval005_HS_OPC) <- c("up","down")
dif_gene_2fold_pval005_HS_OPC[[1]] <- rownames(HS_OPC)[which(HS_OPC[,2] > 1 & HS_OPC[,5] < 0.05)]
dif_gene_2fold_pval005_HS_OPC[[2]] <- rownames(HS_OPC)[which(HS_OPC[,2] < -1 & HS_OPC[,5] < 0.05)]
source("~/R_command/all_GO_heat2.R")
source("~/R_command/all_go_analysis_v2_symbol2.R")
dif_gene_2fold_pval005_HS_OPC_GO <- dif_gene_2fold_pval005_HS_OPC
for(i in 1:2){
 dif_gene_2fold_pval005_HS_OPC_GO[[i]] <- all_go_analysis_v2_symbol2(dif_gene_2fold_pval005_HS_OPC[[i]],"human")
}
all_GO_heat2(dif_gene_2fold_pval005_HS_OPC_GO,1) -> samp
write.table(samp,"dif_gene_2fold_pval005_HS_OPC_GO.txt",sep="\t",quote=F)

summary(factor(OL_OPC_human$Type[which(group == 1)]))
summary(factor(OL_OPC_human$Type[which(group == 2)]))
summary(factor(OL_OPC_human$Type[which(group == 3)]))

#mouse
load("OL_OPC_mouse.dat")
OL_OPC_mouse <- FindClusters(OL_OPC_mouse,resolution=1)

FeaturePlot(OL_OPC_mouse,c("Hspb1","Hspa1b","Olig1","Olig2"),min.cutoff=0,max.cutoff=1,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Pdgfra","Cspg4","Sox10","Olig3"),min.cutoff=0,max.cutoff=1,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Mbp","Mobp","Opalin","Fchsd2"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Slc13a1","Rbfox1","Cntn3","Aff3"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Fth1","Ptgds","Malat1","Nkain2"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Fth1","Ptgds","Malat1","Nkain2"),min.cutoff=2,max.cutoff=4,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Fth1","Ptgds","Malat1","Nkain2"),min.cutoff=4.5,max.cutoff=6,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Pcdh9","Lrp1b","Rpl10","Eef1a1"),min.cutoff=0,max.cutoff=2,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Pcdh9","Lrp1b","Rpl10","Eef1a1"),min.cutoff=1.5,max.cutoff=3.5,cols=c("gray","red"))
FeaturePlot(OL_OPC_mouse,c("Pcdh9","Lrp1b","Rpl10","Eef1a1"),min.cutoff=2,max.cutoff=4,cols=c("gray","red"))

group <- rep(1,ncol(OL_OPC_mouse))
names(group) <- colnames(OL_OPC_mouse)
group[which(OL_OPC_mouse@active.ident == 16)] <- 2
DimPlot(OL_OPC_mouse,label=T) -> p1
plot(p1$data[,1:2],pch=".",col=hue_pal()(3)[group])

OL_OPC_mouse@meta.data$Group <- group

DEG_mouse_c15_c4 <- FindMarkers(OL_OPC_mouse,ident.1=15,ident.2=4)
DEG_mouse_c0_c7 <- FindMarkers(OL_OPC_mouse,ident.1=0,ident.2=7)

OL_rep_mouse_2 <- vector("list",2)
names(OL_rep_mouse_2) <- 1:2
OL_rep_mouse_2[[1]] <- rownames(DEG_mouse_c15_c4[DEG_mouse_c15_c4[,5]<=1e-50 & DEG_mouse_c15_c4[,2] >= 1,])
OL_rep_mouse_2[[2]] <- rownames(DEG_mouse_c15_c4[DEG_mouse_c15_c4[,5]<=1e-50 & DEG_mouse_c15_c4[,2] <= -1,])

OL_rep_mouse_1 <- vector("list",2)
names(OL_rep_mouse_1) <- 1:2
OL_rep_mouse_1[[1]] <- rownames(DEG_mouse_c0_c7[DEG_mouse_c0_c7[,5]<=1e-50 & DEG_mouse_c0_c7[,2] >= 1,])
OL_rep_mouse_1[[2]] <- rownames(DEG_mouse_c0_c7[DEG_mouse_c0_c7[,5]<=1e-50 & DEG_mouse_c0_c7[,2] <= -1,])

exp <- OL_OPC_mouse@assays$RNA@data
OL_score_1 <-  -colMeans(exp[OL_rep_mouse_1[[2]],]) #no up gene
OL_score_2 <- colMeans(exp[OL_rep_mouse_2[[2]],]) - colMeans(exp[OL_rep_mouse_2[[1]],])
OL_score_1 <- scale(OL_score_1)
OL_score_2 <- scale(OL_score_2)
OL_OPC_mouse@meta.data$OL_score_1 <- OL_score_1
OL_OPC_mouse@meta.data$OL_score_2 <- OL_score_2
FeaturePlot(OL_OPC_mouse,"OL_score_1",min.cutoff=-2,max.cutoff=2) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PRGn")))
FeaturePlot(OL_OPC_mouse,"OL_score_2",min.cutoff=-2,max.cutoff=2) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

cor(t(as.matrix(exp)),OL_OPC_mouse$OL_score_1)[,1] -> cor_gene_mouse_1
cor(t(as.matrix(exp)),OL_OPC_mouse$OL_score_2)[,1] -> cor_gene_mouse_2
save(cor_gene_mouse_1,file="cor_gene_mouse_1.dat")
save(cor_gene_mouse_2,file="cor_gene_mouse_2.dat")
write.table(data.frame(sort(cor_gene_mouse_1,decreasing=T)),"cor_gene_OL_OPC_mouse1.rnk",sep="\t",quote=F,col.names=F)
write.table(data.frame(sort(cor_gene_mouse_2,decreasing=T)),"cor_gene_OL_OPC_mouse2.rnk",sep="\t",quote=F,col.names=F)
samp <- sort(cor_gene_mouse_1,decreasing=T)
write.table(data.frame(GENE=toupper(names(samp)),SCORE=samp),"cor_gene_OL_OPC_mouse1_humangene.rnk",sep="\t",quote=F,row.names=F,col.names=F)
samp <- sort(cor_gene_mouse_2,decreasing=T)
write.table(data.frame(GENE=toupper(names(samp)),SCORE=samp),"cor_gene_OL_OPC_mouse2_humangene.rnk",sep="\t",quote=F,row.names=F,col.names=F)

#compare
names(cor_gene_mouse_1) <- toupper(names(cor_gene_mouse_1))
names(cor_gene_mouse_2) <- toupper(names(cor_gene_mouse_2))
gene <- intersect(names(cor_gene_human_1),names(cor_gene_mouse_1))
data.frame(Human=-1*cor_gene_human_1[gene],Mouse=-1*cor_gene_mouse_1[gene]) -> df
write.table(df,"cor_gene_OL_OPC_comp_1.txt",sep="\t",quote=F)

my.col <- rep("black",nrow(df))
names(my.col) <- rownames(df)
my.col[c("FTH1","LRP1B","NKAIN2","PCDH9")] <- "red"
plot(df[,1:2],pch=20,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6),col=my.col)
abline(h=0,v=0,lty=2)

gene <- intersect(names(cor_gene_human_2),names(cor_gene_mouse_2))
data.frame(Human=cor_gene_human_2[gene],Mouse=cor_gene_mouse_2[gene]) -> df
write.table(df,"cor_gene_OL_OPC_comp_2.txt",sep="\t",quote=F)

my.col <- rep("black",nrow(df))
names(my.col) <- rownames(df)
my.col[c("OPALIN","FCHSD2","SLC13A1","RBFOX1","CNTN3","AFF3")] <- "red"
plot(df[,1:2],pch=20,xlim=c(-0.7,0.7),ylim=c(-0.7,0.7),col=my.col)
abline(h=0,v=0,lty=2)

#compare two scores
library(viridis)
library(gplots)
library(fields)

df <- data.frame(Maturation=OL_OPC_human$OL_score_1,Polarization=OL_OPC_human$OL_score_2)
df[df > 2] <- 2
df[df < -2] <- -2
h2d <- hist2d(df,nbin=150,show=F)
samp <- log10(h2d$counts+1)
samp[samp > 1] <- 1
image.plot(h2d$x,h2d$y,samp,col=c("white",viridis(127)))

df <- data.frame(Maturation=OL_OPC_mouse$OL_score_1,Polarization=OL_OPC_mouse$OL_score_2)
df[df > 2] <- 2
df[df < -2] <- -2
h2d <- hist2d(df,nbin=150,show=F)
samp <-	log10(h2d$counts+1)
samp[samp > 1] <- 1
image.plot(h2d$x,h2d$y,samp,col=c("white",viridis(127)))


#OL_OPC_human.cds <- as.cell_data_set(OL_OPC_human)
#OL_OPC_human.cds <- cluster_cells(OL_OPC_human.cds,reduction="UMAP")
#OL_OPC_human.cds <- learn_graph(OL_OPC_human.cds,use_partition=F)
#root_cells <- colnames(OL_OPC_human)[which(OL_OPC_human@active.ident == 8)]
#OL_OPC_human.cds <- order_cells(OL_OPC_human.cds,reduction_method="UMAP",root_cells=root_cells)

source("~/R_command/box_strip_v2.R")
vec <- vector("list",3)
vec[[1]] <- OL_OPC_mouse$OL_score_1[which(OL_OPC_mouse$Type == "Young" & OL_OPC_mouse$Group == 1)]
vec[[2]] <- OL_OPC_mouse$OL_score_1[which(OL_OPC_mouse$Type == "Middle" & OL_OPC_mouse$Group == 1)]
vec[[3]] <- OL_OPC_mouse$OL_score_1[which(OL_OPC_mouse$Type == "Old" & OL_OPC_mouse$Group == 1)]
box_strip_v2(vec,c(-3,3),col=hue_pal()(3)[c(1,3,2)])
abline(h=0,lty=2)

vec <- vector("list",3)
vec[[1]] <- OL_OPC_human$OL_score_2[which(OL_OPC_human$Type == "Young" & OL_OPC_human$Group == 1)]
vec[[2]] <- OL_OPC_human$OL_score_2[which(OL_OPC_human$Type == "Old" & OL_OPC_human$Group == 1)]
box_strip_v2(vec,c(-3,3),col=hue_pal()(3)[c(1,2)])
abline(h=0,lty=2)

#comp polarization
df <- data.frame(Type="a",Score=OL_OPC_human$OL_score_2)
df[which(OL_OPC_human$Type == "Old"),1] <- "b"
df[which(OL_OPC_human$Group == 1),] -> df
ggplot(df,aes(Score,fill=Type)) + geom_density(alpha=0.5) + theme_bw() + xlim(c(-3,3)) + scale_fill_manual(values=hue_pal()(3)[c(1,2)]) + ylim(c(0,0.6))

df <- data.frame(Type="a",Score=OL_OPC_mouse$OL_score_2)
df[which(OL_OPC_mouse$Type == "Old"),1] <- "b"
df[which(OL_OPC_mouse$Type == "Middle"),1] <- "c"
df[which(OL_OPC_human$Group == 1),] -> df
ggplot(df,aes(Score,fill=Type)) + geom_density(alpha=0.5) + theme_bw() + xlim(c(-3,3)) + scale_fill_manual(values=hue_pal()(3)[c(1,2,3)]) + ylim(c(0,0.6))
