library(Seurat)
library(Signac)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(Signac)

load("~/scratch/PD_midbrain_v24/merge.peaks_v2_annote.dat")
human_ATAC <- subset(merge.peaks,subset = Type != "PD")
save(human_ATAC,file="human_ATAC.dat")

load("~/scratch/mouse_midbrain_v16/human.integrated.dat")
human_ATAC@active.ident <- factor(human.integrated@active.ident[colnames(human_ATAC)])
clust <- as.numeric(human_ATAC@active.ident)
ODC_ATAC <- subset(human_ATAC, cells=colnames(human_ATAC)[which(clust >= 2 & clust <= 13)])
ODC_RNA <- subset(human.integrated,cells=colnames(human_ATAC)[which(clust >= 2 & clust <= 13)])

sq <- RSQLite::dbConnect(RSQLite::SQLite(),db(JASPAR2024()))
pfm <- TFBSTools::getMatrixSet(sq,list(species = "Homo sapiens", collection = "CORE"))
ODC_ATAC <- AddMotifs(object = ODC_ATAC,genome=BSgenome.Hsapiens.UCSC.hg38,pfm=pfm)
DCR <- FindMarkers(ODC_ATAC,ident.1=7,test.use="LR",latent.vars='nCount_ATAC',min.pct=0.01)
DEG <- FindMarkers(ODC_RNA,ident.1=7,latent.vars='nCount_RNA')
save(DCR,file="DCR.dat")
save(DEG,file="DEG.dat")
save(ODC_ATAC,file="ODC_ATAC.dat")
save(ODC_RNA,file="ODC_RNA.dat")

source("~/R_command/all_go_analysis_v2_symbol2.R")
source("~/R_command/all_GO_heat2.R")
dif_gene_2fold_pval005 <- vector("list",2)
names(dif_gene_2fold_pval005) <- c("up","down")
dif_gene_2fold_pval005[[1]] <- rownames(DEG)[which(DEG[,2] > 1 & DEG[,5] < 0.05)]
dif_gene_2fold_pval005[[2]] <- rownames(DEG)[which(DEG[,2] < -1 & DEG[,5] < 0.05)]
dif_gene_2fold_pval005 -> dif_gene_2fold_pval005_GO
for(i in 1:2){
all_go_analysis_v2_symbol2(dif_gene_2fold_pval005[[i]],"human") -> dif_gene_2fold_pval005_GO[[i]]
}
all_GO_heat2(dif_gene_2fold_pval005_GO, 1) -> samp
FeaturePlot(ODC_RNA,dif_gene_2fold_pval005[[1]][1:9],min.cutoff=0,max.cutoff=1)

load("~/scratch/PD_midbrain_v24/rmsk_ratio.dat")
exp <- ODC_RNA@assays$RNA@data
ODC_RNA@meta.data$rmsk_ratio <- rmsk_ratio[colnames(ODC_RNA)]
cor_rmsk <- cor(t(as.matrix(exp)),ODC_RNA$rmsk_ratio)
cor_rmsk <- cor_rmsk[,1]
cor_rmsk <- sort(cor_rmsk,decreasing=T)
save(cor_rmsk,file="cor_rmsk.dat")
write.table(data.frame(cor_rmsk),file="cor_rmsk.rnk",sep="\t",quote=F,col.names=F)
FeaturePlot(ODC_RNA,c("DST","AKAP9","CCDC88A","AFMID"),min.cutoff=2,max.cutoff=4,cols=c("gray","red")) #FKBP5
my.col <- rep("black",length(cor_rmsk))
names(my.col) <- names(cor_rmsk)
my.col[c("SLC5A11","SLC25A29","MALAT1","PLP1","PDE4B","MOBP","FCHSD2","OPALIN","RBFOX1","AFF3","FKBP5")] <- "red"
plot(sort(cor_rmsk,decreasing=F),pch=20,ylim=c(-0.6,0.6),col=my.col[names(sort(cor_rmsk,decreasing=F))])
abline(h=0,lty=2)

my.col <- rep("black",nrow(DCR))
my.col[which(DCR[,2] > 1 & DCR[,5] < 0.05)] <- "orange"
my.col[which(DCR[,2] < -1 & DCR[,5] < 0.05)] <- "dodgerblue"
plot(DCR[,2],-log10(DCR[,5]),col=my.col,pch=20,xlim=c(-2,2))

exp <- ODC_RNA@assays$ATAC@data
exp <- exp[rownames(DCR),]
cor_rmsk_ATAC <- cor(t(as.matrix(exp)),rmsk_ratio[colnames(exp)])
cor_rmsk_ATAC <- cor_rmsk_ATAC[,1]
save(cor_rmsk_ATAC,file="cor_rmsk_ATAC.dat")


library(ggplot2)
library(scales)

#link analysis
ODC_RNA@assays$ATAC <- ODC_ATAC@assays$ATAC
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
DefaultAssay(ODC_RNA) <- "ATAC"
Annotation(ODC_RNA) <- annotations
ODC_RNA <- LinkPeaks(ODC_RNA,peak.assay="ATAC", expression.assay="RNA")
save(ODC_RNA,file="ODC_RNA_v2.dat")

ODC_RNA_7 <- subset(ODC_RNA,cells=colnames(ODC_RNA)[which(ODC_RNA@active.ident == 7)])
ODC_RNA_other <- subset(ODC_RNA,cells=colnames(ODC_RNA)[which(ODC_RNA@active.ident != 7)])
ODC_RNA_7 <- LinkPeaks(ODC_RNA,peak.assay="ATAC", expression.assay="RNA")
save(ODC_RNA_7,file="ODC_RNA_7_v2.dat")
ODC_RNA_other <- LinkPeaks(ODC_RNA_other,peak.assay="ATAC", expression.assay="RNA")
save(ODC_RNA_other,file="ODC_RNA_other_v2.dat")

#rmsk analysis
library(venneuler)
link7 <- Links(ODC_RNA_7)
link_other <- Links(ODC_RNA_other)
link7_connect <- paste(link7$peak,":",link7$gene,sep="")
link_other_connect <- paste(link_other$peak,":",link_other$gene,sep="")
length(intersect(link_other_connect,link7_connect))
ven <- venneuler(c(A=3748-3380,B=4128-3380,"A&B"=3380))
ven$labels <- c("","")
plot(ven,col=c('red','green'))

source("~/R_command/exclusive.R")
exclusive(link_other_connect,link7_connect) -> link7_unique

strsplit(link7_unique,":") -> samp
target <- rep("",length(samp))
for(i in 1:length(samp)){
 target[i] <- samp[[i]][2]
}
target <- unique(target)

rmsk <- read.table("combined.peaks_v2_rmsk.txt",sep="\t",header=F)
peak <- unique(rmsk[,4])
rmsk_flag <- rep("non_rmsk",length(cor_rmsk_ATAC))
names(rmsk_flag) <- names(cor_rmsk_ATAC)
rmsk_flag[intersect(names(rmsk_flag),peak)] <- "rmsk"
summary(factor(rmsk_flag))
ggplot(df,aes(Cor,fill=Flag)) + geom_density(alpha=0.5) + theme_bw() + xlim(c(-0.1,0.1)) + scale_fill_manual(values=c("dodgerblue","orange"))

DefaultAssay(ODC_RNA_7) <- "ATAC"

ODC_RNA_7@meta.data$nCount_ATAC <- ODC_ATAC$nCount_ATAC[colnames(ODC_RNA_7)]
ODC_RNA_7@meta.data$nFeature_ATAC <- ODC_ATAC$nFeature_ATAC[colnames(ODC_RNA_7)]

id <- unique(ODC_RNA_7$orig.ident)
list <- read.table("frag_list.txt",sep="\t",header=F)
list <- list[,1]
frag_vec <- vector("list",length(list))
names(frag_vec) <- list
for(i in 1:length(list)){
 frag_vec[[i]] <- CreateFragmentObject(path=list[i],cells=colnames(ODC_RNA_7)[which(ODC_RNA_7$orig.ident == id[i])])
}
Fragments(ODC_RNA_7) <- NULL
Fragments(ODC_RNA_7) <- frag_vec
save(ODC_RNA_7,file="ODC_RNA_7_v2.dat")

ODC_RNA_other@meta.data$nCount_ATAC <- ODC_ATAC$nCount_ATAC[colnames(ODC_RNA_other)]
ODC_RNA_other@meta.data$nFeature_ATAC <- ODC_ATAC$nFeature_ATAC[colnames(ODC_RNA_other)]
Fragments(ODC_RNA_other) <- NULL
Fragments(ODC_RNA_other) <- frag_vec
save(ODC_RNA_other,file="ODC_RNA_other_v2.dat")

type <- rep("a",ncol(ODC_RNA_7))
type[which(ODC_RNA_7$Type == "Old")] <- "b"
L1M <- unique(rmsk[intersect(which(rmsk[,10]!="HAL1M8" & rmsk[,10]!="HAL1ME"),grep("L1M",rmsk[,10])),4])
obj <- ODC_RNA_7
link7_RBFOX1 <- link7[link7$gene=="RBFOX1",]
#link7_RBFOX1[which(link7_RBFOX1$peak == "chr16-6214890-6216693" | link7_RBFOX1$peak == "chr16-6232091-6235488" | link7_RBFOX1$peak == "chr16-6246707-6254003" | link7_RBFOX1$peak == "chr16-5927935-5930313" | link7_RBFOX1$peak == "chr16-6185442-6191509" | link7_RBFOX1$peak == "chr16-6232091-6235488" | link7_RBFOX1$peak == "chr16-6356194-6358100" | link7_RBFOX1$peak == "chr16-5894771-5896746" | link7_RBFOX1$peak == "chr16-6177366-6179208" | link7_RBFOX1$peak == "chr16-6381429-6385436"),] -> link7_RBFOX1_L1

rmsk_flag <- rep("non-rmsk",length(link7_RBFOX1$peak))
names(rmsk_flag) <- link7_RBFOX1$peak
rmsk_peak <- unique(rmsk[,4])
rmsk_flag[intersect(names(rmsk_flag),rmsk_peak)] <- "rmsk"

L1P <- unique(rmsk[intersect(which(rmsk[,10]!="L1PREC2"),grep("L1P",rmsk[,10])),4])
L1M <- unique(rmsk[intersect(which(rmsk[,10]!="HAL1M8" & rmsk[,10]!="HAL1ME"),grep("L1M",rmsk[,10])),4])
L_other <- unique(rmsk[which(rmsk[,10]=="L3" | rmsk[,10]=="L3b" | rmsk[,10]=="L4_C_Mam" | rmsk[,10]=="L4_B_Mam" | rmsk[,10]=="L4_A_Mam" | rmsk[,10]=="L5"),4])
rmsk_flag[intersect(names(rmsk_flag),L1P)] <- "L1P"
rmsk_flag[intersect(names(rmsk_flag),L1M)] <- "L1M"
rmsk_flag[intersect(names(rmsk_flag),L_other)] <- "L_other"

source("~/R_command/box_strip_v2.R")
cor(ODC_RNA@assays$RNA@data["RBFOX1",],t(as.matrix(ODC_RNA@assays$ATAC@data[names(rmsk_flag),])))[1,] -> cor_RBFOX1
vec <- vector("list",5)
#vec[[1]] <- cor_RBFOX1[which(rmsk_flag == "non-rmsk")]
#vec[[2]] <- cor_RBFOX1[which(rmsk_flag == "rmsk" & rmsk_flag != "L_other" & rmsk_flag != "L1M" & rmsk_flag != "L1P")]
#vec[[3]] <- cor_RBFOX1[which(rmsk_flag == "L_other")]
#vec[[4]] <- cor_RBFOX1[which(rmsk_flag == "L1M")]
#vec[[5]] <- cor_RBFOX1[which(rmsk_flag == "L1P")]

vec[[1]] <- link7_RBFOX1$score[which(rmsk_flag == "non-rmsk")]
vec[[2]] <- link7_RBFOX1$score[which(rmsk_flag == "rmsk" & rmsk_flag != "L_other" & rmsk_flag != "L1M" & rmsk_flag != "L1P")]
vec[[3]] <- link7_RBFOX1$score[which(rmsk_flag == "L_other")]
vec[[4]] <- link7_RBFOX1$score[which(rmsk_flag == "L1M")]
vec[[5]] <- link7_RBFOX1$score[which(rmsk_flag == "L1P")]

box_strip_v2(vec,c(0.05,0.3),hue_pal()(5)[5:1])

type <- rep("a",ncol(ODC_RNA_7))
type[which(ODC_RNA_7$Type == "Old")] <- "b"
obj <- ODC_RNA_7
obj@meta.data$Type2 <- type

link7_RBFOX1 <- link7[link7$gene=="RBFOX1",]
link_select <- link7_RBFOX1[which(rmsk_flag == "L1P"),]
link_select$score[link_select$score > 0.2] <- 0.2
#link_select$score <- 1
Links(obj) <- link_select
p <- CoveragePlot(obj,region="chr16-5500000-6500000",assay="ATAC",links=TRUE,tile=F,annotation=TRUE,group.by="Type2",window=5000,min.cutoff=0)
p[[1]][[4]] <- p[[1]][[4]] & scale_color_gradient2(low="gray",mid="purple",high="blue")
p & scale_fill_manual(values = hue_pal()(3)[1:2])

link7_RBFOX1 <- link7[link7$gene=="RBFOX1",]
link_select <- link7_RBFOX1[which(rmsk_flag == "L1M"),]
link_select$score[link_select$score > 0.2] <- 0.2
#link_select$score <- 1
Links(obj) <- link_select
p <- CoveragePlot(obj,region="chr16-5500000-6500000",assay="ATAC",links=TRUE,tile=F,annotation=TRUE,group.by="Type2",window=5000,min.cutoff=0)
p[[1]][[4]] <- p[[1]][[4]] & scale_color_gradient2(low="gray",mid="purple",high="blue")
p & scale_fill_manual(values = hue_pal()(3)[1:2])

link7_RBFOX1 <- link7[link7$gene=="RBFOX1",]
link_select <- link7_RBFOX1[which(rmsk_flag == "L_other"),]
link_select$score[link_select$score > 0.2] <- 0.2
#link_select$score <- 1
Links(obj) <- link_select
p <- CoveragePlot(obj,region="chr16-5500000-6500000",assay="ATAC",links=TRUE,tile=F,annotation=TRUE,group.by="Type2",window=5000,min.cutoff=0)
p[[1]][[4]] <- p[[1]][[4]] & scale_color_gradient2(low="gray",mid="purple",high="blue")
p & scale_fill_manual(values = hue_pal()(3)[1:2])

link7_RBFOX1 <- link7[link7$gene=="RBFOX1",]
link_select <- link7_RBFOX1[which(rmsk_flag == "rmsk"),]
link_select$score[link_select$score > 0.2] <- 0.2
#link_select$score <- 1
Links(obj) <- link_select
p <- CoveragePlot(obj,region="chr16-5500000-6500000",assay="ATAC",links=TRUE,tile=F,annotation=TRUE,group.by="Type2",window=5000,min.cutoff=0)
p[[1]][[4]] <- p[[1]][[4]] & scale_color_gradient2(low="gray",mid="purple",high="blue")
p & scale_fill_manual(values = hue_pal()(3)[1:2])

link7_RBFOX1 <- link7[link7$gene=="RBFOX1",]
link_select <- link7_RBFOX1[which(rmsk_flag == "non-rmsk"),]
link_select$score[link_select$score > 0.1] <- 0.1
#link_select$score <- 1
Links(obj) <- link_select
p <- CoveragePlot(obj,region="chr16-5500000-6500000",assay="ATAC",links=TRUE,tile=F,annotation=TRUE,group.by="Type2",window=5000,min.cutoff=0)
p[[1]][[4]] <- p[[1]][[4]] & scale_color_gradient(low="gray",high="purple")
p & scale_fill_manual(values = hue_pal()(3)[1:2])

df <- data.frame(x=c(0,1),y=c(0,1),z=c(-0.1,0.1))
ggplot(df,aes(x,y)) + geom_point(aes(colour=z)) + scale_color_gradient2(low="gray",mid="purple",high="blue")