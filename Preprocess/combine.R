library(Seurat)
library(Matrix)

load("human_counts.dat")
load("mouse_counts.dat")

homolog <- read.table("/home/ytanaka/projects/def-ytanaka/ytanaka/mart_export.txt",sep="\t",header=T)
homolog <- homolog[which(duplicated(homolog[,3])==F),]
homolog[homolog[,4]!="",] -> homolog
homolog <- homolog[which(duplicated(homolog[,4])==F),]

homolog[,3] -> samp
names(samp) <- homolog[,4]
samp -> homolog

intersect(rownames(human_counts),names(homolog)) -> samp
homolog[samp] -> homolog
samp <- c(rownames(mouse_counts),homolog)
samp[duplicated(samp)] -> homolog

combine <- cbind(human_counts[names(homolog),],mouse_counts[homolog,])
combine <- CreateSeuratObject(counts=combine,project="combined",names.field = 2, names.delim = "-")

org <- rep("mouse",ncol(combine))
org[grep("1_",combine@meta.data$orig.ident)] <- "human"
AddMetaData(combine,metadata=org,col.name="Organism") -> combine

load("human_meta.dat")
load("mouse_meta.dat")
rbind(human_meta[,c("ID","Type","CellType")],mouse_meta[,c("ID","Type","CellType")]) -> samp
AddMetaData(combine,metadata=samp$ID,col.name="ID") -> combine
AddMetaData(combine,metadata=samp$Type,col.name="Type") -> combine
AddMetaData(combine,metadata=samp$CellType,col.name="CellType") -> combine

combine.list <- SplitObject(object=combine,split.by="orig.ident")

for(i in 1:length(combine.list)){
combine.list[[i]] <- NormalizeData(combine.list[[i]],verbose=F)
combine.list[[i]] <- FindVariableFeatures(combine.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

combine.anchors <- FindIntegrationAnchors(combine.list,dims=1:20)
combine.integrated <- IntegrateData(anchorset = combine.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = combine.integrated) <- "integrated"
combine.integrated <- ScaleData(object = combine.integrated, verbose = FALSE)
combine.integrated <- RunPCA(object = combine.integrated, npcs = 20, verbose = FALSE)
combine.integrated <- RunUMAP(object = combine.integrated, reduction = "pca", dims = 1:20)

combine.integrated <- FindNeighbors(combine.integrated,dims=1:20,reduction="pca")
combine.integrated <- FindClusters(combine.integrated)

save(combine.integrated,file="combine.dat")
save(combine.list,file="combine_list.dat")
save(combine.anchors,file="combine_anchors.dat")

library(genefilter)
preclust <- sort(unique(combine.integrated@active.ident))
clust_num <- length(unique(preclust))
ratio <- vector("list",clust_num)
names(ratio) <- preclust
ratio -> pval
ratio -> dif_gene_1p25_pval005
row_count <- nrow(combine.integrated@assays$RNA)
col_count <- ncol(combine.integrated@assays$RNA)
row_split <- 10000

for(i in 1:clust_num){
      print(i)
      select <- numeric(col_count)
      select[which(combine.integrated@active.ident == preclust[i])] <- 1
      rowid <- 1
      while(rowid + row_split -1 < row_count){
      exp <- as.matrix(combine.integrated@assays$RNA[rowid:(rowid+row_split-1),1:col_count])
      ratio[[i]] <- c(ratio[[i]],rowMeans(exp[,which(select==1)]) - rowMeans(exp[,which(select==0)]))
      pval[[i]]  <- rbind(pval[[i]],rowttests(as.matrix(exp),fac=factor(select)))
      rowid <- rowid + row_split
      }
      exp <- as.matrix(combine.integrated@assays$RNA[rowid:row_count,1:col_count])
      ratio[[i]] <- c(ratio[[i]],rowMeans(exp[,which(select==1)]) - rowMeans(exp[,which(select==0)]))
      pval[[i]]  <- rbind(pval[[i]],rowttests(as.matrix(exp),fac=factor(select)))
      dif_gene_1p25_pval005[[i]] <- names(ratio[[i]])[which(ratio[[i]] > log2(1.25) & pval[[i]][,3] < 0.05)]
}
save(ratio,file="ratio.dat")
save(pval,file="pval.dat")
save(dif_gene_1p25_pval005,file="dif_gene_1p25_pval005.dat")

rm(exp)


