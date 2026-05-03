library(Seurat)
library(Matrix)

mouse <- Read10X_h5("merge/outs/count/filtered_feature_bc_matrix.h5")
mouse <- CreateSeuratObject(counts=mouse,project="combined",names.field = 2, names.delim = "-")
mito.gene <-  rownames(mouse)[grep("mt-",rownames(mouse))]
percent.mito <- Matrix::colSums(mouse@assays$RNA[mito.gene, ])/Matrix::colSums(mouse@assays$RNA)
AddMetaData(mouse,metadata=percent.mito,col.name="percent.mito") -> mouse

type <- rep("Middle",ncol(mouse))
type[grep("-1",colnames(mouse))] <- "Old"
type[grep("-2",colnames(mouse))] <- "Old"
type[grep("-3",colnames(mouse))] <- "Old"
type[grep("-4",colnames(mouse))] <- "Old"
type[grep("-9",colnames(mouse))] <- "Young"
type[grep("-10",colnames(mouse))] <- "Young"
type[grep("-11",colnames(mouse))] <- "Young"
type[grep("-12",colnames(mouse))] <- "Young"
AddMetaData(mouse,metadata=type,col.name="Type") -> mouse

ID <- rep("Middle",ncol(mouse))
ID[grep("-1",colnames(mouse))] <- "A1"
ID[grep("-2",colnames(mouse))] <- "A2"
ID[grep("-3",colnames(mouse))] <- "A3"
ID[grep("-4",colnames(mouse))] <- "A4"
ID[grep("-5",colnames(mouse))] <- "M1"
ID[grep("-6",colnames(mouse))] <- "M2"
ID[grep("-7",colnames(mouse))] <- "M3"
ID[grep("-8",colnames(mouse))] <- "M4"
ID[grep("-9",colnames(mouse))] <- "Y1"
ID[grep("-10",colnames(mouse))] <- "Y2"
ID[grep("-11",colnames(mouse))] <- "Y3"
ID[grep("-12",colnames(mouse))] <- "Y4"
AddMetaData(mouse,metadata=ID,col.name="ID") -> mouse


mouse <- subset(mouse,subset= nFeature_RNA > 100 & nFeature_RNA < 7000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.05)

mouse.list <- SplitObject(object=mouse,split.by="orig.ident")

for(i in 1:length(mouse.list)){
mouse.list[[i]] <- NormalizeData(mouse.list[[i]],verbose=F)
mouse.list[[i]] <- FindVariableFeatures(mouse.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

mouse.anchors <- FindIntegrationAnchors(mouse.list,dims=1:20)
mouse.integrated <- IntegrateData(anchorset = mouse.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = mouse.integrated) <- "integrated"
mouse.integrated <- ScaleData(object = mouse.integrated, verbose = FALSE)
mouse.integrated <- RunPCA(object = mouse.integrated, npcs = 20, verbose = FALSE)
mouse.integrated <- RunUMAP(object = mouse.integrated, reduction = "pca", dims = 1:20)

mouse.integrated <- FindNeighbors(mouse.integrated,dims=1:20,reduction="pca")
mouse.integrated <- FindClusters(mouse.integrated)

save(mouse.integrated,file="mouse.dat")
save(mouse.list,file="mouse_list.dat")
save(mouse.anchors,file="mouse_anchors.dat")

library(genefilter)
preclust <- sort(unique(mouse.integrated@active.ident))
clust_num <- length(unique(preclust))
ratio <- vector("list",clust_num)
names(ratio) <- preclust
ratio -> pval
ratio -> dif_gene_1p25_pval005
row_count <- nrow(mouse.integrated@assays$RNA)
col_count <- ncol(mouse.integrated@assays$RNA)
row_split <- 10000

for(i in 1:clust_num){
      print(i)
      select <- numeric(col_count)
      select[which(mouse.integrated@active.ident == preclust[i])] <- 1
      rowid <- 1
      while(rowid + row_split -1 < row_count){
      exp <- as.matrix(mouse.integrated@assays$RNA[rowid:(rowid+row_split-1),1:col_count])
      ratio[[i]] <- c(ratio[[i]],rowMeans(exp[,which(select==1)]) - rowMeans(exp[,which(select==0)]))
      pval[[i]]  <- rbind(pval[[i]],rowttests(as.matrix(exp),fac=factor(select)))
      rowid <- rowid + row_split
      }
      exp <- as.matrix(mouse.integrated@assays$RNA[rowid:row_count,1:col_count])
      ratio[[i]] <- c(ratio[[i]],rowMeans(exp[,which(select==1)]) - rowMeans(exp[,which(select==0)]))
      pval[[i]]  <- rbind(pval[[i]],rowttests(as.matrix(exp),fac=factor(select)))
      dif_gene_1p25_pval005[[i]] <- names(ratio[[i]])[which(ratio[[i]] > log2(1.25) & pval[[i]][,3] < 0.05)]
}
save(ratio,file="ratio.dat")
save(pval,file="pval.dat")
save(dif_gene_1p25_pval005,file="dif_gene_1p25_pval005.dat")

rm(exp)
