library(Seurat)
library(paran)

##### subclustering 
expr <- readRDS("/SATA/rlhua/project/eCCA/seurat_RDSdata/expr.RDS")
celltype <- unique(expr$celltype)
for (i in celltype) {
  sr <- expr[,expr$celltype == i]
  sr <- FindVariableFeatures(sr, selection.method = "vst", nfeatures = 1000, verbose = F)
  sr <- ScaleData(sr,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"))
  mtx <- as.matrix(sr@assays$integrated@scale.data)
  pr <- paran::paran(t(mtx),iterations = 10)
  npca <- sum(pr$Ev/pr$RndEv > 1.5)
  sr <- RunPCA(sr,npcs=npca) # pca
  sr$npca <- npca
  sr <- FindNeighbors(sr, reduction = "pca", dims = 1:npca)
  sr <- FindClusters(sr, resolution = 0.5) # 聚类
  sr <- RunTSNE(sr, dims = 1:npca)
  sr <- RunUMAP(sr, dims = 1:npca) # umap
  sr$tissue <- ifelse(sr$orig.ident %in% c("190065A_P1","190065D_P4","190065E_P5"),"Normal","Tumor")
  saveRDS(sr,paste0("/SATA/rlhua/project/eCCA/seurat_RDSdata/",i,".RDS"))
}






















