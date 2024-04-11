library(scCancer)
library(ggplot2)
library(magrittr)
library(Seurat)
library(cowplot)
library(data.table)

##### scCancer
dir <- list.files("/SATA/rlhua/project/eCCA/cellranger")
savePath <- c("/SATA/rlhua/project/eCCA/scCancer/")
for(fold in dir){
  stat.result <- runScStatistics(
    dataPath = file.path(fold,"outs"),
    savePath = paste0(savePath,"/stat/",fold),
    sampleName = fold,
    species="human",
  )
  
  anno.results <- runScAnnotation(
    dataPath = file.path(fold,"outs"),
    statPath = paste0(savePath,"/stat/",fold),
    savePath = paste0(savePath,"/anno/",fold),
    sampleName = fold,
    geneSet.method = "average",
    pc.use = 50,
    species = "human",
    genome = "hg38"
  )
} 

##### combine samples
single.savePaths <- list.files("/SATA/rlhua/project/eCCA/scCancer/anno",full.names=T)
expr.list <- list()
for(i in 1:length(single.savePaths)){
  sampleName <- gsub(".*/","",single.savePaths[i])
  print(sampleName)
  expr.list[[sampleName]] <- readRDS(paste0(single.savePaths[i], "/expr.RDS"))
} 

##### add scrublet score
fs <- list.files("/SATA/rlhua/project/eCCA/scrublet","doublet.txt",full.names = T)
for(i in 1:length(fs)){
  dt <- fread(fs[i],sep = ",")
  dt$barcode <- gsub("-1$","",dt$barcode)
  cells_singlet <- dt$barcode[dt$doublet_scores < dt$threshold]
  expr.list[[i]] <- expr.list[[i]][,colnames(expr.list[[i]]) %in% cells_singlet]
  scrublet_score <- dt$doublet_scores
  names(scrublet_score) <- dt$barcode
  expr.list[[i]]$scrublet_score <- scrublet_score[colnames(expr.list[[i]])]
}

##### integration, dimensionality reduction and clustering
npc=50
expr.anchors <- Seurat::FindIntegrationAnchors(object.list = expr.list,dims = 1:npc) 
anchors <- expr.anchors@anchors
anchors$cellType1 <- "NULL"
anchors$cellType2 <- "NULL"
anchors$malignType1 <- "NULL"
anchors$malignType2 <- "NULL"
anchors$malignScore1 <- -1
anchors$malignScore2 <- -1
for(oi in expr.anchors@reference.objects){
  cur.ix <- which(anchors$dataset1 == oi)
  anchors$cellType1[cur.ix] <- expr.list[[oi]]@meta.data$Cell.Type[anchors$cell1[cur.ix]] 
  anchors$malignType1[cur.ix] <- expr.list[[oi]]@meta.data$Malign.type[anchors$cell1[cur.ix]] 
  anchors$malignScore1[cur.ix] <- expr.list[[oi]]@meta.data$Malign.score[anchors$cell1[cur.ix]]
  cur.ix <- which(anchors$dataset2 == oi)
  anchors$cellType2[cur.ix] <- expr.list[[oi]]@meta.data$Cell.Type[anchors$cell2[cur.ix]]
  anchors$malignType2[cur.ix] <- expr.list[[oi]]@meta.data$Malign.type[anchors$cell2[cur.ix]]
  anchors$malignScore2[cur.ix] <- expr.list[[oi]]@meta.data$Malign.score[anchors$cell2[cur.ix]]
}
anchors.new <- subset(anchors, cellType1 != "Epithelial" & cellType1 != "Unknown" & cellType2 != "Epithelial" & cellType2 != "Unknown") 
expr.anchors@anchors <- anchors.new 

genes <- c()
for (i in 1:length(expr.list)) {
  genes <- c(genes,rownames(expr.list[[i]]))
}
genes <- unique(genes)
expr <- Seurat::IntegrateData(anchorset = expr.anchors,dims = 1:npc,features.to.integrate =genes) 
expr <- Seurat::ScaleData(expr,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent")) 
mtx <- as.matrix(expr@assays$integrated@scale.data)
pr <- paran::paran(t(mtx),iterations = 10)
npca <- sum(pr$Ev/pr$RndEv > 1.5)
expr <- Seurat::RunPCA(expr,npcs=npca) # pca
expr <- Seurat::FindNeighbors(expr, reduction = "pca", dims = 1:npca)
expr <- Seurat::FindClusters(expr, resolution = 0.8)
expr <- Seurat::FindClusters(expr, resolution = 1.2,n.start = 20,n.iter = 20)
expr <- Seurat::FindClusters(expr, resolution = 0.1)
expr <- Seurat::RunTSNE(expr, perplexity=60, max_iter=2000, stop_lying_iter=1000,dims = 1:npca) # tsne
expr <- Seurat::RunUMAP(expr, dims = 1:npca) # umap

##### cell type annotation
expr$celltype <- ""
expr$celltype[expr$seurat_clusters %in% c(0,1,13)] <- "T cells"
expr$celltype[expr$seurat_clusters %in% c(2,10)] <- "B cells"
expr$celltype[expr$seurat_clusters %in% c(3)] <- "Myeloid"
expr$celltype[expr$seurat_clusters %in% c(4)] <- "Fibroblasts"
expr$celltype[expr$seurat_clusters %in% c(5)] <- "Endothelial"
expr$celltype[expr$seurat_clusters %in% c(6)] <- "Mast cells"
expr$celltype[expr$seurat_clusters %in% c(11)] <- "Schwann"
expr$celltype[expr$seurat_clusters %in% c(7,8,9,12)] <- "Epithelial"
expr$celltype[expr$integrated_snn_res.0.8==15] <- "HPLCs"

saveRDS(expr,file = "/SATA/rlhua/project/eCCA/seurat_RDSdata/expr.RDS")

##### tsne plot
colors <- c("#DB3B32","#CED951","#82BC5E","#476AAE","#8E579B","#EB9B3F","#85CBDB","#DA408F","#A89770")
names(colors) <- sort(unique(expr$celltype))
p1 <- ggplot()+
  geom_point(aes(x=expr@reductions$tsne@cell.embeddings[,"tSNE_1"],
                 y=expr@reductions$tsne@cell.embeddings[,"tSNE_2"],
                 color=expr$celltype),shape=16,size=0.001)+
  scale_color_manual(values = colors)+
  xlab("tSNE_1")+ylab("tSNE_2")+
  theme_classic()+
  theme(aspect.ratio = 1,legend.position = "none",axis.text = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
p2 <- ggplot()+
  geom_point(aes(x=expr@reductions$tsne@cell.embeddings[,"tSNE_1"],
                 y=expr@reductions$tsne@cell.embeddings[,"tSNE_2"],
                 color=expr$orig.ident),shape=16,size=0.001)+
  scale_color_brewer(palette="Set1",labels=c("P1T","P1C","P2T","P3T","P3C","P4T","P4C","P5T"))+
  xlab("tSNE_1")+ylab("tSNE_2")+
  theme_classic()+
  theme(aspect.ratio = 1,legend.position = "none",axis.text = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
expr$tissue <- ifelse(expr$orig.ident %like% "_C","TUMOR","CONTROL")
p3 <- ggplot()+
  geom_point(aes(x=expr@reductions$tsne@cell.embeddings[,"tSNE_1"],
                 y=expr@reductions$tsne@cell.embeddings[,"tSNE_2"],
                 color=expr$tissue),shape=16,size=0.001)+
  scale_color_manual(values = c('palegreen3','dodgerblue4'))+
  xlab("tSNE_1")+ylab("tSNE_2")+
  theme_classic()+
  theme(aspect.ratio = 1,legend.position = "none",axis.text = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
p4 <- ggplot()+
  geom_point(aes(x=expr@reductions$tsne@cell.embeddings[,"tSNE_1"],
                 y=expr@reductions$tsne@cell.embeddings[,"tSNE_2"],
                 color=expr$nCount_RNA),shape=16,size=0.001)+
  scale_color_gradient(low = "grey",high = "blue")+
  xlab("tSNE_1")+ylab("tSNE_2")+
  theme_classic()+
  theme(aspect.ratio = 1,legend.position = "none",axis.text = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
save_plot("/SATA/rlhua/project/eCCA/figure/tSNE_plot_of_origin.png",p1,ncol = 1 ,nrow = 1,base_width = 1.7,base_height = 1.7)
save_plot("/SATA/rlhua/project/eCCA/figure/tSNE_plot_of_sample.png",p2,ncol = 1 ,nrow = 1,base_width = 1.7,base_height = 1.7)
save_plot("/SATA/rlhua/project/eCCA/figure/tSNE_plot_of_celltype.png",p3,ncol = 1 ,nrow = 1,base_width = 1.7,base_height = 1.7)
save_plot("/SATA/rlhua/project/eCCA/figure/tSNE_plot_of_nUMI.png",p4,ncol = 1 ,nrow = 1,base_width = 1.7,base_height = 1.7)

##### violinplot for celltype markers
markers <- c('CD79A','MS4A1','PLVAP','FLT1','EPCAM','KRT18','COL1A1','PDGFRB',"MKI67","TOP2A",'IL1RL1','KIT','CD14','LYZ','MPZ','NCAM1','CD3D','CD3E')
dt <- reshape2::melt(as.matrix(expr@assays$RNA@data[markers,])) %>% as.data.table()
colnames(dt)[1:2] <- c("gene","cell")
dt$celltype <- expr$celltype[dt$cell]
colors <- c("#DB3B32","#CED951","#82BC5E","#476AAE","#8E579B","#EB9B3F","#85CBDB","#DA408F","#A89770")
names(colors) <- sort(unique(expr$celltype))
ggviolin(dt,'celltype','value',scale='width',fill='celltype',color='transparent')+
  facet_grid(gene~.)+theme(strip.background = element_blank(),legend.position='None',axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.3))+
  scale_fill_manual(values = colors)+xlab('')+ylab('')+theme(axis.line = element_line(size = 0.2),strip.text.y = element_text(angle = 180,size=7,hjust=1),axis.text.y = element_text(angle = 90,hjust=0.3),axis.text = element_text(size=7),text = element_text(size=7),axis.ticks=element_line(size=0.2))+scale_y_continuous(limit=c(-0.3,5.3),breaks=c(0,5))
ggsave("/SATA/rlhua/project/eCCA/figure/celltype_markers_violinplot.pdf",height =3.5,width = 2)
