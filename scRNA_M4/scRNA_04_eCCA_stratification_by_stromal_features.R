library(ggplot2)
library(magrittr)
library(Seurat)
library(cowplot)
library(data.table)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(corrplot)

##### clinical information
sample_info <- fread("/SATA/rlhua/project/eCCA/bulkRNAseq/sample_info.tsv",sep = "\t",header = T)
sample_info$status <- ifelse(sample_info$status=="Alive",0,1)
sample_info <- sample_info[sample_name!=""]

##### RNA expression matrix
mtx <- readRDS("/SATA/rlhua/project/eCCA/result/bulkRNAseq/FPKM_gencode.v22.annotation.gtf_AllSampleCountMatrix.RDS") 
geneID2symbol <- fread("/SATA/rlhua/project/eCCA/result/bulkRNAseq/geneID2symbol.tsv",sep = "\t",header = F,col.names = c("ID","symbol")) 
idx <- match(rownames(mtx), geneID2symbol$ID) 
mtx <- mtx[!is.na(idx),]
mtx[["symbol"]] <- geneID2symbol$symbol[idx[!is.na(idx)]] 
mtx <- aggregate(mtx[,-ncol(mtx)],list(mtx$symbol),sum) 
rownames(mtx) <- mtx$Group.1
mtx <- mtx[,-1]
mtx <- mtx[,sample_info$sample]
mtx <- as.matrix(mtx)
mtx <- mtx[rowMeans(mtx)>0,] 
mtx <- log1p(mtx) 

##### markers
dt <- readRDS("/SATA/rlhua/project/eCCA/markers/markers_celltype/clustermarker.RDS")
schwann_markers <- dt$gene[dt$cluster=="Schwann"]
HPLC_markers <- dt$gene[dt$cluster=="HPLCs"]
dt <- readRDS("/SATA/rlhua/project/eCCA/markers/markers_Fibroblasts/clustermarker.RDS")
Fib_C3_markers <- dt$gene[dt$cluster==3]
dt <- readRDS("/SATA/rlhua/project/eCCA/markers/markers_T_cells/clustermarker.RDS")
Treg_markers <- dt$gene[dt$cluster==2]
dt <- readRDS("/SATA/rlhua/project/eCCA/markers/markers_Mast_cells/clustermarker.RDS")
Mast_C2_markers <- dt$gene[dt$cluster==2]
Mast_C3_markers <- dt$gene[dt$cluster==3]
lst_markers <- list(Schwann=schwann_markers,HPLCs=HPLC_markers,Treg=Treg_markers,"Fib3"=Fib_C3_markers,"Mast2"=Mast_C2_markers,"Mast3"=Mast_C3_markers)

##### GSVA
gsva_score <- GSVA::gsva(mtx,lst_markers)

##### clustering
hc <- hclust(dist(t(gsva_score)),method = "ward.D2")
colidx <- colnames(gsva_score)[hc$order]
gsva_score <- gsva_score[,colidx]

##### heatmap
anno_col <- as.data.frame(survival_info[,.(time,status)])
rownames(anno_col) <- survival_info$sample
anno_col$status <- ifelse(anno_col$status==0,"alive","dead")
anno_col$time <- as.character(ceiling(anno_col$time/12))
annotation_colors <- list(time=colorRampPalette(c("#E5FFE5","#006400"))(6),
                          status=c(alive="grey",dead="black"))

p <- pheatmap(gsva_score,border_color = "white",
              color = colorRampPalette(c("blue", "white", "red"))(100),
              cluster_cols =F,
              cluster_rows = T,
              fontsize=5,
              fontsize_row=5,
              scale="none",
              show_colnames=T,
              show_rownames = T,
              cellheight = 5,cellwidth = 5,
              annotation_col=anno_col,
              width=8,height=4,
              legend=T,
              treeheight_row =3,
              annotation_legend = T,
              annotation_colors=annotation_colors,
              clustering_method="ward.D2",
              filename = "/SATA/rlhua/project/eCCA/figure/43eCCA_RNAseq_heatmap.pdf")

##### correlation analysis
dt <- as.data.table(t(combn(rownames(gsva_score),2)))
colnames(dt) <- c("Cell1","Cell2")
lst <- lapply(1:nrow(dt), function(x){
  dt1 <- data.table(cell1=gsva_score[dt$Cell1[x],],cell2=gsva_score[dt$Cell2[x],])
  res <- cor.test(dt1$cell1,dt1$cell2,method = "spearman")
  p <- ggplot(dt1,aes(x=cell1,y=cell2))+geom_point(size=0.2)+
    stat_smooth(method = "loess",level = 0.95)+
    ggtitle(paste0("Cor=",round(res$estimate,3),",p=",round(res$p.value,4)))+
    xlab(dt$Cell1[x])+ylab(dt$Cell2[x])+
    theme(aspect.ratio = 1,text = element_text(size = 5),plot.title = element_text(size = 5))
})
p <- plot_grid(plotlist=lst,ncol = 4,nrow = ceiling(length(lst)/4),align = "hv",axis = "tblr")
fn <- "/SATA/rlhua/project/eCCA/figure/correlation_scatterplot.pdf"
cowplot::save_plot(fn,p,ncol = 4,nrow =ceiling(length(lst)/4),base_height = 1,base_width = 1)

M <- cor(t(gsva_score),method = "spearman")
pdf("/SATA/rlhua/project/eCCA/figure/correlation_corrplot.pdf",height = 1.5,width = 1.5)
corrplot(M,type = "lower",method ="circle", col = colorRampPalette(c("blue", "white", "red"))(100),order = "hclust",hclust.method = "ward.D2",tl.cex = 0.4)
dev.off()
