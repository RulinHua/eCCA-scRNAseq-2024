library(ggplot2)
library(magrittr)
library(Seurat)
library(data.table)
library(survival)
library(survminer)
library(dplyr)

## table of cluster markers
dt_markers <- readRDS("/SATA/rlhua/project/eCCA/markers/markers_Fibroblasts/clustermarker.RDS")

## clinical information
sample_info <- fread("/SATA/rlhua/project/eCCA/bulkRNAseq/sample_info.tsv",sep = "\t",header = T)
sample_info$status <- ifelse(sample_info$status=="Alive",0,1)
sample_info <- sample_info[sample_name!=""]

## RNA expression matrix
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

## survival analysis
clusters <- sort(unique(dt_markers$cluster))
lst <- list()
for (i in 1:length(clusters)) {
  genes <- dt_markers$gene[dt_markers$cluster==clusters[i]]
  gsva_score <- GSVA::gsva(mtx,list(genes))
  dt <- sample_info
  dt$gsva_score <- gsva_score[,dt$sample]
  dt$group <- ifelse(dt$gsva_score > median(dt$gsva_score),"high","low")
  dt$group <- factor(dt$group,levels = c("high","low"))
  fit <- survfit(Surv(time,status)~group,data = dt)
  data.survdiff <- survdiff(Surv(time,status)~group,data = dt)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  ## Kaplan-Meier survival curve
  lst[[clusters[i]]] <- ggsurvplot(fit,
                                   conf.int = F, 
                                   censor = T,
                                   censor.size =2,
                                   size=0.5,
                                   palette = c("#B2182B","#2166AC"),
                                   xlab="Months",
                                   legend=c(0.7,0.95),
                                   title=paste0("cluster ",clusters[i]),
                                   legend.title = "",
                                   font.legend = 6,
                                   legend.labs=c("group=high","group=low"),
                                   pval = paste("CCA(N=43)", paste("P = ",round(p.val,3), sep = ""),sep = "\n"),
                                   pval.size=2,pval.coord=c(0,0.1),
                                   ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size=6),axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6),axis.title.x = element_text(size=6),axis.title.y = element_text(size = 6),axis.line = element_line(size = 0.3),axis.ticks = element_line(size=0.3)))
}
p <- arrange_ggsurvplots(lst, print = T,ncol = 5, nrow = 1)
ggsave("/SATA/rlhua/project/eCCA/figure/Fibroblasts_subtype_survplot.pdf",p,width = 7.5,height = 1.5)
