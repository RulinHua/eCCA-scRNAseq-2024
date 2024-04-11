#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <参考文件>,<差异基因symbol向量RDS文件>,Rscript Hallmarkenrich.R deg.RDS
#功能：超几何分布筛选富集Hallmark通路
#input：<差异基因symbol向量RDS文件>向量
#output：exel文件包含两个表,一个是显著富集的结果,另一个是所有检验结果，可以根据需要重新设置差异阈值

Args<-commandArgs(T)
degfile <- Args[1]
deg=readRDS(degfile)
library(msigdbr)
hallmarkdf=msigdbr(species = "Homo sapiens", category = "H")
hallmarkdf=hallmarkdf[,c("entrez_gene","gs_name")]
#hallmarkdf
hallmarklist=unstack(hallmarkdf)
hallmarklist=lapply(hallmarklist,function(x){sort(unique(as.character(x)))})
#基于超几何分布的富集函数
enrichment<-function(pathIDVsEntrezID,deg.symbol,low=10,high=500){
  IDtrans<-function(symbol){
    library(org.Hs.eg.db)
    #keytypes(org.Hs.eg.db)
    symbol2entrezID.df=AnnotationDbi::select(org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
    symbol2entrezID.df=symbol2entrezID.df[!is.na(symbol2entrezID.df$ENTREZID),]
    symbol2entrezID.df=unique(symbol2entrezID.df)
    return(symbol2entrezID.df)
  }
  
  deg.trans.df<- IDtrans(deg.symbol)
  deg.trans.vec <- deg.trans.df$SYMBOL
  names(deg.trans.vec)<- deg.trans.df$ENTREZID
  queryset<-deg.trans.df$ENTREZID
  AllKEGG.EntrezID<-unique(unlist(pathIDVsEntrezID))
  geneset<-intersect(queryset,AllKEGG.EntrezID)
  pathway.length<-sapply(pathIDVsEntrezID,length)
  pathIDVsEntrezID<-pathIDVsEntrezID[pathway.length>=low & pathway.length<=high]
  statistics=lapply(pathIDVsEntrezID,function(x){
    refset=x
    overlap <- intersect(geneset,refset)
    overlap.entrez <-paste(overlap,collapse = ",")
    k <- length(overlap)# 差异基因中属于hsa pathway的基因个数-1
    M <- length(refset)# 在背景基因下 hsa pathway的基因个数
    N <- length(AllKEGG.EntrezID) # 背景基因个数
    n <- length(geneset)# 差异基因个数
    pvalues <-  phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    fold<- (k/n)/(M/N)
    c(N,M,n,k,fold,pvalues)
  })
  statistics.df<-as.data.frame(do.call("rbind",statistics))
  colnames(statistics.df)<-c("#genes.in.background",
                             "#genes.in.the.pathway",
                             "#genes.in.the.query",
                             "#genes.overlapped",
                             "Enrichment.fold",
                             "P.value")
  statistics.df$BH.FDR<-p.adjust(statistics.df$P.value,method = "BH")
  overlap.ID=lapply(pathIDVsEntrezID,function(x){
    refset=x
    overlapEntrezID <- intersect(geneset,refset)
    overlapSymbol <- deg.trans.vec[overlapEntrezID]
    overlap.entrez <- paste(overlapEntrezID,collapse = ",")
    overlap.symbol <- paste(overlapSymbol,collapse = ",")
    c(overlap.entrez,overlap.symbol)
  })
  overlap.df<-as.data.frame(do.call("rbind",overlap.ID))
  colnames(overlap.df)<-c("EntrezID.overlapped",
                          "Genes.overlapped")
  res.df<-cbind(statistics.df,overlap.df)
  res.df
}

enrich.hallmark<-enrichment(hallmarklist,deg)
enrich.hallmark$Functional.Pathway <- sub("HALLMARK_","",rownames(enrich.hallmark))
significant.hallmark=enrich.hallmark[enrich.hallmark$BH.FDR<=0.05,]
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb,sheetName = 'Hallmark_significant')
addWorksheet(wb,sheetName = 'Hallmark_enrichment')

writeData(wb,sheet = "Hallmark_significant",x = significant.hallmark,rowNames =T)
writeData(wb,sheet = "Hallmark_enrichment",x = enrich.hallmark,rowNames =T)
saveWorkbook(wb, paste0(basename(degfile),"Hallmarkenrichment.xlsx"), overwrite = TRUE)
