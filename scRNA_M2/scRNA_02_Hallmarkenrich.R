#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <gene symbol RDS file>,Rscript Hallmarkenrich.R deg.RDS

Args<-commandArgs(T)
degfile <- Args[1]
deg=readRDS(degfile)
library(msigdbr)
hallmarkdf=msigdbr(species = "Homo sapiens", category = "H")
hallmarkdf=hallmarkdf[,c("entrez_gene","gs_name")]
#hallmarkdf
hallmarklist=unstack(hallmarkdf)
hallmarklist=lapply(hallmarklist,function(x){sort(unique(as.character(x)))})

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
    k <- length(overlap)
    M <- length(refset)
    N <- length(AllKEGG.EntrezID) 
    n <- length(geneset)
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
