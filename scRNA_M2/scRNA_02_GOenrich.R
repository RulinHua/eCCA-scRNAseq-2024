#!/usr/bin/env Rscript
options(warn=-1)
#12/30/20
#Created by Dekang Lv
#USAGE:Rscript cmdfile <gene symbol RDS file>,Rscript GOenrich.R deg.RDS

Args<-commandArgs(T)
degfile <- Args[1]
deg=readRDS(degfile)

###1.get GOIDVsEntrezID which will be used as reference for enrichment in human
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
GOTERM.df <- toTable(GOTERM)
GOTERM.df <- GOTERM.df[, c("go_id", "Term", "Ontology")]
GOTERM.df <- unique(GOTERM.df)
colnames(GOTERM.df)[2]="Functional.Pathway"

#go2gene <- select(org.Hs.eg.db, keys = names(goterms), keytype = "GOALL", columns = c("GOALL", "ENTREZID"))
go2gene <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                                 keys = names(goterms), 
                                 column = "ENTREZID", 
                                 keytype = "GOALL", 
                                 multiVals = "list")
go2gene<-go2gene[sapply(go2gene,function(x){!is.na(x[1])})]
GOIDVsEntrezID<-lapply(go2gene,unique)
#go2gene2=stack(go2gene)
#dim(go2gene2)

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


enrich_GO<-function(GOIDVsEntrezID,deg.symbol,GOTERM){
  enrichment<-enrichment(GOIDVsEntrezID,deg.symbol)
  GO2term.vec <-GOTERM$Functional.Pathway
  names(GO2term.vec) <- GOTERM$go_id
  GO2ontology.vec  <-GOTERM$Ontology
  names(GO2ontology.vec) <- GOTERM$go_id
  enrichment$Functional.Pathway <- GO2term.vec[rownames(enrichment)]
  enrichment$Ontology <- GO2ontology.vec[rownames(enrichment)]
  enrichment<-enrichment[order(enrichment$Ontology),]
  enrichment
}

enrich.GO<-enrich_GO(GOIDVsEntrezID,deg,GOTERM.df)
significant.GO=enrich.GO[enrich.GO$BH.FDR<=0.05,]

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb,sheetName = 'GO_significant')
addWorksheet(wb,sheetName = 'GO_enrichment')

writeData(wb,sheet = "GO_significant",x = significant.GO,rowNames =T)
writeData(wb,sheet = "GO_enrichment",x = enrich.GO,rowNames =T)
saveWorkbook(wb, paste0(basename(degfile),"_GOenrichment.xlsx"), overwrite = TRUE)
