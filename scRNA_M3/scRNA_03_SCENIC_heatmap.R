#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <seuratRDS> <s3.AUCell.txt> <s5.sigdiff.xlsx>

Args<-commandArgs(trailingOnly = T)
print(Args)
expr_RDS=Args[1]
s3_AUCell_txt=Args[2]
s5_sigdiff_xlsx=Args[3]

top=10
#####

##get significantly high-expressed regulons(top-regulon) from s5.sigdiff.xlsx
library(openxlsx)
s5=read.xlsx(s5_sigdiff_xlsx,sheet=2,rowNames=T)
top10s5=data.frame()
for(clust in sort(unique(s5$cluster))){
  sub=s5[s5$cluster==clust & s5$t>0 ,]
  if(nrow(sub)>=top){sub=head(sub,top)}
  top10s5=rbind(top10s5,sub)
}
top10s5$gene=sub("(+)", "", top10s5$gene, fixed = T)
##get significantly high-expressed regulons from s5.sigdiff.xlsx

##calculate the mean expression of top-regulon-corresponding TF for each cluster
library(Seurat)
print("yes")
expr=readRDS(expr_RDS)
cluster_mean_expression=data.frame()
#cluster_mean_expression_org=data.frame()
#cluster_signal2noise=data.frame()
meta=expr@meta.data
#matexpr=as.matrix(expr@assays$integrated@data)
matexpr=as.matrix(expr@assays$RNA@data)
sort(unique(meta$seurat_clusters))
for(clust in sort(unique(meta$seurat_clusters))){
    print(clust)
    cells=rownames(meta)[meta$seurat_clusters==clust]
    othercells=rownames(meta)[meta$seurat_clusters!=clust]
    print(length(cells))
    print(length(othercells))
    exprmat=matexpr[,cells]
    otherexprmat=matexpr[,othercells]
    rowmean=rowMeans(exprmat)
    otherrowmean=rowMeans(otherexprmat)
    rowmeandf=data.frame(value1=rowmean,value2=otherrowmean)
    colnames(rowmeandf)=c(clust,paste0("not_",clust))
    if(ncol(cluster_mean_expression)==0){
      cluster_mean_expression=rowmeandf
    }else{
      cluster_mean_expression=cbind(cluster_mean_expression,rowmeandf)
    }
}
mean_expr=cluster_mean_expression[,2*(1:(ncol(cluster_mean_expression)/2))-1]
top10s5$gene %in% rownames(mean_expr)
mean_expr_top10=mean_expr[top10s5$gene,]
##calculate the mean expression of regulon-corresponding TF for each cluster

##calulate the mean activity of top-regulon for each cluster
##calculate the mean activity of regulon for each cluster
activity=as.data.frame(t(read.table(s3_AUCell_txt,check.names = F)))
dim((activity))
mean_regulon=data.frame()
for(clust in sort(unique(expr$seurat_clusters))){
  cells=colnames(expr)[expr$seurat_clusters==clust]
  cellactivity=apply(activity[,cells],1,mean)
  print(ncol(mean_regulon))
  if(ncol(mean_regulon)==0){
    #print("d")
    mean_regulon=data.frame(cellactivity)
  }else{
    mean_regulon=cbind(mean_regulon,cellactivity)}
}
head(mean_regulon)
colnames(mean_regulon)=sort(unique(expr$seurat_clusters))
head(mean_regulon)
rownames(mean_regulon)=sub("(+)", "", rownames(mean_regulon), fixed = T)
head(mean_regulon)
##calculate the mean activity of regulon for each cluster
top10s5$gene %in% rownames(mean_regulon)
mean_regulon_top10=mean_regulon[top10s5$gene,]
##calulate the mean activity of top-regulon for each cluster

#####
library(pheatmap)
hig=dim(mean_expr_top10)[1]/10
wid=dim(mean_expr_top10)[2]
pdf("regulon_top10.pdf",width=wid,height = hig)
color1=Seurat::CustomPalette(low = "magenta", high = "yellow", mid = "black",k = 50)
pheatmap(mean_regulon_top10,cluster_rows=F,cluster_cols=F,scale="row",color=color1,
         border_color=NA,fontsize=6)
fun_light2saturated<-function(hue=0.65,
                              sfrom=0,sto=1,
                              vfrom=0.952,vto=2/3,
                              length.out=50){
  saturations=seq(from=sfrom, to=sto, length.out=length.out)
  values=seq(from=vfrom, to=vto, length.out=length.out)
  cols_hex=c()
  for(i in 1:length(values)){
    col_hex=hsv(h=hue,s=saturations[i],v=values[i])
    cols_hex=c(cols_hex,col_hex)
  }
  return(cols_hex)
}
blue2red=c(rev(fun_light2saturated()),fun_light2saturated(hue=0.01))
pheatmap(mean_expr_top10,cluster_rows=F,cluster_cols=F,scale="row",color=blue2red,
         border_color=NA,fontsize=6)
dev.off()

