#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <dataframe for markers RDS file>

Args<-commandArgs()
print(Args)
markerdf=readRDS(Args[6])
cmdfile=strsplit( Args[4], split="=" )[[1]][2]
cmddir=dirname(normalizePath(cmdfile))
print(cmddir)
pwd=getwd()
savePath=paste0(pwd,"/M4markerEnrichment_",gsub("\\..*","",basename(Args[6])))
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}
setwd(savePath)

for(clust in unique(markerdf$cluster)){
  markers=markerdf$gene[markerdf$cluster==clust]
  filename=paste0("markers_",clust,".RDS")
  saveRDS(unique(markers),file = filename)
  GOenrichfile=paste0(cmddir,"/human_overrepresentation_analysis/GOenrich.R")
  system(paste("Rscript",GOenrichfile,filename))
  KEGGenrichfile=paste0(cmddir,"/human_overrepresentation_analysis/KEGGenrich.R")
  system(paste("Rscript",KEGGenrichfile,filename))
  Hallmarkenrichfile=paste0(cmddir,"/human_overrepresentation_analysis/Hallmarkenrich.R")
  system(paste("Rscript",Hallmarkenrichfile,filename))
  #Reactomeenrichfile=paste0(cmddir,"/human_overrepresentation_analysis/Reactomeenrich.R")
  #system(paste("Rscript",Reactomeenrichfile,filename))
  system(paste("rm",filename))
}
