library(CellChat)
library(Seurat)
library(patchwork)

##### load seurat object
sr <- readRDS("/SATA/rlhua/project/eCCA/seurat_RDSdata/cell_filtered_RDS/expr.RDS")

##### cellchat
cellchat <- createCellChat(object = sr,  group.by = "celltype",assay="RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean",population.size=T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

##### plotting
pdf("/SATA/rlhua/project/eCCA/figure/cellchat_circleplot.pdf")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

dt_prob <- lapply(dimnames(cellchat@net$prob)[[3]], function(x){
  dt <- reshape2::melt(cellchat@net$prob[,,x])
  dt$lr <- x
  return(dt)
}) %>% rbindlist
dt_pval <- lapply(dimnames(cellchat@net$pval)[[3]], function(x){
  dt <- reshape2::melt(cellchat@net$pval[,,x])
  dt$lr <- x
  return(dt)
}) %>% rbindlist
identical(dt_pval$Var1,dt_prob$Var1)
identical(dt_pval$Var2,dt_prob$Var2)
identical(dt_pval$lr,dt_prob$lr)
dt <- data.table(sender=dt_prob$Var1,receiver=dt_prob$Var2,lr=dt_prob$lr,prob=dt_prob$value,pval=dt_pval$value)
dt <- dt[sender %in% c("Fibroblasts 3","Mast cells 2","Treg") & receiver %in% c("Fibroblasts 3","Mast cells 2","Treg")]
dt$sr <- paste0(dt$sender,"_",dt$receiver)
dt$sr <- gsub("Fibroblasts 3","Fib3",dt$sr) %>% gsub("Mast cells 2","Mast2",.)
dt <- lapply(unique(dt$lr), function(x){
  dt1 <- dt[lr==x]
  if (max(dt1$prob)>0) {
    if (dt1$sender[dt1$prob==max(dt1$prob)] != dt1$receiver[dt1$prob==max(dt1$prob)]) {
      dt1 <- dt1[sender!=receiver]
      dt1$scaled_prob <- scale(dt1$prob)
      return(dt1)
    }
  } else {return(NULL)}
}) %>% rbindlist
dt$pval[dt$pval==0] <- 0.01
dt$lr <- paste(cellchat@DB$interaction[dt$lr,"pathway_name"],dt$lr)
dt$lr <- factor(dt$lr,levels = unique(dt$lr[order(dt$lr)]))

ggplot(dt)+
  geom_tile(aes(x=sr,y=lr),fill="#F5F5F5",colour = "white",size=0.2)+
  geom_point(aes(x=sr,y=lr,color=prob,size=-log10(pval)))+
  coord_equal()+
  scale_size_continuous(range = c(0,1.5))+
  scale_color_gradientn(colors = c("blue","yellow","red"))+
  xlab("Sender_Receiver")+ylab("Ligand_Receptor")+labs(color="mean")+
  # ggtitle("Normal")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 5,angle = 45,vjust = 1,hjust =1 ),axis.text.y = element_text(size=5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),panel.background = element_blank(),
        panel.grid =element_blank(),legend.position = "right")
ggsave("/SATA/rlhua/project/eCCA/figure/cellchat_dotplot.pdf",height = 4,width = 5)

