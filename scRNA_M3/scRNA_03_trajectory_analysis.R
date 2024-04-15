library(Seurat)
library(monocle)

##### load seurat object
sr <- readRDS("/SATA/rlhua/project/eCCA/seurat_RDSdata/Fibroblasts.RDS") 

##### monocle trajectory analysis
data <- as(as.matrix(sr@assays$RNA@counts), 'sparseMatrix')
pdata <- sr@meta.data
pd <- new('AnnotatedDataFrame', data =pdata)
fdata <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fdata)
monocle_cds <- newCellDataSet(data,phenoData = pd,featureData =fd,
                              expressionFamily = negbinomial.size(),
                              lowerDetectionLimit=1)
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1 )

disp_table <- dispersionTable(monocle_cds)
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1*dispersion_fit)$gene_id)
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)

##### plotting
pdf("/SATA/rlhua/project/eCCA/figure/Fibroblasts_trajectory.pdf")
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",cell_size=0.2)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size=0.2)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size=0.2)
p <- cowplot::plot_grid(p1,p2,p3,ncol = 1,nrow = 3,align = "hv",axis = "tblr")
dev.off()

pdf("/SATA/rlhua/project/eCCA/figure/Fibroblasts_trajectory_facet.pdf")
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")+
  facet_wrap(~seurat_clusters)+
  theme(aspect.ratio = 1)
dev.off()