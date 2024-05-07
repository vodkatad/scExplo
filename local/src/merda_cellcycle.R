load('CRC0322_cetux_1_Seurat_prelim.Rdata')
library(Seurat)
cc <- '../../local/share/data/regev_lab_cell_cycle_genes.txt'
cycle <- read.table(cc)
s.genes <- cycle[seq(1,43),]
g2.genes <- cycle[seq(44,97),]
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
pdf('merda.pdf')
PCAPlot(srdata)
dev.off()

res <- 0.2

pdf('merda3.pdf')
srdata <- RunUMAP(srdata, dims=1:npc)
DimPlot(srdata, reduction = "umap")
graphics.off()

pdf('merda01.pdf')
FeaturePlot(srdata, features =  'ATOH1')
dev.off()

srdata <- FindNeighbors(srdata, dims = 1:npc)
srdata <- FindClusters(srdata, resolution = res)

pdf('merda4.pdf')
DimPlot(srdata, reduction = "umap")
graphics.off()
pdf('merda5.pdf')
DimPlot(srdata, reduction="umap", group.by = "Phase")
graphics.off()

cy0 <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters, mito=srdata$percent.mito)

srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:npc)
srdata <- RunTSNE(srdata, dims=1:npc)
srdata <- FindNeighbors(srdata, dims = 1:npc)
srdata <- FindClusters(srdata, resolution = res)

pdf('merda6.pdf')
DimPlot(srdata, reduction = "umap")
graphics.off()
pdf('merda7.pdf')
DimPlot(srdata, reduction="umap", group.by = "Phase")
graphics.off()


pdf('merda8.pdf')
TSNEPlot(srdata, group.by = "Phase")
graphics.off()

pdf('merda9.pdf')
PCAPlot(srdata, group.by = "Phase")
graphics.off()

pdf('merda10.pdf')
FeaturePlot(srdata, features =  'ATOH1')
dev.off()
cy1 <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters, mito=srdata$percent.mito)

srdata <- AddMetaData(
  object = srdata,
  metadata = cy0$srdata.seurat_clusters, col.name= "not_corr_cl")
  
pdf('merda12.pdf')
DimPlot(srdata, reduction="umap", group.by = "not_corr_cl")
graphics.off()


load('CRC0322_cetux_1_Seurat_prelim.Rdata')
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
res <- 0.2
srdata <- RunUMAP(srdata, dims=1:npc)

srdata <- AddMetaData(
  object = srdata,
  metadata = cy1$srdata.seurat_clusters, col.name= "corr_cl")
  

pdf('merda13.pdf')
DimPlot(srdata, reduction="umap", group.by = "corr_cl")
graphics.off()

