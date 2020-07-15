infile <- snakemake@input[["data"]]
cc <- snakemake@input[["cc"]]
markersf <- snakemake@output[["markers"]]
cyclecl <- snakemake@output[["cyclecl"]]
res <- as.numeric(snakemake@wildcards[['res']])
umapcl <- snakemake@output[['umapcl']]
umapcy <- snakemake@output[['umapcy']]
heatmap <- snakemake@output[['heatmap']]
set.seed(42)

load(infile)
library("Seurat")
library(dplyr)
cycle <- read.table(cc)
s.genes <- cycle[seq(1,43),]
g2.genes <- cycle[seq(44,97),]
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))

srdata <- RunUMAP(srdata, dims=1:npc)
srdata <- FindNeighbors(srdata, dims = 1:npc)
srdata <- FindClusters(srdata, resolution = res)

pdf(umapcl)
DimPlot(srdata, reduction = "umap")
graphics.off()
pdf(umapcy)
DimPlot(srdata, reduction="umap", group.by = "Phase")
graphics.off()

markers <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, test.use="negbinom")
write.table(markers, markersf, sep="\t", quote=FALSE)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf(heatmap)
DoHeatmap(srdata, features = top10$gene) + NoLegend()
graphics.off()

cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)

write.table(cy, cyclecl, quote=FALSE, sep="\t")