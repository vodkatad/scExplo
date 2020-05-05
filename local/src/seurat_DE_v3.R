infile <- snakemake@input[["data"]]
outfile <- snakemake@output[["out"]]
outimg <- snakemake@output[["outimg"]]
res <- as.numeric(snakemake@wildcards[["res"]])

load(infile)
library("Seurat")
srdata <- FindClusters(object = srdata, resolution = res)
#srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
pdf(outimg)
DimPlot(srdata, reduction="tsne")
graphics.off()
de <- FindAllMarkers(srdata, logfc.threshold=0)
write.table(de, file=outfile, sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
