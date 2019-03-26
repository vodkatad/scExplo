infile <- snakemake@input[["data"]]
outfile <- snakemake@output[["out"]]
outimg <- snakemake@output[["outimg"]]
res <- as.numeric(snakemake@wildcards[["res"]])

load(infile)
library("Seurat")
srdata <- FindClusters(object = srdata, reduction.type = "pca", dims.use = 1:npc,  resolution = res, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
pdf(outimg)
p1 <- TSNEPlot(srdata, do.return = TRUE, pt.size = 0.5)
plot(p1)
graphics.off()
de <- FindAllMarkers(srdata, return.thresh=1, logfc.threshold=0)
write.table(de, file=outfile, sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
