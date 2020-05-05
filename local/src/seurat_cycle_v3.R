infile <- snakemake@input[["data"]]
cc <- snakemake@input[["cc"]]
outdir <- snakemake@params[["outdir"]]
out <- snakemake@output[["outdata"]]
#if (dir.exists(outdir)) {
#    stop(paste0(outdir, "outdir already exists for seurat_cycle!"))
#}
# do snakemake creates it?
dir.create(outdir)
load(infile)
library("Seurat")
cycle <- read.table(cc)
setwd(outdir)
s.genes <- cycle[seq(1,43),]
g2.genes <- cycle[seq(44,97),]
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
pdf("cycle.pdf")
DimPlot(srdata, reduction="tsne", group.by = "Phase")
graphics.off()
srdata<- RunPCA(object = srdata, features= c(as.character(s.genes), as.character(g2.genes)))
pdf("cycle_pca.pdf")
DimPlot(object = srdata, reduction="pca")
graphics.off()
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata<- RunPCA(object = srdata, features = c(as.character(s.genes), as.character(g2.genes)))
pdf("cycle_pca_scaled.pdf")
DimPlot(object = srdata, reduction="pca")
graphics.off()


pdf("test_TSNEko.pdf")
DimPlot(srdata, reduction='tsne', group.by = "Phase")
graphics.off()

# UPTO
# we need to reset correctly all metrics (PCA and similar) before going on

#srdata <- FindVariableGenes(object = srdata, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = pars$x_low, x.high.cutoff = pars$x_high, y.cutoff = pars$y_cutoff) # params defined...?

#http://satijalab.org/seurat/cell_cycle_vignette.html TODO
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
# long running time: JackStrawPlot(object = srdata, PCs = 1:12)
pdf("elbow_PC.pdf")
ElbowPlot(object = srdata)
graphics.off()

srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
pdf("test_TSNEok.pdf")
DimPlot(srdata, reduction="tsne", group.by = "Phase")
graphics.off()

setwd("..")
save.image(out)
