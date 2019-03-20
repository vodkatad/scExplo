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
srdata <- CellCycleScoring(object = srdata, s.genes = s.genes, g2m.genes = g2.genes, set.ident=TRUE)
p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "Phase")
pdf("cycle.pdf")
plot(p1)
graphics.off()
srdata<- RunPCA(object = srdata, pc.genes = c(as.character(s.genes), as.character(g2.genes)), do.print = FALSE)
pdf("cycle_pca.pdf")
PCAPlot(object = srdata)
graphics.off()
srdata <- ScaleData(object = srdata, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
srdata<- RunPCA(object = srdata, pc.genes = c(as.character(s.genes), as.character(g2.genes)), do.print = FALSE)
pdf("cycle_pca_scaled.pdf")
PCAPlot(object = srdata)
graphics.off()


p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "Phase")
pdf("test_TSNEko.pdf")
plot(p1)
graphics.off()
# we need to reset correctly all metrics (PCA and similar) before going on
srdata <- FindVariableGenes(object = srdata, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = pars$x_low, x.high.cutoff = pars$x_high, y.cutoff = pars$y_cutoff) # params defined...?

#http://satijalab.org/seurat/cell_cycle_vignette.html TODO
srdata <- RunPCA(object = srdata, pc.genes = srdata@var.genes, do.print = FALSE)
# long running time: JackStrawPlot(object = srdata, PCs = 1:12)
pdf("elbow_PC.pdf")
PCElbowPlot(object = srdata)
graphics.off()

srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "Phase")
pdf("test_TSNEok.pdf")
plot(p1)
graphics.off()

setwd("..")
save.image(out)
