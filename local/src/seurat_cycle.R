infile <- snakemake@input[["data"]]
cc <- snakemake@input[["cc"]]
outdir <- snakemake@params[["outdir"]]

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
srdata <- ScaleData(object = srdata, vars.to.regress = c("nGene", "percent.mito","S.Score", "G2M.Score"))
srdata<- RunPCA(object = srdata, pc.genes = c(as.character(s.genes), as.character(g2.genes)), do.print = FALSE)
pdf("cycle_pca_scaled.pdf")
PCAPlot(object = srdata)
graphics.off()
save.image(snakemake@outfile[["outdata"]])
