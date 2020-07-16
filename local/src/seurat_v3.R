
#this seurat analysis:
#(plot) [kmeans_5_clusters]egrassi@hactarlogin$ seff 30481
#Job ID: 30481
#Cluster: hactar
#User/Group: egrassi/egrassi
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:31:09
#CPU Efficiency: 23.45% of 02:12:51 core-walltime
#Job Wall-clock time: 02:12:51
#Memory Utilized: 11.40 GB
#Memory Efficiency: 72.96% of 15.62 GB

dir_10x <- snakemake@params[["dir"]]
#todo
outdir <- snakemake@params[["outdir"]]
outdata <- snakemake@output[["outdata"]]
name <- snakemake@params[["name"]]
### parameters that needs to be defined by hand after looking at the plots
params <- snakemake@input[["pars"]]
pars <- read.table(params, header=TRUE, sep="\t")

print(dir_10x)
dir.create(outdir)
library("Seurat")
library("dplyr")
data <- Read10X(data.dir = dir_10x)
setwd(outdir)
srdata <- CreateSeuratObject(counts = data, min.cells = pars$min_cells, min.features = pars$min_genes, project = name)
#mito.genes <- grep(pattern = "^MT-", x = rownames(x = srdata@data), value = TRUE)
#percent.mito <- Matrix::colSums(srdata@raw.data[mito.genes, ])/Matrix::colSums(srdata@raw.data)
srdata[["percent.mito"]] <- PercentageFeatureSet(srdata, pattern = "^MT-")
#srdata <- AddMetaData(object = srdata, metadata = percent.mito, col.name = "percent.mito")

pdf("violin.pdf")
VlnPlot(object = srdata, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3)
graphics.off()
pdf("mito_numi.pdf")
par(mfrow = c(1, 2))
#GenePlot(object = srdata, gene1 = "nUMI", gene2 = "percent.mito")
#GenePlot(object = srdata, gene1 = "nUMI", gene2 = "nGene")
FeatureScatter(srdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

graphics.off()

#srdata <- FilterCells(object = srdata, subset.names = c("nGene", "percent.mito"), low.thresholds = c(pars$low_g, pars$low_mito), high.thresholds = c(pars$high_g, pars$high_mito))
srdata <- subset(srdata, subset = nFeature_RNA > pars$low_g & nFeature_RNA < pars$high_g & percent.mito < 25) # fuck fuck fuck!
# fuuuUCK TODO RERUN WITH THIS PARAMETER?s
#srdata <- subset(srdata, subset = percent.mito < pars$high_mito)

# pdf("hist_before_norm.pdf")
# hist(colSums(srdata@data),breaks = 100, main = "Total expression before normalisation",xlab = "Sum of expression")
# graphics.off()
# srdata <- NormalizeData(object = srdata, normalization.method = "LogNormalize", scale.factor=10000)
# pdf("hist_after_norm.pdf")
# hist(colSums(srdata@data),breaks = 100, main = "Total expression after normalisation",xlab = "Sum of expression")
# graphics.off()

srdata <- FindVariableFeatures(srdata, selection.method = "vst", nfeatures = 2000)
#srdata <- FindVariableGenes(object = srdata, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = pars$x_low, x.high.cutoff = pars$x_high, y.cutoff = pars$y_cutoff) # params defined...?

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(srdata), 10)

pdf("hist_before_norm.pdf")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(srdata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
graphics.off()

srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito")) # also cycle here! # step taking more time

#http://satijalab.org/seurat/cell_cycle_vignette.html TODO
srdata <- RunPCA(object = srdata, features=VariableFeatures(object=srdata))
# long running time: JackStrawPlot(object = srdata, PCs = 1:12)
pdf("elbow_PC.pdf")
ElbowPlot(object = srdata)
graphics.off()
npc <- pars$npc

srdata <- FindNeighbors(srdata, dims = 1:npc)
srdata <- FindClusters(srdata, resolution = 1.2)
#srdata <- FindClusters(object = srdata, reduction.type = "pca", dims.use = 1:npc,  resolution = 1.2, print.output = 0, save.SNN = TRUE)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
#p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5)

pdf("clusters.pdf")
par(mfrow = c(1, 2))
#plot_grid(p1, p2)
DimPlot(srdata, reduction="tsne")
DimPlot(object = srdata, reduction="pca")
graphics.off()


setwd("..")
save.image(outdata)
