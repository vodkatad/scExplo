osamples <- snakemake@params[["osamples"]]
infile <- snakemake@input[["data"]]
clfile <- snakemake@input[["clfile"]]
outdir <- snakemake@params[["outdir"]]
clusters <- read.csv(clfile) # "/home/egrassi/ba/dataset/cellranger/aggr/outs/analysis/clustering/kmeans_5_clusters/clusters.csv")

load(infile)
library("Seurat")

setwd(outdir)

# obtain sample from barcodes
samples <- strsplit(colnames(srdata@data), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,5), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- srdata@cell.names
all(colnames(srdata@data)==srdata@cell.names)
srdata <- AddMetaData(object = srdata, metadata =  merged)
pdf("samples_plot.pdf")
p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "sample")
p2 <- PCAPlot(object = srdata, dim.1 = 1, dim.2 = 2, group.by = "sample")
plot_grid(p1, p2)
graphics.off()


pdf("tsne_samples.pdf"); 
p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "sample"); 
plot(p1); 
graphics.off()

pdf("pca_samples.pdf"); 
plot(p2); 
graphics.off()

# by hand I tried another resolution
#srdata <- FindClusters(object = srdata, reduction.type = "pca", dims.use = 1:npc,  resolution = 0.6, print.output = 0, save.SNN = TRUE)
#nrclusters
#p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5)
#pdf("tsne_clusters.pdf"); plot(p1); dev.off()

#pdf("pca_clusters.pdf")
#p2 <- PCAPlot(object = srdata, dim.1 = 1, dim.2 = 2)
#plot(p2)
#graphics.off()

#10x clustering is bad?

rownames(clusters) <- clusters$Barcode
clusters$Barcode <- NULL
srdata <- AddMetaData(object = srdata, metadata =  clusters)

pdf("tsne_10x.pdf"); 
p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "Cluster"); 
plot(p1); 
graphics.off()
