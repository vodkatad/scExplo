osamples <- snakemake@params[["osamples"]]
infile <- snakemake@input[["data"]]
clfile <- snakemake@input[["clfile"]]
outdir <- snakemake@params[["outdir"]]
clusters <- read.csv(clfile) # "/home/egrassi/ba/dataset/cellranger/aggr/outs/analysis/clustering/kmeans_5_clusters/clusters.csv")

load(infile)
library("Seurat")

save.image('pippo.Rdata')
setwd(outdir)

# obtain sample from barcodes
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
#all(colnames(srdata)==srdata@cell.names)
save.image('pippo.Rdata')
srdata[['sample']] <- merged[,"sample"]
pdf("samples_plot.pdf")
#p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "sample")
par(mfrow = c(1, 2))
DimPlot(srdata, reduction="tsne", group.by="sample")
DimPlot(srdata, reduction="pca", group.by="sample")
#p2 <- PCAPlot(object = srdata, dim.1 = 1, dim.2 = 2, group.by = "sample")
#plot_grid(p1, p2)
graphics.off()


#pdf("tsne_samples.pdf"); 
#p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "sample"); 
#plot(p1); 
#graphics.off()

#pdf("pca_samples.pdf"); 
#plot(p2); 
#graphics.off()

# by hand I tried another resolution
srdata <- FindClusters(object = srdata,  resolution = 0.1)
#nrclusters
pdf("tsne_clusters.pdf")

DimPlot(srdata, reduction="tsne")
#pdf("pca_clusters.pdf")
#p2 <- PCAPlot(object = srdata, dim.1 = 1, dim.2 = 2)
#plot(p2)
graphics.off()

#10x clustering is bad?

#rownames(clusters) <- clusters$Barcode
#clusters$Barcode <- NULL
#srdata <- AddMetaData(object = srdata, metadata =  clusters)

#pdf("tsne_10x.pdf"); 
#p1 <- TSNEPlot(srdata, do.return = T, pt.size = 0.5, group.by = "Cluster"); 
#plot(p1); 
#graphics.off()
