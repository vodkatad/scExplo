library(Seurat)
in_f <- snakemake@input[[1]]

names <- snakemake@params[['sign']]

d <- read.csv(in_f, row.names=1)

seu <- CreateSeuratObject(counts = as.matrix(d))
seu <- NormalizeData(object = seu)
#seu <- LogNormalize(object = seu)
seu <- FindVariableFeatures(object = seu)
seu <- ScaleData(object = seu)

#seu <- RunPCA(object = seu)
#DimPlot(object = seu, reduction = "pca")
#seu <- FindNeighbors(object = seu, dims = 1:30) # ARBIT
#seu <- FindClusters(object = seu)
#seu <- RunUMAP(object = seu, dims = 1:30)  # ARBIT
#DimPlot(object = seu, reduction = "umap")


slist <- list()
for (i in seq(2, length(snakemake@input))) {
    sign_f <- snakemake@input[[i]]
    sonf <- read.table(sign_f, header=TRUE, sep="\t")
    slist[[i-1]] <- sonf$human_symbol
}
names(slist) <- names
save.image('p.Rdata')

seu <- AddModuleScore(seu, features=slist, name=names(slist))

names_n <- paste0(names, seq(1, length(names)))

#FeaturePlot(object = seu, features = "Onf1", reduction="pca")
#FeaturePlot(object = seu, features = "Stem2", reduction="pca")

res <- FetchData(seu, names_n)
write.table(res, file=snakemake@output[['sign']], sep="\t", quote=FALSE)
#ggplot(data=pdf, aes(x=stem, y=onf))+geom_point()
