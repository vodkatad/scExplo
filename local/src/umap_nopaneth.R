paneth <- c('DEFA6','DLL1','ATOH1', 'GFI1', 'DEFA5')

library(Seurat)

umi_saver_f <- '/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_Ire_cetuxi/CRC0322_NT_1_3000_saverAlivePC_G1.csv.gz'


count <- read.table(gzfile(umi_saver_f), sep=",", header=TRUE, row.names=1)

dim(count)

remove_i <- c()
for (i in seq(1, length(paneth))) {
  fi <- grep(paste0(':', paneth[i], '$'), rownames(count))
  if (length(fi) != 0) {
    print(rownames(count)[fi])
    remove_i <- c(remove_i, fi)
  }
}


counts_nopan <- count[-remove_i,]

seu <- CreateSeuratObject(counts_nopan)
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize",scale.factor = 10000)
seu <- FindVariableFeatures(object = seu, mean.function = ExpMean, dispersion.function = LogVMR)#, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = FALSE)

#seu <- ScaleData(object = seu, vars.to.regress = c("nUMI"))
#nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell.
# > summary(seu@meta.data[,3])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14951   14951   14951   14951   14951   14951 
# hm, does not count 0? or do we have all 0 rows only after saver?
seu <- ScaleData(object = seu, vars.to.regress = c("nCount_RNA"))

npc <- 5
seu <- RunPCA(object = seu, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)
seu <- RunUMAP(object = seu, dims = 1:npc, do.fast = TRUE)
#seu <- ProjectPCA(object = seu, do.print = FALSE) # I never did this

DimPlot(object = seu, reduction="umap")

#rownames(counts_nopan)[grep('CACNA2D2', rownames(counts_nopan))]
FeaturePlot(seu, features = c('ENSG00000196169:KIF19'))
FeaturePlot(seu, features = c('ENSG00000007402:CACNA2D2'))

#SERPINA1
FeaturePlot(seu, features = c('ENSG00000197249:SERPINA1'))

pdf('last_ETS1.pdf')
FeaturePlot(seu, features = c('ENSG00000134954:ETS1'))
dev.off()

#https://github.com/satijalab/seurat/issues/4436
seu <- FindNeighbors(seu, dims = 1:npc)
seu <- FindClusters(object = seu, resolution = 0.1)

markers <- FindAllMarkers(seu, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, test.use="negbinom")
write.table(markers, 'last_322NT_markers01.tsv', sep="\t")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf('last_heat.pdf')
DoHeatmap(seu, features = top10$gene) + NoLegend()
dev.off()

pdf('last_dim.pdf')
DimPlot(object = seu, reduction="umap")
dev.off()


df <- data.frame(row.names=names(seu$orig.ident), cluster=seu$seurat_clusters)
write.table(df, file="last_clu.tsv", sep="\t", quote=F)

####

umi_saver_f <- '/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_Ire_cetuxi/CRC0322_cetux_1_saverAlivePC_G1.csv.gz'


count <- read.table(gzfile(umi_saver_f), sep=",", header=TRUE, row.names=1)
G
dim(count)

remove_i <- c()
for (i in seq(1, length(paneth))) {
  fi <- grep(paste0(':', paneth[i], '$'), rownames(count))
  if (length(fi) != 0) {
    print(rownames(count)[fi])
    remove_i <- c(remove_i, fi)
  }
}


counts_nopan <- count[-remove_i,]

seu <- CreateSeuratObject(counts_nopan)
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize",scale.factor = 10000)
seu <- FindVariableFeatures(object = seu, mean.function = ExpMean, dispersion.function = LogVMR)#, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = FALSE)

#seu <- ScaleData(object = seu, vars.to.regress = c("nUMI"))
#nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell.
# > summary(seu@meta.data[,3])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14951   14951   14951   14951   14951   14951 
# hm, does not count 0? or do we have all 0 rows only after saver?
seu <- ScaleData(object = seu, vars.to.regress = c("nCount_RNA"))

npc <- 5
seu <- RunPCA(object = seu, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)
seu <- RunUMAP(object = seu, dims = 1:npc, do.fast = TRUE)
#seu <- ProjectPCA(object = seu, do.print = FALSE) # I never did this

#DimPlot(object = seu, reduction="umap")

#rownames(counts_nopan)[grep('CACNA2D2', rownames(counts_nopan))]
pdf('lastc_KIF19.pdf')
FeaturePlot(seu, features = c('ENSG00000196169:KIF19'))
dev.off()
pdf('lastc_CACNA2D2.pdf')
FeaturePlot(seu, features = c('ENSG00000007402:CACNA2D2'))
dev.off()

#SERPINA1
pdf('lastc_SERPINA1.pdf')
FeaturePlot(seu, features = c('ENSG00000197249:SERPINA1'))
dev.off()

pdf('lastc_ETS1.pdf')
FeaturePlot(seu, features = c('ENSG00000134954:ETS1'))
dev.off()

#https://github.com/satijalab/seurat/issues/4436
seu <- FindNeighbors(seu, dims = 1:npc)
seu <- FindClusters(object = seu, resolution = 0.1)

markers <- FindAllMarkers(seu, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, test.use="negbinom")
write.table(markers, 'lastc_322NT_markers01.tsv', sep="\t")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf('lastc_heatmap.pdf')
DoHeatmap(seu, features = top10$gene) + NoLegend()
dev.off()


df <- data.frame(row.names=names(seu$orig.ident), cluster=seu$seurat_clusters)
write.table(df, file="lastc_clu.tsv", sep="\t", quote=F)


###
###
