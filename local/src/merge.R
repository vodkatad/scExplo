paneth <- read.table('paneth_cells_cycle.tsv', sep="\t", stringsAsFactors = F)
colnames(paneth) <- c('id','cycle','cl','sample')
paneth$paneth <- "yes"
paneth$barcode <- sapply(samples, function(x) {x[1]})
samples <- strsplit(paneth$id, '-', fixed=T)
paneth$barcode <- sapply(samples, function(x) {x[1]})
table(paneth$sample)
paneth[paneth$sample == "CRC0327_NT_1",'model'] = 1
paneth[paneth$sample == "CRC0327_cetux_1",'model'] = 2
paneth[paneth$sample == "CRC0322_NT_1",'model'] = 3
paneth[paneth$sample == "CRC0322_cetux_1",'model'] = 4
paneth[paneth$sample == "CRC0327_NT_2",'model'] = 5
paneth[paneth$sample == "CRC0327_cetux_2",'model'] = 6
paneth[paneth$sample == "CRC1502_NT_1",'model'] = 7
paneth[paneth$sample == "CRC1502_cetux_1",'model'] = 8
table(paneth$model)
head(paneth)
paneth$cellid <- paste0(paneth$barcode, '-', paneth$model)
head(paneth)
tail(paneth)
mm <- data.frame(id=colnames(srdata))
m <- merge(paneth, mm, by.y="id", by.x="cellid", all.y=T, sort=F)
head(m)
mmm <- m[match(mm$id, m$cellid),]
all(mmm$cellid==colnames(srdata))
srdata[['paneth']] <- mmm[,"ppaneth"]
DimPlot(srdata, reduction="umap", group.by="paneth")
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
samples_id <- data.frame(id=s)
osamples = c("CRC0327_NT_1", "CRC0327_cetux_1", "CRC0322_NT_1", "CRC0322_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2", "CRC1502_NT_1", "CRC1502_cetux_1"])
osamples = c("CRC0327_NT_1", "CRC0327_cetux_1", "CRC0322_NT_1", "CRC0322_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2", "CRC1502_NT_1", "CRC1502_cetux_1")
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
srdata[['model']] <- merged[,"sample"]
DimPlot(srdata, reduction="umap", group.by="model")
DimPlot(srdata, reduction="tsne", group.by="model")
DimPlot(srdata, reduction="tsne", group.by="paneth")
DimPlot(srdata, reduction="tsne", group.by="cycle")
DimPlot(srdata, reduction="tsne", group.by="Phase")
DimPlot(srdata, reduction="tsne")
head(markers_aggr)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10 <- markers_aggr %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(srdata, features = top10$gene) + NoLegend()
DoHeatmap(srdata, features = top10$gene) + NoLegend()
FeaturePlot(srdata, features = c('SOX4','ATOH1') )
VlnPlot(srdata, features = c('SOX4','ATOH1') )
png('heat.png');DoHeatmap(srdata, features = top10$gene) + NoLegend(); dev.off()
top10
top10[top10$cluster=="7","gene"]
markers_aggr[markers_aggr$cluster=="7" & markers_aggr$p_val_adj<0.001,]
ma <- markers_aggr[markers_aggr$cluster=="7" & markers_aggr$p_val_adj<0.001,]
ma[ma$gene =="EREG",]
ma[ma$gene =="AREG",]
ma[ma$gene =="MYC",]
ma[ma$gene =="DLL1",]
ma[ma$gene =="DEFA6",]
ma[ma$gene =="DEFA5",]
ma[order(ma$p_val_adj),]




### 327


dir_10x <- "327_aggr/outs/filtered_feature_bc_matrix"
pars <- read.table("../../local/share/data/seurat_params.txt", header=TRUE, sep="\t")
pars
srdata <- CreateSeuratObject(counts = data, min.cells = pars$min_cells, min.features = 100, project = 'CRC0327')
data <- Read10X(data.dir = dir_10x)
srdata <- CreateSeuratObject(counts = data, min.cells = pars$min_cells, min.features = 100, project = 'CRC0327')
srdata
srdata[["percent.mito"]] <- PercentageFeatureSet(srdata, pattern = "^MT-")
VlnPlot(object = srdata, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3)
srdata <- FindVariableFeatures(srdata, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srdata), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot1 <- VariableFeaturePlot(srdata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
params
pars
srdata_f <- subset(srdata, subset = nFeature_RNA > pars$low_g & nFeature_RNA < pars$high_g & percent.mito < 25)
srdata
srdata_f
srdata_f <- FindVariableFeatures(srdata, selection.method = "vst", nfeatures = 2000)
cycle <- read.table('../../local/share/data//regev_lab_cell_cycle_genes.txt')
s.genes <- cycle[seq(1,43),]
g2.genes <- cycle[seq(44,97),]
set.seed(42)
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
srdata_orig <- srdata
srdata <- srdata_f
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
npc <- 10
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction="umap", group.by = "Phase")
markers327c <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")
top10 <- markers327c %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(srdata, features = top10$gene) + NoLegend()
#DimPlot(srdata, reduction="umap", group.by = "percent.mito")

FeaturePlot(srdata, features = 'percent.mito' )
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
samples_id <- data.frame(id=s)
osamples <- c("CRC0327_NT_1", "CRC0327_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2")
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
head(merged)
srdata[['sample']] <- merged[,"sample"]
DimPlot(srdata, reduction="pca", group.by="sample")
DimPlot(srdata, reduction="umap", group.by="sample")

### new aggr test
load('sacro_lento_aggr.Rdata')
FeaturePlot(srdata, features = 'percent.mito' )
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters, srdata$percent.mito)
library(ggplot)
library(ggplot2)
ggplot(cy, aes(seurat_clusters, percent.mito, fill=seurat_clusters)) +
geom_boxplot() +theme_bw()+ggtitle("Mitochondrial read % in clusters")
head(cy)
ggplot(cy, aes(srdata.seurat_clusters, srdata.percent.mito, fill=srdata.seurat_clusters)) +
geom_boxplot() +theme_bw()+ggtitle("Mitochondrial read % in clusters")
head(colnames(srdata))
head(rownames(srdata))
head(names(srdata))
srdata
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters, srdata$percent.mito,srdata$sample)
paneth <- read.table('paneth_cells_cycle.tsv', sep="\t", stringsAsFactors = F)
colnames(paneth) <- c('id','cycle','cl','sample')
paneth$paneth <- "yes"
paneth$barcode <- sapply(samples, function(x) {x[1]})
samples <- strsplit(paneth$id, '-', fixed=T)
paneth$barcode <- sapply(samples, function(x) {x[1]})
table(paneth$sample)
paneth[paneth$sample == "CRC0327_NT_1",'model'] = 1
paneth[paneth$sample == "CRC0327_cetux_1",'model'] = 2
paneth[paneth$sample == "CRC0322_NT_1",'model'] = 3
paneth[paneth$sample == "CRC0322_cetux_1",'model'] = 4
paneth[paneth$sample == "CRC0327_NT_2",'model'] = 5
paneth[paneth$sample == "CRC0327_cetux_2",'model'] = 6
paneth[paneth$sample == "CRC1502_NT_1",'model'] = 7
paneth[paneth$sample == "CRC1502_cetux_1",'model'] = 8
table(paneth$model)
head(paneth)
paneth$cellid <- paste0(paneth$barcode, '-', paneth$model)
head(paneth)
tail(paneth)
m <- merge(paneth, mm, by.y="id", by.x="cellid", all.y=T, sort=F)
head(m)
mmm <- m[match(mm$id, m$cellid),]
all(mmm$cellid==colnames(srdata))
m <- data.frame(id=colnames(srdata))
mm <- data.frame(id=colnames(srdata))
m <- merge(paneth, mm, by.y="id", by.x="cellid", all.y=T, sort=F)
mmm <- m[match(mm$id, m$cellid),]
all(mmm$cellid==colnames(srdata))
srdata[['paneth']] <- mmm[,"ppaneth"]
head(mmm)
mmm$ppaneth <- "no"
mmm[!is.na(mmm$cl) ]$ppaneth <- "no"
mmm[!is.na(mmm$cl), ]$ppaneth <- "no"
table(mmm$ppaneth)
head(paneth)
tail(mmm)
mmm[!is.na(mmm$paneth), ]$ppaneth <- "yes"
table(mmm$ppaneth)
osamples
osamples <- c("CRC0327_NT_1", "CRC0327_cetux_1", "CRC0322_NT_1", "CRC0322_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2", "CRC1502_NT_1", "CRC1502_cetux_1")
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
srdata[['sample']] <- merged[,"sample"]
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters, srdata$percent.mito,srdata$sample)
ggplot(cy, aes(srdata.sample, srdata.percent.mito, fill=srdata.sample)) +
geom_boxplot() +theme_bw()+ggtitle("Mitochondrial read % in clusters")
ggplot(cy, aes(srdata.sample, srdata.percent.mito, fill=srdata.sample)) +
geom_boxplot() +theme_bw()+ggtitle("Mitochondrial read % in samples")


