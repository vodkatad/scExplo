setwd('/mnt/trcanmed/snaketree/prj/scRNA/dataset/scRNA_Ire_cetuxi')
library(Seurat)
cycle <- read.table('../../local/share/data//regev_lab_cell_cycle_genes.txt')
library(ggplot2)
library(dplyr)
s.genes <- cycle[seq(1,43),]
g2.genes <- cycle[seq(44,97),]
set.seed(42)

load("Sample_S26093_499_seurat/prelim.Rdata")
title <- 'CRC0327_NT_1'


srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
nt1 <- c2
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)


load("Sample_S26369_503_seurat/prelim.Rdata")
title <- 'CRC0327_NT_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
nt2 <- c2
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)

load("Sample_S26370_504_seurat/prelim.Rdata")
title <- 'CRC0327_Cetuxi_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
cetuxi2 <- c2
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)

load('Sample_S26094_500_seurat/prelim.Rdata')
title <- 'CRC0327_Cetuxi_1'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)
cetuxi1 <- c2

### 1502

load("Sample_S26371_505_seurat/prelim.Rdata")
title <- 'CRC1502_NT_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
nt2 <- c2
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)

load("Sample_S26371_506_seurat/prelim.Rdata")
title <- 'CRC1502_Cetuxi_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
cetuxi2 <- c2
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)

nt2$treat <- 'nt2'
cetuxi2$treat <- 'cetuxi2'
d <- rbind(nt2, cetuxi2)
d$x <- paste0(d$Var1, d$treat)

ggplot(d, aes(x=x, y=Freq, fill=Var1))+ geom_bar(stat = "identity")+theme_bw()+scale_x_discrete(labels=c('cetuxi','nt','cetuxi','nt','cetuxi','nt'))+theme(text = element_text(size=20))
mat <- matrix(d$Freq, byrow=F, ncol=2)
fisher.test(mat)$p.value


### 322 but gene squilibria here

load("Sample_S26095_501_seurat/prelim.Rdata")
title <- 'CRC0322_NT_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
nt2 <- c2
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)

load("Sample_S26096_502_seurat/prelim.Rdata")
title <- 'CRC0322_Cetuxi_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))
ElbowPlot(object = srdata)
srdata <- RunTSNE(object = srdata, dims.use = 1:npc, do.fast = TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
c2 <- as.data.frame(table(srdata$Phase))
cetuxi2 <- c2
ggplot(c2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+theme_bw()+ggtitle(title)

## NT_1 done "by hand"
load("Sample_S26369_503_seurat/prelim.Rdata")
title <- 'CRC0327_NT_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")

DimPlot(srdata, reduction = "umap",group.by="Phase")
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)
table(cy[cy$srdata.seurat_clusters==5,'cycle'])
table(cy[cy$srdata.seurat_clusters!=5,'cycle'])


markers <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")

markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
head(markers)
dim(markers[markers$p_val_adj < 0.05,])
m <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
m$gene
FeaturePlot(srdata, features = m$gene )
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(srdata, features = top10$gene) + NoLegend()
##
load("Sample_S26371_506_seurat/prelim.Rdata")
title <- 'CRC1502_Cetuxi_2'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")

DimPlot(srdata, reduction = "umap",group.by="Phase")
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)
table(cy[cy$srdata.seurat_clusters==5,'cycle'])
table(cy[cy$srdata.seurat_clusters!=5,'cycle'])


markers <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")
library(dplyr)
markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
head(markers)
dim(markers[markers$p_val_adj < 0.05,])
m <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
m$gene
FeaturePlot(srdata, features = m$gene )
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(srdata, features = top10$gene) + NoLegend()

####start here
###327 replicate 1
load("Sample_S26093_499_seurat/prelim.Rdata")
title <- 'CRC0327_NT_1'

srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction = "umap",group.by="Phase")
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)
table(cy[cy$srdata.seurat_clusters==5,'cycle'])
table(cy[cy$srdata.seurat_clusters!=5,'cycle'])

markers327nt1<- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")


load("Sample_S26094_500_seurat/prelim.Rdata")
title <- 'CRC0327_Cetuxi_1'

srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction = "umap",group.by="Phase")
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)
table(cy[cy$srdata.seurat_clusters==5,'cycle'])
table(cy[cy$srdata.seurat_clusters!=5,'cycle'])

markers327c1 <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")

### 327 replicate 2

load("Sample_S26369_503_seurat/prelim.Rdata")
title <- 'CRC0327_NT_2'

srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction = "umap",group.by="Phase")
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)
table(cy[cy$srdata.seurat_clusters==5,'cycle'])
table(cy[cy$srdata.seurat_clusters!=5,'cycle'])

markers327nt2 <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")


load("Sample_S26370_504_seurat/prelim.Rdata")
title <- 'CRC0327_Cetuxi_2'

srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction = "umap",group.by="Phase")
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)
table(cy[cy$srdata.seurat_clusters==5,'cycle'])
table(cy[cy$srdata.seurat_clusters!=5,'cycle'])

markers327c2 <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")


##

# "Sample_S26095_501","Sample_S26096_502"
#"CRC0322_NT_1", "CRC0322_cetux_1",
load("Sample_S26095_501_seurat/prelim.Rdata")
title <- 'CRC0322_NT_1'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction="umap", group.by = "Phase")
markers322nt <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")
top10 <- markers322nt %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(srdata, features = top10$gene) + NoLegend()
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)

dim(cy[cy$srdata.seurat_clusters==3,])
dim(cy)

table(cy[cy$srdata.seurat_clusters==3,'cycle'])
table(cy[cy$srdata.seurat_clusters!=3,'cycle'])
###

load("Sample_S26096_502_seurat/prelim.Rdata")
title <- 'CRC0322_CETUXI_1'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction="umap", group.by = "Phase")


markers322c <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")

top10 <- markers322c %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(srdata, features = top10$gene) + NoLegend()
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)

dim(cy[cy$srdata.seurat_clusters==3,])
dim(cy)

table(cy[cy$srdata.seurat_clusters==3,'cycle'])
table(cy[cy$srdata.seurat_clusters!=3,'cycle'])

set.seed(42)
### 1502
#"Sample_S26371_505","Sample_S26371_506"
# "CRC1502_NT_1", "CRC1502_cetux_1"
#defensina/ATOH1 su e AREG/EREG 



load("Sample_S26371_505_seurat/prelim.Rdata")
title <- 'CRC1502_NT_1'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction="umap", group.by = "Phase")
markers <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")
markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
head(markers)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(srdata, features = top10$gene) + NoLegend()
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)

dim(cy[cy$srdata.seurat_clusters==4,])
dim(cy)

table(cy[cy$srdata.seurat_clusters==4,'cycle'])
table(cy[cy$srdata.seurat_clusters!=4,'cycle'])
###


load("Sample_S26371_506_seurat/prelim.Rdata")
title <- 'CRC1502_Cetuxi_1'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction="umap", group.by = "Phase")


markers1502c <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")
top10 <- markers1502c %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(srdata, features = top10$gene) + NoLegend()
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)

dim(cy[cy$srdata.seurat_clusters==3,])
dim(cy)

table(cy[cy$srdata.seurat_clusters==3,'cycle'])
table(cy[cy$srdata.seurat_clusters!=3,'cycle'])

###

markers[markers$cluster==3 & markers$gene == 'ATOH1',]
markers[markers$cluster==3 & markers$gene == 'DEFA6',]
markers[markers$cluster==3 & markers$gene == 'AREG',]
markers[markers$cluster==3 & markers$gene == 'EREG',]
markers[markers$cluster==3 & markers$gene == 'MYC',]
markers[markers$cluster==3 & markers$gene == 'DLL1',]
markers[markers$cluster==3 & markers$gene == 'SOX4',]

markers[markers$cluster==4 & markers$gene == 'ATOH1',]
markers[markers$cluster==4 & markers$gene == 'DEFA6',]
markers[markers$cluster==4 & markers$gene == 'AREG',]
markers[markers$cluster==4 & markers$gene == 'EREG',]
markers[markers$cluster==4 & markers$gene == 'MYC',]
markers[markers$cluster==4 & markers$gene == 'DLL1',]
markers[markers$cluster==4 & markers$gene == 'SOX4',]

## 327 bad



load("Sample_S26094_500_seurat/prelim.Rdata")
title <- 'CRC0327_cetuxi_1'
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
DimPlot(srdata, reduction="tsne", group.by = "Phase")
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)
srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)
DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction="umap", group.by = "Phase")
markers327c <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")
top10 <- markers327c %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(srdata, features = top10$gene) + NoLegend()
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)

dim(cy[cy$srdata.seurat_clusters==4,])
dim(cy)

table(cy[cy$srdata.seurat_clusters==4,'cycle'])
table(cy[cy$srdata.seurat_clusters!=4,'cycle'])

####

load('aggr_seurat/prelim.Rdata')
srdata <- CellCycleScoring(object = srdata, s.features = s.genes, g2m.features = g2.genes, set.ident=TRUE)
srdata <- ScaleData(object = srdata, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
srdata <- RunPCA(object = srdata, features=VariableFeatures(srdata))

srdata <- RunUMAP(srdata, dims=1:15)

srdata <- FindNeighbors(srdata, dims = 1:15)
srdata <- FindClusters(srdata, resolution = 0.2)

DimPlot(srdata, reduction = "umap")
DimPlot(srdata, reduction="umap", group.by = "Phase")
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)
write.table(cy, file='aggr_clusters_0.2.tsv', sep="\t", quote=FALSE)
srdata <- FindClusters(srdata, resolution = 0.08)
cy <- data.frame(cycle=srdata$Phase, srdata$seurat_clusters)

write.table(cy, file='aggr_clusters_0.08.tsv', sep="\t", quote=FALSE)
markers <- FindAllMarkers(srdata, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use="negbinom")
m <- as.data.frame(markers)
write.table(m, file="aggr_cluster_markerks008.tsv", quote=F, sep='\t')

