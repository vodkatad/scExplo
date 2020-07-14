load('/mnt/trcanmed/snaketree/prj/scRNA/dataset/scRNA_Ire_cetuxi/aggr-cc_seurat/prelim.Rdata')
srdata
library(seurat)
library(Seurat)
FeaturePlot(srdata, features=c('SOX9'));
FeaturePlot(srdata, features=c('TP53INP2'));
FeaturePlot(srdata, features=c('CORO1A'));
FeaturePlot(srdata, features=c('FN1'));
str(srdata)
head(srdata[['nFeature_RNA']])
genes <- srdata[['nFeature_RNA']]
genes$sample <- sapply(rownames(genes), function(x) {strsplit(x)[[1]][2]})
genes$sample <- sapply(rownames(genes), function(x) {strsplit(x,'-')[[1]][2]})
head(genes)
dim(genes)
library(ggplot2)
ggplot(genes, aes(x=nFeature_RNA, y=sample)) +
geom_boxplot()
ggplot(genes, aes(y=nFeature_RNA, x=sample)) +
geom_boxplot()
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
osamples <- c("CRC0327_NT_1", "CRC0327_cetux_1", "CRC0322_NT_1", "CRC0322_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2", "CRC1502_NT_1", "CRC1502_cetux_1")
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
head(merged)
m <- merge(genes, merged)
head(m)
head(genes)
m <- merge(genes, merged, by="row.names")
head(m)
ggplot(genes, aes(y=nFeature_RNA, x=sample.y)) +
geom_boxplot()
ggplot(m, aes(y=nFeature_RNA, x=sample.y)) +
geom_boxplot()
nerged
merged
cells_327_cetuxi1 <- cells[grepl('-2', cells)]
cells_327_cetuxi2 <- cells[grepl('-6', cells)]
nerged
transkate
translate
cells_327_nt1 <- cells[grepl('-1', cells)]
cells_327_nt2 <- cells[grepl('-5', cells)]
mah <- GetAssayData(object = srdata[['RNA']], slot='data')
de_replicates_cetuxi_327 <- FindMarkers(mah, cells.1 = cells_327_cetuxi1, cells.2 = cells_327_cetuxi2)
cells_327_cetuxi1
cells <- colnames(srdata)
cells_327_cetuxi1 <- cells[grepl('-2', cells)]
cells_327_cetuxi2 <- cells[grepl('-6', cells)]
de_replicates_cetuxi_327 <- FindMarkers(mah, cells.1 = cells_327_cetuxi1, cells.2 = cells_327_cetuxi2)
de_replicates_nt_327 <- FindMarkers(mah, cells.1 = cells_327_nt1, cells.2 = cells_327_nt2)
cells_327_nt1 <- cells[grepl('-1', cells)]
cells_327_nt2 <- cells[grepl('-5', cells)]
de_replicates_nt_327 <- FindMarkers(mah, cells.1 = cells_327_nt1, cells.2 = cells_327_nt2)
head(de_replicates_cetuxi_327)
nrow(de_replicates_cetuxi_327[de_replicates_cetuxi_327$p_val_adj <0.05,])
nrow(de_replicates_nt_327[de_replicates_cetuxi_327$p_val_adj <0.05,])
de_replicates_cent1_327 <- FindMarkers(mah, cells.1 = cells_327_cetuxi1, cells.2 = cells_327_nt1)
de_replicates_cent2_327 <- FindMarkers(mah, cells.1 = cells_327_cetuxi2, cells.2 = cells_327_nt2)
nrow(de_replicates_cent1_327[de_replicates_cent1_327$p_val_adj <0.05,])
nrow(de_replicates_cent2_327[de_replicates_cent2_327$p_val_adj <0.05,])
cent2 <- nrow(de_replicates_cent2_327[de_replicates_cent2_327$p_val_adj <0.05,])
cent1 <- nrow(de_replicates_cent1_327[de_replicates_cent1_327$p_val_adj <0.05,])
samecells <- nrow(de_replicates_cetuxi_327[de_replicates_cetuxi_327$p_val_adj<0.05,])
diffcells <- nrow(de_replicates_cetuxi_327[de_replicates_cetuxi_327$p_val_adj<0.05,])
diffcells <- nrow(de_replicates_nt_327[de_replicates_nt_327$p_val_adj<0.05,])
df <- data.frame(de_cetuxint1 = cent1, de_cetuxint2=cent2, de_replicates_diff_cetuxi=diffcells, de_replicates_same_nt= nrow(de_replicates_nt_327[de_replicates_nt_327$p_val_adj<0.05,]))
df
df2 <- t(df)
df2
colnames(df2)<-'ndeg'
df2$comparison <- rownames(df2)
df2 <- t(df)
colnames(df2)<-'ndeg'
names(df2)
df2
df2 <- as.data.frame(df2)
df2$comparison <- rownames(df2)
ggplot(genes, aes(y=ndeg, x=comparison)) +geom_bar()
ggplot(df2, aes(y=ndeg, x=comparison)) +geom_bar()
ggplot(df2, aes(y=ndeg, x=comparison)) +geom_points()
ggplot(df2, aes(y=ndeg, x=comparison)) +geom_point()
df2
head(de_replicates_cent1_327)
cent1 <- rownames(de_replicates_cent1_327[de_replicates_cent1_327$p_val_adj <0.05,])
cent2 <- rownames(de_replicates_cent2_327[de_replicates_cent2_327$p_val_adj <0.05,])
length(intersect(cent1, cent2))
cet <- rownames(de_replicates_cetuxi_327[de_replicates_cetuxi_327$p_val_adj <0.05,])
nt <- rownames(de_replicates_nt_327[de_replicates_nt_327$p_val_adj <0.05,])
length(intersect(cet, nt))
intersect(cent1, cent2)
rescetnt <- data.frame(genes=intersect(cet, nt)))
rescetnt <- data.frame(genes=intersect(cet, nt))
resde <- data.frame(genes=intersect(cent1,cent2))
universe <- rownames(de_replicates_cent1_327)
dim(de_replicates_cent1_327)
dim(de_replicates_cent2_327)
dim(de_replicates_cetuxi_327)
universe <- union(rownames(de_replicates_cent1_327),rownames(de_replicates_cent2_327),rownames(de_replicates_cetuxi_327), rownames(de_replicates_nt_327))
universe <- union(c(rownames(de_replicates_cent1_327),rownames(de_replicates_cent2_327),rownames(de_replicates_cetuxi_327), rownames(de_replicates_nt_327)))
universe <- unique(c(rownames(de_replicates_cent1_327),rownames(de_replicates_cent2_327),rownames(de_replicates_cetuxi_327), rownames(de_replicates_nt_327)))
head(universe)
length(universe)
head(rescetnt)
universe <- data.frame(genes=universe))
universe <- data.frame(genes=universe)
head(universe)
verse <- read.table(universef, header=FALSE, sep="\t")
goenrich <- function(ontology, namedgenes, id) {
GOdata <- new("topGOdata",
ontology = ontology,
allGenes = namedgenes,
geneSel = function(x) {x!=0},
nodeSize = 10,
annot = annFUN.org, mapping = "org.Hs.eg.db", ID = id)
allGO = usedGO(object = GOdata)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(allGO), numChar=1000)
allRes$pnom <- as.numeric(allRes$classicFisher)
allRes$classicFisher <- NULL
allRes
}
goenrichall <- function(genes, universe, ontologies, id) {
namedgenes <- rep(0, nrow(universe))
names(namedgenes) <- universe[,1]
namedgenes[names(namedgenes) %in% genes] <- 1
allont <- lapply(ontologies, goenrich, namedgenes, id)
garbage <- lapply(seq(1,length(ontologies)), function(i) {allont[[i]]$ont <<- ontologies[[i]]})
allontdf <- do.call(rbind, allont)
allontdf
}
library("org.Hs.eg.db",  quietly=TRUE)
library(topGO,  quietly=TRUE)
goenrichall(resde, universe, list('BP'),'symbol')
dego <- goenrichall(resde, universe, list('BP'),'symbol')
derepli <- goenrichall(rescetnt, universe, list('BP'),'symbol')
head(dego)
head(derepli)
head(rescetnt)
derepli <- goenrichall(rescetnt$genes, universe, list('BP'),'symbol')
dego <- goenrichall(resde$genes, universe, list('BP'),'symbol')
head(derepli)
derepli$padj <- p.adjust(derepli$pnom, method='FDR')
derepli$padj <- p.adjust(derepli$pnom, method='fdr')
dego$padj <- p.adjust(dego$pnom, method='fdr')
dego[dego$padj < 0.05,]
nrow(dego[dego$padj < 0.05,])
nrow(dereply[dereply$padj < 0.05,])
nrow(derepl[derepl$padj < 0.05,])
nrow(derepli[derepli$padj < 0.05,])
head(derepli[order(dereply$padj),]
)
head(derepli[order(derepli$padj),])
head(dego[order(dego$padj),])
length(intersect(rownames(de_replicates_cent2_327), rownames(de_replicates_cetuxi_327)))
length(intersect(rownames(de_replicates_cent2_327), rownames(de_replicates_nt_327)))
l1 <- intersect(rownames(de_replicates_cent2_327), rownames(de_replicates_nt_327))
l2 <- intersect(rownames(de_replicates_cent2_327), rownames(de_replicates_cetuxi_327))
length(intersect(l1,l2))
head(derepli[order(derepli$padj),],n=30)
head(dego[order(dego$padj),],n=30)
resde$genes
head(de_replicates_cent2_327)
de_replicates_cent1_327[grepl(rownames(de_replicates_cent2_327),'DEFA5'),]
de_replicates_cent1_327[grepl('DEFA5',rownames(de_replicates_cent2_327)),]
de_replicates_cent1_327[grepl('TP53INP2',rownames(de_replicates_cent2_327)),]
de_replicates_cent1_327[grepl('GRIN2C',rownames(de_replicates_cent2_327)),]
de_replicates_cent1_327[grepl('GRIN2C',rownames(de_replicates_cent1_327)),]
de_replicates_cent1_327[grepl('MALAT1',rownames(de_replicates_cent1_327)),]
de_replicates_cent1_327[grepl('TNS4',rownames(de_replicates_cent1_327)),]
