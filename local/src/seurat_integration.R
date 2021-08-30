#!/usr/bin/env Rscript
library(getopt)
library(Seurat)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'input', 'i', 1, 'character',
  'kind', 'k', 1, 'character',
  'res', 'r', 1, 'numeric',
  'prefix', 'p', 1, 'character',
  'genes', 'g', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$prefix) | is.null(opt$genes) | is.null(opt$res) | is.null(opt$kind) | is.null(opt$input) | !is.null(opt$help)) {
  cat(getopt(opts, usage=TRUE))
  stop('-p, -g, -r, -k and - i are all mandatory')
}


input_files <- opt$input
prefix <- opt$prefix
genes_f <- opt$genes
res <- opt$res
kind <- opt$kind

KINDS <- c('all','G0')#,'correct')
if (!kind %in% KINDS) { 
  stop('-k can only be all, G0 or correct')
  stop('correct still to be implemented!')
}

#CRC0322_res0.1_G1_clu_cycle.tsv, CRC0322_res0.1_G1_UMAP.png, CRC0322_res0.1_G1_markers.tsv, CRC0322_res0.1_G1_markersall.tsv, CRC0322_res0.1_G1_markers, CRC0322_res0.1_G1_violins
# prefix == CRC0322_res0.1_G1_
output_clucy_f <- paste0(prefix, 'clu_cycle.tsv')
output_umap <- paste0(prefix, 'UMAP.png')
output_markers <- paste0(prefix, 'markers.tsv')
output_markersall <- paste0(prefix, 'markersall.tsv')
outputd_markers <- paste0(prefix, '_markers')
outputd_violins <- paste0(prefix, '_violins')

input_files_l <- unlist(strsplit(input_files ','))

set.seed(42)
nt1df <- read.table('saver_CRC0327_NT_1.csv',sep=',',header=TRUE,row.names=1)
nt2df <- read.table('saver_CRC0327_NT_2.csv',sep=',',header=TRUE,row.names=1)
#ctx1df <- read.table('CRC0327_cetux_1_4000_saver_correctdowns.csv',sep=',',header=TRUE,row.names=1)
ctx1df <- read.table('CRC0327_cetux_1_4000_saver_correctdowns.csv',sep=',',header=TRUE,row.names=1)
ctx2df <- read.table('saver_CRC0327_cetux_2.csv',sep=',',header=TRUE,row.names=1)
nt1 <- CreateSeuratObject(nt1df)
nt2 <- CreateSeuratObject(nt2df)
ctx1 <- CreateSeuratObject(ctx1df)
ctx2 <- CreateSeuratObject(ctx2df)
data_list <- list(nt1=nt1, nt2=nt2, ctx1=ctx1, ctx2=ctx2)
data_list <- lapply(X = data_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
features <- SelectIntegrationFeatures(object.list = data_list) # nfeatures 5000?
anchors <- FindIntegrationAnchors(object.list = data_list, anchor.features = features)
CRC0327_saver <- IntegrateData(anchorset = anchors) 
DefaultAssay(CRC0327_saver) <- "integrated"
CRC0327_saver <- ScaleData(CRC0327_saver, verbose = FALSE)
npc <- 30 # chose wisely?
CRC0327_saver <- RunPCA(CRC0327_saver, npcs = npc, verbose = FALSE) # chose?
CRC0327_saver <- RunUMAP(CRC0327_saver, reduction = "pca", dims = 1:npc)
CRC0327_saver <- FindNeighbors(CRC0327_saver, reduction = "pca", dims = 1:npc)
cells <- colnames(CRC0327_saver)
samples <- substring(cells, 20, 21)
samples <- ifelse(samples == 1, 'NT_1', ifelse(samples == 2, "NT_2", ifelse(samples == 3, "CTX_1", "CTX_2")))
CRC0327_saver$sample <- samples

cycle <- read.table('regev_lab_cell_cycle_genes_onlyensg.txt')
s.genes <- cycle[seq(1,44),]
#+1 -2
g2m.genes <- cycle[seq(45,97),]



CRC0327_saver <- CellCycleScoring(object = CRC0327_saver, s.features = s.genes, g2m.features = g2m.genes, set.ident=TRUE)

CRC0327_saver <- ScaleData(CRC0327_saver, verbose = FALSE)
####CRC0327_saver <- ScaleData(CRC0327_saver,  vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"), features = rownames(CRC0327_saver))
npc <- 30 # chose wisely?
CRC0327_saver <- RunPCA(CRC0327_saver, npcs = npc, verbose = FALSE) # chose?
CRC0327_saver <- RunUMAP(CRC0327_saver, reduction = "pca", dims = 1:npc)
CRC0327_saver <- FindNeighbors(CRC0327_saver, reduction = "pca", dims = 1:npc)

pdf <- data.frame(id=colnames(CRC0327_saver), phase= CRC0327_saver$Phase)
write.table(pdf, file="cycle.tsv", sep="\t", quote=F)
table(CRC0327_saver$Phase, CRC0327_saver$sample)
table(CRC0327_saver$Phase)

toRemove <- pdf[pdf$phase != 'G1', 'id']
CRC0327_saver <- CRC0327_saver[,!colnames(CRC0327_saver) %in% toRemove]


CRC0327_saver <- RunPCA(CRC0327_saver, npcs = npc, verbose = FALSE) # chose?
CRC0327_saver <- RunUMAP(CRC0327_saver, reduction = "pca", dims = 1:npc)
CRC0327_saver <- FindNeighbors(CRC0327_saver, reduction = "pca", dims = 1:npc)

### wwwaffancul### wwwaffancul### wwwaffanculo
resolution <- 0.5
CRC0327_saver <- FindClusters(CRC0327_saver, resolution = resolution)
pdf(paste0('CRC0327_saver', resolution, '.pdf'))
DimPlot(CRC0327_saver, reduction = "umap")
graphics.off()

resolution <- 0.2
CRC0327_saver <- FindClusters(CRC0327_saver, resolution = resolution)
pdf(paste0('CRC0327_saver', resolution, '.pdf'))
DimPlot(CRC0327_saver, reduction = "umap")
graphics.off()

resolution <- 0.1
CRC0327_saver <- FindClusters(CRC0327_saver, resolution = resolution)
pdf(paste0('CRC0327_saver', resolution, '.pdf'))
DimPlot(CRC0327_saver, reduction = "umap")
graphics.off()

resolution <- 0.05
CRC0327_saver <- FindClusters(CRC0327_saver, resolution = resolution)
pdf(paste0('CRC0327_saver', resolution, '.pdf'))
DimPlot(CRC0327_saver, reduction = "umap")
graphics.off()

resolution <- 0.01
CRC0327_saver <- FindClusters(CRC0327_saver, resolution = resolution)
pdf(paste0('CRC0327_saver', resolution, '.pdf'))
DimPlot(CRC0327_saver, reduction = "umap")
graphics.off()

# bad habits and global variables my friend with R live sessions
clusterres <- function(resolution) {
	CRC0327_saver <- FindClusters(CRC0327_saver, resolution = resolution)
	#pdf(paste0('CRC0327_saver', resolution, '.pdf'))
	#DimPlot(CRC0327_saver, reduction = "umap")
	#graphics.off()
	
	df <- data.frame(row.names=names(CRC0327_saver$orig.ident), sample=CRC0327_saver$sample, cluster=CRC0327_saver$seurat_clusters)
	write.table(df, file=paste0("CRC0327_saver_seuint", resolution, ".tsv"), sep="\t", quote=F)
	n_cl <- length(unique(df$cluster))
	n_markers_vs_all <- data.frame(cl=seq(0, (n_cl-1)), n=rep(NA, n_cl))
	DefaultAssay(CRC0327_saver) <- "RNA" # stupit stupid!!
	for (i in seq(0,(n_cl-1))) {
		n <- 0
		tryCatch(
		{
			cl <- FindConservedMarkers(CRC0327_saver, ident.1 = i, verbose = FALSE, grouping.var="sample")
			write.table(cl, file=paste0("res_", resolution, "_cl_",i,".tsv"), sep="\t", quote=F)
			n <- nrow(cl[cl$max_pval < 0.05,])
		}, error= 
		{
			function(cond) {
				return(0)
			}
		})
		n_markers_vs_all[i+1, 'n'] <- n                
	}
	write.table(n_markers_vs_all, file=paste0("markersall_",resolution, ".tsv"), sep="\t", quote=F)
	
	n_markers <- matrix(rep(NA, n_cl*n_cl), nrow=n_cl, ncol=n_cl)
	for (i in seq(0,(n_cl-1))) {
		j <- i+1
		while (j <= (n_cl-1)) {
		#for (j in seq(i+1, (n_cl-1))) {
			n <- 0
			tryCatch(
			{
				cl <- FindConservedMarkers(CRC0327_saver, ident.1 = i, ident.2 = j,verbose = FALSE, grouping.var="sample")
				write.table(cl, file=paste0("res_", resolution, "_cl_", i,"vs",j,".tsv"), sep="\t", quote=F)
				n <- nrow(cl[cl$max_pval < 0.05,])
			}, error =
			{	
				function(cond) {
					return(0)
				}
			})
			n_markers[i+1, j+1] <- n                
			j <- j + 1
		}
	}
	write.table(n_markers, file=paste0("markers_",resolution, ".tsv"), sep="\t", quote=F)
}

clusterres(0.5)
clusterres(0.2)
clusterres(0.1)
clusterres(0.05)
clusterres(0.01)

