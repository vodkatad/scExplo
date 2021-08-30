#!/usr/bin/env Rscript
library(getopt)
library(Seurat)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'input', 'i', 1, 'character',
  'kind', 'k', 1, 'character',
  'res', 'r', 1, 'numeric',
  'npc', 'n', 1, 'numeric',
  'cycle', 'c', 1, 'character',
  'prefix', 'p', 1, 'character',
  'genes', 'g', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$prefix) | is.null(opt$genes) | is.null(opt$res) | is.null(opt$kind) | is.null(opt$input) | is.null(opt$ccgenes) | is.null(opt$npc) | !is.null(opt$help)) {
  cat(getopt(opts, usage=TRUE))
  stop('-p, -g, -r, -k, -c -n and - i are all mandatory')
}


input_files <- opt$input
prefix <- opt$prefix
genes_f <- opt$genes
resolution <- opt$res
kind <- opt$kind
ccgenes_f <- opt$cycle
npc <- opt$npc

KINDS <- c('all','G1')#,'correct')
if (!kind %in% KINDS) { 
  stop('-k can only be all, G1 or correct')
  stop('correct still to be implemented!')
}

#CRC0322_res0.1_G1_clu_cycle.tsv, CRC0322_res0.1_G1_UMAP.png, CRC0322_res0.1_G1_markers.tsv, CRC0322_res0.1_G1_markersall.tsv, CRC0322_res0.1_G1_markers, CRC0322_res0.1_G1_violins
# prefix == CRC0322_res0.1_G1_
output_clucy_f <- paste0(prefix, 'clu_cycle.tsv')
output_umap_f <- paste0(prefix, 'UMAP.pdf')
output_PC_f <- paste0(prefix, 'PC.pdf')
output_markers_f <- paste0(prefix, 'markers.tsv')
output_markersall_f <- paste0(prefix, 'markersall.tsv')
outputd_markers <- paste0(prefix, '_markers')
outputd_violins <- paste0(prefix, '_violins')

input_files_l <- unlist(strsplit(input_files ','))

set.seed(42)

createSeurat <- function(file) {
  df <- read.table(file, sep=',', header=TRUE, row.names=1)
  res <- CreateSeuratObject(df)
  return(res)
}

data_list <- createSeurat(input_files_l)

data_list <- lapply(X = data_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

features <- SelectIntegrationFeatures(object.list = data_list) # nfeatures 5000?
anchors <- FindIntegrationAnchors(object.list = data_list, anchor.features = features)
sdata <- IntegrateData(anchorset = anchors) 
DefaultAssay(sdata) <- "integrated"
sdata <- ScaleData(sdata, verbose = FALSE)

pdf(output_PC_f)
ElbowPlot(sdata)
dev.off()

# Is this step needed if we repeat it after scaling on cell cycle?
sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE)
sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
cells <- colnames(sdata)
samples <- substring(cells, 20, 21) # the number assigned after cell barcode to each sample
# do we want to do something similar? or map to input names?
# FIXME TODO
#samples <- ifelse(samples == 1, 'NT_1', ifelse(samples == 2, "NT_2", ifelse(samples == 3, "CTX_1", "CTX_2")))
#CRC0327_saver$sample <- samples

cycle <- read.table(ccgenes_f)
s.genes <- cycle[seq(1,44),]
#+1 -2
g2m.genes <- cycle[seq(45,97),]

sdata <- CellCycleScoring(object = sdata, s.features = s.genes, g2m.features = g2m.genes, set.ident=TRUE)

# if kind == "correct"
#sdata <- ScaleData(sdata, verbose = FALSE)
#sdata <- ScaleData(sdata,  vars.to.regress = c("S.Score", "G2M.Score")) #, features = rownames(CRC0327_saver))
#sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE) # chose?
#sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
#sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)

#write.table(cc_df, file="cycle.tsv", sep="\t", quote=F)
#table(CRC0327_saver$Phase, CRC0327_saver$sample)
#table(CRC0327_saver$Phase)
if (kind == 'G1') {
  cc_df <- data.frame(id=colnames(sdata), phase=sdata$Phase)
  toRemove <- cc_df[cc_df$phase != 'G1', 'id']
  sdata <- sdata[,!colnames(sdata) %in% toRemove]
  sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE) # chose?
  sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
  sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
} #else if (kind == "correct") {
#sdata <- ScaleData(sdata, verbose = FALSE)
#sdata <- ScaleData(sdata,  vars.to.regress = c("S.Score", "G2M.Score")) #, features = rownames(CRC0327_saver))
#sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE) # chose?
#sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
#sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
#}

sdata <- FindClusters(sdata, resolution = resolution)
pdf(output_umap_f)
DimPlot(sdata, reduction = "umap")
graphics.off()

	
df <- data.frame(row.names=names(sdata$orig.ident), sample=sdata$sample, cluster=sdata$seurat_clusters, cycle=sdata$Phase)
write.table(df, file=output_clucy_f, sep="\t", quote=F)

n_cl <- length(unique(df$cluster))
n_markers_vs_all <- data.frame(cl=seq(0, (n_cl-1)), n=rep(NA, n_cl))
DefaultAssay(sdata) <- "RNA"

setwd(outputd_markers)
for (i in seq(0,(n_cl-1))) {
	n <- 0
	tryCatch(
	{
		cl <- FindConservedMarkers(sdata, ident.1 = i, verbose = FALSE, grouping.var="sample")
		write.table(cl, file=paste0("res_", resolution, "_cl_",i,".tsv"), sep="\t", quote=F)
		n <- nrow(cl[cl$max_pval < 0.05,])
	}, error= function(cond) {} # we swallow and keep 0
	)
	n_markers_vs_all[i+1, 'n'] <- n                
}
setwd(..)
write.table(n_markers_vs_all, file=output_markers_f, sep="\t", quote=F)
setwd(outputd_markers)
	
n_markers <- matrix(rep(NA, n_cl*n_cl), nrow=n_cl, ncol=n_cl)
for (i in seq(0,(n_cl-1))) {
	j <- i+1
	while (j <= (n_cl-1)) {
	#for (j in seq(i+1, (n_cl-1))) {
		n <- 0
		tryCatch(
		{
			cl <- FindConservedMarkers(sdata, ident.1 = i, ident.2 = j,verbose = FALSE, grouping.var="sample")
			write.table(cl, file=paste0("res_", resolution, "_cl_", i,"vs",j,".tsv"), sep="\t", quote=F)
			n <- nrow(cl[cl$max_pval < 0.05,])
		}, error = function(cond) {}
    )
		n_markers[i+1, j+1] <- n                
		j <- j + 1
	}
}
write.table(n_markers, file=output_markersall_f, sep="\t", quote=F)
# plot of numbers? I'd say separate rule here

genes <- read.table(genes_f, sep="\t", header=F)
colnames(genes) <- c('symbol', 'ensg')

genes_expl <- do.call(rbind, lapply(genes$symbol, function(x) { y <- genes[genes$symbol==x,, drop=F]; 
                                                                ensgs <- strsplit(as.character( y$ensg),','); 
                                                                res=unlist(unique(ensgs)); 
                                                                data.frame(sym=rep(x,length(res)), ens=res)}))

setwd(..)
setwd(outputd_violins)
for (i in seq(1, nrow(genes))) {
 	tryCatch({
	  pdf(paste0(genes[i,'sym'], "_", genes[i,'ens'], '_res', resolution ,'.pdf')); 
	  print(VlnPlot(sdata, genes[i,'ens']));
  	dev.off()
  },
	error=function(cond) {}
  ) # do we need to swallow here too? TODO FIXME
}