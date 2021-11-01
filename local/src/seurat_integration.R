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
  'base', 'b', 1, 'character',
  'prefix', 'x', 1, 'character',
  'genes', 'g', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

save.image('pippo.Rdata')
if (is.null(opt$prefix) | is.null(opt$genes) | is.null(opt$res) | is.null(opt$kind) | is.null(opt$input) | is.null(opt$cycle) | is.null(opt$npc) | is.null(opt$base) | !is.null(opt$help)) {
  cat(getopt(opts, usage=TRUE))
  stop('-b -x, -g, -r, -k, -c -n and - i are all mandatory')
}


input_files <- opt$input
prefix <- opt$prefix
genes_f <- opt$genes
resolution <- opt$res
kind <- opt$kind
ccgenes_f <- opt$cycle
npc <- opt$npc
based <- opt$base

KINDS <- c('all','G1','correct')
if (!kind %in% KINDS) { 
  stop('-k can only be all, G1 or correct')
  stop('correct still to be implemented!')
}


#CRC0322_res0.1_G1_clu_cycle.tsv, CRC0322_res0.1_G1_UMAP.png, CRC0322_res0.1_G1_markers.tsv, CRC0322_res0.1_G1_markersall.tsv, CRC0322_res0.1_G1_markers, CRC0322_res0.1_G1_violins
# prefix == CRC0322_res0.1_G1_

output_clucy_f <- paste0(prefix, 'clu_cycle.tsv')
output_cy_f <- paste0(prefix, 'cycle.tsv')
output_umap_f <- paste0(prefix, 'UMAP.pdf')
output_PC_f <- paste0(prefix, 'PC.pdf')
output_markers_f <- paste0(prefix, 'markers.tsv')
output_markersall_f <- paste0(prefix, 'markersall.tsv')
outputd_markers <- paste0(based, prefix, 'markers')
outputd_violins <- paste0(based, prefix, 'violins')
cdir <- getwd()

print(outputd_violins)
print(outputd_markers)
print(cdir)
print(genes_f)

input_files_l <- unlist(strsplit(input_files, ','))

samples <- basename(input_files_l)
samples <- gsub('_saverAlive.csv.gz','', samples, fixed=TRUE)

set.seed(42)

save.image('pluto.Rdata')

createSeurat <- function(file) {
  df <- read.table(gzfile(file), sep=',', header=TRUE, row.names=1)
  res <- CreateSeuratObject(df)
  return(res)
}

data_list <- lapply(input_files_l, createSeurat)

data_list <- lapply(X = data_list, FUN = function(x) {
    ## mito removal
    #x[["percent.mito"]] <- PercentageFeatureSet(x, pattern = "^MT-") # need to work on ENSG :(
    #x <- subset(x, subset = percent.mito < 25)
    #Raf:
    #Per ribosomali tutto quello sotto il 50% delle conte totali è buono sopra sono i doppietti
    #Per rna umano sano in genere la soglia del mitocondriali è al disotto del 10%
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

features <- SelectIntegrationFeatures(object.list = data_list) # nfeatures 5000?
anchors <- FindIntegrationAnchors(object.list = data_list, anchor.features = features)
sdata <- IntegrateData(anchorset = anchors) 
DefaultAssay(sdata) <- "integrated"

print('integrated')
sdata <- ScaleData(sdata, verbose = FALSE)

print('scaled')

# Is this step needed if we repeat it after scaling on cell cycle?
sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE)

print('PCA1')
pdf(output_PC_f)
ElbowPlot(sdata)
dev.off()

sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
cells <- colnames(sdata)

# add info on samples name to the object
samples_n <- sapply(strsplit(cells, "_"), function(x){x[[2]]})
sdata$sample <- samples[as.numeric(samples_n)]

cycle <- read.table(ccgenes_f)
s.genes <- cycle[seq(1,44),]
#+1 -2
g2m.genes <- cycle[seq(45,97),]

sdata <- CellCycleScoring(object = sdata, s.features = s.genes, g2m.features = g2m.genes, set.ident=TRUE)

print('cc')

df <- data.frame(row.names=names(sdata$orig.ident), sample=sdata$sample, cycle=sdata$Phase)
write.table(df, file=output_cy_f, sep="\t", quote=F)

if (kind == 'G1') {
  cc_df <- data.frame(id=colnames(sdata), phase=sdata$Phase)
  toRemove <- cc_df[cc_df$phase != 'G1', 'id']
  sdata <- sdata[,!colnames(sdata) %in% toRemove]
  sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE) # chose?
  sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
  sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
  print('PCA2')
} else if (kind == "correct") {
  #sdata <- ScaleData(sdata, verbose = FALSE)
  sdata <- ScaleData(sdata,  vars.to.regress = c("S.Score", "G2M.Score")) #, features = rownames(CRC0327_saver))
  sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE) # chose?
  sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
  sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
}

sdata <- FindClusters(sdata, resolution = resolution)
pdf(output_umap_f)
DimPlot(sdata, reduction = "umap")
graphics.off()


print('clustered')
df <- data.frame(row.names=names(sdata$orig.ident), sample=sdata$sample, cluster=sdata$seurat_clusters, cycle=sdata$Phase)
write.table(df, file=output_clucy_f, sep="\t", quote=F)

n_cl <- length(unique(df$cluster))
n_markers_vs_all <- data.frame(cl=seq(0, (n_cl-1)), n=rep(NA, n_cl))
DefaultAssay(sdata) <- "RNA"

print(outputd_markers)
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
setwd(cdir)
write.table(n_markers_vs_all, file=output_markers_f, sep="\t", quote=F)
setwd(outputd_markers)
	
print('markers1')
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
setwd(cdir)
write.table(n_markers, file=output_markersall_f, sep="\t", quote=F)
# plot of numbers? I'd say separate rule here
print('markers2')

setwd(cdir)
genes <- read.table(genes_f, sep="\t", header=F)
colnames(genes) <- c('symbol', 'ensg')

print('genes')
genes_expl <- do.call(rbind, lapply(genes$symbol, function(x) {y <- genes[genes$symbol==x,, drop=F]; ensgs <- strsplit(as.character( y$ensg),','); res=unlist(unique(ensgs)); data.frame(sym=rep(x,length(res)), ens=res)}))

setwd(outputd_violins)
for (i in seq(1, nrow(genes_expl))) {
 	tryCatch({
	  pdf(paste0(genes_expl[i,'sym'], "_", genes_expl[i,'ens'], '_res', resolution ,'.pdf'))
	  print(VlnPlot(sdata, genes_expl[i,'ens']))
  	  dev.off()
  },
	error=function(cond) {}
  ) # do we need to swallow here too? TODO FIXME
}

# https://www.biostars.org/p/9478172/
norm_data <- GetAssayData(object = sdata, slot = "data")
write.table(norm_data, file=gzfile('norm.tsv.gz'), sep="\t", quote=F)
if (any(Assays(sdata) == 'scale_data')) {
  scaled_data <- GetAssayData(object = sdata, slot = "scale_data")
  write.table(scaled_data, file=gzfile('scaled.tsv.gz'), sep="\t", quote=F)
}