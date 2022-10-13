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
  'prefix', 'x', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$prefix) | is.null(opt$res) | is.null(opt$kind) | is.null(opt$input) | is.null(opt$cycle) | is.null(opt$npc) | is.null(opt$base) | !is.null(opt$help)) {
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
}

output_clucy_f <- paste0(prefix, 'clu_cycle.tsv')
output_PC_f <- paste0(prefix, 'PC.pdf')
output_scaled_f <- paste0(prefix, 'scaled.tsv.gz')

cdir <- getwd()

input_files_l <- unlist(strsplit(input_files, ','))

samples <- basename(input_files_l)
samples <- gsub('_saverAlive.csv.gz','', samples, fixed=TRUE)

set.seed(42)

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
sdata <- RunPCA(sdata, verbose = FALSE)

print('PCA1')
pdf(output_PC_f)
ElbowPlot(sdata)
dev.off()

#sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
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

if (kind == 'G1') {
  cc_df <- data.frame(id=colnames(sdata), phase=sdata$Phase)
  toRemove <- cc_df[cc_df$phase != 'G1', 'id']
  sdata <- sdata[,!colnames(sdata) %in% toRemove]
  sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE) # chose?
  #sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
  sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
  print('PCA2')
} else if (kind == "correct") {
  #sdata <- ScaleData(sdata, verbose = FALSE)
  sdata <- ScaleData(sdata,  vars.to.regress = c("S.Score", "G2M.Score")) #, features = rownames(CRC0327_saver))
  sdata <- RunPCA(sdata, npcs = npc, verbose = FALSE) # chose?
  #sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:npc)
  sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:npc)
}

sdata <- FindClusters(sdata, resolution = resolution)
#pdf(output_umap_f)
#DimPlot(sdata, reduction = "umap")
#graphics.off()


print('clustered')
df <- data.frame(row.names=names(sdata$orig.ident), sample=sdata$sample, cluster=sdata$seurat_clusters, cycle=sdata$Phase)
write.table(df, file=output_clucy_f, sep="\t", quote=F)

DefaultAssay(sdata) <- "RNA"

# https://www.biostars.org/p/9478172/
norm_data <- GetAssayData(object = sdata, slot = "data")
write.table(norm_data, file=gzfile('norm.tsv.gz'), sep="\t", quote=F)
if (any(Assays(sdata) == 'scale_data')) {
  scaled_data <- GetAssayData(object = sdata, slot = "scale_data")
  write.table(scaled_data, file=gzfile(output_scaled_f), sep="\t", quote=F)
}
