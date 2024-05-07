#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'vande', 'v', 1, 'character',
  'scratch', 's', 1, 'character',
  'res', 'r', 1, 'numeric',
  'mod_rcasc', 'm', 1, 'character',
  'normalized', 'n', 0, 'logical',
  'variableFeat', 'f', 1, 'character',
  'pca', 'p', 1, 'numeric'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$vande) | !is.null(opt$help) | is.null(opt$pca) | is.null(opt$res) | is.null(opt$mod_rcasc)) {
    cat(getopt(opts, usage=TRUE))
    stop('-v, -s, -r, -m and -p are mandatory')
}

if (is.null(opt$normalized)) {
  opt$normalized <- FALSE
}
if (opt$normalized & is.null(opt$variableFeat)) {
  stop('If working with normalized/scaled data (-n) I need also the variableFeatures given by Seurat (-f smt.rds)!')
}

source(opt$mod_rcasc)

SCRATCH <- opt$scratch
SEPARATOR <- ','
setwd(dirname(opt$vande))
# seed has been 157 for all res 0.2, set to 173 for cellcyclecorr
# 143
seuratBootstrap(group="docker", scratch.folder=SCRATCH, file=opt$vande, nPerm=1, permAtTime=1, percent=10, separator=SEPARATOR, pcaDimensions=opt$pca, seed = 173, resolution=opt$res, logTen=0, isNormalized=opt$normalized, variableFeatures=opt$variableFeat)
# output is..the clustering file
