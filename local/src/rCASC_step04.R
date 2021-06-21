#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'vande', 'v', 1, 'character',
  'scratch', 's', 1, 'character',
  'pca', 'p', 1, 'numeric'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$vande) | !is.null(opt$help) | is.null(opt$pca)) {
    cat(getopt(opts, usage=TRUE))
    stop('-v, -p -s and -o are mandatory')
}

SCRATCH <- opt$scratch
SEPARATOR <- ','
#save.image('pluto.Rdata')
setwd(dirname(opt$vande))
seuratBootstrap(group="docker", scratch.folder=SCRATCH, file=opt$vande, nPerm=90, permAtTime=10, percent=10, separator=SEPARATOR, pcaDimensions=opt$pca, seed = 143)
# output is..the clustering file
