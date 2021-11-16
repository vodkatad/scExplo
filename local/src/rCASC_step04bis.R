#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'vande', 'v', 1, 'character',
  'scratch', 's', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$vande) | is.null(opt$scratch) | !is.null(opt$help)) {
    cat(getopt(opts, usage=TRUE))
    stop('-v and -s are mandatory')
}

SCRATCH <- opt$scratch
SEPARATOR <- ','
#save.image('pluto.Rdata')
setwd(dirname(opt$vande))
seurat_ccycle(group="docker", scratch.folder=SCRATCH, file=opt$vande, separator=SEPARATOR, seed = 157)
