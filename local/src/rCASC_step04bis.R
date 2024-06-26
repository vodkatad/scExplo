#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'vande', 'v', 1, 'character',
  'scratch', 's', 1, 'character', 
  'source', 'c', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$vande) | is.null(opt$scratch) | !is.null(opt$help) | is.null(opt$source)) {
    cat(getopt(opts, usage=TRUE))
    stop('-v and -s are mandatory')
}

source(opt$source)

SCRATCH <- opt$scratch
SEPARATOR <- ','
SEED <- 157
#save.image('pluto.Rdata')
setwd(dirname(opt$vande))
seurat_ccycle(group="docker", scratch.folder=SCRATCH, file=opt$vande, separator=SEPARATOR, seed = SEED)
#seuratPCAEval(group = "docker", scratch.folder=SCRATCH, file=opt$vande, separator=",", logTen = 0, seed = SEED)
# this generates 
#{sample}/Results/filtered_annotated_saver_ribomito_{sample}/PCE_bowPlot.pdf

