#!/usr/bin/env Rscript

library("getopt")
library("GSEABase")

opts <- matrix(c(
        'help', 'h', 0, 'logical',
        'input' , 'i', 1, 'character',
        'rds'  , 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$input) || is.null(opt$rds)) {
    usage <- getopt(opts, usage=TRUE)
    stop(usage)
}

l <- unlist(strsplit(opt$input,','))
os <- lapply(l, function(x) readRDS(file=x))
genesetcoll <- GeneSetCollection(os)
saveRDS(genesetcoll, file = opt$rds)

