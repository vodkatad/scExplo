#!/usr/bin/env Rscript

library("singscore")
library('getopt')
opts <- matrix(c(
        'help', 'h', 0, 'logical',
        'expression' , 'e', 1, 'character',
        'signature'  , 's', 1, 'character',
        'outfile','o',1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$signature) || is.null(opt$expression) || is.null(opt$outfile)) {
    usage <- getopt(opts, usage=TRUE)
    stop(usage)
}

expr_file <- opt$expression
rds_sign <- opt$signature

geneset <- readRDS(rds_sign)
expr_data <- read.table(gzfile(expr_file), sep="\t", header=TRUE, row.names=1)
rankData <- rankGenes(expr_data)
# TODO warning about missing genes?
scorelist <- lapply(geneset, function(x) { simpleScore(rankData, upSet = x) } )
n <- ncol(expr_data)
scoredf <- as.data.frame(sapply(scorelist, function(x) { if (length(x$TotalScore)!=0) {x$TotalScore} else {rep(NA, n)} }))
rownames(scoredf) <- colnames(expr_data)
colnames(scoredf) <- names(geneset)

res <- t(scoredf)
write.table(res, file=opt$outfile, quote=FALSE, sep="\t")
# do we want also pvalues?
