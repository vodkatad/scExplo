#!/usr/bin/env Rscript
library(getopt)
library(rCASC)
opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'clustering', 'c', 1, 'character',
  'expr', 'e', 1, 'character',
  'genes', 'g', 1, 'character',
  'pdf', 'p', 1, 'character',
  'source', 't', 1, 'character',
  'scratch', 's', 1, 'character'
  ), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$clustering) | !is.null(opt$help) | is.null(opt$pdf) | is.null(opt$genes) | is.null(opt$source) | is.null(opt$expr)) {
    cat(getopt(opts, usage=TRUE))
    stop('-v, -p -s and -o are mandatory')
}

source(opt$source)

SCRATCH <- opt$scratch
SEPARATOR <- ','
setwd(dirname(opt$clustering))

expr <- read.table(opt$expr, sep=SEPARATOR, header=TRUE, row.names=1)
col_sum <- apply(expr, 2, sum)
tmp1 <- t(expr)/col_sum
tmp1 <- t(tmp1)
tmp1 <- tmp1 * 1000000
cpm_expr <- log2(tmp1+1)

fn_cpm <- substr(opt$expr, 1, nchar(opt$expr)-4)
fn_cpm <- paste0(fn_cpm, '_log2_pc1_cpm.csv')
write.csv(cpm_expr, file=fn_cpm, row.names=TRUE)
geneVisualization(group="docker", scratch.folder=SCRATCH, file=fn_cpm,  clustering.output=opt$clustering, geneList=opt$genes, separator=SEPARATOR, finalName=opt$pdf)
