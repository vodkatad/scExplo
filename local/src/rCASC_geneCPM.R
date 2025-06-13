#!/usr/bin/env Rscript
library(getopt)
library(rCASC)
opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'expr', 'e', 1, 'character'
  ), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$expr)) {
    cat(getopt(opts, usage=TRUE))
    stop('-e is mandatory')
}

SEPARATOR <- ','
setwd(dirname(opt$expr))

expr <- read.table(opt$expr, sep=SEPARATOR, header=TRUE, row.names=1)
row_sum <- apply(expr, 1, sum)
expr <- expr[row_sum !=0, ]
col_sum <- apply(expr, 2, sum)
tmp1 <- t(expr)/col_sum
tmp1 <- t(tmp1)
tmp1 <- tmp1 * 1000000
cpm_expr <- log2(tmp1+1)

fn_cpm <- substr(opt$expr, 1, nchar(opt$expr)-4)
fn_cpm <- paste0(fn_cpm, '_log2_pc1_cpm.csv')
write.csv(cpm_expr, file=fn_cpm, row.names=TRUE)
