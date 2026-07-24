library(CytoTRACE)
in_f <- snakemake@input[[1]]

mat <- read.csv(in_f)
rownames(mat) <- mat$X
mat$X <- NULL
mat <- as.matrix(mat)
mode(mat) <- 'integer'
ctt <- CytoTRACE(mat, ncores=3)
res <- ctt$CytoTRACE
write.table(res, file=snakemake@output[['sign']], sep="\t", quote=FALSE)
