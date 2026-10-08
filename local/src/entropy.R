in_f2 <- snakemake@input[['expr']]
d2 <- read.csv(in_f2)
rownames(d2) <- d2$X
d2$X <- NULL
entropy <- apply(d2, 2, function(x) { - sum(x *log(x, base=2))} )
write.table(as.data.frame(entropy), file=snakemake@output[['en']], sep="\t", quote=F)