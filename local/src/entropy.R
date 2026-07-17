in_f2 <- snakemake@input[['expr']]
d2 <- read.csv(in_f2)
rownames(d2) <- d2$X
d2$X <- NULL

# added after the run on public datasets which had no all 0 rows (or they would have NaN entropy) - no zeros in log2cpm after imputation
d3 <- apply(d2, 1, function(x) {all(x!=0)})
d4 <- d2[d3,]

entropy <- apply(d4, 2, function(x) { - sum(x *log(x, base=2))} )
write.table(as.data.frame(entropy), file=snakemake@output[['en']], sep="\t", quote=F)
