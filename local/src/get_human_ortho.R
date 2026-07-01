library(babelgene)

input_f <- snakemake@input[['genes']]
output_f <- snakemake@output[['outf']]
log_f <- snakemake@log[['log']]

d <- read.table(input_f, sep="\t", header=TRUE)
hs <- orthologs(genes = d$mm, human=FALSE, species='mouse')

sink(log_f)
dim(d)
dim(hs)
sink()

write.table(hs[, 'human_symbol', drop=F], file=output_f, sep="\t", quote=FALSE, row.names=FALSE)
