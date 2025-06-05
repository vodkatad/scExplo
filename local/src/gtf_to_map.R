library(rtracklayer)
gtf_f <- snakemake@input[['gtf']]
map_f <- snakemake@output[['map']]

x <- data.frame(rtracklayer::import(gtf_f))
dd <- unique(data.frame(gs=x$gene_name, ens=x$gene_id))
write.table(dd, file=map_f, sep="\t", quote=FALSE, row.names=FALSE)

