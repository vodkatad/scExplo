library(rtracklayer)
gtf_f <- snakemake@input[['gtf']]
map_f <- snakemake@output[['map']]
map_f_disamb <- snakemake@output[['map_disamb']]

x <- data.frame(rtracklayer::import(gtf_f))
dd <- unique(data.frame(gs=x$gene_name, ens=x$gene_id))
save.image('pippo.Rdata')

tto_disamb <- as.data.frame(table(dd$gs))
to_disamb <- as.character(tto_disamb[tto_disamb$Freq > 1, 1] )

ddd <- dd[!dd$gs %in% to_disamb,]
dd_disamb <- dd[dd$gs %in% to_disamb,]

choose_ensg <- function(gene, data) {
	myd <- data[data$gs == gene, ]
	myd[1,, drop=FALSE]
}
dd_disamb <- dd_disamb[order(dd_disamb$gs, dd_disamb$ens),]
dd_disambiguato <- lapply(to_disamb, choose_ensg, dd_disamb)
dd_disambiguato <- do.call(rbind, dd_disambiguato)

dd_disambiguato <- rbind(dd_disambiguato, ddd)
write.table(dd, file=map_f, sep="\t", quote=FALSE, row.names=FALSE)
write.table(dd_disambiguato, file=map_f_disamb, sep="\t", quote=FALSE, row.names=FALSE)

