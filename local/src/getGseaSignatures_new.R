#!/usr/bin/env Rscript

library("GSEABase")


gmt_file <- snakemake@input[["gtm"]]
outputfile <- snakemake@output[["rds"]]


genesetcollection <- getGmt(gmt_file, geneIdType=SymbolIdentifier()) # is cat needed? FIXME
saveRDS(genesetcollection, file = outputfile)

