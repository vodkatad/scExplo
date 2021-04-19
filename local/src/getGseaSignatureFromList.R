#!/usr/bin/env Rscript

library("getopt")
library("GSEABase")

opts <- matrix(c(
        'help', 'h', 0, 'logical',
        'genes' , 'g', 1, 'character',
        'name' , 'n', 1, 'character',
        'rds'  , 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$genes) || is.null(opt$rds) || is.null(opt$name)) {
    usage <- getopt(opts, usage=TRUE)
    stop(usage)
}
genes_file <- opt$genes
outputfile <- opt$rds

genes <- read.table(genes_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)

geneset <- GeneSet(genes[,1], geneIdType=SymbolIdentifier()) # does not work? Y?
setName(geneset) <- opt$name
#manual boh
# > details(GeneSet(geneIds(c('TP53')), geneIdType=SymbolIdentifier()))
# Error in (function (classes, fdef, mtable)  : 
#   unable to find an inherited method for function ‘geneIds’ for signature ‘"character"’
# > details(GeneSet('TP53', geneIdType=SymbolIdentifier()))
# setName: NA 
# geneIds: TP53 (total: 1)
# geneIdType: Symbol
# collectionType: Null 
# setIdentifier: godot:44786:Mon Apr 19 17:43:47 2021:5
# description: 
# organism: 
# pubMedIds: 
# urls: 
# contributor: 
# setVersion: 0.0.1
# creationDate: 
# > 
saveRDS(geneset, file = outputfile)

