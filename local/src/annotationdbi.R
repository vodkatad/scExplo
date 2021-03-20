#!/usr/bin/env Rscript
# Translate a column using AnnotationDbi and the given library (defualt or.Hs.eg.db)
# can work with headerless or headered files (-d), if -a is given the translation is added as a column after the translated one (-n) otherwise
# the translated column takes its place.


library("getopt")

opts <- matrix(c(
        'help', 'h', 0, 'logical',
        'inputfile' , 'i', 1, 'character',
        'outputfile'  , 'o', 1, 'character',
        'to' , 't', 1, 'character',
        'from'  , 'f', 1, 'character',
        'colnumber'  , 'n', 1, 'integer',
        'library', 'l', 1, 'character',
        'hasheader', 'd', 0, 'logical',
        'newname', 'a', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

# TODO: probably this if can be skipped if one has the time to fully read getopt manual and set arguments as mandatory
if (is.null(opt$inputfile) || is.null(opt$outputfile) || is.null(opt$colnumber) || !is.null(opt$help) || is.null(opt$to) || is.null(opt$from)) {
    usage <- getopt(opts, usage=TRUE)
    stop(usage)
}
library(AnnotationDbi)

if (is.null(opt$hasheader)) {
 opt$hasheader <- FALSE
}
    
if (is.null(opt$library))  {
    opt$library <- "org.Hs.eg.db"
}
library(opt$library, character.only=TRUE)

if (!opt$hasheader) {
    IDtsv <- read.table(gzfile(opt$inputfile), stringsAsFactors = FALSE, sep = "\t", quote = "", na.strings=c("", " ", "NA"))
} else {
    IDtsv <- read.table(gzfile(opt$inputfile), stringsAsFactors = FALSE, sep = "\t", quote = "", na.strings=c("", " ", "NA"), header=TRUE)
}
keys <- as.character(IDtsv[,opt$colnumber])
merge <- function(x) { paste(x, collapse=',')}
geneIDSymbols <- mapIds(get(opt$library), keys=keys, column=opt$to, keytype=opt$from, multiVals=merge)

if (!is.null(opt$newname)) {
    IDtsv[,opt$newname] <- as.character(geneIDSymbols)
    # 1 ... ncol(IDtsv)
    # 1 ... ncol(IDtsv)+1
    ncols <- ncol(IDtsv)
    neworder <- unique(c(seq(1, opt$colnumber), ncols, seq((opt$colnumber+1), (ncols-1))))
    IDtsv <- IDtsv[, neworder]
} else {
    IDtsv[,opt$colnumber] <- as.character(geneIDSymbols)
}

if (!opt$hasheader) {
    write.table(IDtsv, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file = gzfile(opt$outputfile))
} else {
    write.table(IDtsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, file = gzfile(opt$outputfile))
}
