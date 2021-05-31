#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'counts', 'c', 1, 'character',
  'gtf', 'g', 1, 'character',
  'output', 'o', 1, 'character',
  'scratch', 's', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$counts) | !is.null(opt$help) | is.null(opt$gtf) | is.null(opt$output) | is.null(opt$scratch)) {
    cat(getopt(opts, usage=TRUE))
    stop('-c, -g, -s and -o are mandatory')
}

INPUT <- basename(opt$counts)

#WD="/path_to_raw_and_saver_data/" # need to be in the same directory?
#SCRATCH <- "/home/rcalogero/scratch/"  
SCRATCH <- opt$scratch
SEPARATOR <- ","
GTF <- opt$gtf
ls
#setwd(WD)
mitoRiboUmi(group = "docker", scratch.folder=SCRATCH, file=INPUT, separator=SEPARATOR, gtf.name=GTF, bio.type="protein_coding", umiXgene=3) # TODO set 3 as param
# cannot be run in parallel if they write files with the same names.
file.rename(from="Ribo_mito.pdf", to=opt$output)
