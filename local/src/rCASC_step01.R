#!/usr/bin/env Rscript
library(getopt)
library(SAVER)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'counts', 'c', 1, 'character',
  'saver', 's', 1, 'character',
  'cores', 'p', 2, 'numeric'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$counts) | !is.null(opt$help) | is.null(opt$saver)) {
    cat(getopt(opts, usage=TRUE))
    stop('-s and -c are mandatory')
}

cores <- 1
if (!is.null(opt$cores)) {
    cores <- opt$cores
}
WD <- dirname(opt$counts)
INPUT <- basename(opt$counts)
SEPARATOR <- ","
print(WD)
print(INPUT)
cwd <- getwd()
setwd(WD)
raw.data <- read.table(INPUT, sep=SEPARATOR, header = TRUE, row.names=1)
#dataset <- as.matrix(raw.data)
dim(raw.data)
setwd(cwd)
save.image('pippo.Rdata')
dataset.saver <- saver(raw.data, ncores = cores, estimates.only = TRUE)

write.table(dataset.saver, gzfile(opt$saver), sep=SEPARATOR, col.names=NA)
