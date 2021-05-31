library(getopt)
library(SAVER)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'counts', 'c', 1, 'character',
  'saver', 's', 1, 'character'), ncol=4, byrow=TRUE)
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

cwd <- getwd()
setwd(WD)
raw.data <- read.table(INPUT, sep=SEPARATOR, header = TRUE, row.names=1)
dataset <- as.matrix(raw.data)
dim(dataset)
setwd(cwd)
dataset.saver <- saver(dataset, ncores = cores, estimates.only = TRUE)

write.table(dataset.saver, gzfile(opt$saver), sep=SEPARATOR, col.names=NA)
