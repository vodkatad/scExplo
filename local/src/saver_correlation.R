#!/usr/bin/env Rscript
library(getopt)
library(SAVER)

cores <- 1
INPUT <- snakemake@input[['data']]
SEPARATOR <- ","
raw.data <- read.table(INPUT, sep=SEPARATOR, header = TRUE, row.names=1)
dataset <- as.matrix(raw.data)
saver1 <- saver(dataset, ncores = cores, estimates.only = FALSE)
corr<-cor.genes(saver1)

write.table(corr, snakemake@output[['corr']], sep=SEPARATOR)
