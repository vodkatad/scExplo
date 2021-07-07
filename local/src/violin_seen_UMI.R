#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'vande', 'v', 1, 'character',
  'clusters', 'c', 1, 'character',
  'data', 'd', 1, 'character',
  'output', 'o', 1, 'character',
  'thr','t', 1, 'numeric'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$vande) | is.null(opt$clusters) | is.null(opt$data) | is.null(opt$thr) | is.null(opt$output) | !is.null(opt$help)) {
  cat(getopt(opts, usage=TRUE))
  stop('-v, -c, -d, -t and - o are mandatory')
}


library(ggplot2)
print(opt$thr)
cl <- read.table(opt$clusters, sep=",", header=TRUE)
data <- read.table(opt$data, sep=",", header=TRUE, row.names=1)
data2 <- read.table(opt$vande, sep=",", header=TRUE)
genes <- strsplit(as.character(data2$X),":")
gg <- unlist(lapply(genes, function(x) {x[1]}))
data3 <- data[rownames(data) %in% gg,]
save.image('pippo.Rdata')
less <- apply(data3, 2, function(x) {sum(x<opt$thr)})
clu <- unique(cl$Belonging_Cluster)
#clus <- lapply(clu, function(x) { cells <- cl[cl$Belonging_Cluster ==  x, 'cellName']; less[names(less) %in% cells] })
#res <- data.frame(cl= cl$Belonging_Cluster[order(cl$Belonging_Cluster)], genes_lower=unlist(clus))
#clu <- unique(cl$Belonging_Cluster)
clus <- lapply(clu, function(x) { cells <- cl[cl$Belonging_Cluster ==  x, 'cellName']; less[names(less) %in% cells] })
res <- data.frame(cl= cl$Belonging_Cluster[order(cl$Belonging_Cluster)], genes_lower=unlist(clus))
res$cl <- as.factor(res$cl)
ggplot(data=res, aes(x=cl, y=genes_lower))+geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape = NA)+geom_point(alpha=0.5, size=0.7)+theme_bw()

ggsave(opt$output)
