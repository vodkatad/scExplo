#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'vande', 'v', 1, 'character',
  'scratch', 's', 1, 'character',
  'output', 'o', 1, 'character',
  'name', 'n', 1, 'character',
  'cycle', 'y', 1, 'character',
  'output_c', 't', 1, 'character',
  'cls', 'c', 1, 'numeric'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$vande)  | is.null(opt$scratch) | !is.null(opt$help) | is.null(opt$output) | is.null(opt$cls) | is.null(opt$name) | is.null(opt$cycle) | is.null(opt$output_c)) {
    cat(getopt(opts, usage=TRUE))
    stop('-v, -n -c -s -y -t and -o are mandatory')
}

SCRATCH <- opt$scratch
SEPARATOR <- ','
vande_f <- opt$vande
output_d <- opt$output # this is the clustering file not output!!!
cycle_f <- opt$cycle
CLS <- opt$cls
PROJECTNAME <- opt$name
SCRATCH <- opt$scratch
#$CLS = 5
#CC="/path_to_cell_cycle_estimation"
#CLSOUT=paste(WD,"path_to_Results",CLS,"path_to_clustering.output_file",sep="/")
#PROJECTNAME="XX"

# we determine CLS using the clustering file
cc <- read.table(cycle_f, SEPARATOR, header=F, stringsAsFactors=F)
cls <- read.table(output_d, SEPARATOR, header=T, stringsAsFactors=F, row.names=1)
cls.u <- unique(cls$Belonging_Cluster)
CLS <- length(cls.u)
print(paste0('N.clusters= ', CLS))
#pbulk
bulkClusters(group="docker", scratch.folder=SCRATCH, file=vande_f, separator=SEPARATOR, cl=output_d)

#cometsc
# n of threads == n of clusters is recommended?
cometsc(group="docker", file=vande_f, scratch.folder=SCRATCH,  threads=CLS, X=0.15, K=1, counts="True", skipvis="False", nCluster=CLS, separator=SEPARATOR)

#pblkAE
autoencoder4pseudoBulk(group="docker", scratch.folder=SCRATCH, file=vande_f, separator=SEPARATOR, permutation=50, nEpochs=1000, projectName=PROJECTNAME, bN=output_d)

save.image('pippo.Rdata')
#associating cell cycle to clusters

cc.table <- matrix(data= rep(0, (3 * CLS)), ncol=3)
cc.table <- as.data.frame(cc.table)
names(cc.table) <- c("G1","G2M","S")
rownames(cc.table) <- paste("cl", cls.u, sep="")

for(i in cls.u){
  cls.tmp <- cls[which(cls$Belonging_Cluster == i),]
  cc.tmp <- cc[which(cc[,1] %in% rownames(cls.tmp)),]
  cc.v <- c(0,0,0)
  names(cc.v) <- c("G1", "G2M", "S")
  cc.df <- as.data.frame.matrix(table(cc.tmp))
  cc.df.sum <- apply(cc.df, 2, sum)
  cc.table[i,which(names(cc.table) %in% names(cc.df.sum))] <- as.numeric(cc.df.sum)
}

cc.table <- cc.table[order(rownames(cc.table)),]
write.table(cc.table, opt$output_c, sep="\t", col.names=NA)