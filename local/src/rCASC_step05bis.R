#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'output', 'o', 1, 'character',
  'cycle', 'y', 1, 'character',
  'output_c', 't', 1, 'character'
), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (!is.null(opt$help) | is.null(opt$output) | is.null(opt$cycle) | is.null(opt$output_c)) {
    cat(getopt(opts, usage=TRUE))
    stop('-t, -y and -o are mandatory')
}

SEPARATOR <- ','
output_d <- opt$output # this is the clustering file not output!!!
cycle_f <- opt$cycle
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

#associating cell cycle to clusters
save.image('pippo.Rdata')
cc.table <- matrix(data= rep(0, (3 * CLS)), ncol=3)
cc.table <- as.data.frame(cc.table)
names(cc.table) <- c("G1","G2M","S")
rownames(cc.table) <- cls.u
#rownames(cc.table) <- paste("cl", cls.u, sep="")


for(i in cls.u){
  cls.tmp <- cls[which(cls$Belonging_Cluster == i),]
  cc.tmp <- cc[which(cc[,1] %in% rownames(cls.tmp)),]
  cc.df <- as.data.frame.matrix(table(cc.tmp))
  cc.df.sum <- apply(cc.df, 2, sum)
  #cc.table[i,which(names(cc.table) %in% names(cc.df.sum))] <- as.numeric(cc.df.sum)
  cc.table[which(rownames(cc.table)==i),which(names(cc.table) %in% names(cc.df.sum))] <- as.numeric(cc.df.sum)
}
cc.table <- cc.table[order(as.numeric(rownames(cc.table))),]
rownames(cc.table) <- paste("cl", rownames(cc.table), sep="")

t1 <- table(cls$Belonging_Cluster)
t2 <- rowSums(cc.table)
if (!all(as.numeric(t1)==as.numeric(t2))) {
	stop('Qualquadra non cosa in cc-cl!')
}

write.table(cc.table, opt$output_c, sep="\t", col.names=NA)
