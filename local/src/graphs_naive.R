library(pheatmap)
d <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/Graph_NN/notebook/NT_GNN_vanilla/filtered_CRC0322_NT_1_3000.csv', row.names=1, sep=",", header=T, stringsAsFactors = F)

count_expr <- function(x, thr=1) {
  sum(x>thr)
}

ll <- apply(d, 1, count_expr, thr=0)
ll2 <- apply(d, 1, count_expr, thr=1)
ll3 <- apply(d, 1, count_expr, thr=7)

summary(apply(d, 1, count_expr, thr=20))

> summary(alll)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.000   0.000   0.000   1.711   0.000  16.710 
> summary(alll[alll!=0])
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
4.787   7.006   7.683   7.927   8.633  16.710 

nt <- list.files(path = "/mnt/cold1/snaketree/prj/scRNA/dataset/Graph_NN/notebook/NT_GNN_vanilla/")
setwd('/mnt/cold1/snaketree/prj/scRNA/dataset/Graph_NN/notebook/NT_GNN_vanilla/')
load <- function(path) {
  read.table(path, row.names=1, sep=",", header=T, stringsAsFactors = F)
}

tables <- lapply(nt, load)

lapply(tables, dim)
cg <- Reduce(intersect, lapply(tables, colnames)) #5342 genes common to all samples

genes_n_cell <- function(df, thr=0) {
  colSums(df[,cg] > thr)
}

ngenes <- sapply(tables, genes_n_cell)
colnames(ngenes) <- nt

cc <- colSums(ngenes)

dfp <- t(t(ngenes)/cc)


pheatmap(dfp, show_rownames = F)

d <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/Graph_NN/notebook/GSE132465_GNN_vanilla/filtered_SMC01-T.csv', row.names=1, sep=",", header=T, stringsAsFactors = F)

ll <- apply(d, 1, count_expr, thr=0)
ll2 <- apply(d, 1, count_expr, thr=1)
ll3 <- apply(d, 1, count_expr, thr=7)


alll <- unlist(d)
