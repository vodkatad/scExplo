library(mixtools) 
data <- read.table('/tmp/filtered_annotated_saver_ribomito_CRC0322LMO_CTX72h_1_singleron0722_log2_pc1_cpm.csv', sep=",", header=T)
rownames(data) <- data$X
data$X <- NULL
f5 <- c('ATOH1$', 'DLL1$', 'GFI1$','DEFA5$', 'DEFA6$')
ll <- lapply(f5, function(x) { grepl(x, rownames(data))})
wanted <- ll[[1]]
for (i in seq(2, length(ll))) { wanted <- wanted | ll[[i]] }
geneListM <- rownames(data)[wanted]
geneListM
mainMatrix <- as.matrix(data)
mainMatrix2=colSums(mainMatrix[geneListM,])
head(mainMatrix2)
summary(data$AAACATCGAACCGAGAACGTATCA)
sum(data$AAACATCGAACCGAGAACGTATCA[c(1,210,300,4000,5000)])
glmone <- function(logcpm, maxit=10000) {
    mixmdltemp=try(mixmdl <- normalmixEM(logcpm,k = 2,maxit=maxit,maxrestarts=1000000000000000000))
    return(list(means=mixmdltemp$mu, post=as.data.frame(mixmdltemp$posterior)))
}
get_gene <- function(g, data) {
    data[g,]
}
extract_fit <- function(g, data) {
    mat <- get_gene(g, data)
    glmone(mat)
}
listadilista <- lapply(geneListM, extract_fit)
listadilista <- lapply(geneListM, extract_fit, mainMatrix)
head(listadilista)
names(listadilista)
listadilista[[1]]
length(listadilista)
lapply(listadilista, function(x) {x$means})
for (i in seq(1, length(geneListM))) { # cbind rather than print + add cell ids
    write.table(listadilista[[i]]$post, paste0(geneListM[i], '.csv'), sep=",", quote=F)
}
# then compare with dd <- read.table('/data/reanalysis_on_AIsc/paneth_only/all/log2cpmfisher/singctx_log2cpmfisher.csv', sep= ",", header=T, row.names=1)
# then repeat for lmx singleron
history(n=100)
savehistory('/tmp/wip.R')
