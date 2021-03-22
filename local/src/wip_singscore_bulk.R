library(pheatmap)

h1 <- read.table('CRC0327_NT_1_h.scores.tsv', sep="\t")
h2 <- read.table('CRC0327_NT_2_h.scores.tsv', sep="\t")
colnames(h1) <- paste0('NT1_', colnames(h1))
colnames(h2) <- paste0('NT2_', colnames(h2))
hallm <- cbind(h1, h2)
pheatmap(hallm)
pheatmap(hallm[complete.cases(hallm),], clustering_distance_rows='correlation', clustering_distance_cols='correlation')

h1 <- read.table('CRC0327_NT_1_c2.scores.tsv', sep="\t")
h2 <- read.table('CRC0327_NT_2_c2.scores.tsv', sep="\t")
colnames(h1) <- paste0('NT1_', colnames(h1))
colnames(h2) <- paste0('NT2_', colnames(h2))
hallm <- cbind(h1, h2)
pheatmap(hallm[complete.cases(hallm),], show_rownames = FALSE)
pheatmap(cor(hallm[complete.cases(hallm),]))

require(TxDb.Hsapiens.UCSC.hg38.knownGene) 

geneLengths <- function(egs)
{
  #require(org.Hs.eg.db)
  exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by='gene')    
  #egs    = unlist(  mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
  sapply(egs,function(eg)
  {
    exons = exons.db[[eg]]
    if (! is.null(exons)) {
      exons = reduce(exons)
      sum( width(exons) )
    } else {
       0
    } 
  })
}


d <- read.table(gzfile('CRC0327_NT_1_entrez_shaved.tsv.gz'), sep="\t", header=T, row.names=1)

l <- geneLengths(as.character(rownames(d)))
l <- l[l!=0]
d <- d[rownames(d) %in% names(l),]
d1 <- t(sapply(rownames(d), function(x) { d[rownames(d)==x,]/l[names(l)==x]} ))
write.table(d1, 'test_1NT.tsv.gz', sep="\t", quote=F, row.names=TRUE, col.names=T)


d <- read.table(gzfile('CRC0327_NT_2_entrez_shaved.tsv.gz'), sep="\t", header=T, row.names=1)

l <- geneLengths(as.character(rownames(d)))
l <- l[l!=0]
d <- d[rownames(d) %in% names(l),]
d2 <- t(sapply(rownames(d), function(x) { d[rownames(d)==x,]/l[names(l)==x]} ))
write.table(d2, 'test_2NT.tsv.gz', sep="\t", quote=F, row.names=TRUE, col.names=T)

#egrassi@godot:/mnt/trcanmed/snaketree/prj/scRNA/dataset/CRC0327_pseudobulks$ /mnt/trcanmed/snaketree/prj/scRNA/local/bin/singscore -s h.rds -o test_1_h.scores.tsv -e test_1NT.tsv.gz 
#test_1_h.scores.tsv


h1 <- read.table('test_1_h.scores.tsv', sep="\t")
h2 <- read.table('test_2_h.scores.tsv', sep="\t")
colnames(h1) <- paste0('NT1_', colnames(h1))
colnames(h2) <- paste0('NT2_', colnames(h2))
hallm <- cbind(h1, h2)
pheatmap(hallm)
pheatmap(cor(hallm))

############ pippo
ncores <- 10
permuteResult <-
  generateNull(
    upSet = geneset[[1]],
    rankData = rankData,
    centerScore = TRUE,
    knownDirection = TRUE,
    B = 1000,
    ncores = ncores,
    seed = 1,
    useBPPARAM = NULL
  )

pvals <- getPvals(permuteResult, scoredf, subSamples = 1:5)
