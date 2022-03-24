library(pheatmap)
library(reshape)

d <- read.table('/scratch/trcanmed/scExplo/dataset/mixture_models_ae_pb/H.tsv', sep="\t", header=F)
colnames(d) <- c('sample', 'pathway', 'nes', 'fdr')

d[d$fdr > 0.1, 'nes'] <- 0

wide <- as.data.frame(cast(d, formula="sample~pathway", value="nes", add.missing=TRUE, fill=0))
rownames(wide) <- wide$sample
wide$sample <- NULL
range <- max(abs(wide))
pheatmap(t(wide), breaks = seq(-range, range, length.out = 100))


d <- read.table('/scratch/trcanmed/scExplo/dataset/mixture_models_ae_pb/C2.tsv', sep="\t", header=F)
colnames(d) <- c('sample', 'pathway', 'nes', 'fdr')

d[d$fdr > 0.05, 'nes'] <- 0

wide <- as.data.frame(cast(d, formula="sample~pathway", value="nes", add.missing=TRUE, fill=0))
rownames(wide) <- wide$sample
wide$sample <- NULL
range <- max(abs(wide))
pheatmap(t(wide), breaks = seq(-range, range, length.out = 100), show_rownames = F)


#sel <- wide[,apply(wide, 2, function(x){any(x != 0)})]
sel <- wide[,apply(wide, 2, function(x){sum(x != 0)==4})]
range <- max(abs(sel))
pheatmap(t(sel), breaks = seq(-range, range, length.out = 100), show_rownames = T)
