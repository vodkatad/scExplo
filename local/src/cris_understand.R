library(cellrangerRkit, quietly=TRUE)
library(pheatmap)


data <- snakemake@input[["data"]]
heatmap <- snakemake@output[["heatmap"]]
outable <- snakemake@output[["tsv"]]
load(data)
#save.image('p.Rdata')
# test <- wanted_mtx[,c(1,10)]
cris_data <- merge(ensg_cris, crisGenes, by.x="ALIAS",by.y="Gene.ID")
classes <- apply(wanted_mtx, 2, function(x) { expr <- x[x >= 1]; d <- cris_data[cris_data$ENSEMBL %in% names(expr),]; table(d$Gene.ID.1);})
write.table(classes, file=outable, sep="\t", quote=FALSE)
tclasses <- t(classes)
df <- as.data.frame( table(cris_data[,"Gene.ID.1"]))
garbage <- lapply(levels(df$Var1), function(x) {tclasses[,x] <<- tclasses[,x]/df[df$Var1==x,"Freq"]})
pheatmap(as.matrix(tclasses), scale = "none", show_rownames=F, show_colnames = T, cluster_cols=F, filename=heatmap)
#sink(file=snakemake@output[["simple"]])
#sinkf <- snakemake@output[["simple"]]
cat('called cris', nrow(wanted_mtx))
cat('expressed cris', nrow(wanted_mtx_expr))
tt <- table(cris_data[cris_data$ENSEMBL %in% n,"Gene.ID.1"])
df2 <- as.data.frame(tt)
m <- merge(df, df2, by="Var1")
m$p <- m$Freq.y/m$Freq.x
write.table(as.data.frame(m), sep="\t", quote=F)
