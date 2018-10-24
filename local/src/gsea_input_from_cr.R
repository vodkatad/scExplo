clustersFile <- snakemake@input[["clusters"]]
debug <- snakemake@params[["debug"]]
out <- snakemake@output[[1]]

if (debug == "yes") {
  save.image(file=paste0(out,'.debug','.RData'))
}

data <- read.table(clustersFile, sep=",", header=TRUE)

cols <- colnames(data)

# wrong iterative idea
#clusters <- cols[grep("Cluster", cols, fixed=TRUE)]
#totcl <- max(as.numeric(sapply(clusters, function(x) {strsplit(x, '.', fixed=TRUE)[[1]][2]})))
#res <- data.frame(matrix(nrow=totcl*, ncol=4))
# then foreach cluster flatten...

library(reshape)
toremove <- grep("UMI.Count", cols, fixed=TRUE)
toremove2 <- grep("Adjusted", cols, fixed=TRUE)
toremove3 <- grep("Gene.ID", cols, fixed=TRUE)
wide_data <- data[,-c(toremove, toremove2, toremove3)]
data_long <- melt(wide_data, id.vars=c("Gene.Name"))
colnames(data_long) <- c("geneid","name","sort")
write.table(data_long, file=out, sep="\t", quote=FALSE, col.names = TRUE, row.names=FALSE)
