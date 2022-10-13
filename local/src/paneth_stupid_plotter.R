library(ggplot2)

tsne_f <- snakemake@input[["tsne"]]
raf_f <- snakemake@input[["raf"]]
outplot_f <- snakemake@output[["plot"]]

pval_thr <- as.numeric(snakemake@wildcards[["pval"]])

dir.create(file.path(snakemake@params[["outdir"]]), showWarnings = FALSE)

paneth_pvals <- read.table(raf_f, sep=",", header=TRUE, stringsAsFactors=FALSE)
tsne_coords <- read.table(tsne_f, sep=",", header=TRUE, stringsAsFactors=FALSE)
colnames(paneth_pvals)[1] <- 'cellName_phase'
paneth_pvals$cellName <- sapply(paneth_pvals$cellName_phase, function(x) {strsplit(x, '_')[[1]][1]})
save.image('p.Rdata')
merged_data <- merge(paneth_pvals, tsne_coords, by="cellName")
# TODO where is the lost single cell? add check here and ask Raf

merged_data$ispaneth <- ifelse(merged_data$V1 < pval_thr, 'Paneth', 'Not Paneth')

p <- ggplot(merged_data, aes(x=xChoord, y=yChoord, color=ispaneth)) +
        geom_point() +theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
        ) + scale_color_manual(values=c("black","darkgoldenrod")) +
        theme(legend.title=element_blank()) +
        theme_bw() + ggtitle(snakemake@wildcards[["sample"]])


ggsave(outplot_f, plot=p)