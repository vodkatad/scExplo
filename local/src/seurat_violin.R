infile <- snakemake@input[["data"]]
my_outdir <- snakemake@output[["outdir"]]
genes_comma <- snakemake@params[['genes']]
set.seed(42)

genes <- unlist(strsplit(genes_comma,','))
dir.create(my_outdir)

load(infile)
library("Seurat")
library(ggplot2)

violin <- function(gene, srdata) {
    tryCatch({
        VlnPlot(srdata, features=gene)
        ggsave(paste0(my_outdir, '/', gene, '_vln.pdf'))
    },
    error=function(cond) {
        message(cond)
    })
}
garbage <- lapply(as.list(genes), violin, srdata)