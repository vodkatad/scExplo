library(openxlsx)

de_xlsx <- snakemake@input[['xlsx']]
dir_dir <- snakemake@output[['dir']] # we always will hail dir dir as our lord and master

data <- read.xlsx(de_xlsx, sheet="all_de", rowNames=TRUE)

cols <- colnames(data)
lfc_cols <- cols[grepl('log2FoldChange', cols, fixed=TRUE)]

if (dir.exists(dir_dir)) {
   unlink(dir_dir)
}
dir.create(dir_dir)
setwd(dir_dir)
genes <- sapply(strsplit(rownames(data), '.', fixed=TRUE),  '[[', 2) 
for (c in lfc_cols) {
  name <- strsplit(c, '.', fixed=TRUE)[[1]][1]
  res <- data.frame(gene=genes, Freq=data[, c]) # colnames are fixed in GSEA_analysis.R from the freqs of outliers in biobanca
  res <- res[order(res$Freq),]
  fn <- paste0(name, '.rnk')
  write.table(res, fn, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
