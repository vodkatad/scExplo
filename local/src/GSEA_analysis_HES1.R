### GSEA enrichment analysis

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- snakemake@input[["gene_res_freq"]]
GSEA_r <- snakemake@output[["GSEA_r"]]
GSEA_ridgeplot <- snakemake@output[["GSEA_ridgeplot"]]
type <- snakemake@wildcards[["msign"]]


gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)

###order

print(head(gene_res_df))
geneList <- gene_res_df[,'values']
names(geneList) <- as.character(gene_res_df[,'genes'])
geneList <- sort(geneList, decreasing = TRUE)



m_t2g <- msigdbr(species = "Homo sapiens", category = type) %>% 
  dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici



em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1,nPerm=10000)
ciao<-em[(em$enrichmentScore>0.5 | em$enrichmentScore< -0.5)& em$p.adjust<0.05,asis=TRUE]

write.table(em@result, file = GSEA_r, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

ridgeplot(ciao, showCategory = 20,decreasing=TRUE )
graphics.off()
ggsave(GSEA_ridgeplot,width = 29.7, height = 21.0, units = "cm")


