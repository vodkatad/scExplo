### GO enrichment analysis
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(dplyr)
library(xtable)

gene_res_f <- snakemake@input[["intersection"]]
gene_univ_f <- snakemake@input[["universe"]]
out_dir <- snakemake@output[["out_dir"]]
print('=============================')
print(gene_res_f)
print(gene_univ_f)
gene_res_df <- read.table(gene_res_f, quote = "", sep = ",", header = TRUE,row.names = 1)

#gene_univ_f <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/genes_residuals_universe.tsv"
gene_univ_df <- read.table(gene_univ_f, quote = "", sep = ",", header = TRUE,row.names = 1)

###order ##scegliere un treshold per discriminare quali geni tenere in base alle frequenze
#gene_freq_10 <- subset(gene_res_df, gene_res_df$Freq > threshold)
#geneList <- gene_freq_10$gene
geneList <- as.character(gene_res_df$gene)

### as.character universe
geneUni <- gene_univ_df$gene
geneUni <- as.character(geneUni)

#m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  #dplyr::select(gs_name, human_gene_symbol)

egocc <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "CC",
                pAdjustMethod = "BH",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

egomf <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "MF",
                pAdjustMethod = "BH",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

egobp <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)


dir.create(out_dir)

setwd(out_dir)
pdf('results.pdf',width=12,height=12)
barplot(egocc, showCategory = 20,fontsize_row = 10,fontsize_col = 10)+ggtitle("CC")
#ggsave("CC.pdf")
barplot(egomf, showCategory = 20,fontsize_row = 10,fontsize_col = 10)+ggtitle("MF")
#ggsave("MF.pdf")
barplot(egobp, showCategory = 20,fontsize_row = 10,fontsize_col = 10)+ggtitle("BP")
#ggsave("BP.pdf")
graphics.off()
# call also for MF, BP, produce three separated plots (either use a directory for output or three separate filenames, as you prefer :))
# rbind the ego@result in a single dataframe after having addead a column 'ontology' with CC, MF or BP , use p.adjust to correct the nominal pvalue in a single run (https://www.biostars.org/p/12182/)

egocc@result$ontology <- "CC"
egobp@result$ontology <- "BP"
egomf@result$ontology <- "MF"
egoall_df <- rbind(egocc@result, egobp@result, egomf@result)
egoall_df$p.adjust <- p.adjust(egoall_df$pvalue, method='BH')
write.table(egoall_df, file = 'GO_results.tsv', quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

#save.image('GO.Rdata')
# then print a single dataframe with all the information together
# we should understand how to do the same for GSEA...

#KEGG and REACTOME
print('converto amichetti ')
univdf <- bitr(geneUni, fromType = "SYMBOL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db) # not all maps

targetdf <- bitr(geneList, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db) # not all maps
print('convertiti')
print(targetdf)


rea <- enrichPathway(gene          = targetdf$ENTREZID,
                universe      = univdf$ENTREZID, 
                pvalueCutoff = 1, readable=FALSE)
df4 <- as.data.frame(rea)
pdf('results_KEGG_reactome.pdf',width=12,height=12)
#dotplot(kk, showCategory=20,fontsize_row = 10,fontsize_col = 10)+ggtitle("KEGG")
write.table(df4, file = 'reactome_results.tsv', quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)
dotplot(rea, showCategory=20)
graphics.off()
