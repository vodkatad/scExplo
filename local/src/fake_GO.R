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
#library(ReactomePA)
library(msigdbr)
library(dplyr)
#library(xtable)

gene_res_f <- '/home/sborgato/Sofia_Rstudio/GO_list/cas9_solo_327_UP.tsv' #selezionare lista GENI ATTENZIONEEEEEE!
gene_univ_f <- '/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/ko_atoh1/EGF_vs_CET_327only/fpkm.tsv.gz'

universe<-read.table(gene_univ_f, quote = "", sep = "\t", header = TRUE)
fpkm<-as.data.frame(universe)
universe$gene<-row.names(universe)


gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE,row.names = 1)


geneList <- as.character(gene_res_df$gene)


geneUni <- universe$gene
geneUni <- as.character(geneUni)


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


pdf("fake_GO/GO_results.pdf")
barplot(egocc, showCategory = 30,fontsize_row = 10,fontsize_col = 10)+ggtitle("CC")

barplot(egomf, showCategory = 30,fontsize_row = 10,fontsize_col = 10)+ggtitle("MF")

barplot(egobp, showCategory = 30,fontsize_row = 10,fontsize_col = 10)+ggtitle("BP")
graphics.off()
egocc@result$ontology <- "CC"
egobp@result$ontology <- "BP"
egomf@result$ontology <- "MF"
egoall_df <- rbind(egocc@result, egobp@result, egomf@result)
egoall_df$p.adjust <- p.adjust(egoall_df$pvalue, method='BH')
write.table(egoall_df, file = 'fake_GO/GO_results.tsv', quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)
#cazzo<-as.data.frame(egomf)
#geni_wnt<-cazzo[cazzo$ID=='GO:0090090','geneID']
#geni_wnt<-strsplit(geni_wnt,'/')
#geni_wnt<-unlist(geni_wnt)
#fpkm[geni_wnt,]
#ATTENZIONE SCEGLI tra MF CC o BP
geni<-egomf[egomf$p.adjust<0.05,'geneID']
geni_vec<-c()
conto<-0
for (el in geni){
  conto=conto+length(strsplit(el,'/')[[1]])
  tmp<-strsplit(el,'/')
  for (e in tmp){
    geni_vec<-c(geni_vec,e)
  }
  
}

print(unique(geni_vec))
#save.image('GO.Rdata')
# then print a single dataframe with all the information together
# we should understand how to do the same for GSEA...

#KEGG and REACTOME

