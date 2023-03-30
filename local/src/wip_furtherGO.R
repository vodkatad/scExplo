library(clusterProfiler)
library(org.Hs.eg.db)
#library(ggnewscale)
#library(ReactomePA)
library(msigdbr)
library(dplyr)


universe <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/Alberto/CRC0327_cetux_2_universe.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)
genes <- read.table('/mnt/cold1/snaketree/prj/scRNA/local/share/data/Alberto_CRC0327_cetux_2.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)

univdf <- bitr(universe$gene, fromType = "SYMBOL",
               toType = c("ENTREZID"),
               OrgDb = org.Hs.eg.db) # not all maps

targetdf <- bitr(genes$gene, fromType = "SYMBOL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db) # not all maps

lost_target <- setdiff(genes$gene, targetdf$SYMBOL)
lost_univ <- setdiff(universe$gene, univdf$SYMBOL)

lost_target

msig_C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

C2 <- enricher(gene          = targetdf$ENTREZID,
              universe      = univdf$ENTREZID, TERM2GENE = msig_C2,  pvalueCutoff  = 0.05)

df5 <- as.data.frame(C2)
dotplot(C2, showCategory=nrow(df5))

