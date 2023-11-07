# library(GEOquery)
# gset <- getGEO("GSE108277", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# ex <- exprs(gset)
# # log2 transform
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# ex <- log2(ex) }
library(lumi)
load('/scratch/trcanmed/scExplo/local/share/data/xeno_loess_final.rda')
#xeno <-  load('/scratch/trcanmed/scExplo/local/share/data/xeno_loess_full.rda')
xeno_final <- xeno_loess # filtrato per probe non crosshyb + detected e probe con varianza maggiore scelta per ogni gene

#load('/scratch/trcanmed/scExplo/local/share/data/xeno_loess_full.rda')
#xeno <- xeno_loess # dato totale
metadata <-  pData(xeno_final)
metadata_cronici <- metadata[metadata$TRATTAMENTO %in% c('CRONICO', 'NT'),]

# extract only cronici/nt
expr_df <- exprs(xeno_final)

trad_gene <- pData(featureData(xeno_final))[, 'TargetID', drop=FALSE]
# already selected 1 probe per gene?
stopifnot(nrow(trad_gene) == length(unique(trad_gene$TargetID)))

m_expr_df <- merge(expr_df, trad_gene, by="row.names")
m_expr_df$Row.names <- NULL
rownames(m_expr_df) <- m_expr_df$TargetID
m_expr_df$TargetID <- NULL

#cronici_df <- expr_df[, colnames(expr_df) %in% gsub('X', '', rownames(metadata_cronici), fixed=TRUE)]


metadata$id <- gsub('X', '', rownames(metadata), fixed=TRUE)
write.table(metadata, file="/scratch/trcanmed/scExplo/local/share/data/metadata_GSE108277.tsv", sep="\t", quote=FALSE, row.names=FALSE)

write.table(m_expr_df, file=gzfile("/scratch/trcanmed/scExplo/local/share/data/expr_GSE108277.tsv.gz"), sep="\t", quote=FALSE, row.names=FALSE)
