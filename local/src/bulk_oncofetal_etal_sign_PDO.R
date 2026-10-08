##/mnt/trcanmed/snaketree/stash/
library("GSEABase")
library('GSVA')
library(dplyr)
library(tidyr)
library(reshape)
library(ggplot2)
##### load Moorman signature as R object
sign_f <- '/mnt/trcanmed/snaketree/stash/collection_signature_def.tsv'
sign <- read.table(sign_f, sep="\t", header=T)
sign <- sign[sign$Kept, ]
sign_l <- unique(sign$Annotation)

make_gene_set <- function(signature, data) {
  data <- data[data$Annotation == signature,]
  #geneset <- GeneSet(data$Gene, geneIdType=SymbolIdentifier()) # does not work? Y?
  #setName(geneset) <- signature
  #return(geneset)
  return(data$Gene)
}

signatures <- lapply(sign_l, make_gene_set, sign)
names(signatures) <- sign_l
#signatures_s <- GeneSetCollection(signatures)

signatures_s <- geneIdsToGeneSetCollection(signatures)
##################################################################
##### load tmm whole
##### subset to the wanted P-M (or not)
#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ grep X_BASALE /mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe | grep -f /tmp/pm > ~/pm_basale_RNASeq
#samples_f <- '~/pm_basale_RNASeq'

meta <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe', sep="\t", header=T)
meta <- meta[grepl('LMO_BASALE', meta$type),]
w <- meta$sample_id_R

tmm_f <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz' # H_
expr_data <- read.table(gzfile(tmm_f), sep="\t", header=TRUE, row.names=1)
rownames(expr_data) <- gsub("H_", "", rownames(expr_data))

#samples <- read.table(samples_f, sep="\t", header=F)

expr_w <- expr_data[, w]
#expr_w <- expr_data[, grepl('LMX', colnames(expr_data)) | grepl('PRX', colnames(expr_data))] # ridondante
#expr_w <- expr_datam$stem_index <- m$RSC3 - m$CBC5
##### score
#expr_data <- log(expr_data+1, base=2) # train test cesta has run with this also CMP
#ssgsea.norm
#Barbie  et  al.   (2009)  normalizing  the  scores  by  the  absolute  difference
#between the minimum and the maximum,  as described in their paper.   Whenssgsea.norm=FALSEthis last normalization step is skipped

gsva_par <- gsvaParam(exprData=as.matrix(expr_w), geneSets=signatures_s,  kcdf="Gaussian")
res <- gsva(gsva_par, verbose=T)

models <- unique(substr(colnames(res), 0, 7))

tres <- as.data.frame(t(res))
tres$stem_index <- tres$RSC - tres$CBC

write.table(tres, file='/mnt/trcanmed/snaketree/stash/onf_etal_scores_all_pdo.tsv', sep="\t", quote=F)

tres$gen <- rownames(tres)
long <- melt(tres, id.vars='gen')
long$smodel <- substr(long$gen,0,7)
long$model <- substr(long$gen ,0,10)

long$type <- ifelse(grepl('LMX', long$model), 'Met', 'Pri')
ggplot(data=long[long$variable=='onf', ], aes(x=reorder(gen, value), y=value, fill=type))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                                             axis.text.x=element_blank(),
                                                                                                             axis.ticks.x=element_blank())

long$type <- factor(long$type, levels=c('Pri', 'Met'))
ggplot(data=long, aes(x=type, y=value))+ geom_violin()+
  geom_boxplot(width=0.3) + facet_wrap(~variable)

ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s,]
  print(s)
  print(wilcox.test(myd[myd$type=='Pri', 'value'], myd[myd$type=='Met', 'value']))
}

## TODO Marti chiedere media repliche e rifare stat
# todo provare pairs (vedere bulk_moorman_sign_PriMets.R)
# TODO KRAS da Fra
# correlazioni tra coppie signature

ggplot(data=long, aes(x=reorder(gen, value), y=value, fill=type))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                     axis.text.x=element_blank(),
                                                                                     axis.ticks.x=element_blank())+facet_wrap(~variable)

d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_update0625.tsv', sep='\t', header=T)
krasmut <- unique(d[d$KRAS!='wt', 'CASE'])
long$KRAS <- ifelse(long$smodel %in% krasmut, 'Mut', 'WT')

long$KRAS <- factor(long$KRAS, levels=c('WT', 'Mut'))

ggplot(data=long[long$variable=='onf', ], aes(x=reorder(gen, value), y=value, fill=KRAS))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                                             axis.text.x=element_blank(),
                                                                                                             axis.ticks.x=element_blank())
ggplot(data=long, aes(x=KRAS, y=value))+ geom_violin()+
  geom_boxplot(width=0.3) + facet_wrap(~variable)


ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s,]
  print(s)
  print(wilcox.test(myd[myd$KRAS=='WT', 'value'], myd[myd$KRAS=='Mut', 'value']))
}
##################

# Ma evil

tmm_f <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/MA_evil/tmm.tsv.gz' # H_
expr_data <- read.table(gzfile(tmm_f), sep="\t", header=TRUE, row.names=1)
#rownames(expr_data) <- gsub("H_", "", rownames(expr_data))
gsva_par <- gsvaParam(exprData=as.matrix(expr_data), geneSets=signatures_s,  kcdf="Gaussian")
res <- gsva(gsva_par, verbose=T)

models <- unique(substr(colnames(res), 0, 7))

tres <- as.data.frame(t(res))
tres$stem_index <- tres$RSC - tres$CBC

###
