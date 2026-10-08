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

# for ulisse's version 1.60.0 of GSEAbase
gene_sets_list <- signatures 
gsc <- GeneSetCollection(
  lapply(names(gene_sets_list), function(setName) {
    GeneSet(
      geneIds = unique(gene_sets_list[[setName]]),
      geneIdType = SymbolIdentifier(),
      setName = setName
    )
  })
)
#signatures_s <- GeneSetCollection(signatures)

#signatures_s <- geneIdsToGeneSetCollection(signatures) # godot version
##### load tmm whole
##### subset to the wanted P-M (or not)
#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ grep X_BASALE /mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe | grep -f /tmp/pm > ~/pm_basale_RNASeq
#samples_f <- '~/pm_basale_RNASeq'

#meta <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_samples_data', sep="\t", header=T)
#meta <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_samples_data', sep="\t", header=T)
meta <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_bcat.vs.scr/bcat_WT_samples_data', sep="\t", header=T)

tmm_f <- '/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_tmm.tsv.gz' # H_
expr_data <- read.table(gzfile(tmm_f), sep="\t", header=TRUE, row.names=1)

#gsva_par <- gsvaParam(exprData=as.matrix(expr_w), geneSets=signatures_s,  kcdf="Gaussian") # godot
#res <- gsva(gsva_par, verbose=T)
res <- gsva(as.matrix(expr_data), gsc, kcdf="Gaussian", 'gsva') # ulisse
## switch to jupy.sh to have a gsea that works (with param)

tres <- read.table('/mnt/trcanmed/snaketree/stash/tcf7l2_pla.tsv', sep="\t", header=T)
#tres <- read.table('/mnt/trcanmed/snaketree/stash/tcf7l2_pla_mut.tsv', sep="\t", header=T) #'/mnt/trcanmed/snaketree/stash/tcf7l2_pla_mut.tsv'
m <- merge(tres, meta, by="row.names")
m$model_r <- paste0(m$model, m$replicates)
ggplot(data=m, aes(x=geno, y=CBC))+geom_boxplot(outlier.shape=NULL)+geom_jitter(height=0, aes(color=model, group=model_r))+theme_bw()
m$geno <- factor(m$geno, levels=c('NE','N2'))
ggplot(data=m, aes(x=geno, y=onf))+geom_boxplot(outlier.shape=NA)+geom_point(aes(color=model),position = position_jitter(width = 0.3, seed = 3922))+
  geom_path(
    aes(group = model_r, color=model), 
    size = 0.4, alpha = 0.5, 
    position = position_jitter(width = 0.3, seed = 3922)
  ) +
  theme_bw()
