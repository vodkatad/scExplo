library("GSEABase")
library('GSVA')
library(dplyr)
library(tidyr)
library(reshape)
library(ggplot2)
##### load Moorman signature as R object
sign_f <- '~/2023-07-12062C-Supplementary_Table_4.txt'
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
##### load tmm whole
##### subset to the wanted P-M (or not)
#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ grep X_BASALE /mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe | grep -f /tmp/pm > ~/pm_basale_RNASeq
samples_f <- '~/pm_basale_RNASeq'
tmm_f <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz' # H_
expr_data <- read.table(gzfile(tmm_f), sep="\t", header=TRUE, row.names=1)
rownames(expr_data) <- gsub("H_", "", rownames(expr_data))

samples <- read.table(samples_f, sep="\t", header=F)

#expr_w <- expr_data[, samples$V1]
expr_w <- expr_data[, grepl('LMX', colnames(expr_data)) | grepl('PRX', colnames(expr_data))]
#expr_w <- expr_data
##### score
#expr_data <- log(expr_data+1, base=2) # train test cesta has run with this also CMP
#ssgsea.norm
#Barbie  et  al.   (2009)  normalizing  the  scores  by  the  absolute  difference
#between the minimum and the maximum,  as described in their paper.   Whenssgsea.norm=FALSEthis last normalization step is skipped

gsva_par <- gsvaParam(exprData=as.matrix(expr_w), geneSets=signatures_s,  kcdf="Gaussian")
res <- gsva(gsva_par, verbose=T)

models <- unique(substr(colnames(res), 0, 7))

tres <- as.data.frame(t(res))


tres$gen <- rownames(tres)
long <- melt(tres, id.vars='gen')
long$smodel <- substr(long$gen,0,7)
long$model <- substr(long$gen ,0,10)
long <- long[long$gen %in% samples$V1,] #####

ave_long <- long %>% group_by(model, variable) %>%
  summarise(
    score = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

ave_long <- as.data.frame(ave_long)
ave_long$pm <- substr(ave_long$model, 9, 9)
ave_long$smodel <- substr(ave_long$model, 0, 7)
ave_long$pmr <- ifelse(ave_long$pm == 'R', 'P', 'M')
ave_long$pmr <- factor(ave_long$pmr, levels=c('P','M'))
ggplot(data=ave_long, aes(x=pmr, y=score))+geom_point(size=1)+geom_line(aes(group=smodel), size=0.5)+facet_wrap(~variable)+theme_bw(base_size = 20)

# plot
pairedtest <- function(sign, data) {
  myd <- data[data$variable==sign,]
  p <- myd[myd$pmr == 'P',]
  m <- myd[myd$pmr == 'M',]
  mod <- intersect(unique(p$smodel), unique(m$smodel))
  p <- p[p$smodel %in% mod,]
  m <- m[m$smodel %in% mod,]
  p <- p[order(p$smodel),]
  m <- m[order(m$smodel),]
  wilcox.test(x=m$score, y=p$score, alternative='greater', paired=T)
}

wres <- lapply(sign_l, pairedtest, ave_long)
names(wres) <- sign_l

noncan <- c('EMT', 'Injury Repair', 'Squamous', 'Osteoblast', 'Neuroendocrine', 'Endoderm Development')
ave_long$canonical <- ifelse(ave_long$variable %in% noncan, 'No', 'Yes')

for (sm in unique(ave_long$smodel)) {
  test <- ave_long[ave_long$smodel==sm,]
  print(ggplot(data=test, aes(x=pmr, y=score, color=canonical))+geom_point(size=1)+geom_line(aes(group=variable), size=0.5)+theme_bw(base_size = 20)+ggtitle(sm))
}

largest_deltas <- function(smodel, data) {
  myd <- data[data$smodel == smodel,]
  if (length(unique(myd$pmr)) == 1) {
    return(c(NA, NA))
  } else {
    deltas <- c()
    mydb <- myd
    myd <- myd[myd$canonical=='Yes',]
    for (s in unique(myd$variable)) {
      delta <- myd[myd$variable==s & myd$pmr=='M', 'score']-myd[myd$variable==s & myd$pmr=='P', 'score']
      deltas <- c(deltas, delta)
    }
    max_c <- max(deltas)
    myd <- mydb[mydb$canonical=='No',]
    deltas <- c()
    for (s in unique(myd$variable)) {
      delta <- myd[myd$variable==s & myd$pmr=='M', 'score']-myd[myd$variable==s & myd$pmr=='P', 'score']
      deltas <- c(deltas, delta)
    }
    max_nc <- max(deltas)
    return(c(max_c, max_nc))
  }
}

lll <- as.data.frame(t(sapply(unique(ave_long$smodel), largest_deltas, ave_long)))
lll <- lll[!is.na(lll$V1),]
colnames(lll) <- c('Canonical', 'Non Canonical')
library(pheatmap)
pheatmap(lll)
tt <- test[test$canonical=='Yes',]

mr <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/SourceData/bestbet_0.12_0.24.tsv',sep="\t", header=T)

mr$smodel <- substr(rownames(mr), 0, 7)
deltas <- c()
for (m in unique(mr$smodel)) {
  mm <- mr[mr$smodel == m,]
  deltas <- c(deltas, mm[1, 'intercept']-mm[2, 'intercept'])
}
dfm <- data.frame(smodel=unique(mr$smodel), deltamr=deltas)

mm <- merge(lll, dfm, by.x="row.names", by.y='smodel')
rownames(mm) <- mm$Row.names
mm$Row.names <- NULL


pheatmap(mm, scale='column')



d <- as.matrix(lll)

minv <- min(d)
maxv <- max(d)
#d[d < -4] <- -4
#d[d > 4] <- 4

neutral_value <- 0
#bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
#bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#e1e1e1", "#e1e1e1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))
#pheatmap(matrix, breaks = seq(-rg, rg, length.out = 100))
pheatmap(d, cluster_rows = F, cluster_cols=F,
         breaks = bk, color=my_palette, na_col = "#FFFFFF")



ggplot(data=mm, aes(x=deltamr, y=`Non Canonical`))+geom_point(size=1)+geom_smooth(method='lm', size=1)+theme_bw(base_size = 20)
ggplot(data=mm, aes(x=deltamr, y=`Canonical`))+geom_point(size=1)+geom_smooth(method='lm', size=1)+theme_bw(base_size = 20)
cor.test(mm$deltamr, mm$`Non Canonical`)
cor.test(mm$deltamr, mm$`Canonical`)

tt %>%  group_by(variable) %>%
  summarise(
    score = mean(value, na.rm = TRUE),
    .groups = "drop"
  )