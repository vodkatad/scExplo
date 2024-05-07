library(ggplot2)
library(ggrepel)
LGR5_f <- '/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMO_BASALE-LGR5_tmm.tsv'
KRT20_f <- '/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMO_BASALE-KRT20_tmm.tsv'

# TODO remember to look at replicates

LGR5 <- read.table(LGR5_f, sep="\t", header=F)
KRT20 <- read.table(KRT20_f, sep="\t", header=F)


LGR5_f <- '/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMO_BASALE-LGR5_tmm_ave.tsv'
KRT20_f <- '/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMO_BASALE-KRT20_tmm_ave.tsv'

# TODO remember to look at replicates

LGR5 <- read.table(LGR5_f, sep="\t", header=T)
KRT20 <- read.table(KRT20_f, sep="\t", header=T)
colnames(LGR5)[1] <- 'LGR5'
colnames(KRT20)[1] <- 'KRT20'
m <- merge(LGR5, KRT20, by="model")

w <- c('CRC0322','CRC0327', 'CRC0069','CRC0542')
m$sel <- ifelse(m$model %in% w, 'yes','no')
m$label <- ifelse(m$model %in% w, as.character(m$model), '')
ggplot(data=m, aes(x=log2(LGR5), y=log2(KRT20), color=sel, label=label))+geom_point()+scale_color_manual(values=c('black','red'))+geom_text_repel()
