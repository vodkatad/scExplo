library(ggplot2)
d <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_sign.tsv.gz'), sep="\t", header=T)

d$plasticity <- d$esophagus4-d$colon3

ggplot(data=d, aes(x=onf2, y=plasticity))+geom_point(size=0.5)+geom_smooth(method='lm')+theme_bw(base_size = 20)

d$type <- ifelse(grepl(".T", rownames(d), fixed=T), 'tumor', 'normal')
ggplot(data=d, aes(x=onf2, y=stem1))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

do <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_sign_others.tsv.gz'), sep="\t", header=T)
m <- merge(d, do, by="row.names")

ggplot(data=m, aes(x=onf2, y=coreHRC1))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=m, aes(x=onf2, y=ganesh_moorman_fetal2))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=m, aes(x=coreHRC1, y=ganesh_moorman_fetal2))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


d2 <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_entropy.tsv.gz'), sep="\t", header=T)
m <- merge(d, d2, by="row.names")

ggplot(data=m, aes(x=plasticity, y=entropy))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


ggplot(data=m, aes(x=onf2, y=entropy))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=m, aes(x=stem1, y=entropy))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


ggplot(data=d, aes(x=type, y=onf2))+geom_boxplot()+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=d, aes(x=type, y=stem1))+geom_boxplot()+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

#     217 esophagus_mucosa

#?

dd <- d[grepl('SMC', rownames(d)),]
ggplot(data=dd, aes(x=onf2, y=plasticity))+geom_point(size=0.5)+geom_smooth(method='lm')+theme_bw(base_size = 20)

ggplot(data=dd, aes(x=onf2, y=stem1))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

