d <- read.table('~/CytoTRACE_results.txt', sep="\t")
d2 <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/SMC25-T_entropy.tsv', sep="\t")

m <- merge(d, d2, by="row.names")

plot(m$entropy, m$CytoTRACE)

d3 <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/SMC25-T_sign.tsv', sep="\t")
m <- merge(m, d3, by.y="row.names", by.x="Row.names")

plot(m$onf2, m$CytoTRACE)
plot(m$onf2, m$entropy)
  

dd <- read.csv('~/test.csv')
m2 <- merge(dd, d, by='row.names')
cor.test(m2$x, m2$CytoTRACE)
plot(m2$x, m2$CytoTRACE)
# cool


dd2 <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores//SMC25-T_ctt.tsv', sep= "\t")
m2 <- merge(dd2, d, by='row.names')
cor.test(m2$x, m2$CytoTRACE)
plot(m2$x, m2$CytoTRACE)
# cool2

m <- merge(d, d2, by="row.names")
plot(m$entropy, m$CytoTRACE)

m3 <- merge(dd2, d2, by="row.names")
plot(m$entropy, m$x)

d3 <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/SMC25-T_sign.tsv', sep="\t")
m <- merge(m, d3, by.y="row.names", by.x="Row.names")
plot(m$onf2, m$CytoTRACE)

m <- merge(m3, d3, by.y="row.names", by.x="Row.names")
plot(m$onf2, m$x)