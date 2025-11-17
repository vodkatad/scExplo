library(CRISclassifier)
setwd('~')
d <- read.csv('/mnt/cold2/snaketree/prj/PPH/local/share/data/saver_mat/filtered_annotated_saver_ribomito_CRC1502_NT_1_log2_pc1_cpm.csv')
colnames(d)[1] <- 'Symbol2'
d$Symbol <- sapply(strsplit(d$Symbol2, ':'), function(x){x[[2]][1]})
d$Symbol2 <- NULL

dd <- d[ ,c(ncol(d),seq(2, ncol(d)-1))]

# remove genes with all zeroes
sta <- apply(dd[, -1], 1, function(x) {all(x==0)})
ddd <- dd[!sta,]
# N= 3563 genes removed

# then we remove cells with poor support
summary(colMeans(ddd[,-1]))
sta_c <- apply(ddd[, -1], 2, function(x) {sum(x>3.667)})

summary(sta_c)
dddd <- ddd[, c(T, sta_c > 8382)]

write.table(dddd, file=gzfile('test2.txt.gz'), sep="\t", quote=F, row.names=F) # -> no row.names :(
cris_classifier('test2.txt.gz', output.name='test', nresmpl=1)

############################################

# Try from counts: /mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_Ire_cetuxi/CRC1502_NT_1_dir/filtered_annotated_CRC1502_NT_1.tsv
#Then, we filtered out low-quality cells with less than 1000 genes supported by at least 4 reads. In total, 4616 cells passed all the described criteria
# CPM

d <- read.csv('/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_Ire_cetuxi/CRC1502_NT_1_dir/filtered_annotated_CRC1502_NT_1.tsv')
colnames(d)[1] <- 'Symbol2'
d$Symbol <- sapply(strsplit(d$Symbol2, ':'), function(x){x[[2]][1]})
d$Symbol2 <- NULL

dd <- d[ ,c(ncol(d),seq(2, ncol(d)-1))]

# remove genes with all zeroes
sta <- apply(dd[, -1], 1, function(x) {all(x==0)})
ddd <- dd[!sta,]
# N= 3755 genes removed

# then we remove cells with poor support
sta_c <- apply(ddd[, -1], 2, function(x) {sum(x>=4)})
#> summary(sta_c)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   118.0   482.0   529.7   791.0  4324.0 

summary(sta_c)
dddd <- ddd[, c(T, sta_c > 1500)]


sta <- apply(dddd[, -1], 1, function(x) {all(x==0)})
dddd <- dddd[!sta,]
# other 648...

# CPM
gs <- dddd[,1, drop=F]
expr <- dddd[,-1]
col_sum <- apply(expr, 2, sum)
tmp1 <- t(expr)/col_sum
tmp1 <- t(tmp1)
tmp1 <- tmp1 * 1000000
cpm_expr <- log2(tmp1+1)

ddddd <- cbind(gs, tmp1)
write.table(ddddd, file=gzfile('test_nosaver.txt.gz'), sep="\t", quote=F, row.names=F) # -> no row.names :(
cris_classifier('test_nosaver.txt.gz', output.name='test', nresmpl=1)

# also tried without cpm to no avail..

#CRIS classifier uses a dataset-wide computation of the expression Z-score value for each gene in each sample: these Z-score normalized profiles are then compared with the five CRIS class templates (i.e., centroids) 

#For single-cell data, all genes with at least one read per gene in a sample were included, for a total of 95% of the CRIS gene signatures.

#Benjamini–Hochberg false discovery rate (BH.FDR) < 0.2, as previously reported


d <- read.csv('/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_Ire_cetuxi/CRC1502_NT_1_dir/filtered_annotated_CRC1502_NT_1.tsv')
colnames(d)[1] <- 'Symbol2'
d$Symbol <- sapply(strsplit(d$Symbol2, ':'), function(x){x[[2]][1]})
d$Symbol2 <- NULL

dd <- d[ ,c(ncol(d),seq(2, ncol(d)-1))]

# remove genes with all zeroes
sta <- apply(dd[, -1], 1, function(x) {all(x==0)})
ddd <- dd[!sta,]
# N= 3755 genes removed

# then we remove cells with poor support
sta_c <- apply(ddd[, -1], 2, function(x) {sum(x>=4)})
#> summary(sta_c)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   118.0   482.0   529.7   791.0  4324.0 

summary(sta_c)
dddd <- ddd[, c(T, sta_c > 1000)]


sta <- apply(dddd[, -1], 1, function(x) {all(x==0)})
dddd <- dddd[!sta,]
# other 648...

# CPM
gs <- dddd[,1, drop=F]
expr <- dddd[,-1]
col_sum <- apply(expr, 2, sum)
tmp1 <- t(expr)/col_sum
tmp1 <- t(tmp1)
tmp1 <- tmp1 * 1000000
cpm_expr <- log2(tmp1+1)

ddddd <- cbind(gs, tmp1)
write.table(ddddd, file=gzfile('test_nosaver.txt.gz'), sep="\t", quote=F, row.names=F) # -> no row.names :(
cris_classifier('test_nosaver.txt.gz', output.name='test', nresmpl=1)

# unire più campioni? TODO

# provare con quelli di Ba:
# Fig 2, il più het https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01176-5#Sec2
d <- read.csv('/mnt/trcanmed/snaketree/stash/CRC0022.csv')
d1 <- read.csv('/mnt/trcanmed/snaketree/stash/CRC0066.csv')
d2 <- read.csv('/mnt/trcanmed/snaketree/stash/CRC0076.csv')
d3 <- read.csv('/mnt/trcanmed/snaketree/stash/CRC0177.csv')
#d4 <- read.csv('/mnt/trcanmed/snaketree/stash/CRC0475.csv')
d5 <- read.csv('/mnt/trcanmed/snaketree/stash/CRC0515.csv')

dorig <- d
d <- cbind(d, d1[,-1])
d <- cbind(d, d2[,-1])
d <- cbind(d, d3[,-1])
#d <- cbind(d, d4[,-1])
d <- cbind(d, d5[,-1]) # 15766 cells, matches with Claudio
# ensg to gs
map <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_GSE132465/ensg_gs_map.tsv', header=T, sep="\t")

colnames(d)[1] <- 'Symbol2'
m <- merge(d, map, by.x="Symbol2", by.y="ens") # warning just renames sames barcode I think
m$Symbol2 <- NULL

dd <- m[ ,c(ncol(m),seq(2, ncol(m)-1))]
colnames(dd)[1] <- 'Symbol'

# remove genes with all zeroes
sta <- apply(dd[, -1], 1, function(x) {all(x==0)})
ddd <- dd[!sta,]

# then we remove cells with poor support
sta_c <- apply(ddd[, -1], 2, function(x) {sum(x>=4)})

summary(sta_c)
dddd <- ddd[, c(T, sta_c > 1000 & sta_c < 1500)] # < 1500 to badly remove doublets? nopeg TODO try doublet finder
#dddd <- ddd[, c(T, sta_c > 1000)]

sta2 <- apply(dddd[, -1], 1, function(x) {all(x==0)})
dddd <- dddd[!sta2,]

# CPM
gs <- dddd[,1, drop=F]
expr <- dddd[,-1]
col_sum <- apply(expr, 2, sum)
tmp1 <- t(expr)/col_sum
tmp1 <- t(tmp1)
tmp1 <- tmp1 * 1000000
cpm_expr <- log2(tmp1+1)

ddddd <- cbind(gs, tmp1)
#ddddd <- cbind(gs, cpm_expr)
write.table(ddddd, file=gzfile('test_nosaver.txt.gz'), sep="\t", quote=F, row.names=F) # -> no row.names :(
cris_classifier('test_nosaver.txt.gz', output.name='test', nresmpl=1)

cris_multilabel('test_nosaver.txt.gz', output.name='testmulti', nresmpl=1)
# In total, 4616 cells passed all the described criteria.  these are 5488 - doublets are an issue?
