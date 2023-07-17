library(pheatmap)
data <- read.table(gzfile('/mnt/trcanmed/temp/scTimeCourse_fpkm.tsv.gz'), sep="\t", header=TRUE)
meta <- read.table('/mnt/trcanmed/temp/scTimeCourse_samples', sep="\t", header=TRUE)

#ogenes <- genes
genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_322_CTX_WNT.txt', header=FALSE)
#genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_327_2_CTX_WNT.txt', header=FALSE)
#genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_1502_CTX_WNT.txt', header=FALSE)
#genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_542_CTX_EE.txt', header=FALSE)
length(intersect(genes$V1, ogenes$V1))
nrow(genes)
nrow(ogenes)


sdata <- data[rownames(data) %in% genes$V1,]

nrow(genes)==nrow(sdata)
setdiff(genes$V1, rownames(sdata)) # TODO look for synomyms
#egrassi@godot:/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops$ grep ACP3 /mnt/trcanmed/bioinfotree/prj/paired_scRNA/dataset/v1/gene_len 
# why it's not there?

# WNT 322 proverei ctx (short time e long timecourse) normalizzato su EGF long (8 o 12g)

meta <- meta[grepl('CRC0322', meta$id),]
rownames(meta) <- meta$id
meta$id <- NULL

sdata <- sdata[, grepl('CRC0322', colnames(sdata)),]

pheatmap(log(sdata+1), annot_cols=meta)

lfpkm <- log(sdata+1)

norm <- data.frame(CRC0322_3d_CTX=lfpkm[, 'CRC0322_3d_CTX'] - lfpkm[, 'CRC0322_8d_EGF'],
                   CRC0322_7d_CTX=lfpkm[, 'CRC0322_7d_CTX'] - lfpkm[, 'CRC0322_8d_EGF'],
                   CRC0322_8d_EGF=lfpkm[, 'CRC0322_8d_EGF'] - lfpkm[, 'CRC0322_8d_EGF'],
                   CRC0322_8d_NOEGF=lfpkm[, 'CRC0322_8d_NOEGF'] - lfpkm[, 'CRC0322_8d_EGF'], row.names=rownames(sdata))

pheatmap(norm, annot_cols=meta)

minv <- min(norm)
maxv <- max(norm)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("#FFCC99", "darkred"))(n = length(bk2)-1)))

pheatmap(norm, annot_cols=meta,  breaks = bk, color=my_palette, show_rownames = FALSE)

# metagene
metag <- colMeans(lfpkm)
dp <- data.frame(metag)
dp$sample <- rownames(dp)
dp <- dp[c('CRC0322_4d_EGF', 'CRC0322_8d_EGF', 'CRC0322_12d_EGF', 'CRC0322_4d_NOEGF', 'CRC0322_8d_NOEGF', 
           'CRC0322_12d_NOEGF','CRC0322_3d_CTX' ,'CRC0322_7d_CTX'),]
dp$x <- seq(1, nrow(dp))
dp$g <- sapply(strsplit(dp$sample, "_"), function(x){x[3]})
ggplot(data=dp, aes(x=reorder(sample,x), y=metag, group=g, color=g))+
  geom_point()+geom_line()+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("metagene WNT")+xlab('')

genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_322_CTX_paneth.txt', header=FALSE)

sdata <- data[rownames(data) %in% genes$V1,]

nrow(genes)==nrow(sdata)
setdiff(genes$V1, rownames(sdata)) # TODO look for synomyms
#egrassi@godot:/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops$ grep ACP3 /mnt/trcanmed/bioinfotree/prj/paired_scRNA/dataset/v1/gene_len 
# why it's not there?

# WNT 322 proverei ctx (short time e long timecourse) normalizzato su EGF long (8 o 12g)

meta <- meta[grepl('CRC0322', meta$id),]
rownames(meta) <- meta$id
meta$id <- NULL

sdata <- sdata[, grepl('CRC0322', colnames(sdata)),]

pheatmap(log(sdata+1), annot_cols=meta)

lfpkm <- log(sdata+1)

norm <- data.frame(CRC0322_3d_CTX=lfpkm[, 'CRC0322_3d_CTX'] - lfpkm[, 'CRC0322_8d_EGF'],
                   CRC0322_7d_CTX=lfpkm[, 'CRC0322_7d_CTX'] - lfpkm[, 'CRC0322_8d_EGF'],
                   CRC0322_8d_EGF=lfpkm[, 'CRC0322_8d_EGF'] - lfpkm[, 'CRC0322_8d_EGF'],
                   CRC0322_8d_NOEGF=lfpkm[, 'CRC0322_8d_NOEGF'] - lfpkm[, 'CRC0322_8d_EGF'], row.names=rownames(sdata))

minv <- min(norm)
maxv <- max(norm)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(norm, annot_cols=meta,  breaks = bk, color=my_palette, show_rownames = FALSE)

# metagene
metag <- colMeans(lfpkm)
dp <- data.frame(metag)
dp$sample <- rownames(dp)
dp <- dp[c('CRC0322_4d_EGF', 'CRC0322_8d_EGF', 'CRC0322_12d_EGF', 'CRC0322_4d_NOEGF', 'CRC0322_8d_NOEGF', 
           'CRC0322_12d_NOEGF','CRC0322_3d_CTX' ,'CRC0322_7d_CTX'),]
dp$x <- seq(1, nrow(dp))
dp$g <- sapply(strsplit(dp$sample, "_"), function(x){x[3]})
ggplot(data=dp, aes(x=reorder(sample,x), y=metag, group=g, color=g))+
  geom_point()+geom_line()+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("metagene Paneth")+xlab('')
# cesso
#TFF3        -0.87827960   -2.392220792              0     -0.680776438
# TODO try with Sofia's cliques?

####33 ssgsea on timecourse
data <- read.table('~/Alberto_Alberto-scores.txt', sep="\t")
data <- data[, grepl('CRC0322', colnames(data))]
data <- as.data.frame(t(data))

data <- data[c('CRC0322_4d_EGF', 'CRC0322_8d_EGF', 'CRC0322_12d_EGF', 'CRC0322_4d_NOEGF', 'CRC0322_8d_NOEGF', 
           'CRC0322_12d_NOEGF','CRC0322_3d_CTX' ,'CRC0322_7d_CTX'),]
data$x <- seq(1, nrow(data))
data$g <- sapply(strsplit(rownames(data), "_"), function(x){x[3]})
data$sample <- row.names(data)
ggplot(data=data, aes(x=reorder(sample,x), y=`322_CTX_WNT`, group=g, color=g))+
  geom_point()+geom_line()+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("322_CTX_WNT")+xlab('')

# LGR5 322 proverei EGF long (8 e 12g) su CTX timecourse long 
####################genes_pop_322_EGF8d_LGR5.txt
genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_322_EGF8d_LGR5.txt', header=FALSE)

sdata <- data[rownames(data) %in% genes$V1,]

nrow(genes)==nrow(sdata)
setdiff(genes$V1, rownames(sdata)) # TODO look for synomyms
#egrassi@godot:/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops$ grep ACP3 /mnt/trcanmed/bioinfotree/prj/paired_scRNA/dataset/v1/gene_len 
# why it's not there?

meta <- meta[grepl('CRC0322', meta$id),]
rownames(meta) <- meta$id
meta$id <- NULL

sdata <- sdata[, grepl('CRC0322', colnames(sdata)),]

pheatmap(log(sdata+1), annot_cols=meta)

lfpkm <- log(sdata+1)

norm <- data.frame(CRC0322_8d_EGF=lfpkm[, 'CRC0322_8d_EGF']-lfpkm[, 'CRC0322_7d_CTX'],
                   CRC0322_12d_EGF=lfpkm[, 'CRC0322_12d_EGF'] - lfpkm[, 'CRC0322_7d_CTX'],
                   CRC0322_7d_CTX=lfpkm[, 'CRC0322_7d_CTX'] - lfpkm[, 'CRC0322_7d_CTX'],
                   CRC0322_8d_NOEGF=lfpkm[, 'CRC0322_8d_NOEGF'] - lfpkm[, 'CRC0322_7d_CTX'], row.names=rownames(sdata))

minv <- min(norm)
maxv <- max(norm)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(norm, annot_cols=meta,  breaks = bk, color=my_palette, show_rownames = FALSE)

# metagene
metag <- colMeans(lfpkm)
dp <- data.frame(metag)
dp$sample <- rownames(dp)
dp <- dp[c('CRC0322_4d_EGF', 'CRC0322_8d_EGF', 'CRC0322_12d_EGF', 'CRC0322_4d_NOEGF', 'CRC0322_8d_NOEGF', 
           'CRC0322_12d_NOEGF','CRC0322_3d_CTX' ,'CRC0322_7d_CTX'),]
dp$x <- seq(1, nrow(dp))
dp$g <- sapply(strsplit(dp$sample, "_"), function(x){x[3]})
ggplot(data=dp, aes(x=reorder(sample,x), y=metag, group=g, color=g))+
  geom_point()+geom_line()+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("metagene LGR5")+xlab('')
# KRT20 322 proverei EGF long (8 e 12g) su CTX timecourse long 
####################genes_pop_322_EGF8d_LGR5.txt
genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_322_EGF8d_KRT20.txt', header=FALSE)

sdata <- data[rownames(data) %in% genes$V1,]

nrow(genes)==nrow(sdata)
setdiff(genes$V1, rownames(sdata)) # TODO look for synomyms
#egrassi@godot:/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops$ grep ACP3 /mnt/trcanmed/bioinfotree/prj/paired_scRNA/dataset/v1/gene_len 
# why it's not there?

meta <- meta[grepl('CRC0322', meta$id),]
rownames(meta) <- meta$id
meta$id <- NULL

sdata <- sdata[, grepl('CRC0322', colnames(sdata)),]

pheatmap(log(sdata+1), annot_cols=meta)

lfpkm <- log(sdata+1)

norm <- data.frame(CRC0322_8d_EGF=lfpkm[, 'CRC0322_8d_EGF']-lfpkm[, 'CRC0322_7d_CTX'],
                   CRC0322_12d_EGF=lfpkm[, 'CRC0322_12d_EGF'] - lfpkm[, 'CRC0322_7d_CTX'],
                   CRC0322_7d_CTX=lfpkm[, 'CRC0322_7d_CTX'] - lfpkm[, 'CRC0322_7d_CTX'],
                   CRC0322_8d_NOEGF=lfpkm[, 'CRC0322_8d_NOEGF'] - lfpkm[, 'CRC0322_7d_CTX'], 
                   CRC0322_8d_NOEGF=lfpkm[, 'CRC0322_12d_NOEGF'] - lfpkm[, 'CRC0322_7d_CTX'], row.names=rownames(sdata))


minv <- min(norm)
maxv <- max(norm)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(norm, annot_cols=meta,  breaks = bk, color=my_palette, show_rownames = FALSE)

# metagene
metag <- colMeans(lfpkm)
dp <- data.frame(metag)
dp$sample <- rownames(dp)
dp <- dp[c('CRC0322_4d_EGF', 'CRC0322_8d_EGF', 'CRC0322_12d_EGF', 'CRC0322_4d_NOEGF', 'CRC0322_8d_NOEGF', 
           'CRC0322_12d_NOEGF','CRC0322_3d_CTX' ,'CRC0322_7d_CTX'),]
dp$x <- seq(1, nrow(dp))
dp$g <- sapply(strsplit(dp$sample, "_"), function(x){x[3]})
ggplot(data=dp, aes(x=reorder(sample,x), y=metag, group=g, color=g))+
  geom_point()+geom_line()+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("metagene KRT20")+xlab('')

### EE 
#542 EE si puo fare ctc vs nt perche nel nt non c'e' popolazione EE
data <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/scRNA_paired/fpkm.tsv.gz'), sep="\t", header=TRUE)
meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/scRNA_paired/samples_data', sep="\t", header=TRUE)

genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_542_CTX_EE.txt', header=FALSE)

sdata <- data[rownames(data) %in% genes$V1,]

nrow(genes)==nrow(sdata)
setdiff(genes$V1, rownames(sdata)) # TODO look for synomyms
#egrassi@godot:/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops$ grep ACP3 /mnt/trcanmed/bioinfotree/prj/paired_scRNA/dataset/v1/gene_len 
# why it's not there?

# WNT 322 proverei ctx (short time e long timecourse) normalizzato su EGF long (8 o 12g)

pheatmap(log(sdata+1), annot_cols=meta)

lfpkm <- log(sdata+1)


norm <- data.frame(CRC0542_Cetux_72h_presorting=lfpkm[, 'CRC0542_Cetux_72h_presorting']-lfpkm[, 'CRC0542_NT_72h_presorting'],
                   CRC0542_Cetux_72h_postsorting=lfpkm[, 'CRC0542_Cetux_72h_postsorting'] - lfpkm[, 'CRC0542_NT_72h_postsorting'],
                   CRC0542_Cetux_1w_presorting=lfpkm[, 'CRC0542_Cetux_1w_presorting'] - lfpkm[, 'CRC0542_NT_1w_postsorting'],
                   CRC0542_Cetux_1w_postsorting=lfpkm[, 'CRC0542_Cetux_1w_postsorting'] - lfpkm[, 'CRC0542_NT_1w_postsorting'],
                   row.names=rownames(sdata))


minv <- min(norm)
maxv <- max(norm)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(norm, annot_cols=meta,  breaks = bk, color=my_palette, show_rownames = FALSE)

########################### treated ssgsea scores Alberto
d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/Alberto_Alberto-scores.tsv', sep="\t")
ann <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/samples_data', sep="\t", header=TRUE)
rownames(ann) <- ann$id
ann$id <- NULL

ann <- ann[order(ann$treat),]
d <- d[, match(rownames(ann), colnames(d))]


minv <- min(d)
maxv <- max(d)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette)

ann <- ann[order(ann$sample, ann$treat),]
d <- d[, match(rownames(ann), colnames(d))]
pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette)

## PIS on this ssGSEA
fc <- function(model, data, sign, ann) {
  d <- data[rownames(data)==sign,]
  dt <- t(d)
  m <- merge(dt, ann, by="row.names")
  d <- m[m$sample==model, ]
  d <- d[order(d$Row.names),] # we resort to using gen id to identify treat/nt, but we check
  if (nrow(d) == 4) {
    fc <- mean(c(d[2, sign] - d[1, sign], d[4, sign] - d[3, sign])) 
    stopifnot(d[1, 'treat'] == "NT" && d[2, 'treat'] == "cetuxi" && d[3, 'treat'] == "NT" && d[4, 'treat'] == "cetuxi")
  } else {
    fc <- d[2, sign] - d[1, sign]
    stopifnot(d[1, 'treat'] == "NT" && d[2, 'treat'] == "cetuxi")
  }
  return(fc)
}

panethIndScore <- sapply(unique(ann$sample), fc, d, '322_CTX_paneth', ann)
scores <- data.frame(row.names=unique(ann$sample), PIS=panethIndScore)
scores$lPIS <- scores$PIS
scores$model <- rownames(scores)

ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS)))+geom_col(fill='blue')+ylab("delta SSGSEA CRC0322 Paneth")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))


panethIndScore <- sapply(unique(ann$sample), fc, d, '322_CTX_WNT', ann)
scores <- data.frame(row.names=unique(ann$sample), PIS=panethIndScore)
scores$lPIS <- scores$PIS
scores$model <- rownames(scores)

ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS)))+geom_col(fill='blue')+ylab("delta SSGSEA CRC0322 WNT")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))

## test on R

d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/Alberto_Alberto-scores.tsv', sep="\t")
ann <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/samples_data', sep="\t", header=TRUE)
rownames(ann) <- ann$id
ann$id <- NULL

ann <- ann[order(ann$treat),]
d <- d[, match(rownames(ann), colnames(d))]


minv <- min(d)
maxv <- max(d)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette)

ann <- ann[order(ann$sample, ann$treat),]
d <- d[, match(rownames(ann), colnames(d))]
pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette)

## PIS on this ssGSEA
fc <- function(model, data, sign, ann) {
  d <- data[rownames(data)==sign,]
  dt <- t(d)
  m <- merge(dt, ann, by="row.names")
  d <- m[m$sample==model, ]
  d <- d[order(d$Row.names),] # we resort to using gen id to identify treat/nt, but we check
  if (nrow(d) == 4) {
    fc <- mean(c(d[2, sign] - d[1, sign], d[4, sign] - d[3, sign])) 
    stopifnot(d[1, 'treat'] == "NT" && d[2, 'treat'] == "cetuxi" && d[3, 'treat'] == "NT" && d[4, 'treat'] == "cetuxi")
  } else {
    fc <- d[2, sign] - d[1, sign]
    stopifnot(d[1, 'treat'] == "NT" && d[2, 'treat'] == "cetuxi")
  }
  return(fc)
}

panethIndScore <- sapply(unique(ann$sample), fc, d, '322_CTX_paneth', ann)
scores <- data.frame(row.names=unique(ann$sample), PIS=panethIndScore)
scores$lPIS <- scores$PIS
scores$model <- rownames(scores)

ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS)))+geom_col(fill='red')+ylab("delta SSGSEA CRC0322 Paneth")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))


panethIndScore <- sapply(unique(ann$sample), fc, d, '322_CTX_WNT', ann)
scores <- data.frame(row.names=unique(ann$sample), PIS=panethIndScore)
scores$lPIS <- scores$PIS
scores$model <- rownames(scores)

ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS)))+geom_col(fill='red')+ylab("delta SSGSEA CRC0322 WNT")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))


########################### whole ssgsea scores Alberto

d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/Alberto_Alberto-scores.tsv', sep="\t")
ann <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', sep="\t", header=TRUE)
rownames(ann) <- ann$id
ann$id <- NULL

ann <- ann[order(ann$type),]
d <- d[, match(rownames(ann), colnames(d))]
ann$irino <- NULL

ann$class <- substr(rownames(ann), 8, 10)
ann$type_redux <- ifelse(grepl('BASALE', ann$type), 'BASAL', 
                                                   ifelse(grepl('cetux', ann$type),
                                                          ifelse(grepl('NT', ann$type), 'NT', 
                                                                 'CTX'), 'BASAL'))

ann2 <- ann
ann$type <- NULL
minv <- min(d)
maxv <- max(d)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette, show_colnames = FALSE)

ann <- ann[order(ann$type_redux),]
d <- d[, match(rownames(ann), colnames(d))]

pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette, show_colnames = FALSE)

ann3 <- ann[ann$class== "LMO" & ann$type_redux == "BASAL",]
d1 <- d[,colnames(d) %in%  rownames(ann3)]
pheatmap(d1, annotation_col=ann, cluster_rows = F, cluster_cols=T, breaks = bk, color=my_palette, show_colnames = FALSE)

ann3 <- ann[ann$class== "LMX" & ann$type_redux == "BASAL",]
d3 <- d[,colnames(d) %in%  rownames(ann3)]
pheatmap(d3, annotation_col=ann, cluster_rows = F, cluster_cols=T, breaks = bk, color=my_palette, show_colnames = FALSE)

ann3 <- ann[ann$class== "LMH" & ann$type_redux == "BASAL",]
d4 <- d[,colnames(d) %in%  rownames(ann3)]
pheatmap(d4, annotation_col=ann, cluster_rows = F, cluster_cols=T, breaks = bk, color=my_palette, show_colnames = FALSE)

td3 <- as.data.frame(t(d3[rownames(d3)== "322_CTX_WNT",]))
td4 <- as.data.frame(t(d4[rownames(d4)== "322_CTX_WNT",]))
td3$model <- substr(rownames(td3), 0, 7)
td4$model <- substr(rownames(td4), 0, 7)
m <- merge(td3, td4, by="model")
cor.test(m$`322_CTX_WNT.x`,m$`322_CTX_WNT.y`)

td3 <- as.data.frame(t(d3[rownames(d3)== "322_CTX_paneth",]))
td4 <- as.data.frame(t(d4[rownames(d4)== "322_CTX_paneth",]))
td3$model <- substr(rownames(td3), 0, 7)
td4$model <- substr(rownames(td4), 0, 7)
m <- merge(td3, td4, by="model")
cor.test(m$`322_CTX_paneth.x`,m$`322_CTX_paneth.y`)

td3 <- as.data.frame(t(d3[rownames(d3)== "322_CTX_paneth",]))
td1 <- as.data.frame(t(d4[rownames(d1)== "322_CTX_paneth",]))
td3$model <- substr(rownames(td3), 0, 7)
td1$model <- substr(rownames(td1), 0, 7)
m <- merge(td3, td1, by="model")
cor.test(m$`322_CTX_paneth.x`,m$`322_CTX_paneth.y`)

td3 <- as.data.frame(t(d3[rownames(d3)== "322_CTX_WNT",]))
td1 <- as.data.frame(t(d1[rownames(d1)== "322_CTX_WNT",]))
td3$model <- substr(rownames(td3), 0, 7)
td1$model <- substr(rownames(td1), 0, 7)
m <- merge(td3, td1, by="model")
cor.test(m$`322_CTX_WNT.x`,m$`322_CTX_WNT.y`)

### GSEA on PDO CEtuxI
library(fgsea)

pathways <- readRDS('/mnt/trcanmed/snaketree/prj/scRNA/dataset/CRC0327_pseudobulks/coll_genes_pop_Alberto.rds')
gene_res_df <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/treat_cutoff0.05-cetuxi.vs.NT.gseain.tsv', quote = "", sep = "\t", header = TRUE)
###order
geneList <- gene_res_df$sort
names(geneList) <- as.character(gene_res_df$geneid)
geneList <- sort(geneList, decreasing = TRUE)

p <- list()
for (i in 1:length(pathways)) {
  p[[pathways[[i]]@setName]] <- pathways[[i]]@geneIds
} 

df <- fgsea(p, geneList, minSize=15, maxSize=1000, nperm=10000, nproc=6) 
plotEnrichment(p[['322_CTX_paneth']][[as.character(row[1])]],geneList) +labs(title='322_CTX_paneth')

