library(pheatmap)
data <- read.table(gzfile('/mnt/trcanmed/temp/scTimeCourse_fpkm.tsv.gz'), sep="\t", header=TRUE)
meta <- read.table('/mnt/trcanmed/temp/scTimeCourse_samples', sep="\t", header=TRUE)

#ogenes <- genes
#genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_322_CTX_WNT.txt', header=FALSE)
#genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_327_2_CTX_WNT.txt', header=FALSE)
#genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_1502_CTX_WNT.txt', header=FALSE)
genes <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops/genes_pop_542_CTX_EE.txt', header=FALSE)
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

