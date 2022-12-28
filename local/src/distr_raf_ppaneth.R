library(ggplot2)


setwd('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/')

data <- read.table('CRC0322LMX_NT1w_1_singleron0722_dir/filtered_annotated_saver_ribomito_CRC0322LMX_NT1w_1_singleron0722_log2_pc1_cpm.csv', sep=",", header=TRUE)

paneth <- read.table('/mnt/cold1/calogero/reanalysis_on_AIsc/paneth_only/all/log2cpmfisher/singpdxnt_log2cpmfisher.csv', sep=",", header=TRUE, stringsAsFactors = FALSE)



paneth$paneth <- ifelse(paneth$V1 < 0.001, 'paneth', 'nopaneth')
paneth$panethName <- sapply(paneth$X, function(x) {strsplit(x, '_')[[1]][1]})

rownames(data) <- data$X
data$X <- NULL
tdata <- t(data)
tdata  <- as.data.frame(tdata)
tdata$paneth <- ifelse(rownames(tdata) %in% paneth[paneth$paneth == "paneth", 'panethName'], 'paneth', 'nonpaneth')

w <- 'DEFA6'
i <- grep(w, colnames(tdata))
colnames(tdata)[i] <- w
ggplot(data=tdata, aes_string(x=w, color='paneth'))+geom_density()

for (w in c('DEFA5', 'ATOH1', 'DLL1', 'GFI1')) {
  i <- grep(w, colnames(tdata), fixed=TRUE)
  print(i)
  colnames(tdata)[i] <- w
  print(ggplot(data=tdata, aes_string(x=w, color='paneth'))+geom_density())
}


setwd('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/')

data <- read.table('CRC0322LMO_CTX72h_1_singleron0722_dir/filtered_annotated_saver_ribomito_CRC0322LMO_CTX72h_1_singleron0722_log2_pc1_cpm.csv', sep=",", header=TRUE)

paneth <- read.table('/mnt/cold1/calogero/reanalysis_on_AIsc/paneth_only/all/log2cpmfisher/singctx_log2cpmfisher.csv', sep=",", header=TRUE, stringsAsFactors = FALSE)



paneth$paneth <- ifelse(paneth$V1 < 0.001, 'paneth', 'nopaneth')
paneth$panethName <- sapply(paneth$X, function(x) {strsplit(x, '_')[[1]][1]})

rownames(data) <- data$X
data$X <- NULL
tdata <- t(data)
tdata  <- as.data.frame(tdata)
tdata$paneth <- ifelse(rownames(tdata) %in% paneth[paneth$paneth == "paneth", 'panethName'], 'paneth', 'nonpaneth')

w <- 'DEFA6'
i <- grep(w, colnames(tdata))
colnames(tdata)[i] <- w
ggplot(data=tdata, aes_string(x=w, color='paneth'))+geom_density()

for (w in c('DEFA5', 'ATOH1', 'DLL1', 'GFI1')) {
  i <- grep(w, colnames(tdata), fixed=TRUE)
  print(i)
  colnames(tdata)[i] <- w
  print(ggplot(data=tdata, aes_string(x=w, color='paneth'))+geom_density())
}

############
base <- '/data/reanalysis_on_AIsc/paneth_only/all/log2cpmfisher/'
files <- c("crc322ctx1_log2cpmfisher.csv",
           "crc322ctx_log2cpmfisher.csv",
           "crc322nt1_log2cpmfisher.csv",
           "crc322nt_log2cpmfisher.csv",
           "crc327ctx1_log2cpmfisher.csv",
           "crc327ctx2_log2cpmfisher.csv",
           "crc327nt1_log2cpmfisher.csv",
           "crc327nt2_log2cpmfisher.csv",
           "crc542ctx_log2cpmfisher.csv",
           "crc542nt_log2cpmfisher.csv",
           "crc69ctx_log2cpmfisher.csv",
           "crc69nt_log2cpmfisher.csv",
           "singctx_log2cpmfisher.csv",
           "singnt_log2cpmfisher.csv",
           "singpdxctx_log2cpmfisher.csv",
           "singpdxnt_log2cpmfisher.csv")
data <- as.list(paste0(base, files))

names(data) <- c('CRC0322_cetux_1',
                 'CRC0322_3d_CTX',
                 'CRC0322_NT_1_3000',
                 'CRC0322_8d_NOEGF',
                 "CRC0327_cetux_1_4000",
                 "CRC0327_cetux_2",
                 "CRC0327_NT_1",
                 "CRC0327_NT_2",
                 "CRC0542_CTX72h_1",
                 "CRC0542_NT72h_1",
                 "CRC0069_CTX72h_1",
                 "CRC0069_NT72h_1",
                 'CRC0322LMO_CTX72h_1_singleron0722', 
                 'CRC0322LMO_NT72h_1_singleron0722', 
                 'CRC0322LMX_CTX1w_1_singleron0722', 
                 'CRC0322LMX_NT1w_1_singleron072')

plot_p <- function(name, list) {
  pvals <- read.table(list[[name]], sep=",", header=TRUE, stringsAsFactors = FALSE)
  colnames(pvals) <- c('cell_phase', 'np', 'p')
  pvals$padj_np <- p.adjust(pvals$np, method="BH")
  pvals$padj_p <- p.adjust(pvals$p, method="BH")
  
  print(name)
  print(nrow(pvals))
  print(table(pvals$np < 0.01))
  print(table(pvals$padj_np < 0.01))
  print(table(pvals$padj_np < 0.01 & pvals$padj_p > 0.05))
  print(cor.test(log10(pvals$padj_np+1e-100), log10(pvals$padj_p+1e-100)))
  
  #print(plot(log10(pvals$np+1e-100), log10(pvals$p + 1e-100)))
  pvals$isPaneth_np <- ifelse(pvals$padj_np < 0.01, 'yes', 'no')

  print(ggplot(data=pvals, aes(x=padj_p, fill=isPaneth_np))+geom_histogram(position='dodge')+
    scale_x_log10()+theme_bw()+theme(text=element_text(size=15))+geom_vline(xintercept=0.05)+
    ggtitle(name))
}

gg <- sapply(names(data), plot_p, data)

name <- 'CRC0322LMX_CTX1w_1_singleron0722'

pvals <- read.table(data[[name]], sep=",", header=TRUE, stringsAsFactors = FALSE)
colnames(pvals) <- c('cell_phase', 'np', 'p')
pvals$padj_np <- p.adjust(pvals$np, method="BH")
pvals$padj_p <- p.adjust(pvals$p, method="BH")


ggplot(data=pvals, aes(x=padj_p))+geom_histogram(position='dodge')+
        scale_x_log10()+theme_bw()+theme(text=element_text(size=15))

pvals[,'padj_np'] <- pvals[,'padj_np']+1e-100
ggplot(data=pvals, aes(x=padj_np))+geom_histogram(position='dodge')+
  scale_x_log10()+theme_bw()+theme(text=element_text(size=15))