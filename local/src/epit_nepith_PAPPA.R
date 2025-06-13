library(reshape)
library(ggplot2)

#setwd('/mnt/cold1/snaketree/prj/scRNA/dataset/Scanpy_prova/count_matrix/GSE144735/')
setwd('/mnt/cold1/snaketree/prj/scRNA/dataset/Scanpy_prova/count_matrix/GSE132465/')

wanted <- 'PAPPA'
epi <- list.files(path = ".", pattern = '*epithelial*')
nepi <- list.files(path = ".", pattern = '*notEpithelial*')
thr <- 1

length(epi)==length(nepi)
data <- data.frame()
datap <- data.frame()
for (i in seq(1, length(epi))) {
  #sid <- stringr::str_extract(string = epi[i], pattern = "KUL\\d+-.")
  sid <- stringr::str_extract(string = epi[i], pattern = "SMC\\d+-.")
  expr_e <- read.csv(epi[i], row.names=1)
  expr_ne <- read.csv(nepi[i], row.names=1)
  data_e <- data.frame(tissue=rep("epithelial", ncol(expr_e)), sample=rep(sid, ncol(expr_e)), expr= unlist(expr_e[rownames(expr_e)==wanted,]))
  data_ne <- data.frame(tissue=rep("notEpithelial", ncol(expr_ne)), sample=rep(sid, ncol(expr_ne)), expr= unlist(expr_ne[rownames(expr_ne)==wanted,]))
  data <- rbind(data, data_e)
  data <- rbind(data, data_ne)
  perc_e <- sum(data_e$expr >= thr) / nrow(data_e)
  perc_ne <- sum(data_ne$expr >= thr) / nrow(data_ne)
  datap <- rbind(datap, c(sid, perc_e, perc_ne))
}

data$x <- paste0(data$sample, '_', data$tissue)
data <- data[order(data$x),]
data$lexpr <- log2(data$expr+1)

lab <- c()
labs <- c(unique(data$sample), unique(data$sample))
labs <- labs[order(labs)]
sk <- T
for (i in seq(1, length(labs))) {
  if (sk) {
    lab <- c(lab, labs[i])
    sk <- F
  } else {
    lab <- c(lab, ' ')
    sk <- T
  }
}

ggplot(data=data, aes(x=x, y=lexpr))+geom_violin(aes(color=tissue))+
  ylab('log2_1p UMI')+theme_bw(base_size=15)+
  scale_x_discrete(labels=lab)+ theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

colnames(datap) <- c('sample', 'frac_epithelial', 'frac_nepithelial')

long <- melt(datap, id="sample")
long$value <- as.numeric(long$value)
ggplot(data=long, aes(x=sample, y=value, fill=variable))+geom_col(position="dodge")+
  ylab('frac of cells PAPPA >= 3 UMI')+theme_bw(base_size=15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
