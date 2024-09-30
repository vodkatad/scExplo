library('pheatmap')
library(ggplot2)
#data_path<-'/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/HES1_sorted_spearman/CRC0322_cetux_1_HES1_correlated.tsv'


data_path <- snakemake@input[["data"]]
out <- snakemake@output[["plot"]]
data<-read.table(data_path,stringsAsFactors = FALSE,header = TRUE,sep = '\t')
title <- snakemake@wildcards[["sample"]]


p<-ggplot(data, aes(values)) +
  geom_density()+ggtitle(title)

ggsave(filename=out,p,width = 29.7, height = 21.0, units = "cm")