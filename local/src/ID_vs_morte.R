library(data.table)
library(mclust)
library(cluster)
library(reshape2)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(resample)
library(dplyr)
library(viridis)
library(igraph)
library(patchwork)
library(diptest)
#data_path<-'/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/ff_data/CRC0322_NT_1_3000_log2cpm.csv'
#meta_path<-'/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/kmeans/metageni/CRC0322_NT_1_3000_metagene.csv'
#cluster_path<-'/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/kmeans/CRC0322_NT_1_3000/CRC0322_NT_1_3000_kmeans_2comp.csv'

id_path<-snakemake@input[['id']]
cc_path<-snakemake@input[['cc']]

output_plot<-snakemake@output[['plot_out']]
id<- read.table(file = id_path,sep=",",header = TRUE,stringsAsFactors = FALSE)
id<-as.data.frame(id)
id$cell_id <- gsub("-", ".", id$cell_id)
row.names(id)<-id$cell_id
id$cell_id<-NULL
#sono arrivata qua 

cc_data<-read.table(file = cc_path,row.names = 1,sep=",",header = FALSE,stringsAsFactors = FALSE)
cc_data<-as.data.frame(cc_data)
colnames(cc_data)<-c('cc')
cc_data$cell_id<-rownames(cc_data)
id$cell_id<-rownames(id)

id <- id %>%
  left_join(cc_data, by = "cell_id") %>%
  mutate(cc = ifelse(is.na(cc), "filtered", cc))

rownames(id)<-id$cell_id
id$cell_id<-NULL
id$cc<-as.factor(id$cc)
pdf(output_plot)

y_min <- 0
y_max <- ceiling(max(id$id))
#capo<-max(x_max,y_max)

y_breaks <- pretty(c(0, y_max), n = 5) 
y_max<-max(y_breaks)


p<-ggplot(data=id,aes(x=cc,y=id,color=cc))+geom_boxplot()+geom_jitter(alpha=0.5,height = 0)+
xlab('cc')+ylab('ID')+
theme_classic()+theme(axis.ticks.x = element_blank() )+
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max), breaks = y_breaks)

print(p)


graphics.off()

