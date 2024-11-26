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
data_path<-snakemake@input[['data']]
meta_path<-snakemake@input[['meta']]
cluster_path<-snakemake@input[['cluster']]
output_plot<-snakemake@output[['plot_out']]
meta<- read.table(file = meta_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta<-as.data.frame(meta)
data<-read.table(file = data_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
data<-as.data.frame(data)
cluster<-read.table(file = cluster_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
cluster<-as.data.frame(cluster)
data<-merge(data,meta, by='row.names')
row.names(data)<-data$Row.names
data$Row.names<-NULL
data<-merge(data,cluster,by='row.names')
data[data$isPaneth=='filtered','isPaneth']<-'Others'
data[data$isPaneth=='nPaneth','isPaneth']<-'Others'
colori<-c(rainbow(10)[2],rainbow(10)[8])
x_min <- 0
x_max <- ceiling(max(data$x))
y_min <- 0
y_max <- ceiling(max(data$DLL1))
capo<-max(x_max,y_max)
print(max(data$x))


breaks <- pretty(c(0, capo), n = 5) 
re<-max(breaks)


p<-ggplot(data=data,aes(x=x,y=DLL1,color=isPaneth))+geom_point()+scale_color_manual(values = colori)+
  xlab('metagene')+ylab('DLL1')+
  theme_classic()+theme()+
  scale_x_continuous(expand = c(0, 0), limits = c(0, re), breaks = breaks)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, re), breaks = breaks)
#pdf(output_plot)
print(p)
graphics.off()
#median(data$DLL1)
#corr.test(data$DLL1,data$HES1)
