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
id_path<-snakemake@input[['id']]
cc_path<-snakemake@input[['cc']]

output_plot<-snakemake@output[['plot_out']]
meta<- read.table(file = id_path,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta<-as.data.frame(meta)
meta$cell_id <- gsub("-", ".", meta$cell_id)
row.names(meta)<-meta$cell_id
meta$cell_id<-NULL

data<-read.table(file = data_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
data<-as.data.frame(data)

cc_data<-read.table(file = cc_path,row.names = 1,sep=",",header = FALSE,stringsAsFactors = FALSE)
cc_data<-as.data.frame(cc_data)
colnames(cc_data)<-c('cc')

data<-merge(data,meta, by='row.names')
row.names(data)<-data$Row.names
data$Row.names<-NULL

data<-merge(data,cc_data, by='row.names')
row.names(data)<-data$Row.names
data$Row.names<-NULL
print(head(cc_data))
data$cc<-as.factor(data$cc)
#colori<-c(rainbow(10)[2],rainbow(10)[8])
pdf(output_plot)
for(gene in colnames(data)){
  if (gene!='id' & gene!='cc'){
    print(gene)
    x_min <- 0
    x_max <- ceiling(max(data$id))
    y_min <- 0
    y_max <- ceiling(max(data[[gene]]))
    #capo<-max(x_max,y_max)
    x_breaks <- pretty(c(0, x_max), n = 5) 
    y_breaks <- pretty(c(0, y_max), n = 5) 
    print(x_breaks)
    x_max<-max(x_breaks)
    y_max<-max(y_breaks)
    
    
    p<-ggplot(data=data,aes(x=id,y=data[[gene]],color=cc,shape = cc))+geom_point(alpha=0.5)+
    xlab('ID')+ylab(gene)+
    theme_classic()+theme()+
    scale_x_continuous(expand = c(0, 0), limits = c(0, x_max), breaks = x_breaks)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_max), breaks = y_breaks)

    print(p)}

}

graphics.off()

