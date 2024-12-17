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

id_nt_path<-snakemake@input[['id_nt']]
id_cet_path<-snakemake@input[['id_cet']]
cc_path_nt<-snakemake@input[['cc']][1]
print(cc_path_nt)
cc_path_cet<-snakemake@input[['cc']][2]

output_plot<-snakemake@output[['out']]
id_nt<- read.table(file = id_nt_path,sep=",",header = TRUE,stringsAsFactors = FALSE)
id_cet<- read.table(file = id_cet_path,sep=",",header = TRUE,stringsAsFactors = FALSE)
id_nt<-as.data.frame(id_nt)
id_cet<-as.data.frame(id_cet)
id_nt$cell_id <- gsub("-", ".", id_nt$cell_id)
id_cet$cell_id <- gsub("-", ".", id_cet$cell_id)
row.names(id_nt)<-id_nt$cell_id
row.names(id_cet)<-id_cet$cell_id
id_nt$cell_id<-NULL
id_cet$cell_id<-NULL

cc_data_nt<-read.table(file = cc_path_nt,row.names = 1,sep=",",header = FALSE,stringsAsFactors = FALSE)
cc_data_cet<-read.table(file = cc_path_cet,row.names = 1,sep=",",header = FALSE,stringsAsFactors = FALSE)
cc_data_nt<-as.data.frame(cc_data_nt)
cc_data_cet<-as.data.frame(cc_data_cet)
colnames(cc_data_nt)<-c('cc')
colnames(cc_data_cet)<-c('cc')

data_nt<-merge(id_nt,cc_data_nt, by='row.names')
data_cet<-merge(id_cet,cc_data_cet, by='row.names')
row.names(data_nt)<-data_nt$Row.name
row.names(data_cet)<-data_cet$Row.name
data_nt$Row.name<-NULL
data_cet$Row.name<-NULL

data_cet$trattamento<-'cet'
data_nt$trattamento<-'nt'
data<-rbind(data_nt,data_cet)
data$Row.names<-NULL
print(head(data))

#sono arrivata qua 




data$cc<-as.factor(data$cc)
data$trattamento<-as.factor(data$trattamento)
pdf(output_plot)

y_min <- 0
y_max <- ceiling(max(data$id))
#capo<-max(x_max,y_max)

y_breaks <- pretty(c(0, y_max), n = 5) 
y_max<-max(y_breaks)


p<-ggplot(data=data,aes(x=trattamento,y=id))+geom_boxplot()+geom_jitter(alpha=0.5,height = 0,aes(color=cc))+
xlab('trattamento')+ylab('ID')+
theme_classic()+theme(axis.ticks.x = element_blank() )+
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max), breaks = y_breaks)

print(p)


graphics.off()

