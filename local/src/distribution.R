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
library(ggpubr)

#output_plot<-snakemake@output[[1]]
clust<-snakemake@input$c[1]
print(clust)
metagene<-snakemake@input$metagene[1]
output_plot<-snakemake@output$out[1]
kmeans<- read.table(file = clust,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta<- read.table(file = metagene,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)

a<-ggplot(data=meta, aes(x=x)) +geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("Metagene Distribution")
print(a)
 

#dato<-rbind(dato_cet,dato_nt)
kmeans$isPaneth<-as.character(kmeans$isPaneth)
kmeans$isPaneth[kmeans$isPaneth=='filtered']<-'nPaneth'
kmeans$isPaneth<-as.factor(kmeans$isPaneth)
kmeans<-merge(kmeans,meta,by='row.names')

b<-ggplot(data=kmeans, aes(x=x,fill=isPaneth)) +geom_density(alpha=0.8) +
  ggtitle("Metagene Distribution")
print(b)

library(plyr)
mu <- ddply(kmeans, "isPaneth", summarise, grp.mean=mean(x))
b<-ggplot(data=kmeans, aes(x=x,fill=isPaneth,color=isPaneth)) +geom_density(alpha=0.5) +
  ggtitle("Metagene Distribution ")+theme_bw()
print(b)

g<-ggarrange(a,b, ncol = 2, nrow = 1)
ggsave(file=output_plot, plot=g, width=16, height=8)
graphics.off()