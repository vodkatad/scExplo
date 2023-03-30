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
input<-snakemake@input[['data']]

kmeans_out<-snakemake@output[['out']]
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)

stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
cinque<-c("ATOH1","GFI1","DLL1","DEFA5","DEFA6")

cinque_df<-dato_t[,cinque]

cl <- kmeans(cinque_df, 2)

meta_mu<-apply(cl$centers,1,geometric.mean)
ordine<-order(unlist(meta_mu))
ordinate<-c('nPaneth','Paneth')
#ordinate

final_ordinate<-c()
for (i in seq(1,2)){
  final_ordinate[ordine[i]] <-ordinate[i] 
}
cluster_id<-cl$cluster

isPaneth<-c()
for (el in cluster_id){
  isPaneth<-append(isPaneth,final_ordinate[el])
}
posteriors<-cbind(cluster_id,isPaneth)


posteriors<-data.frame(posteriors)
posteriors$cluster_id <- as.numeric(posteriors$cluster_id) 
library (vegan)
sil<-silhouette(posteriors$cluster_id,dist(cinque_df))

silu<-as.data.frame.matrix(sil)
rownames(silu)<-rownames(posteriors)

silu$preSilh<-posteriors$isPaneth

silu$postSilh<-posteriors$isPaneth
silu$preSilh<-as.character(silu$preSilh)
silu$postSilh<-as.character(silu$postSilh)
silu$postSilh[silu$sil_width<0 & silu$preSilh=='Paneth']<-'filtered'

silu$isPaneth<-silu$postSilh
silu$postSilh[silu$postSilh=='filtered']<-0
silu$postSilh[silu$postSilh=='nPaneth']<-0
silu$postSilh[silu$postSilh=='Paneth']<-1


silu$isPaneth<-as.factor(silu$isPaneth)

write.table(silu,file=kmeans_out,sep=',',quote=FALSE)
