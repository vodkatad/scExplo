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
umap<-snakemake@input[['umap']]

kmeans_out<-snakemake@output[['out']]
plot<-snakemake@output[['out_plot']]
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
#dato_t<- transpose(dato)
#rownames(dato_t) <- colnames(dato)
#colnames(dato_t)<-rownames(dato)

#stringa_split<-function(stringa){
  #res<-strsplit(stringa , split = ":")[[1]][2]
  #return(res)
#}
#colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
cinque<-c("HES1")

cinque_df<-dato[,cinque]



cl <- kmeans(cinque_df, 2)

meta_mu<-apply(cl$centers,1,geometric.mean)
ordine<-order(unlist(meta_mu))
ordinate<-c('HES_neg','HES_pos')
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
print(rownames(head(cinque_df)))
rownames(silu)<-rownames(dato)

silu$preSilh<-posteriors$isPaneth

silu$postSilh<-posteriors$isPaneth
silu$preSilh<-as.character(silu$preSilh)
silu$postSilh<-as.character(silu$postSilh)
silu$postSilh[silu$sil_width<0 & silu$preSilh=='HES_pos']<-'filtered'

silu$isPaneth<-silu$postSilh
silu$postSilh[silu$postSilh=='filtered']<-0
silu$postSilh[silu$postSilh=='HES_neg']<-0
silu$postSilh[silu$postSilh=='HES_pos']<-1


silu$isPaneth<-as.factor(silu$isPaneth)
print('silu')
print(head(silu))
chords<-read.table(file = umap,row.names = 1,sep=",",header = TRUE)
data_merged<-merge(silu,chords,by='row.names')


write.table(silu,file=kmeans_out,sep=',',quote=FALSE)

chords<-read.table(file = umap,row.names = 1,sep=",",header = TRUE)
print('cose')
print(rownames(head(chords)))
print(rownames(head(silu)))
data_merged<-merge(silu,chords,by='row.names')
pdf(plot)
ggplot(data_merged, aes(x=X_umap1, y=X_umap2,color=isPaneth)) + geom_point() +labs(
  title = "K-means2_results")
graphics.off()

