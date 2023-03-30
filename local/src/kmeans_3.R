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
tsne<-snakemake@input[['tsne']]
kmeans_out<-snakemake@output[['out']]
plot_out<-snakemake@output[['tsne_out']]
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
metagene_cinque<-apply(cinque_df,1,geometric.mean)
cl <- kmeans(cinque_df, 3)

meta_mu<-apply(cl$centers,1,geometric.mean)
ordine<-order(unlist(meta_mu))
ordinate<-c('Low','Medium','High')
#ordinate
final_ordinate<-c()
for (i in seq(1,3)){
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
sil<-silhouette(posteriors$cluster_id,dist(dato_t))

silu<-as.data.frame.matrix(sil)
rownames(silu)<-rownames(posteriors)

silu$preSilh<-posteriors$isPaneth
silu$postSilh<-posteriors$isPaneth
silu$preSilh<-as.character(silu$preSilh)
silu$postSilh<-as.character(silu$postSilh)
silu$postSilh[silu$sil_width<0 & silu$preSilh=='High']<-'filtered_High'
silu$postSilh[silu$sil_width<0 & silu$preSilh=='Medium']<-'filtered_Medium'
silu$postSilh<-as.factor(silu$postSilh)
#take chords
chords<-read.table(file = tsne,row.names = 1,sep=",",header = TRUE)
data_merged<-merge(silu,chords,by='row.names')

pdf(plot_out)
ggplot(data_merged, aes(x=xChoord, y=yChoord,color=metagene_cinque)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  title = "Metagene FF")
ggplot(data_merged, aes(x=xChoord, y=yChoord,color=postSilh)) + geom_point() +labs(
  title = "K-means3_results")
#ggplot(data_merged, aes(x=xChoord, y=yChoord,color=preSilh)) + geom_point()

p <- ggplot(data_merged, aes(x=preSilh, y=silhouetteValue)) + 
  geom_violin()+labs(
  title = "3_comp_silhouette")
p<-p + geom_jitter(shape=16, position=position_jitter(0.2))
p + stat_summary(fun=mean, geom="point", size=8, color="red")+labs(x = "cluster")

d <- ggplot(data_merged, aes(x=postSilh, y=silhouetteValue)) + 
  geom_violin()+labs(
  title = "3_comp_after silhouette filtering")
d<-d + geom_jitter(shape=16, position=position_jitter(0.2))
d + stat_summary(fun=mean, geom="point", size=8, color="red")+labs(x = "cluster")
graphics.off()
write.table(posteriors,file=kmeans_out,sep=',',quote=FALSE)
