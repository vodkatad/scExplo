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

#solo per stronzo LMX_1
id_cell<-snakemake@input[['cell_id']] 
dato_cell<-read.table(file = id_cell,row.names = 1,sep=",",header = TRUE)
dato_t<-dato_t[row.names(dato_cell),]
print(length(row.names(dato_t)))


stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
cinque<-c("ATOH1","GFI1","DLL1","DEFA5","DEFA6")
wnt<-c('WIF1' ,'APCDD1'  ,'FGF20' ,'WNT6' ,'LEF1', 'BAMBI')
cinque_df<-dato_t[,cinque]
#wnt_df<-dato_t[,wnt]
metagene_cinque<-apply(cinque_df,1,mean)
#metagene_wnt<-apply(wnt_df,1,mean)
cl <- kmeans(cinque_df, 2)
print(cl)
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
sil<-silhouette(posteriors$cluster_id,dist(dato_t))

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
#take chords
chords<-read.table(file = tsne,row.names = 1,sep=",",header = TRUE)
data_merged<-merge(silu,chords,by='row.names')

pdf(plot_out)
ggplot(data_merged, aes(x=xChoord, y=yChoord,color=metagene_cinque)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  title = "Metagene FF")
#ggplot(data_merged, aes(x=xChoord, y=yChoord,color=metagene_wnt)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
 # title = "Metagene WNT")
ggplot(data_merged, aes(x=xChoord, y=yChoord,color=isPaneth)) + geom_point() +labs(
  title = "K-means2_results")
#ggplot(data_merged, aes(x=xChoord, y=yChoord,color=preSilh)) + geom_point()

p <- ggplot(data_merged, aes(x=preSilh, y=sil_width)) + 
  geom_violin()+labs(
  title = "2_comp_silhouette")
p<-p + geom_jitter(shape=16, position=position_jitter(0.2))
p + stat_summary(fun=mean, geom="point", size=8, color="red")+labs(x = "cluster")

d <- ggplot(data_merged, aes(x=isPaneth, y=sil_width)) + 
  geom_violin()+labs(
  title = "2_comp_after silhouette filtering")
d<-d + geom_jitter(shape=16, position=position_jitter(0.2))
d + stat_summary(fun=mean, geom="point", size=8, color="red")+labs(x = "cluster")
graphics.off()
write.table(silu,file=kmeans_out,sep=',',quote=FALSE)
