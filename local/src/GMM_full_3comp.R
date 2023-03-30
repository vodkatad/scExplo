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
summary<-snakemake@output[['plot']]
posterior_out<-snakemake@output[['out']]
silh_out<-snakemake@output[['sil']]
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
pdf(summary)
set.seed(42)

GMM <- densityMclust(data=cinque_df, G=3, modelNames='VVV',plot=TRUE)
graphics.off()
posteriors<-GMM$z
mu<-data.frame(GMM$parameters$mean)
mu
meta_mu<-list()
for(el in mu){
  meta_mu<-append(meta_mu,mean(el))
}
meta_mu
ordine<-order(unlist(meta_mu))
ordine
ordinate<-c('Low','Medium','High')
#ordinate
final_ordinate<-c()
for (i in seq(1,3)){
  final_ordinate[ordine[i]] <-ordinate[i] 
}

print(meta_mu)
colnames(posteriors)<-as.list(final_ordinate)
max_id<-argmax(posteriors,rows = TRUE)
isPaneth<-c()

for (el in max_id){
  isPaneth<-append(isPaneth,final_ordinate[el])
}
posteriors<-cbind(posteriors,isPaneth)
posteriors<-cbind(posteriors,max_id)

posteriors<-data.frame(posteriors)
posteriors$max_id <- as.numeric(posteriors$max_id) 
library (vegan)
sil<-silhouette(posteriors$max_id,dist(dato_t))

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
  title = "GMM_3_comp_results")
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
write.table(posteriors,file=posterior_out,sep=',',quote=FALSE)
write.table(silu,file=silh_out,sep=',',quote=FALSE)