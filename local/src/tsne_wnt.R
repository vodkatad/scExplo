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
wnt<-c('WIF1' ,'APCDD1'  ,'FGF20' ,'WNT6' ,'LEF1', 'BAMBI')

cinque_df<-dato_t[,cinque]
wnt_df<-dato_t[,wnt]
metagene_cinque<-apply(cinque_df,1,mean)
metagene_wnt<-apply(wnt_df,1,mean)

chords<-read.table(file = tsne,row.names = 1,sep=",",header = TRUE)


pdf(plot_out)
ggplot(chords, aes(x=xChoord, y=yChoord,color=metagene_cinque)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  title = "Metagene FF")
ggplot(chords, aes(x=xChoord, y=yChoord,color=metagene_wnt)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  title = "Metagene WNT")
  for (gene in cinque){
     color<-dato_t[,gene]
    print(ggplot(chords, aes(x=xChoord, y=yChoord,color=color)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  title = gene))
  }



graphics.off()

