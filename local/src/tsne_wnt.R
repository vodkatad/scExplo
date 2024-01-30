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
#tsne<-snakemake@input[['tsne']]
#meta_ou<-snakemake@output[['metagene_out']]
data_dll4<-snakemake@output[['data_dll4']]
#meta_wnt_out<-snakemake@output[['wnt_out']]
#plot_out<-snakemake@output[['tsne_out']]
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)

stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
#cinque<-c("ATOH1","GFI1","DLL1","DEFA5","DEFA6")
#wnt<-c('WIF1' ,'APCDD1' ,'LEF1', 'BAMBI')#,'WNT6' 'FGF20' 
#manno_<-c('COLCA2', 'MDK', 'FABP2', 'TFF3', 'ATOH1', 'ANXA13', 'DLL1', 'SMAD9', 'HES6', 'RETNLB', 'LFNG', 'RNASE1', 'TPM1', 'COL4A4', 'AMIGO2', 'HEPACAM2', 'ELAPOR1', 'ODF2L', 'IL13RA1', 'SOX4', 'KLK11', 'FOXA2', 'TGFBI')
manno_<-c('DLL4','NOTCH1','NOTCH2')
#cinque_df<-dato_t[,cinque]
#wnt_df<-dato_t[,wnt]
#metagene_cinque<-apply(cinque_df,1,mean)
#write.csv(metagene_cinque, meta_ou, row.names=TRUE)
#metagene_wnt<-apply(wnt_df,1,mean)
#write.csv(metagene_wnt, meta_wnt_out, row.names=TRUE)
data_manno<-dato_t[,manno_]
write.csv(data_manno, data_dll4, row.names=TRUE)
#metagene_clique<-apply(data_manno,1,mean)

#chords<-read.table(file = tsne,row.names = 1,sep=",",header = TRUE)


#pdf(plot_out)
#ggplot(chords, aes(x=xChoord, y=yChoord,color=metagene_cinque)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
#title = "Metagene FF")
#ggplot(chords, aes(x=xChoord, y=yChoord,color=metagene_wnt)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  #title = "Metagene WNT")

#ggplot(chords, aes(x=xChoord, y=yChoord,color=metagene_clique)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  #title = "Metagene Clique")
  #for (gene in cinque){
     #color<-dato_t[,gene]
    #print(ggplot(chords, aes(x=xChoord, y=yChoord,color=color)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  #title = gene))
  #}
  #for (gene in manno_){
     #color<-dato_t[,gene]
    #print(ggplot(chords, aes(x=xChoord, y=yChoord,color=color)) + geom_point()+scale_color_gradientn(colours = rainbow(5))+labs(
  #title = gene))
  #}



#graphics.off()

