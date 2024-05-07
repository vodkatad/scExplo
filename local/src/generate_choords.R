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
library(Rtsne)
options(expressions = 5e5)
input<-snakemake@input[['log']]
output<-snakemake@output[['data_out']]
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)

stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
duplicate<-!duplicated(colnames(dato_t))
dato_t<-dato_t[,duplicate]
variance<-colVars(as.matrix(dato_t[sapply(dato_t, is.numeric)]))
keep<-names(variance[variance>0.1])
dato_t<-dato_t[,c(keep)]

choord_generated<-Rtsne(dato_t, is_distance = FALSE,verbose=TRUE,normalize = FALSE)
choord_<-as.data.frame(as.matrix(choord_generated$Y))
colnames(choord_)<-c("xChoord","yChoord")
row.names(choord_)<-row.names(dato_t)
write.table(choord_,file=output,sep=',',quote=FALSE)
