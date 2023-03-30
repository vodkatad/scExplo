library(data.table)
library(mclust)
library(cluster)
library(dplyr)
library(ggplot2)
input<-snakemake@input[['tsne']]
silh<-snakemake@input[['sil']]
output<- snakemake@output[['res']]

silu<-read.table(file = silh,row.names = 1,sep=",",header = TRUE)
silu$paneth<-silu$postSil

silu$paneth[silu$postSil==1]<-'paneth'
silu$paneth[silu$postSil==0]<-'nPaneth'
silu$paneth[silu$preSilh==1 & silu$postSil==0]<-'filtered'
#silu$colori <- as.factor(silu$colori)
silu$preSilh <- as.factor(silu$preSilh)
silu$postSil <- as.factor(silu$postSil)
chords<- read.table(file = input,row.names = 1,sep=",",header = TRUE)

data_merged<-merge(silu,chords,by='row.names')

pdf(output)
ggplot(data_merged, aes(x=xChoord, y=yChoord,color=postSil)) + geom_point()
ggplot(data_merged, aes(x=xChoord, y=yChoord,color=preSilh)) + geom_point()
ggplot(data_merged, aes(x=xChoord, y=yChoord,color=paneth)) + geom_point()
graphics.off()