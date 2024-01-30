library(data.table)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(dplyr)
library(viridis)
library(ggpubr)
library(wesanderson)

input<-snakemake@input[['data']]
kmeans_i<-snakemake@input[['k_out']]
tsne<-snakemake@input[['tsne']]
plot_out<-snakemake@output[['plot_out']]
pdf_<-snakemake@output[['pdf']]
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
kmeans<-read.table(file = kmeans_i,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
metagene_cinque<-apply(dato,1,mean)
metagene_cinque[metagene_cinque>10]<-10
#write.csv(metagene_cinque, meta_ou, row.names=TRUE)

chords<-read.table(file = tsne,row.names = 1,sep=",",header = TRUE)

prova<-merge(chords,kmeans,by='row.names')
prova$isPaneth<-factor(prova$isPaneth,levels = c("nPaneth" ,"filtered","Paneth"))# 

chords<-chords[row.names(dato),]
p<-ggplot(chords, aes(x=xChoord, y=yChoord,color=metagene_cinque)) + 
  geom_point(size=0.4)+scale_color_gradientn(colours = rainbow(10),limits=c(0,10))+#+scale_color_viridis(limits=c(0,10),direction = -1)
labs(title = "", color="metagene")+xlab('tsne1')+ylab('tsne2')+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
colori<-c(rainbow(10)[2],rainbow(10)[5],rainbow(10)[8])#

k<-ggplot(prova, aes(x=xChoord, y=yChoord,color=isPaneth)) + geom_point(size=0.4)+
 labs(title = "", color="")+scale_color_manual(values = colori)+#scale_color_viridis(discrete=TRUE,direction = -1)
xlab('tsne1')+ylab('tsne2')+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())


fine <- ggarrange( p,k, ncol = 2, common.legend = FALSE)
ggsave(plot_out, plot=fine, width=180, height=89, units="mm")
pdf(pdf_)
for (i in colnames(dato)){
  print(i)
  p<-ggplot(chords, aes(x=xChoord, y=yChoord,color=dato[[i]])) + 
  geom_point(size=0.4)+scale_color_gradientn(colours = rainbow(10),limits=c(0,10))+#+scale_color_viridis(limits=c(0,10),direction = -1)
labs(title =i, color="metagene")+xlab('tsne1')+ylab('tsne2')+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
  print(p)
}
graphics.off()