library(data.table)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(dplyr)
library(viridis)
library(ggpubr)
library(wesanderson)


input<-snakemake@input[['count']]
tsne<-snakemake@input[['umap']]
#plot_out<-snakemake@output[['plot_out']]
pdf_<-snakemake@output[['pdf']]
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(input)

data<- read.table(file = input,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
righe<-colnames(data)
colonne<-rownames(data)
data[is.na(data)]<-0
data<-data+1
cpm <-apply(data,2, function(x) (x/sum(x))*1000000)
log_cpm <- log2(cpm)
dato<- transpose(as.data.frame(log_cpm))
rownames(dato) <- righe
colnames(dato)<-colonne


stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
#colnames(dato) <- sapply(colnames(dato),FUN=stringa_split)
#print(ncol(dato))
#print(colnames(dato))
geni<-c('ENSG00000154856','ENSG00000142541','ENSG00000172238')#ENSG00000142541-RPL.., ATOH1:ENSG00000172238 , APCDD1:ENSG00000154856
dato<-dato[,geni]

colnames(dato)<-c('APCDD1','RPL13A','ATOH1')

#print(head(log_cpm))
#row.names(cpm)<-rownames(dato)
#colnames(cpm)<-colnames(dato)



#metagene_cinque<-apply(dato[row.names(kmeans),],1,mean)
#metagene_cinque[metagene_cinque>10]<-10
#dato[dato>10]<-10
#write.csv(metagene_cinque, meta_ou, row.names=TRUE)

chords<-read.table(file = tsne,row.names = 1,sep=",",header = TRUE)

dato<-merge(dato,chords,by='row.names')
rownames(dato)<-dato$Row.names
dato$Row.names<-NULL
print(head(dato))


#chords<-merge(chords,metagene_cinque,by='row.names')
#print(head(chords))
#p<-ggplot(chords, aes(x=X_umap1, y=X_umap2,color=y)) + 
  #geom_point(size=0.4)+scale_colour_gradientn(colours = rainbow(10),limits=c(0,10),guide = guide_colourbar())+#+scale_color_viridis(limits=c(0,10),direction = -1)
#labs(title = "", color="metagene")+xlab('umap1')+ylab('umap2')+guides(color = guide_colourbar(barwidth = 0.5, barheight = 6))+
#theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.title = element_blank())
#colori<-c(rainbow(10)[2],rainbow(10)[8])#,rainbow(10)[5]

#k<-ggplot(prova, aes(x=X_umap1, y=X_umap2,color=isPaneth)) + geom_point(size=0.4)+
 #labs(title = "", color="")+scale_color_manual(values = colori)+#scale_color_viridis(discrete=TRUE,direction = -1)
#xlab('umap1')+ylab('umap2')+guides(fill=guide_legend(title=" "))+guides(color = guide_legend(override.aes = list(size = 2)))+
#theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.title = element_blank())


#fine <- ggarrange( p,k, ncol = 2, common.legend = FALSE)
#ggsave(plot_out, plot=fine, width=180, height=89, units="mm")

pdf(pdf_)
for (i in c('APCDD1','RPL13A','ATOH1')){
  print(i)
  p<-ggplot(dato, aes(x=X_umap1, y=X_umap2,color=dato[,i])) + 
  geom_point(size=1)+scale_colour_gradientn(colours = rainbow(10))+#+scale_color_viridis(limits=c(0,10),direction = -1)
labs(title =i, color=i)+xlab('umap1')+ylab('umap2')+guides(color = guide_colourbar(barwidth = 0.5, barheight = 6))+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
  print(p)
}
graphics.off()









