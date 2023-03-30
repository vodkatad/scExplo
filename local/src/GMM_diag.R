library(data.table)
library(mclust)
library(cluster)
input<-snakemake@input[[1]]
output<- snakemake@output[['out']]
plot<- snakemake@output[['plot']]
output_silh<- snakemake@output[['sil']]
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
pdf(plot)
GMM <- densityMclust(data=cinque_df, G=2, modelNames='VVI',plot=TRUE)
graphics.off()
posteriors<-data.frame(GMM$z)
mu<-data.frame(GMM$parameters$mean)
meta_mu<-list()
for(el in mu){
  meta_mu<-append(meta_mu,exp(mean(log(el))))
}
if(meta_mu[[1]][1]>meta_mu[[2]][1]){
  colnames(posteriors)<-list("pPaneth","pNPpaneth")
}else{
  colnames(posteriors)<-list("pNPaneth","pPaneth")
}
posteriors$isPaneth<- with(posteriors, ifelse(pPaneth>0.99,as.integer(1),as.integer(0)))

sil<-silhouette(posteriors$isPaneth,dist(dato_t))
silu<-as.data.frame.matrix(sil)
rownames(silu)<-rownames(posteriors)

silu$preSilh<-posteriors$isPaneth
silu$postSil<-posteriors$isPaneth
silu$postSil[silu$sil_width<0 & silu$preSilh==1]<-0

write.table(posteriors,file=output,sep=',',quote=FALSE)
write.table(silu,file=output_silh,sep=',',quote=FALSE)

