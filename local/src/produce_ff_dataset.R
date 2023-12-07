library(data.table)

input<-snakemake@input[['data']]

out<-snakemake@output[['out']]
#leggo dato
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

write.table(cinque_df,file=out,sep=',',quote=FALSE)
