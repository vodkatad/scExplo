library(data.table)
library(dplyr)


input<-snakemake@input[['data']]
variance_path<-snakemake@input[['variance']]

out<-snakemake@output[['out']]
#leggo dato
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)
variance<- read.table(file = variance_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
geni<-rownames(variance)
stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
#cinque<-c("ATOH1","GFI1","DLL1","DEFA5","DEFA6")


subset_df <- dato_t %>% select(all_of(intersect(geni, colnames(dato_t))))

write.table(subset_df,file=out,sep=',',quote=FALSE)
