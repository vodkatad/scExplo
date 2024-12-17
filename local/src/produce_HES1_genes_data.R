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
#cinque<-c()
cinque<-c("ATOH1","GFI1","DLL1","DEFA5","DEFA6","OLFM4","NOTCH1","NOTCH2","HES1","Z83851.1",
"MCOLN2","LRTOMT","LINC00501","GNA14","FCMR","EPHB3","DSG4","DGUOK-AS1","C5orf52","ACTL8","AC007249.2","ABCB11","VDR",       
"TMC8","SFXN2","SAMHD1","PARP12","MREG","LRRC37A3","LINC02441" ,"HES4", "GHR"  ,"FCAMR","ESAM","ERG" ,"DEPTOR" ,"CPM" ,"CNN3",  
"C15orf62","BRCC3","B4GALNT3","ASL","ALG1","AL645608.8","AL355987.4","AL139246.5","AC105460.1")

cinque<-intersect(colnames(dato_t),cinque)

cinque_df<-dato_t[,cinque]

write.table(cinque_df,file=out,sep=',',quote=FALSE)
