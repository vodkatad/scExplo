library(data.table)
library(ggplot2)
library(pheatmap)
input<-snakemake@input[['data']]
input_s<-snakemake@input[['annot']]

data<-read.table(input,sep=',',header=TRUE,row.names = 1,stringsAsFactors = FALSE)
sign<-read.table(input_s,sep=',',header=TRUE,row.names=1,stringsAsFactors = FALSE)
print(head(sign))
#print(input.str.contains('_C2_'))
if (grepl("_C2_", input)) {
  print('C2')
  data <- data[grepl('REACTOME', rownames(data)), ]
  sign <- sign[grepl('REACTOME', rownames(sign)), ]
} 
print(nrow(data))
print(colnames(data))
ordine<-c("CRC0322_cetux_1","CRC0322_NT_1_3000","CRC0327_cetux_2","CRC0327_NT_2","CRC0542_CTX72h_1","CRC0542_NT72h_1","CRC1502_cetux_1" ,
"CRC1502_NT_1")
data<-data[,ordine]
sing<-sign[,ordine]
#row.names(sign)<-row.names(data)
#colnames(sign)<-colnames(data)
if(nrow(data)==0){
    print('empty')
  pdf(snakemake@output[['plot']])
  graphics.off()
}else{
pdf(snakemake@output[['plot']],width=12,height=12)
min<-floor(min(data,na.rm=TRUE))
max <- ceiling(max(data,na.rm=TRUE))
v<-max(abs(min),abs(max))
minv<- -v
maxv<- v
print(minv)
print(maxv)
rg <- max(abs(data),na.rm=TRUE);
#pdf(snakemake@output[['plot']],width=12,height=12)
neutral_value <- 0
#bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
#bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("blue",
                                            "white"))(n = length(bk1)-1),
                "white", "white",
                c(colorRampPalette(colors = c("white", "red"))(n
                                                                   = length(bk2)-1)))


pheatmap(data,breaks = bk, color=my_palette,fontsize_col =10,fontsize_row=8,cluster_rows=FALSE,cluster_cols=FALSE, show_rownames=TRUE,display_numbers = sign)

#fontsize_row = 10
graphics.off()
}