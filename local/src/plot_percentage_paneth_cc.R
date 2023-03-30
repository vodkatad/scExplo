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
print(snakemake@input)
cetux<-snakemake@input$clust[1]
nt<-snakemake@input$clust[2]
print('&&&&&&&&&&&&&&&&&&&&&&')
print(cetux)
print(class(cetux))
dato_cet<- read.table(file = cetux,row.names = 1,sep=",",header = TRUE)
dato_nt<- read.table(file = nt,row.names = 1,sep=",",header = TRUE)
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dato_cet$sample_name<- strsplit(cetux , split = "/")[[1]][10]
output_plot<-snakemake@output[[1]]
dato_nt$sample_name<-strsplit(nt , split = "/")[[1]][10]
dato<-rbind(dato_cet,dato_nt)
dato$isPaneth<-as.character(dato$isPaneth)
dato$isPaneth[dato$isPaneth=='filtered']<-'nPaneth'
dato$isPaneth<-as.factor(dato$isPaneth)
agg_df <- aggregate(dato$isPaneth, by=list(dato$sample_name,dato$isPaneth), FUN=length)
agg_df<-agg_df %>%
  group_by(Group.1) %>%
  mutate(percent = x/sum(x)*100)
#count paneth
agg_df<-as.data.frame(agg_df)
agg_df$percent<-round(agg_df$percent,2)

cc_cetux<-snakemake@input$cc[1]
dato_cc_cet<- read.table(file = cc_cetux,row.names = 1,sep=",",header = FALSE)
dato_cc_cet<-merge(dato_cet,dato_cc_cet,by='row.names')

cc_nt<-snakemake@input$cc[2]
dato_cc_nt<- read.table(file = cc_nt,row.names = 1,sep=",",header = FALSE)
dato_cc_nt<-merge(dato_nt,dato_cc_nt,by='row.names')

dato_cc<-rbind(dato_cc_nt,dato_cc_cet)
dato_cc$isPaneth<-as.character(dato_cc$isPaneth)
# dato_cc$isPaneth<-dato_cc$postSilh
dato_cc$isPaneth[dato_cc$isPaneth=='filtered']<-'nPaneth'
dato_cc$isPaneth<-as.factor(dato_cc$isPaneth)


agg_data<-aggregate(dato_cc$isPaneth, by=list(dato_cc$sample_name,dato_cc$isPaneth,dato_cc$V2), FUN=length)
agg_data<-agg_data %>%
  group_by(Group.1) %>%
  mutate(percent = x/sum(x)*100)

agg_data<-as.data.frame(agg_data)
agg_data$percent<-round(agg_data$percent,2)
pdf(output_plot)
ggplot(data=agg_df, aes(x = Group.1,y=percent,fill=Group.2)) +geom_bar(stat="identity")+geom_text(aes(label=percent), vjust=1.6,color="black", size=3.5)+ xlab("Sample") + ylab("%") +labs(fill='')

ggplot(data=agg_data, aes(x = Group.3,y=percent,fill=Group.2)) +geom_bar(stat="identity")+
  geom_text(aes(label=percent), vjust=1.6,color="black", size=3.5) +
  facet_wrap(~Group.1) +xlab("Cell Cycle") + ylab("%")  +labs(fill='')         
graphics.off()