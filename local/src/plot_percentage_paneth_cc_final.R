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
library("ggpubr")


cetux<-snakemake@input$clust[1]
nt<-snakemake@input$clust[2]
output_plot<-snakemake@output[[1]]
cc_cetux<-snakemake@input$cc[1]
cc_nt<-snakemake@input$cc[2]

dato_cet<- read.table(file = cetux,row.names = 1,sep=",",header = TRUE)
dato_nt<- read.table(file = nt,row.names = 1,sep=",",header = TRUE)
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dato_cet$sample_name<- strsplit(cetux , split = "/")[[1]][10]
dato_nt$sample_name<- strsplit(cetux , split = "/")[[1]][10]
sample_name<-strsplit(cetux , split = "/")[[1]][10]


#dato<-rbind(dato_cet,dato_nt)
dato_cet$isPaneth<-as.character(dato_cet$isPaneth)
dato_cet$isPaneth[dato_cet$isPaneth=='filtered']<-'nPaneth'
dato_cet$isPaneth<-as.factor(dato_cet$isPaneth)

dato_nt$isPaneth<-as.character(dato_nt$isPaneth)
dato_nt$isPaneth[dato_nt$isPaneth=='filtered']<-'nPaneth'
dato_nt$isPaneth<-as.factor(dato_nt$isPaneth)


dato_cc_cet<- read.table(file = cc_cetux,row.names = 1,sep=",",header = FALSE)
dato_cc_nt<- read.table(file = cc_nt,row.names = 1,sep=",",header = FALSE)
dato_cet<-merge(dato_cet,dato_cc_cet,by='row.names')
dato_nt<-merge(dato_nt,dato_cc_nt,by='row.names')

dato_cet.summary <- dato_cet %>% group_by(isPaneth,V2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  #group_by(V2) %>% 
  mutate(percent =total_count/sum(total_count))
         #pos =0.8-( cumsum(percent) - 0.5*percent))

cet<-ggplot(dato_cet.summary,  aes(x=isPaneth, y=percent, fill=V2)) +
  geom_bar(stat='identity',  width = .7, colour="black", lwd=0.1) +
  geom_text(aes(label=ifelse(percent >= 0.02, paste0(sprintf("%.0f", percent*100),"%"),"")), position=position_stack(vjust=0.5),colour="black")  +
  scale_y_continuous(labels = scales::percent,limits = c(0, 1))+
  labs(y="", x="",fill="",title =' Cetuximab')
cet

dato_nt.summary <- dato_nt %>% group_by(isPaneth,V2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  #group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count))
         #pos =1-( cumsum(percent) - 0.5*percent))

                                                   

nt<-ggplot(dato_nt.summary,  aes(x=isPaneth, y=percent, fill=V2)) +
  geom_bar(stat='identity',  width = .7, colour="black", lwd=0.1) +
  geom_text(aes(label=ifelse(percent >= 0.02, paste0(sprintf("%.0f", percent*100),"%"),"")), position=position_stack(vjust=0.5),colour="black")+
  scale_y_continuous(labels = scales::percent,limits = c(0, 1))+
  labs(y="", x="",fill="",title =' Not Treated')


c<-ggarrange(nt,cet,ncol=2,common.legend = TRUE, legend = 'right')
c<-annotate_figure(c, top = text_grob(sample_name, color = "red", face = "bold", size = 14))

ggsave(file=output_plot, plot=c, width=10, height=8)