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
library(wesanderson)
library(ggpubr)
library(reshape2) 
library(reshape) 
library(tidyr)




cetux<-snakemake@input$clust[1]
nt<-snakemake@input$clust[2]
output_plot<-snakemake@output[[1]]
cc_cetux<-snakemake@input$cc[1]
cc_nt<-snakemake@input$cc[2]
log<-snakemake@output[['log']]
print(log)


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

dato_cet.summary <- dato_cet %>% group_by(isPaneth) %>% #V2
  summarise(total_count=n(),.groups = 'drop') %>% 
  #group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count))
         #pos =1-( cumsum(percent) - 0.5*percent))

dato_cet.summary_paneth <-dato_cet %>% group_by(isPaneth,V2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count),mucca=1-( cumsum(percent) - 0.5*percent))
colori<-c(rainbow(10)[2],rainbow(10)[8])
cet<-ggplot(dato_cet.summary,  aes(x=isPaneth, y=percent,fill=isPaneth)) +#, fill=V
  geom_bar(stat='identity',  width = .7, fill=colori, lwd=0.1) +
  geom_text(aes(label=ifelse(percent >= 0.03, paste0(sprintf("%.0f", percent*100),"%"),"")), position=position_stack(vjust=0.5),colour="white")  +
  scale_y_continuous(labels = scales::percent,limits = c(0,1))+
  labs(y="", x="")
#cet



#dato_cet.summary_paneth$isPaneth<-as.factor(dato_cet.summary_paneth)
colori_ciclo<-c("#5BBCD6","#FF0000", "#00A08A")  
ciclo<-ggplot(dato_cet.summary_paneth,aes(x=" ",y=percent, fill=V2)) +
  geom_bar(width = 1, stat = "identity")+scale_fill_manual(values=colori_ciclo) +geom_text(aes(y=mucca,label=ifelse(percent >= 0.03, paste0(sprintf("%.0f", percent*100),"%"),"")),colour="black")+
  coord_polar("y", start=0) +
  facet_grid(.~ isPaneth) +theme_void()
#ciclo
design<-"
  1
  3
"
l<-ciclo+cet+ plot_layout(design = design)+plot_annotation(title = 'Cetuximab')
#print(l)

dato_nt.summary <- dato_nt %>% group_by(isPaneth) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  #group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count))
         #pos =1-( cumsum(percent) - 0.5*percent))

dato_nt.summary_paneth <-dato_nt %>% group_by(isPaneth,V2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count),mucca=1-( cumsum(percent) - 0.5*percent))


nt<-ggplot(dato_nt.summary,  aes(x=isPaneth, y=percent,fill=isPaneth)) +#, fill=V
  geom_bar(stat='identity',  width = .7, fill=colori, lwd=0.1) +
  geom_text(aes(label=ifelse(percent >= 0.03, paste0(sprintf("%.0f", percent*100),"%"),"")), position=position_stack(vjust=0.5),colour="white")  +
  scale_y_continuous(labels = scales::percent,limits = c(0,1))+
  labs(y="", x="")

#nt

nt_ciclo<-ggplot(dato_nt.summary_paneth,aes(x=" ",y=percent, fill=V2)) +
  geom_bar(width = 1, stat = "identity")+scale_fill_manual(values=colori_ciclo) +geom_text(aes(y=mucca,label=ifelse(percent >= 0.03, paste0(sprintf("%.0f", percent*100),"%"),"")),colour="black")+
  coord_polar("y", start=0) +
  facet_grid(.~ isPaneth) +theme_void()

#nt_ciclo

m<-nt_ciclo+nt+ plot_layout(design = design)+plot_annotation(title = 'Not treated')

#print(m)
c<-ggarrange(m,l,ncol=2,common.legend = TRUE, legend = 'right')
c<-annotate_figure(c, top = text_grob(sample_name, color = "red", face = "bold", size = 14))

ggsave(output_plot, plot=c,, width=180, height=180, units="mm")


cont_cet<- dato_cet.summary_paneth[,c('isPaneth','V2','total_count')]
spread_df <- spread(cont_cet, key = V2, value = total_count)
spread_df<-as.data.frame(spread_df)
row.names(spread_df)<-spread_df$isPaneth
spread_df$isPaneth<-NULL

dato_cet.summary_ciclo <-dato_cet %>% group_by(V2) %>% 
  summarise(count_cet=n(),.groups = 'drop') 
dato_nt.summary_ciclo <-dato_nt %>% group_by(V2) %>% 
  summarise(count_nt=n(),.groups = 'drop') 

paneth_cc_cet<-as.data.frame(dato_cet.summary_paneth[dato_cet.summary_paneth$isPaneth=='Paneth',])
row.names(paneth_cc_cet)<-paneth_cc_cet$V2
or_cet<-(paneth_cc_cet['S','total_count']+paneth_cc_cet['G2M','total_count'])/sum(paneth_cc_cet$total_count)
print('odds_ratio cetuximab')
print(or_cet)
paneth_cc_nt<-as.data.frame(dato_nt.summary_paneth[dato_nt.summary_paneth$isPaneth=='Paneth',])
row.names(paneth_cc_nt)<-paneth_cc_nt$V2
or_nt<-(paneth_cc_nt['S','total_count']+paneth_cc_nt['G2M','total_count'])/sum(paneth_cc_nt$total_count)
print('odds_ratio_nt')
print(or_nt)

sink(log)
print('Cetuximab Paneth vs Not Paneth')
cont_cet<- dato_cet.summary_paneth[,c('isPaneth','V2','total_count')]
spread_df <- spread(cont_cet, key = V2, value = total_count)
spread_df<-as.data.frame(spread_df)
row.names(spread_df)<-spread_df$isPaneth
spread_df$isPaneth<-NULL
print('contigency table cetuximab Paneth vs Not Paneth')
print(t(spread_df))
chisq <- chisq.test(t(spread_df))
print(chisq)

print('Not treated Paneth vs Not Paneth')
cont_nt<- dato_nt.summary_paneth[,c('isPaneth','V2','total_count')]
spread_df <- spread(cont_nt, key = V2, value = total_count)
spread_df<-as.data.frame(spread_df)
row.names(spread_df)<-spread_df$isPaneth
spread_df$isPaneth<-NULL
print('contigency table not treated Paneth vs Not Paneth')
print(t(spread_df))
chisq <- chisq.test(t(spread_df))
print(chisq)

dato_cet.summary_ciclo<-as.data.frame(dato_cet.summary_ciclo)
row.names(dato_cet.summary_ciclo)<- dato_cet.summary_ciclo$V2
dato_cet.summary_ciclo$V2<-NULL

dato_nt.summary_ciclo<-as.data.frame(dato_nt.summary_ciclo)
row.names(dato_nt.summary_ciclo)<- dato_nt.summary_ciclo$V2
dato_nt.summary_ciclo$V2<-NULL
cet_nt<-rbind(t(dato_cet.summary_ciclo),t(dato_nt.summary_ciclo))
chisq <- chisq.test(t(cet_nt))
print('contigency table NT vs Treated')
print(t(cet_nt))
print(chisq)
#Paneth vs Paneth
treated_paneth<-dato_cet.summary_paneth[dato_cet.summary_paneth$isPaneth=='Paneth',]
treated_paneth<-treated_paneth[,c('V2',"total_count")]
setnames(treated_paneth, old = c('total_count'),
         new = c('total_count_treated'))
treated_paneth<-as.data.frame(treated_paneth)
row.names(treated_paneth)<-treated_paneth$V2
treated_paneth$V2<-NULL


nt_paneth<-dato_nt.summary_paneth[dato_nt.summary_paneth$isPaneth=='Paneth',]
nt_paneth<-nt_paneth[,c('V2',"total_count")]
setnames(nt_paneth, old = c('total_count'),
         new = c('total_count_not_treated'))
nt_paneth<-as.data.frame(nt_paneth)
row.names(nt_paneth)<-nt_paneth$V2
nt_paneth$V2<-NULL

paneth<-rbind(t(nt_paneth),t(treated_paneth))
chisq <- chisq.test(paneth)
print('chi-squared test Paneth not treated vs paneth treated')
print('contigency table')
print(paneth)
print(chisq)


#Not Paneth vs Not Paneth 

treated_paneth<-dato_cet.summary_paneth[dato_cet.summary_paneth$isPaneth=='nPaneth',]
treated_paneth<-treated_paneth[,c('V2',"total_count")]
setnames(treated_paneth, old = c('total_count'),
         new = c('total_count_treated'))
treated_paneth<-as.data.frame(treated_paneth)
row.names(treated_paneth)<-treated_paneth$V2
treated_paneth$V2<-NULL


nt_paneth<-dato_nt.summary_paneth[dato_nt.summary_paneth$isPaneth=='nPaneth',]
nt_paneth<-nt_paneth[,c('V2',"total_count")]
setnames(nt_paneth, old = c('total_count'),
         new = c('total_count_not_treated'))
nt_paneth<-as.data.frame(nt_paneth)
row.names(nt_paneth)<-nt_paneth$V2
nt_paneth$V2<-NULL

paneth<-rbind(t(nt_paneth),t(treated_paneth))
chisq <- chisq.test(paneth)
print('chi-squared test nPaneth not treated vs npaneth treated')
print('contigency table')
print(paneth)
print(chisq)


sink()
