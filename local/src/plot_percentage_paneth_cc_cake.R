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

library(reshape2) 
library(reshape) 
library(tidyr)


clust<-snakemake@input[['clust']]
cc<-snakemake@input[['cc']]
plot_out<-snakemake@output[['plot_out']]
log_out<-snakemake@output[['log']]


dato<- read.table(file = clust,row.names = 1,sep=",",header = TRUE)

#dato<-rbind(dato_cet,dato_nt)
dato$isPaneth<-as.character(dato$isPaneth)
dato$isPaneth[dato$isPaneth=='filtered']<-'Others'
dato$isPaneth[dato$isPaneth=='nPaneth']<-'Others'
dato$isPaneth<-as.factor(dato$isPaneth)



dato_cc<- read.table(file = cc,row.names = 1,sep=",",header = FALSE)
dato<-merge(dato,dato_cc,by='row.names')

dato.summary <- dato %>% group_by(isPaneth) %>% #V2
  summarise(total_count=n(),.groups = 'drop') %>% 
  #group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count))
         #pos =1-( cumsum(percent) - 0.5*percent))

dato.summary_paneth <-dato %>% group_by(isPaneth,V2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count),mucca=1-( cumsum(percent) - 0.5*percent))
colori<-c('blue','red')
paneth<-ggplot(dato.summary,  aes(x=isPaneth, y=percent,fill=isPaneth)) +#, fill=V
  geom_bar(stat='identity',  width = .7, fill=colori, lwd=0.1) +
  geom_text(aes(label=ifelse(percent >= 0.03, paste0(sprintf("%.0f", percent*100),"%"),"")), position=position_stack(vjust=0.5),colour="white")  +
  scale_y_continuous(labels = scales::percent,limits = c(0,1))+
  labs(y="", x="")
#cet
cont<- dato.summary_paneth[,c('isPaneth','V2','total_count')]
spread_df <- spread(cont, key = V2, value = total_count)
spread_df<-as.data.frame(spread_df)
row.names(spread_df)<-spread_df$isPaneth
spread_df$isPaneth<-NULL

chisq <- chisq.test(t(spread_df))
#dato_cet.summary_paneth$isPaneth<-as.factor(dato_cet.summary_paneth)
colori_ciclo<-c("#5BBCD6","#FF0000", "#00A08A")  
ciclo<-ggplot(dato.summary_paneth,aes(x=" ",y=percent, fill=V2)) +
  geom_bar(width = 1, stat = "identity")+scale_fill_manual(values=colori_ciclo) +geom_text(aes(y=mucca,label=ifelse(percent >= 0.03, paste0(sprintf("%.0f", percent*100),"%"),"")),colour="black")+
  coord_polar("y", start=0) +
  facet_grid(.~ isPaneth) +theme_void()+ggtitle(paste('chi_squared: pvalue:',round(chisq$p.value,15)))+theme(plot.title = element_text(hjust = 0.5))

#ciclo
design<-"
  1
  3
"
l<-ciclo+paneth+ plot_layout(design = design)
#print(l)



#print(m)

pdf(plot_out)
print(l)
graphics.off()



