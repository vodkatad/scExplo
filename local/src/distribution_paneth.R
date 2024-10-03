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
#library(ggpubr)

#output_plot<-snakemake@output[[1]]
meta_path<-snakemake@input$metagene
clu_path<-snakemake@input$cluster



output_plot<-snakemake@output$out[1]
#kmeans<- read.table(file = clust,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta<- read.table(file = meta_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta<-as.data.frame(meta)

cluster<-read.table(file = clu_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
cluster<-as.data.frame(cluster)

df_combined<-merge(meta,cluster,by='row.names')
print(head(df_combined))
rownames(df_combined)<-df_combined$Row.names
df_combined$Row.names<-NULL
df_combined<-df_combined[df_combined$isPaneth!='filtered',]

# Combina i due dataframe


# Crea il grafico
x_min <- 0
x_max <- ceiling(max(df_combined$x))

# Crea una sequenza di tick breaks, inclusi il valore minimo e massimo
breaks_x <- pretty(c(x_min, x_max), n = 5) 
x_max<-max(breaks_x)
print(breaks_x) # "pretty" genera una serie di break esteticamente piacevoli
print(x_max)
# Crea il grafico con tick specificati manualmente
df_combined$isPaneth <- factor(df_combined$isPaneth, levels = c("Paneth", "nPaneth"), labels = c("Paneth", "Others"))
colori<-c(rainbow(10)[2],rainbow(10)[5],rainbow(10)[8])
# Crea il grafico con il nuovo nome nella legenda
a <- ggplot(df_combined, aes(x = x, color = isPaneth)) +
  geom_density(aes(y = after_stat(count * 0.05)), position = "identity", bw = 0.05, size = 1) +
  scale_color_manual(name = "Cluster", 
                     values = c("Others" = colori[1], "Paneth" = colori[3])) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) +  # Specifica i tick manualmente
  ggtitle("Metagene Distribution") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Sfondo bianco
        panel.grid.major = element_blank(),                          # Rimuovi griglie maggiori
        panel.grid.minor = element_blank(),                          # Rimuovi griglie minori
        axis.line = element_line(color = "black"),                   # Colore nero per gli assi
        axis.ticks = element_line(color = "black"),                  # Tick marks neri
        axis.ticks.length = unit(0.2, "cm"),                         # Lunghezza dei tick
        axis.title.x = element_text(size = 12),                      # Etichetta asse X
        axis.title.y = element_text(size = 12))                      # Etichetta asse Y
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))


p <- ggplot(data=df_combined, aes(x=x)) + 
  geom_histogram(aes(fill='white'), alpha=0.3, binwidth=0.05, position = 'identity')

ggp <- ggplot_build(a)
maxy <- ceiling(max(ggp$data[[1]]$count * 0.05) )# [[2]] if histogram is removed
breaks_y <- pretty(c(0, maxy), n = 5) 
maxy<-max(breaks_y)
print(maxy)
print(breaks_y)
a<-a+scale_y_continuous(expand = c(0, 0), limits = c(0, maxy),breaks=breaks_y) +
# a<-ggplot() +#+geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
# #geom_histogram(fill='blue', binwidth=0.05, position = 'identity')+#+
# geom_density(data=meta_cet,position="identity", aes(x=x,y=after_stat(count*0.05)), bw = 0.05,color = "red")+
# geom_density(data=meta_nt,position="identity", aes(x=x,y=after_stat(count*0.05)), bw = 0.05,color = "black")+
# scale_color_manual(name = "Treatment", 
#                      values = c("black" = "Not treated", "red" = "cetuximab")) +
# ggtitle("Metagene Distribution")
# #scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
# #scale_x_continuous(breaks=x_breaks,limits=c(0, max(x_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
# #unmute_theme#+theme(legend.position="none", axis.text.x = element_blank(), 
#      #axis.ticks.x = element_blank(),
#     #legend.spacing.y = uni#t(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))      
pdf(output_plot)
print(a)
graphics.off()
 

#dato<-rbind(dato_cet,dato_nt)
# kmeans$isPaneth<-as.character(kmeans$isPaneth)
# kmeans$isPaneth[kmeans$isPaneth=='filtered']<-'nPaneth'
# kmeans$isPaneth<-as.factor(kmeans$isPaneth)
# kmeans<-merge(kmeans,meta,by='row.names')

# b<-ggplot(data=kmeans, aes(x=x,fill=isPaneth)) +geom_density(alpha=0.8) +
#   ggtitle("Metagene Distribution")
# print(b)

# library(plyr)
# mu <- ddply(kmeans, "isPaneth", summarise, grp.mean=mean(x))
# b<-ggplot(data=kmeans, aes(x=x,fill=isPaneth,color=isPaneth)) +geom_density(alpha=0.5) +
#   ggtitle("Metagene Distribution ")+theme_bw()
# print(b)

# g<-ggarrange(a,b, ncol = 2, nrow = 1)
# ggsave(file=output_plot, plot=g, width=16, height=8)
# graphics.off()