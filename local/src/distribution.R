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
cet<-snakemake@input$c[1]
nt<-snakemake@input$c[2]


output_plot<-snakemake@output$out[1]
#kmeans<- read.table(file = clust,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta_cet<- read.table(file = cet,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta_cet<-as.data.frame(meta_cet)
meta_nt<-read.table(file = nt,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta_nt<-as.data.frame(meta_nt)
print(length(meta_cet$x))
print(length(meta_nt$x))
meta_nt$Treatment <- "Not Treated"
meta_cet$Treatment<- "Cetuximab"


# Combina i due dataframe
df_combined <- rbind(meta_cet, meta_nt)
print(head(df_combined))

# Crea il grafico
x_min <- 0
x_max <- ceiling(max(df_combined$x))

# Crea una sequenza di tick breaks, inclusi il valore minimo e massimo
breaks_x <- pretty(c(x_min, x_max), n = 5) 
x_max<-max(breaks_x)
print(breaks_x) # "pretty" genera una serie di break esteticamente piacevoli
print(x_max)


a <- ggplot(df_combined, aes(x = x, color = Treatment)) +
  #geom_histogram(aes(fill='white'), alpha=0.3, binwidth=0.05, position = 'identity')+
  
  geom_density(position = "identity", bw = 0.05, size = 0.8) +#(aes(y=after_stat(scaled))),
  scale_color_manual(name = "Treatment", 
                     values = c("Cetuximab" = "red", "Not Treated" = "black")) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) +  # Specifica i tick manualmente
    # Mantiene l'asse Y gestito automaticamente
  ggtitle("Metagene Distribution") + labs(y = "Fraction of total cells(A.U.)", x = "Metagene")+
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        panel.grid.major = element_blank(),                          
        panel.grid.minor = element_blank(),                         
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.ticks.length = unit(0.2, "cm"),                         
        axis.title.x = element_text(size = 12),                      
        axis.title.y = element_text(size = 12))                      
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))


# p <- ggplot(data=df_combined, aes(x=x)) + 
#   geom_histogram(aes(fill='white'), alpha=0.3, binwidth=0.05, position = 'identity')
density_data <- ggplot_build(a)$data[[1]]

#ggp <- ggplot_build(a)
#maxy <- ceiling(max(ggp$data[[1]]$count * 0.05) )
maxy <- ceiling(max(density_data$y))# [[2]] if histogram is removed
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
