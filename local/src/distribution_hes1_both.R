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
x_max <- ceiling(max(df_combined$HES1))
print('SUKA')
# Crea una sequenza di tick breaks, inclusi il valore minimo e massimo
breaks_x <- pretty(c(x_min, x_max), n = 5) 
x_max<-max(breaks_x)
print(breaks_x) # "pretty" genera una serie di break esteticamente piacevoli
print(x_max)

df_combined$isPaneth <- factor(df_combined$isPaneth, levels = c("Paneth", "nPaneth"), labels = c("Paneth", "Others"))
colori<-c(rainbow(10)[2],rainbow(10)[5],rainbow(10)[8])

a <- ggplot(df_combined, aes(x = HES1, color = isPaneth)) +
  geom_density(aes(y = after_stat(count * 0.05)), position = "identity", bw = 0.05, size = 1) +
  scale_color_manual(name = "Cluster", 
                     values = c("Others" = colori[1], "Paneth" = colori[3])) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) +  # Specifica i tick manualmente
  ggtitle("HES1 Distribution") +
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


ggp <- ggplot_build(a)
data <- ggp$data[[1]]

data$y_orig <- data$y
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
data$y_scaled <- range01(data$y_orig)

maxy <- ceiling(max(data$y_scaled))
#maxy <- ceiling(max(ggp$data[[1]]$count * 0.05) )# [[2]] if histogram is removed
breaks_y <- pretty(c(0, maxy), n = 5) 
maxy<-max(breaks_y)
#labels_y <- sprintf("%.2f", breaks_y / maxy)
#print(maxy)
#print(breaks_y)
#a<-a+scale_y_continuous(expand = c(0, 0), limits = c(0, maxy),breaks=breaks_y,labels=labels_y) +
#  labs(y = "Relative cells number (A.U.)", x = "Metagene")

data$coloripaneth <- ifelse(data$colour == colori[1], 'Others', 'Paneth')

a <- ggplot(data=data, aes(x=x, y=y_scaled, color=coloripaneth))+geom_line()+
  scale_color_manual(name = "Cluster", 
                     values = c("Others" = colori[1], "Paneth" = colori[3])) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) +  # Specifica i tick manualmente
  ggtitle("HES1 Distribution") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Sfondo bianco
        panel.grid.major = element_blank(),                          # Rimuovi griglie maggiori
        panel.grid.minor = element_blank(),                          # Rimuovi griglie minori
        axis.line = element_line(color = "black"),                   # Colore nero per gli assi
        axis.ticks = element_line(color = "black"),                  # Tick marks neri
        axis.ticks.length = unit(0.2, "cm"),                         # Lunghezza dei tick
        axis.title.x = element_text(size = 12),                      # Etichetta asse X
        axis.title.y = element_text(size = 12))+
        scale_y_continuous(expand = c(0, 0), limits = c(0, maxy),breaks=breaks_y) +
  labs(y = "Relative cells number (A.U.)", x = "HES1")                      # Etichetta asse Y
  #guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))


pdf(output_plot)
print(a)
graphics.off()
 

