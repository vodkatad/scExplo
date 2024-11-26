library(data.table)

library(ggplot2)

library(uwot)
#input<-'/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_Ire_cetuxi/CRC0327_cetux_2_dir/filtered_annotated_saver_ribomito_CRC0327_cetux_2_log2_pc1_cpm.csv'
umap<-snakemake@input[['umap']]
kmeans<-snakemake@input[['kmeans']]
plot_out<-snakemake@output[['plot_out']]


umap<- read.table(file = umap,row.names = 1,sep=",",header = TRUE,,stringsAsFactors = FALSE)
kmeans<-read.table(file = kmeans,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
kmeans$isPaneth[kmeans$isPaneth=='filtered']<-'Others'
kmeans$isPaneth[kmeans$isPaneth=='nPaneth']<-'Others'

data<-merge(umap,kmeans,by='row.names')
print(head(data))



colori<-c('blue','red')
pdf(plot_out)
ggplot(data, aes(x=x, y=y,color=isPaneth)) + geom_point(size=0.4)+
  labs(title = "Kmeans Results", color="")+scale_color_manual(values = colori)+#scale_color_viridis(discrete=TRUE,direction = -1)
  xlab('umap1')+ylab('umap2')+
  theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
graphics.off()
