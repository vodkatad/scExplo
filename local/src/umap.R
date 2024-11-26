library(data.table)

library(ggplot2)

library(uwot)
#input<-'/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_Ire_cetuxi/CRC0327_cetux_2_dir/filtered_annotated_saver_ribomito_CRC0327_cetux_2_log2_pc1_cpm.csv'
input<-snakemake@input[['data']]
data_out<-snakemake@output[['out']]
plot_out<-snakemake@output[['umap_out']]
n_neighbors<-as.numeric(snakemake@wildcards[['n_neighbors']])
min_dist<-as.numeric(snakemake@wildcards[['min_dist']])
print(n_neighbors)
print(min_dist)

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)
#parti da qui 
umap_result <- umap(dato_t,min_dist = min_dist,n_neighbors = n_neighbors)
umap_result<-as.data.frame(umap_result)
colnames(umap_result)<-c('x','y')
write.csv(umap_result, data_out, row.names = TRUE)
umap_result$ATOH1<-dato_t[,grep('ATOH1',colnames(dato_t))]
colori<-c(rainbow(10)[2],rainbow(10)[5],rainbow(10)[8])
pdf(plot_out)
ggplot(umap_result, aes(x=x, y=y,color=ATOH1)) + geom_point(size=0.4)+
  labs(title = "Kmeans", color="")+scale_color_gradientn(colours = colorRampPalette(c("blue","cyan",'green',"yellow","red"))(100))+#scale_color_viridis(discrete=TRUE,direction = -1)
  xlab('tsne1')+ylab('tsne2')+
  theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
graphics.off()
