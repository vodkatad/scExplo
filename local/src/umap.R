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

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)


