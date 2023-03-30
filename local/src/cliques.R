library(data.table)
library(mclust)
library(reshape2)
library('psych')
library(pheatmap)
library(resample)
library(umap)
library(ggplot2)
library(dplyr)
library(viridis)
library(igraph)
library(patchwork)
library('stats')
if(!require('corrplot')) {
  install.packages('corrplot')
  library('corrplot')
}
stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
input<-input<-snakemake@input[['data']]
label_input<-snakemake@input[['label']]
out_universe<-snakemake@output[['universe']]
out_cliques<-snakemake@output[['clique']]
out_plot<-snakemake@output[['plot']]
corr_th<-as.double(snakemake@wildcards[['corr']])
print(out_cliques)
print(out_universe)

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
label<-read.table(file=label_input,row.names = 1,sep=',',header=TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)
label<-label[c('postSilh')]

colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
duplicate<-!duplicated(colnames(dato_t))
dato_t<-dato_t[,duplicate]

variance<-colVars(as.matrix(dato_t[sapply(dato_t, is.numeric)]))
keep<-names(variance[variance>0.1])
dato_<-dato_t[,c(keep)]
dato_<- merge(dato_, label, by = 'row.names')
dato_<-as.data.frame(sapply(dato_, as.numeric))
na_count <-sapply(dato_, function(y) sum(length(which(is.na(y)))))
a<-cor(dato_, method = c("pearson"))

#label clique
dato_1<-a[,c('postSilh')]
dato_1<-dato_1[abs(dato_1)>corr_th]
dato_1<-na.omit(dato_1)
adj<-a[c(names(dato_1)),c(names(dato_1))]
adj[adj>=corr_th]<-1
adj[adj<corr_th]<-0
graph<-graph_from_adjacency_matrix(adjmatrix = adj)
max<-largest_cliques(graph)
clique_graph <- induced_subgraph(graph, max[[1]])
cliccano<-names(max[[1]])
cliccano<-cliccano[cliccano != "postSilh"]
pdf(out_plot)
plot(clique_graph,layout=layout.kamada.kawai, vertex.color="yellow",vertex.label.cex=0.6,vertex.size=23,main='Label clique')
graphics.off()
clique_df<-as.data.frame(cliccano)
colnames(clique_df)<-'gene'
write.table(clique_df, file = out_cliques, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)
universe_df<-as.data.frame(keep)
colnames(universe_df)<-'gene'
write.table(universe_df, file = out_universe, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)