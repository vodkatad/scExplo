# provare con metagene

d <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/HER/CRC0322_NT_1_3000_log2cpm.csv', row.names=1, sep=",")

library(ggplot2)
ggplot(data=d, aes(x=DLL1, y=HES1))+geom_point()+theme_bw()
cet <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/HER/CRC0322_cetux_1_log2cpm.csv', row.names=1, sep=",")
ggplot(data=cet, aes(x=DLL1, y=HES1))+geom_point()+theme_bw()

paneth <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/kmeans/CRC0322_cetux_1/CRC0322_cetux_1_kmeans_2comp_cinque.csv', sep=",", row.names=1)
m <- merge(paneth, cet, by="row.names")
cet <- ggplot(data=m, aes(x=DLL1, y=HES1, color=isPaneth))+geom_point()+theme_bw(base_size=15)+ggtitle('CRC0322 Cetuximab')
panethnt <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/kmeans/CRC0322_NT_1_3000/CRC0322_NT_1_3000_kmeans_2comp_cinque.csv', sep=",", row.names=1)
mnt <- merge(panethnt, d, by="row.names")
nt <- ggplot(data=mnt, aes(x=DLL1, y=HES1, color=isPaneth))+geom_point()+theme_bw(base_size = 15)+ggtitle('CRC0322 NT')


ggsave('~/CRC0322_CET.png', plot=cet)
ggsave('~/CRC0322_NT.png', plot=nt)