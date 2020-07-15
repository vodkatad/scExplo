paneth <- read.table('paneth_cells_cycle.tsv', sep="\t", stringsAsFactors = F)
colnames(paneth) <- c('id','cycle','cl','sample')
paneth$paneth <- "yes"
paneth$barcode <- sapply(samples, function(x) {x[1]})
samples <- strsplit(paneth$id, '-', fixed=T)
paneth$barcode <- sapply(samples, function(x) {x[1]})
table(paneth$sample)
paneth[paneth$sample == "CRC0327_NT_1",'model'] = 1
paneth[paneth$sample == "CRC0327_cetux_1",'model'] = 2
paneth[paneth$sample == "CRC0322_NT_1",'model'] = 3
paneth[paneth$sample == "CRC0322_cetux_1",'model'] = 4
paneth[paneth$sample == "CRC0327_NT_2",'model'] = 5
paneth[paneth$sample == "CRC0327_cetux_2",'model'] = 6
paneth[paneth$sample == "CRC1502_NT_1",'model'] = 7
paneth[paneth$sample == "CRC1502_cetux_1",'model'] = 8
table(paneth$model)
head(paneth)
paneth$cellid <- paste0(paneth$barcode, '-', paneth$model)
head(paneth)
tail(paneth)
m <- merge(paneth, mm, by.y="id", by.x="cellid", all.y=T, sort=F)
head(m)
mmm <- m[match(mm$id, m$cellid),]
all(mmm$cellid==colnames(srdata))
srdata[['paneth']] <- mmm[,"ppaneth"]
DimPlot(srdata, reduction="umap", group.by="paneth")
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
samples_id <- data.frame(id=s)
osamples = c("CRC0327_NT_1", "CRC0327_cetux_1", "CRC0322_NT_1", "CRC0322_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2", "CRC1502_NT_1", "CRC1502_cetux_1"])
osamples = c("CRC0327_NT_1", "CRC0327_cetux_1", "CRC0322_NT_1", "CRC0322_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2", "CRC1502_NT_1", "CRC1502_cetux_1")
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
srdata[['model']] <- merged[,"sample"]
DimPlot(srdata, reduction="umap", group.by="model")
DimPlot(srdata, reduction="tsne", group.by="model")
DimPlot(srdata, reduction="tsne", group.by="paneth")
DimPlot(srdata, reduction="tsne", group.by="cycle")
DimPlot(srdata, reduction="tsne", group.by="Phase")
DimPlot(srdata, reduction="tsne")
head(markers_aggr)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10 <- markers_aggr %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(srdata, features = top10$gene) + NoLegend()
DoHeatmap(srdata, features = top10$gene) + NoLegend()
FeaturePlot(srdata, features = c('SOX4','ATOH1') )
VlnPlot(srdata, features = c('SOX4','ATOH1') )
png('heat.png');DoHeatmap(srdata, features = top10$gene) + NoLegend(); dev.off()
top10
top10[top10$cluster=="7","gene"]
markers_aggr[markers_aggr$cluster=="7" & markers_aggr$p_val_adj<0.001,]
ma <- markers_aggr[markers_aggr$cluster=="7" & markers_aggr$p_val_adj<0.001,]
ma[ma$gene =="EREG",]
ma[ma$gene =="AREG",]
ma[ma$gene =="MYC",]
ma[ma$gene =="DLL1",]
ma[ma$gene =="DEFA6",]
ma[ma$gene =="DEFA5",]
ma[order(ma$p_val_adj),]
