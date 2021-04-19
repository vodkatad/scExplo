load('/mnt/trcanmed/snaketree/prj/scRNA/dataset/scRNA_Ire_cetuxi_reseq/aggr_0.2.umap_cl.pdf.Rdata')
DimPlot(srdata, reduction = "umap")
ggsave("~/tesire/aggr_basic_umap02.svg")
#Saving 14.1 x 9.11 in image
DimPlot(srdata, reduction = "umap",group.by="Phase")
ggsave("~/tesire/aggr_cycle.svg")
#DimPlot(srdata, reduction = "umap", group.by="percent.mito")


paneth <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/dataset/scRNA_Ire_cetuxi/paneth_cells_cycle.tsv', sep="\t", stringsAsFactors = F)
colnames(paneth) <- c('id','cycle','cl','sample')
paneth$paneth <- "yes"
paneth$barcode <- sapply(samples, function(x) {x[1]})
samples <- strsplit(paneth$id, '-', fixed=T)
paneth$barcode <- sapply(samples, function(x) {x[1]})
table(paneth$sample)
paneth[paneth$sample == "CRC0327_NT_1",'model'] = 1
paneth[paneth$sample == "CRC0327_cetux_1",'model'] = 7
paneth[paneth$sample == "CRC0322_NT_1",'model'] = 8
paneth[paneth$sample == "CRC0322_cetux_1",'model'] = 2
paneth[paneth$sample == "CRC0327_NT_2",'model'] = 3
paneth[paneth$sample == "CRC0327_cetux_2",'model'] = 4
paneth[paneth$sample == "CRC1502_NT_1",'model'] = 5
paneth[paneth$sample == "CRC1502_cetux_1",'model'] = 6
table(paneth$model)
head(paneth)
paneth$cellid <- paste0(paneth$barcode, '-', paneth$model)
head(paneth)
tail(paneth)
mm <- data.frame(id=colnames(srdata))
m <- merge(paneth, mm, by.y="id", by.x="cellid", all.y=T, sort=F)
head(m)
#m$id <- as.character(m$id)
#m$cellid <- as.character(m$cellid)
#mm$id <- as.character(m$id)
#mmm <- m[match(colnames(srdata), m$cellid),]
mmm <- m[match(mm$id, m$cellid),]
mmm[is.na(mmm$paneth),'paneth'] <- 'no'
all(mmm$cellid==colnames(srdata))
srdata[['paneth']] <- mmm[,"paneth"]
DimPlot(srdata, reduction="umap", group.by="paneth")
ggsave("~/tesire/aggr_panethss.svg")
samples <- strsplit(colnames(srdata), '-', fixed=T)
s <- sapply(samples, function(x) {x[2]})
cells <- sapply(samples, function(x) {x[1]})
#samples_id <- data.frame(id=s)
osamples = c("CRC0327_NT_1", "CRC0322_cetux_1", "CRC0327_NT_2", "CRC0327_cetux_2", "CRC1502_NT_1", "CRC1502_cetux_1","CRC0327_cetux_1", "CRC0322_NT_1")
#merged <- merge(samples_id, translate, by="id")
#rownames(merged) <- colnames(srdata)
samples_id <- data.frame(id=s)
translate <- data.frame(id=seq(1,length(osamples)), sample=osamples)
merged <- merge(samples_id, translate, by="id")
rownames(merged) <- colnames(srdata)
srdata[['model']] <- merged[,"sample"]
DimPlot(srdata, reduction="umap", group.by="model")
ggsave("~/tesire/aggr_model.svg")

# less clusters
srdata <- FindClusters(srdata, resolution = 0.1)
DimPlot(srdata, reduction="umap")
ggsave("~/tesire/aggr_basic_umap01.svg")

FeaturePlot(srdata, reduction="umap", features = "percent.mito")
ggsave("~/tesire/aggr_percmito.svg")
