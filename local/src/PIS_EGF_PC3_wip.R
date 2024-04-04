
#egrassi@godot:/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/Alberto_pops$ cut -f 2 /tmp/genelist_pc3_th0.1.tsv  | sed 1d  | tr -d '"'  > genes_pop_PC3.txt

d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/Alberto_PC3-scores.tsv', sep="\t", header=TRUE)

meta <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe', sep="\t", header=TRUE, stringsAsFactors = FALSE)

meta_o <- meta[grepl('LMO', meta$type),]
keep <- c('LMO_BASALE', 'LMO_cetuxi_72h_NT')
meta_egf <- meta_o[grepl(keep[1], meta_o$type) | grepl(keep[2], meta_o$type),]

scores <- d[, colnames(d) %in% meta_egf$sample_id_R]

ncol(scores) == nrow(meta_egf)

meta_egf$treat <- ifelse(grepl('BASALE', meta_egf$type), 'EGF', 'noEGF')
meta_egf <- meta_egf[order(meta_egf$treat),]
scores <- scores[, match(meta_egf$sample_id_R, colnames(scores))]

rownames(meta_egf) <- meta_egf$sample_id_R
meta_egf <- meta_egf[, c('treat', 'w3_cetuxi'), drop=F]

minv <- min(unlist(scores))
maxv <- max(unlist(scores))
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

#pheatmap(scores, annotation_col=meta_egf, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette)

pheatmap(scores, annotation_col=meta_egf, cluster_rows = F, cluster_cols=F, scale="row", show_colnames = FALSE, color=colorRampPalette(c("navy", "white", "red"))(50))



ann <- meta_egf
ann$sample <- substr(rownames(ann), 0, 7)
bscores <- scores
## PIS on this ssGSEA
fc <- function(model, data, sign, ann) {
  d <- data[rownames(data)==sign,]
  dt <- t(d)
  m <- merge(dt, ann, by="row.names")
  d <- m[m$sample==model, ]
  d <- d[order(d$Row.names),] # we resort to using gen id to identify treat/nt, but we check
  if (length(intersect(unique(d$treat), c('noEGF', 'EGF')))) {
    fc <- mean(d[d$treat=="EGF", sign])-mean(d[d$treat=="noEGF", sign])
  } else {
    return(NA)
  }
  return(fc)
}

panethIndScore <- sapply(unique(ann$sample), fc, bscores, 'PC3up', ann)
scores <- data.frame(row.names=unique(ann$sample), PIS=panethIndScore)
scores$lPIS <- scores$PIS
scores$model <- rownames(scores)
scores <- scores[!is.na(scores$PIS),]
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS)))+geom_col(fill='blue')+ylab("delta PC3")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))


panethIndScore <- sapply(unique(ann$sample), fc, bscores, 'PC3down', ann)
scores <- data.frame(row.names=unique(ann$sample), PIS=panethIndScore)
scores$lPIS <- scores$PIS
scores$model <- rownames(scores)
scores <- scores[!is.na(scores$PIS),]
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS)))+geom_col(fill='blue')+ylab("delta PC3")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))


# pheatmap on logFPKM?
d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/Alberto_Alberto-scores.tsv', sep="\t")
ann <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/samples_data', sep="\t", header=TRUE)
rownames(ann) <- ann$id
ann$id <- NULL

ann <- ann[order(ann$treat),]
d <- d[, match(rownames(ann), colnames(d))]


minv <- min(d)
maxv <- max(d)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette)

ann <- ann[order(ann$sample, ann$treat),]
d <- d[, match(rownames(ann), colnames(d))]
pheatmap(d, annotation_col=ann, cluster_rows = F, cluster_cols=F, breaks = bk, color=my_palette)