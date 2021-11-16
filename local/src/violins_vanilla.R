library(ggplot2)
data <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_longer/allres0.05_G1_violins/norm.tsv.gz'), sep='\t', header=T)

cellsinfo <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_longer/allres0.05_G1_clu_cycle.tsv', sep="\t", header=T)

my_vio <- function(gene, umi, annot, annot_col_x, annot_col_points, th) {
  df <- t(umi[rownames(umi)==gene,, drop=FALSE])
  #colnames(df) <- gene
  mdf <- merge(df, annot, by="row.names")
  mdf$Row.names <- NULL
  mdf[,annot_col_x] <- as.factor(mdf[, annot_col_x])
  mdf[,annot_col_points] <- as.factor(mdf[, annot_col_points])
  ggplot(data=mdf, aes_string(x=annot_col_x, y=gene))+geom_violin()+geom_jitter(aes_string(color=annot_col_points), alpha=0.2)+th
}

cellsinfo$sample <- as.character(cellsinfo$sample)
cellsinfo$sample <- gsub('_saverAlivePC.csv.gz','', cellsinfo$sample, fixed=TRUE)
my_vio('ENSG00000198719', data, cellsinfo, 'cluster', 'sample', current_theme)
my_vio('ENSG00000198719', data, cellsinfo, 'sample', 'cluster', current_theme)

g <- 'ENSG00000117707'
my_vio(g, data, cellsinfo, 'cluster', 'sample', current_theme)
my_vio('ENSG00000198719', data, cellsinfo, 'sample', 'cluster', current_theme)

g <- 'ENSG00000163586'
my_vio(g, data, cellsinfo, 'cluster', 'sample', current_theme)
my_vio(g, data, cellsinfo, 'sample', 'cluster', current_theme)

data <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_longer/all2res0.05_G1_violins/norm.tsv.gz'), sep='\t', header=T)
cellsinfo <- read.table('/mnt/cold1/snaketree/prj/scRNA/dataset/rCASC_longer/all2res0.05_G1_clu_cycle.tsv', sep="\t", header=T)


g <- 'ENSG00000117707'
my_vio(g, data, cellsinfo, 'cluster', 'sample', current_theme)
#my_vio(g, data, cellsinfo, 'sample', 'cluster', current_theme)



############3

textSize <- 2
current_theme <-
  theme_bw() +
  theme(
    strip.text = element_text(size = rel(textSize)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = rel(textSize)),
    axis.text.x = element_text(size = rel(textSize), angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0,
                               size = rel(textSize)),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    plot.title = element_text(
      face = "bold",
      size = rel(textSize),
      hjust = 0.5
    )
  )

