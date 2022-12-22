library(ggplot2)
library(ggsankey)
library(psych)


cc_us <- read.table('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/CRC0322LMO_NT72h_1_singleron0722_us_dir/filtered_annotated_CRC0322LMO_NT72h_1_singleron0722_us_cellCycle.csv', sep= ",", header=FALSE)
cc_orig <- read.table('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/CRC0322LMO_NT72h_1_singleron0722_dir/filtered_annotated_CRC0322LMO_NT72h_1_singleron0722_cellCycle.csv', sep= ",", header=FALSE)
m <- merge(cc_us, cc_orig, by="V1")

colnames(m) <- c('cellid', 'cycle_us', 'cycle_orig')
rownames(m) <- m$cellid
m$cellid <- NULL
df <- m %>%
make_long(cycle_us, cycle_orig)


ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)

# useful metrics
k <- cohen.kappa(m)
k$kappa
dim(m)
table(m[,1] == m[,2])



