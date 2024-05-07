abe <- read.table('~/pnp_abe.csv', sep=',', header=F)

fitR1 <- read.table('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/CRC0322LMO_CTX72h_1_singleron0722_dir/ENSG00000172238:ATOH1.csv', sep=",", header=T)

apply(abe, 2, function(x) {cor.test(x, fitR$comp.1)$estimate})

fitR2 <- read.table('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/CRC0322LMO_CTX72h_1_singleron0722_dir/ENSG00000162676:GFI1.csv', sep=",", header=T)

apply(abe, 2, function(x) {cor.test(x, fitR$comp.1)$estimate})

fitR3 <- read.table('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/CRC0322LMO_CTX72h_1_singleron0722_dir/ENSG00000164816:DEFA5.csv', sep=",", header=T)

apply(abe, 2, function(x) {cor.test(x, fitR$comp.1)$estimate})


fitR4 <- read.table('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/CRC0322LMO_CTX72h_1_singleron0722_dir/ENSG00000164822:DEFA6.csv', sep=",", header=T)

apply(abe, 2, function(x) {cor.test(x, fitR$comp.1)$estimate})



fitR5 <- read.table('/mnt/cold2/snaketree/prj/scRNA/dataset/invivo_singleron_07_2022_rCASC/CRC0322LMO_CTX72h_1_singleron0722_dir/ENSG00000198719:DLL1.csv', sep=",", header=T)

apply(abe, 2, function(x) {cor.test(x, fitR$comp.1)$estimate})


dd <- read.table('/data/reanalysis_on_AIsc/paneth_only/all/log2cpmfisher/singctx_log2cpmfisher.csv', sep= ",", header=T, row.names=1)


library(poolr)
  
myp <- cbind(fitR1$comp.1,fitR2$comp.1,fitR3$comp.1,fitR4$comp.1,fitR5$comp.1)
fi <- apply(myp, 1, function(x) {fisher(x)$p})


