#!/usr/bin/env Rscript
library(getopt)
library(rCASC)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'counts', 'c', 1, 'character',
  'gtf', 'g', 1, 'character',
  'scratch', 's', 1, 'character',
  'saver', 'a', 1, 'character',
  'saverla', 'f', 1, 'character',
  'saverfarm', 'r', 1, 'character',
  'saverfarmva', 'v', 1, 'character'
  ), ncol=4, byrow=TRUE)
opt <- getopt(opts)

# TODO put the lots of hardcoded thresholds in a configuration file (per sample)
if (is.null(opt$counts) | !is.null(opt$help) | is.null(opt$gtf) | is.null(opt$scratch) | is.null(opt$saver) | is.null(opt$saverla)
   |is.null(opt$saverfarm) | is.null(opt$saverfarmva)) {
    cat(getopt(opts, usage=TRUE))
    stop('everything is mandatory :P')
}

#WD = "/path_to_saver_and_raw_data"
#INPUT="CRC0322_NT_1_bis.csv"
SEPARATOR=","
#GTF="Homo_sapiens.GRCh38.101.gtf"
#SCRATCH="/home/rcalogero/scratch"
INPUT <- opt$counts
SCRATCH <- opt$scratch
SEPARATOR <- ","
GTF <- opt$gtf
#odir <- dirname(opt$output)
saver_f <- opt$saver
saver_farm_f <- opt$saverfarm
saver_farmva_f <- opt$saverfarmva
saver_last <- opt$saverla

#setwd(odir)
#filtering on raw data
print('############################# scanno1')
scannobyGtf(group = "docker", file=INPUT, gtf.name=GTF, biotype="protein_coding", mt = FALSE, ribo.proteins = FALSE, umiXgene = 3, riboStart.percentage = 1, riboEnd.percentage = 50, mitoStart.percentage = 1, mitoEnd.percentage = 40, thresholdGenes = 250)
# thir produces (in the directory where INPUT is)
# filtered_annotated_CRC0542_CTX72h_1.tsv  annotated_CRC0542_CTX72h_1.tsv CRC0542_CTX72h_1_annotated_genes.pdf  (who uses them then?) filteredStatistics.txt
# so we want another directory for saver???

# calculating cell cycle, here not on inputed data!
# its input is filtered_annotated_CRC0542_CTX72h_1.tsv
fa_name <- file.path(dirname(INPUT), paste0('filtered_annotated_', basename(INPUT)))
print('############################# ccycle')
seurat_ccycle(group = "docker", scratch.folder=SCRATCH, file=fa_name, separator=SEPARATOR, seed=111)
# this generates filtered_annotated_CRC0542_CTX72h_1_cellCycle.tsv in the fa_name dir - but who uses THIS?

#annotating saver data
print('############################# scanno2')
scannobyGtf(group = "docker", file=saver_f, gtf.name=GTF, biotype="protein_coding", mt = FALSE, ribo.proteins = FALSE, umiXgene = 3, riboStart.percentage = 0, riboEnd.percentage = 100, mitoStart.percentage = 0, mitoEnd.percentage = 100, thresholdGenes = 1)
# this generates annotated_saver_{sample}.tsv in dir

print('############################# saver read/writetables')
raw <- read.table(fa_name, sep=SEPARATOR, header=T, row.names=1)

saver_fa_name <- file.path(dirname(saver_f), paste0('annotated_', basename(saver_f)))
saver <- read.table(saver_fa_name, sep=SEPARATOR, header=T, row.names=1)
# we read annotated_saver_{sample}.tsv and filter the cells (keep only those in raw, which comes from the first scannobyGtf, no too high ribo/mito)
saver_ribomito <- saver[,which(names(saver)%in%names(raw))]
write.table(saver_ribomito, saver_farm_f, sep=",", col.names=NA)

#filtering
print('############################# topx1')
topx(group = "docker", file=saver_farm_f, threshold=10000, separator=SEPARATOR, type="variance")
print('############################# topx2')
topx(group = "docker", file=saver_farmva_f, threshold=5000, separator=SEPARATOR, type="expression")
#system("mkdir VandE")
#file.rename(from=paste("filtered_expression_filtered_variance_filtered_annotated_saver_ribomito", INPUT, sep="_"), to="VandE.csv")
#system("mv VandE.csv VandE")

print('############################# pcaeval')
seuratPCAEval(group = "docker", scratch.folder=SCRATCH, file=saver_last, separator=",", logTen = 0, seed = 111)
###seuratPCAEval(group = "docker", scratch.folder=SCRATCH, file=saver_farm_f, separator=",", logTen = 0, seed = 111)
# this generates 
#CRC0069_CTX72h_1_dir/Results/filtered_expression_filtered_variance_filtered_annotated_saver_ribomito_CRC0069_CTX72h_1/PCE_bowPlot.pdf
