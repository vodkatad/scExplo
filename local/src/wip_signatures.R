library(ggplot2)
d <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_sign.tsv.gz'), sep="\t", header=T)

d$plasticity <- d$esophagus4-d$colon3


ggplot(data=d, aes(x=onf2, y=plasticity))+geom_point(size=0.5)+geom_smooth(method='lm')+theme_bw(base_size = 20)

d$type <- ifelse(grepl(".T", rownames(d), fixed=T), 'tumor', 'normal')
ggplot(data=d, aes(x=onf2, y=stem1))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

do <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_sign_others.tsv.gz'), sep="\t", header=T)
m <- merge(d, do, by="row.names")
m$stem_index <- m$RSC3 - m$CBC5

ggplot(data=m, aes(x=onf2, y=coreHRC1))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=m, aes(x=onf2, y=ganesh_moorman_fetal2))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=m, aes(x=coreHRC1, y=ganesh_moorman_fetal2))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


ggplot(data=m, aes(x=stem1, y=CBC5))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


ggplot(data=m, aes(x=stem_index, y=onf2))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


d2 <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_entropy.tsv.gz'), sep="\t", header=T)
m <- merge(m, d2, by.x="Row.names", by.y='row.names')

### ERC GPS
tum <- m[m$type!='normal',]
library(MASS)
library(ggrastr)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#https://slowkow.com/notes/ggplot2-color-by-density/
compare <- function(x, y, log, nx, ny, ggt, legend) {
  zx <- length(x[x==0.000000001])
  zy <- length(y[y==0.000000001])
  zz <- x==0.000000001 & y==0.000000001
  zz <- sum(zz)
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=scale(x), y=scale(y), density=get_density(x,y, n=100))
  y_breaks <- guess_ticks(d$y)
  x_breaks <- guess_ticks(d$x)
  if (legend) {
    p <- ggplot(d, aes(x=x, y=y, color=density)) +rasterize(geom_point(size=0.1), dpi=300)+theme_bw()+theme(text = element_text(size=20))+xlab(nx)+ylab(ny)+scale_color_viridis_c()+ggt+
        scale_y_continuous(breaks=y_breaks,limits=c(min(y_breaks), max(y_breaks)),expand = c(0, 0))+
        scale_x_continuous(breaks=x_breaks,limits=c(min(x_breaks), max(x_breaks)),expand = c(0, 0))+
          theme(legend.position = 'none')
  } else {
    p <- ggplot(d, aes(x=x, y=y, color=density)) +rasterize(geom_point(size=0.1), dpi=300)+theme_bw()+theme(text = element_text(size=20))+xlab(nx)+ylab(ny)+scale_color_viridis_c()+ggt+
      scale_y_continuous(breaks=y_breaks,limits=c(min(y_breaks), max(y_breaks)),expand = c(0, 0))+
      scale_x_continuous(breaks=x_breaks,limits=c(min(x_breaks), max(x_breaks)),expand = c(0, 0))
      
  } 
  #+labs(caption=paste0(pe$estimate, ', pval=', pe$p.value)
  return(list(plot=p, num=c(zx, zy, zz, pe$estimate, pe$p.value)))
}

textSize <- 5
largerSize <- textSize + 2

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
unmute_theme <- theme(
  text = element_text(size = textSize, family='sans'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(color = "black"),
  panel.background = element_blank()
)

# function that given values to be plotted on an axis will return:
# vector of breaks, trying to guess which max will be the best one
# this will be used as scale_y_continuous(breaks=  and as ylim(min, max) to have the - also limits-c()
# last tick at the extremity of the axis.
# other parameter is n. of ticks
guess_ticks <- function(values, nticks=5, fixed_max=NULL) {
  vmax <- max(values)
  if (is.null(fixed_max)) { 
    round_max <- ceiling(round(vmax, digits=3))
  } else {
    round_max <- fixed_max
  }
  fixed_min <- floor(round(min(values), digits=3))
  my_breaks <- seq(fixed_min, round_max, length.out=nticks)
  return(my_breaks)
}

svg('~/1.svg', width=1.744, height=1.744)
compare(x=tum$onf2, y=tum$coreHRC1, log=FALSE, nx='OnF', ny='coreHRC', ggt=unmute_theme, legend=T)
dev.off()
svg('~/2.svg', width=1.744, height=1.744)
compare(x=tum$onf2, y=tum$coreHRC1, log=FALSE, nx='OnF', ny='coreHRC', ggt=unmute_theme, legend=F)
dev.off()
svg('~/3.svg', width=1.744, height=1.744)
compare(x=tum$onf2, y=tum$stem_index, log=FALSE, nx='OnF', ny='stem index', ggt=unmute_theme, legend=T)
dev.off()
svg('~/4.svg', width=1.744, height=1.744)
compare(x=tum$onf2, y=tum$stem_index, log=FALSE, nx='OnF', ny='stem index', ggt=unmute_theme, legend=F)
dev.off()
svg('~/5.svg', width=1.744, height=1.744)
compare(x=tum$onf2, y=tum$entropy, log=FALSE, nx='OnF', ny='entropy', ggt=unmute_theme, legend=T)
dev.off()
svg('~/6.svg', width=1.744, height=1.744)
compare(x=tum$onf2, y=tum$entropy, log=FALSE, nx='OnF', ny='entropy', ggt=unmute_theme, legend=F)
dev.off()

ck <- m[m$type != 'normal',]

ggplot(data=m, aes(x=plasticity, y=entropy))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


ggplot(data=m, aes(x=onf2, y=entropy))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

#### cytotrace
ctt <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_ctt.tsv.gz'), sep="\t", header=T)
rownames(m) <- m$Row.names
m$Row.names <- NULL
m2 <- merge(m, ctt, by="row.names")
ggplot(data=m2, aes(x=plasticity, y=x))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))
ggplot(data=m2, aes(x=entropy, y=x))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=m, aes(x=stem1, y=entropy))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


ggplot(data=d, aes(x=type, y=onf2))+geom_boxplot()+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=d, aes(x=type, y=stem1))+geom_boxplot()+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

#     217 esophagus_mucosa

#?

dd <- d[grepl('SMC', rownames(d)),]
ggplot(data=dd, aes(x=onf2, y=plasticity))+geom_point(size=0.5)+geom_smooth(method='lm')+theme_bw(base_size = 20)

ggplot(data=dd, aes(x=onf2, y=stem1))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

############ us
d <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_ussign.tsv.gz'), sep="\t", header=T)

d$plasticity <- d$esophagus4-d$colon3


ggplot(data=d, aes(x=onf2, y=plasticity))+geom_point(size=0.5)+geom_smooth(method='lm')+theme_bw(base_size = 20)

d$type <- ifelse(grepl("cetux", rownames(d), fixed=T), 'cetux', 'NT')
ggplot(data=d, aes(y=onf2, x=type, color=type))+geom_violin(size=0.5)+geom_boxplot(width=0.3, size=0.5)+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=d, aes(y=CBC5, x=type, color=type))+geom_violin(size=0.5)+geom_boxplot(width=0.3, size=0.5)+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))


ggplot(data=d, aes(y=stem1, x=type, color=type))+geom_violin(size=0.5)+geom_boxplot(width=0.3, size=0.5)+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))
ggplot(data=d, aes(y=plasticity, x=type, color=type))+geom_violin(size=0.5)+geom_boxplot(width=0.3, size=0.5)+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

d2 <- read.table(gzfile('/mnt/cold1/snaketree/prj/scRNA/dataset/plasticy_scores/all_usentropy.tsv.gz'), sep="\t", header=T)

rownames(d) <- gsub('_ussign.tsv', '', rownames(d), fixed=T)
rownames(d2) <- gsub('_usentropy.tsv', '', rownames(d2), fixed=T)
m <- merge(d, d2, by="row.names")

ggplot(data=m, aes(x=plasticity, y=entropy))+geom_point(size=0.5, aes(color=type))+geom_smooth(method='lm')+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

ggplot(data=m, aes(y=entropy, x=type, color=type))+geom_violin(size=0.5)+geom_boxplot(width=0.3, size=0.5)+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))

m$sample <- sapply(strsplit(x=m$Row.names, split='_'), function(x) {paste0(c(x[1], x[3]), collapse="-")})
ggplot(data=m, aes(y=entropy, x=type, color=type))+geom_violin(size=0.5)+geom_boxplot(width=0.3, size=0.5)+
  theme_bw(base_size = 20)+scale_color_manual(values=c('blue', 'red'))+facet_wrap(~sample)
