 
glmone <- function(logcpm, maxit=10000) {
    mixmdltemp=try(mixmdl <- normalmixEM(logcpm,k = 2,maxit=maxit,maxrestarts=1000000000000000000))
    return(list(means=mixmdltemp$mu, post=as.data.frame(mixmdltemp$posterior)))
}

get_gene <- function(g, data) {
    data[g,]
}

extract_fit <- function(g, data) {
    mat <- get_gene(g, data)
    glmone(mat)
}

listadilista <- lapply(geneListM, extract_fit)

for (i in seq(1, length(geneListM))) {
    write.table(listadilista[[i]]$post, paste0(geneListM[i], '.csv'), sep=",", quote=F)
}