library(data.table)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)


input<-snakemake@input[['data']]
variance<-snakemake@input[['variance']]


pdf_<-snakemake@output[['pdf']]
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
gene_variance <- apply(dato, 1, var)

# Step 2: Ordinare i geni in base alla varianza (in ordine decrescente)
ordered_genes <- order(gene_variance, decreasing = TRUE)

# Step 3: Selezionare i primi 1000 geni più variabili
top <- dato[ordered_genes[1:1000], ]

top <- as.matrix(top)

#top[top > 7] <- 7

top[is.na(top)] <- 0
top[is.nan(top)] <- 0
top[is.infinite(top)] <- 0
top <- top[apply(top, 1, var) > 0, ]


library(pheatmap)

# Crea una palette di colori che varia tra viola e verde
my_palette <- colorRampPalette(c("purple", "white", "green"))(256)

#print(head(rownames(dato)))
# Genera la heatmap senza nomi di righe, colonne e dendrogrammi
top <- scale(top)
p<-pheatmap(top, 
                  # Scala i valori per ogni gene (riga)
         cluster_rows = TRUE,   # Non mostrare il dendrogramma per le righe
         cluster_cols = TRUE,   # Non mostrare il dendrogramma per le colonne
         show_rownames = FALSE,  # Nascondi i nomi delle righe
         show_colnames = FALSE,  # Nascondi i nomi delle colonne
         col = my_palette, 
         breaks=seq(min(top), max(top),length.out = 257)     # Usa la palette viola-verde
        )

pdf(pdf_)
print(p)
graphics.off()








