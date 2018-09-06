#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))

rm(list=ls())

options(warn=-1)

args = commandArgs(trailingOnly=TRUE)

#infile <- "/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/test/output/tmp/RSCU_files/final_results_RSCU_minus_stop_monocodonic.txt"

infile <- args[1]

df <- read.table(file=infile, sep="\t", header=TRUE)

outdir <- args[2]

outfile <- paste0(outdir, "RSCU_heatmap.pdf")

#outfile <- c("RSCU_heatmap.pdf")

pdf(outfile, width = 10, height = 10)

x_3 <- df

x_3$NORM <- scale(x_3$RSCU, center = TRUE, scale = TRUE)

x <- acast(x_3, GENOME~CODON, value.var="NORM")

my_palette <- colorRampPalette(c("blue","white","red"))(n = 51)

(jColors <-
   with(x_3,
        data.frame(LABEL = levels(LABEL),
                   COLOR = I(brewer.pal(nlevels(LABEL), name = 'Accent')))))

x_3$COLOR <- jColors$COLOR[match(x_3$LABEL, jColors$LABEL)]

species2color <- unique((x_3[,c(1,7)]))

la <- species2color[with(species2color, order(GENOME)),]

distance1 = dist(as.matrix(x), method = "euclidean")
distance2 = dist(as.matrix(t(x)), method = "euclidean")
cluster1 = hclust(distance1, method = "ward.D2")
cluster2 = hclust(distance2, method = "ward.D2")

heatmap.2(x,
          density.info="histogram",
          RowSideColors=la$COLOR,
          col = my_palette,
          cexRow=0.6,
          cexCol=1.1,
          margins =c(5,5),
          Rowv=as.dendrogram(cluster1),
          Colv=as.dendrogram(cluster2),
          labRow=FALSE,
          trace=c("none"))
          par(lend = 1)           # square line ends for the color legend

legend(y=1.1, x=0.85, xpd=TRUE,      # location of the legend on the heatmap plot
       legend = jColors$LABEL, # category labels
       col = jColors$COLOR,# color key
       lty= 1,             # line style
       lwd = 10,           # line width
       cex=1.3,
       box.lwd = 1, box.col = "white",bg = "transparent"
)

dev.off()

