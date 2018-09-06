#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pvclust))
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(gridExtra))

rm(list=ls())

options(warn=-1)

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]

df <- read.table(file=infile, sep="\t", header=TRUE)

outdir <- args[2]

outfile <- paste0(outdir, "dinuc_heatmap.pdf")

x_3 <- df[,c(1,2,3,6)]

x_3$NORM <- scale(x_3$ODDS, center = TRUE, scale = TRUE)

x <- acast(x_3, GENOME~DINUC, value.var="NORM")

my_palette <- colorRampPalette(c("blue","white","red"))(n = 21)

(jColors <-
   with(x_3,
        data.frame(LABEL = levels(LABEL),
                   COLOR = I(brewer.pal(nlevels(LABEL), name = 'Accent')))))

x_3$COLOR <- jColors$COLOR[match(x_3$LABEL, jColors$LABEL)]

species2color <- unique((x_3[,c(1,6)]))

la <- species2color[with(species2color, order(GENOME)),]

distance1 = dist(as.matrix(x), method = "euclidean")
distance2 = dist(as.matrix(t(x)), method = "euclidean")

cluster1 = hclust(distance1, method = "ward.D2")
cluster2 = hclust(distance2, method = "ward.D2")

pdf(file=outfile, width=10, height = 10)

heatmap.2(x, density.info="histogram", RowSideColors=la$COLOR, col = my_palette, cexRow=0.6, cexCol=2, margins =c(6,6), Rowv=as.dendrogram(cluster1), Colv=as.dendrogram(cluster2), labRow=FALSE, trace=c("none"))
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