#suppressPackageStartupMessages(library(corrplot))
#suppressPackageStartupMessages(library(reshape2))

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pvclust))
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(gridExtra))


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

infile <- "/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/final_results_dinuc_hosts.txt"
#infile <- "/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/final_results_dinuc_virus.txt"
#infile <- "/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/lala.txt"

df <- read.table(file=infile, sep="\t", header=TRUE)

x_3 <- df[,c(1,2,3,6)]

x_3$NORM <- scale(x_3$ODDS, center = TRUE, scale = TRUE)
x <- acast(x_3, GENOME~DINUC, value.var="NORM")

(jColors <-
  with(x_3,
       data.frame(LABEL = levels(LABEL),
                  COLOR = I(brewer.pal(nlevels(LABEL), name = 'Accent')))))

x_3$COLOR <- jColors$COLOR[match(x_3$LABEL, jColors$LABEL)]
x_3$LABEL <- jColors$LABEL[match(x_3$LABEL, jColors$LABEL)]

species2color <- unique((x_3[,c(1,6)]))

species2label <- unique((x_3[,c(1,2)]))

la <- species2label[with(species2label, order(GENOME)),]

x_2 <- as.data.frame(x)

x_2$LABEL <- la$LABEL

corr <- cor(x)
p.mat <- cor.mtest(corr)
opa = "Hosts"
#  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
#  col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
#                             "cyan", "#007FFF", "blue", "#00007F"))

pdf("Virus_corrplot.pdf", width = 10, height = 10)

corrplot(corr,
         method = "circle",
#         title=opa,
         tl.col = "black",
         tl.cex = 1,
         type="full",
         #           sig.level = 0.0001,
         sig.level = 0.05/(sum(15:1)), #Bonferroni
         hclust.method = c("ward.D2"),
         #           order = c("alphabet"),
         order = c("hclust"),
         p.mat = p.mat,
         mar=c(0,0,1,0),
         #           addrect=10,
         col=col(21),
         addCoef.col = "black", # Add coefficient of correlation
         number.cex = 0.8,
         insig = "blank")

dev.off()

#species2color <- unique((x_3[,c(1,2)]))

#df[ order(as.numeric(row.names(df))),]

#la <- species2color[with(species2color, order(GENOME)),]

#my_palette <- colorRampPalette(c("blue","white","red"))(n = 21)

#distance = dist(as.matrix(corr), method = "euclidean")

#cluster = hclust(distance, method = "ward.D2")

#heatmap.2(corr, density.info="histogram", col = my_palette, cexRow=1, cexCol=1, margins =c(6,6), Rowv=as.dendrogram(cluster), Colv=as.dendrogram(cluster), trace=c("none"))

