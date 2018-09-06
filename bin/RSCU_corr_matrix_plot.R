suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(reshape2))

rm(list=ls())

options(warn=-1)

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

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]

#infile <- "/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/test/output/tmp/RSCU_files/final_results_RSCU_minus_stop_monocodonic.txt"

df <- read.table(file=infile, sep="\t", header=TRUE)

outdir <- args[2]

clust_mode <- args[3]

if (clust_mode == "alphabet_aa") {
  mode = "alphabet"
} else if (clust_mode == "alphabet_nt") {
  mode = "alphabet"
}else if (clust_mode == "hclust") {
  mode = "hclust"
} else {
  stop("Mode must be either alphabet_aa, alphabet_nt or hclust, please check!")
}

outfile <- paste0(outdir, clust_mode, "_", "RSCU_corrplot.pdf")

#outfile <- c("RSCU_corrplot.pdf")

x_3 <- df

x_3$NORM <- scale(x_3$RSCU, center = TRUE, scale = TRUE)

if (clust_mode == "alphabet_nt") {
  x_3$CODON_PLUS_AA <- paste(x_3$CODON, "-", x_3$AA, sep="")
} else if (clust_mode == "alphabet_aa") {
  x_3$CODON_PLUS_AA <- paste(x_3$AA, "-", x_3$CODON, sep="")
} else {
  x_3$CODON_PLUS_AA <- paste(x_3$CODON, "-", x_3$AA, sep="")
}

pdf(file=outfile, width=10, height=10)

class_list <- as.list(as.matrix(unique(df$LABEL)))

la <- lapply(class_list, function(x)
             x_3[x_3$LABEL == x,])

names(la) <- class_list

ncols <- (round((length(class_list)/2)))

par(mfrow=c(2,ncols))
for (i in 1:length(la)) {
  x <- acast(la[[i]], GENOME~CODON_PLUS_AA, value.var="NORM")
  opa <- class_list[i]
  corr <- cor(x)
  p.mat <- cor.mtest(corr)
  #  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
#    col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
#                               "cyan", "#007FFF", "blue", "#00007F"))
  
  corrplot(corr,
           title=opa,
           tl.col = "black",
           tl.cex = 0.5,
           type="full",
           sig.level = 0.01,
           hclust.method = c("ward.D2"),
           order = mode,
           p.mat = p.mat,
           mar=c(0,0,1,0),
           col=col(21),
           insig = "blank")
}

dev.off()
