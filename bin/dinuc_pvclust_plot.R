library(reshape2)
library(pvclust)

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]

#infile <- ("/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/test/output/tmp/dinuc_files/final_results_dinuc.txt")

df <- read.table(file=infile, sep="\t", header=TRUE)

outdir <- args[2]

bootstrap <- args[3]

outfile1 <- paste0(outdir, "dinuc_pvclust_features.pdf")

outfile2 <- paste0(outdir, "dinuc_pvclust_genomes.pdf")

x_3 <- df[,c(1,2,3,6)]

x_3$NORM <- scale(x_3$ODDS, center = TRUE, scale = TRUE)

x <- acast(x_3, GENOME~DINUC, value.var="NORM")

fit <- pvclust(x, method.hclust="ward.D2",method.dist="euclidean", parallel=TRUE)

fit2 <- pvclust(t(x), method.hclust="ward.D2", method.dist="euclidean", parallel=TRUE)

pdf(outfile1, width = 10, height = 10)

plot(fit, print.num=FALSE, cex.pv=1.5, cex=1.5, col.pv=c(2,0,0), main="Agrupamento de dinucleotídeos", mar=c(5, 8, 4, 1), hang = -1, lty=1, lwd=3) # dendogram with p values

dev.off()

pdf(outfile2, width = 20, height = 10)

plot(fit2, print.num=FALSE, cex.pv=1, cex=0.4, col.pv=c(2,0,0), main="Agrupamento de dinucleotídeos", mar=c(5, 8, 4, 1), hang = -1, lty=1, lwd=2) # dendogram with p values

dev.off()
