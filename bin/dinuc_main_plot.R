#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

rm(list=ls())

options(warn=-1)

is.even <- function(x) x %% 2 == 0

args = commandArgs(trailingOnly=TRUE)

#infile <- ("/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/test/output/tmp/dinuc_files/final_results_dinuc.txt")

infile <- args[1]

df <- read.table(file=infile, sep="\t", header=TRUE)

class_list <- as.list(as.matrix(unique(df$LABEL)))

outdir <- args[2]

#outdir <- c("/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/test/output/plots/")

setwd(outdir)

outfile <- c("dinuc_main_plot.pdf")

pdf(file=outfile, onefile=FALSE)

theme_update(plot.title = element_text(hjust = 0.5))
la <- lapply(class_list, function(x)
     ggplot(data=df[df$LABEL == x,], aes(x=DINUC,y=ODDS), xlab="", ylab="dinucleotide odds ratio", label=GENOME) + 
     xlab("dinucleotide") +
     ylab("odds ratio") +
  #ggplot(data= subset(data_frame, FAMILY != "Ifla" | FAMILY != "Dici"), aes(x=DINUC,y=ODDS), xlab="", ylab="dinucleotide odds ratio", label=GENOME) + 
  #geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size=0.8, position = position_jitter(width = 0.2)) + 
  #geom_jitter(data = subset(data_frame, FAMILY == "Ifla"), color="red") +
  #geom_jitter(data = subset(data_frame, FAMILY == "Dici"), color="blue") +
  #geom_errorbar(stat = "hline", yintercept = "median", col="red", width=0.8, size=1.3, aes(ymax=..y..,ymin=..y..)) + 
  geom_hline(yintercept=1) +
  geom_hline(yintercept=0.78, linetype="dashed") +
  geom_hline(yintercept=1.24, linetype="dashed") +
  stat_summary(fun.y=median, geom="errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.75, size=1, color="red") +
  #geom_crossbar(aes(ymin = 0, ymax = Sepal.Width), size=1,col="red", width = .5)
  #geom_text(data = data_frame, aes(x = names(round_sd) , y = 2, label = paste(round_sd)), size=3) +
  #geom_text(label = data_frame$GENOME, size=2) +
  theme(text = element_text(size=10), axis.text.x = element_text(size=7)) +
  ylim(0,2) +
  annotate("text",  x=Inf, y = Inf, label = x, vjust=1.5, hjust=1.5, size=4)
)

ncols <- round((length(class_list))/2)
ml <- marrangeGrob(la, nrow=2, ncol=ncols, top='DINUCLEOTIDE ODDS RATIO')
print(ml)
#ggsave(outfile, ml, width=10, height=10)
dev.off()
