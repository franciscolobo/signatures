library(ggplot2)
library(ggrepel)
library(tidyr)
library(gridExtra)

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]

#infile <- "/Users/chico/projects/doutorado_tarcisio_final/bin/SIGNATURES/test/output/tmp/ENC_files/final_results_ENC.txt"

df <- read.table(file=infile, sep="\t", header=TRUE)

class_list <- as.list(as.matrix(unique(df$LABEL)))

df$CG_freq <- df$C_freq + df$G_freq

outdir <- args[2]

flag_names <- args[3]

flag_NCP <- args[4]

if (flag_NCP == "TRUE") {
  outfile <- c("ENC_Ncp_plot_per_group.pdf")
}

if (flag_NCP == "FALSE") {
  outfile <- c("ENC_Nc_plot_per_group.pdf")
}
  
outfile_path <- paste0(outdir, outfile)

#pdf(file = outfile_path, width = 8, height = 5)

#pdf(file = outfile, width = 8, height = 5)

pdf(file=outfile_path, width = 8, height = 5, onefile=FALSE)

#function to compute theoretical ENC from GC content only
theoretical_ENC <- function(GC_content, ... ){
  expected_ENC <- 2 + GC_content + (29/((GC_content)^2+(1-GC_content)^2))
  return(expected_ENC)
}

la <- lapply(class_list, function(x)
  if (flag_NCP == "TRUE") {
    ggplot(df[df$LABEL == x, ], aes(x=df[df$LABEL == x, ]$CG_freq, y=df[df$LABEL == x, ]$Ncp)) +
    geom_point() +
    labs(x = "GC content", y="frequency") +
#    labs(colour = "Groups") +
    stat_function(fun = theoretical_ENC, colour = "black") +
    annotate("text",  x=Inf, y = Inf, label = x, vjust=1.5, hjust=1.5, size=4) +
    xlim(0,1)
  } else {
    ggplot(df[df$LABEL == x, ], aes(x=df[df$LABEL == x, ]$CG_freq, y=df[df$LABEL == x, ]$Nc)) +
    geom_point() +
    labs(x = "GC content", y="frequency") +
#    labs(colour = "Groups") +
    stat_function(fun = theoretical_ENC, colour = "black") +
    annotate("text",  x=Inf, y = Inf, label = x, vjust=1.5, hjust=1.5, size=4) +
    xlim(0,1)
  }
)

if (flag_NCP == "TRUE") {
  ncols <- round((length(class_list))/2)
  ml <- marrangeGrob(la, nrow=2, ncol=ncols, top='Ncp plot')
  print(ml)
}

if (flag_NCP == "FALSE") {
  ncols <- round((length(class_list))/2)
  ml <- marrangeGrob(la, nrow=2, ncol=ncols, top='Nc plot')
  print(ml)
}

#if (flag_names == "TRUE") {
#  plot(p + geom_label_repel(label.size=NA, fill=NA, max.iter = 3e3, force = 1, fontface = 'bold', color = 'black', box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"), segment.color = '#000000') +  ggtitle("NC plot") + labs(x="GC content", y="Ncp") + theme(panel.spacing = unit(10, "lines")))
#} else {
#  plot(p)
#}

dev.off()