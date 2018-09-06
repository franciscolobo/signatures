library(ggplot2)
library(ggrepel)
library(tidyr)

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]

df <- read.table(file=infile, sep="\t", header=TRUE)

df$CG_freq <- df$C_freq + df$G_freq

outdir <- args[2]

flag_names <- args[3]

flag_NCP <- args[4]

print(flag_NCP)

if (flag_NCP == "TRUE") {
  outfile <- c("ENC_Ncp_main_plot.pdf")
} else if (flag_NCP == "FALSE") {
  outfile <- c("ENC_Nc_main_plot.pdf")
}
  
outfile_path <- paste0(outdir, outfile)

pdf(file = outfile_path, width = 8, height = 5)
#pdf(file = outfile, width = 8, height = 5)

#function to compute theoretical ENC from GC content only
theoretical_ENC <- function(GC_content, ... ){
  expected_ENC <- 2 + GC_content + (29/((GC_content)^2+(1-GC_content)^2))
  return(expected_ENC)
}
  
  #dummy function for ggplot2
  #p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
if (flag_NCP == "TRUE") {
  p <- ggplot(df, aes(x=df$CG_freq, y=df$Ncp, colour=df$LABEL)) +
  geom_point() +
  labs(x = "GC content", y="frequency") +
  labs(colour = "Groups") +
  stat_function(fun = theoretical_ENC, , colour = "black") + xlim(0,1)
} else {
  p <- ggplot(df, aes(x=df$CG_freq, y=df$Nc, colour=df$LABEL)) +
  geom_point() +
  labs(x = "GC content", y="frequency") +
  labs(colour = "Groups") +
  stat_function(fun = theoretical_ENC, colour = "black") + xlim(0,1)
}

if (flag_names == "TRUE") {
  plot(p + geom_label_repel(label.size=NA, fill=NA, max.iter = 3e3, force = 1, fontface = 'bold', color = 'black', box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"), segment.color = '#000000') +  ggtitle("NC plot") + labs(x="GC content", y="Ncp") + theme(panel.spacing = unit(10, "lines")))
} else {
  plot(p)
}

dev.off()