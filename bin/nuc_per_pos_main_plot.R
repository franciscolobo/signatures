suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

rm(list=ls())

options(warn=-1)

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]

df <- read.table(file=infile, sep="\t", header=TRUE)

class_list <- as.list(as.matrix(unique(df$LABEL)))

outdir <- args[2]

outfile <- c("nuc_per_pos_main_plot.pdf")

outfile_path <- paste0(outdir, outfile)

pdf(file=outfile_path, onefile=FALSE)

theme_update(plot.title = element_text(hjust = 0.5))
la <- lapply(class_list, function(x)
  p <- ggplot(df[df$LABEL == x,], aes(x=NUCLEOTIDE, y=VALUE), xlab="Nucleotide", ylab="Frequency") +
#  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.2, h = 0), size = 1.5, stroke = 0, shape = 16) +
  stat_summary(fun.y=median, geom="errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.75, size=1, color="red") +
  labs(x = "nt per codon position", y="frequency") +
  ggtitle(label=x) +
  facet_grid(. ~ POSITION)
)

ncols <- round((length(class_list))/2)
ml <- marrangeGrob(la, nrow=2, ncol=ncols, top='Nucleotide frequency per codon position')
#plot(ml)
print(ml)

dev.off()

