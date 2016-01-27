library(ggplot2)
library(scales)

# Read in counts
dat <- read.delim("/data/joung/CIRCLE-Seq/complete_analysis/160122_937aa31/output/U2OS_EMX1_counts.txt", comment.char = "", header=TRUE)
dat$ratioPos = (1+dat$Nuclease_Position_Coverage) / (1+dat$Control_Position_Coverage)

# Density plot
ggplot(dat, aes(log2(ratioPos))) + geom_histogram(binwidth=1) + theme_bw() + scale_y_continuous(labels=comma) + ggtitle("log2((1+nuclease_pos_reads)/(1+control_pos_reads))")

# Estimate empirical background CDF
bg <- log2(dat$ratioPos[dat$ratioPos<1])
background_cdf <- ecdf(-bg)

# Calculate p-values
dat$pvalue <- 1 - background_cdf(log2(dat$ratioPos))
minPval <- 1/length(bg)
dat$pvalue <- pmax(dat$pvalue, minPval)

