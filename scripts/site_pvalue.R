#!/apps/lab/aryee/R/R-3.2.3/bin/Rscript --vanilla

# Usage example using test data from the circleseq repository:
#  ./site_pvalue.R ../test/U2OS_EMX1_counts.txt ../test/U2OS_EMX1_counts_pval.txt 

# Usage example using a larger test dataset on erisone:
#  ./site_pvalue.R /data/joung/CIRCLE-Seq/complete_analysis/160122_937aa31/output/U2OS_EMX1_counts.txt U2OS_EMX1_counts_pval.txt 

library("ggplot2")
library("scales")

args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
#infile <- "../test/U2OS_EMX1_counts.txt"
# Read in counts
message("Reading ", infile)
dat <- read.delim(infile, comment.char = "", header=TRUE)

# Condition on having observed at least one read
bg <- dat$Control_Position_Reads
bg <- bg[bg>0] 

# Model control distribution as exponential
message("Calculating p-values")
lambda <- mean(bg)
pval <- 1 - pexp(dat$Nuclease_Position_Reads, rate=1/lambda)
dat$pvalue <- pval

# Model control distribution empirically
background_cdf <- ecdf(bg)
pval_empirical <- 1 - background_cdf(dat$Nuclease_Position_Reads)

message("Saving diagnostic plots to pvalue_diagnostics.pdf")
# Diagnostic plots
pdf(file="pvalue_diagnostics.pdf", width=6, height=2.5)
p <- ggplot(dat, aes(1+Control_Position_Reads)) + scale_x_continuous(limits=c(0,100)) + scale_y_log10(labels=comma) + geom_histogram(binwidth=2, na.rm=TRUE) + theme_bw() + ggtitle("Control_Position_Reads")
suppressWarnings(print(p))
p <- ggplot(dat, aes(1+Nuclease_Position_Reads)) + scale_x_continuous(limits=c(0,100)) + scale_y_log10(labels=comma) + geom_histogram(binwidth=2, na.rm=TRUE) + theme_bw() + ggtitle("Nuclease_Position_Reads")
suppressWarnings(print(p))
idx <- sample(length(pval), min(length(pval), 10000))
plot(pval_empirical[idx], pval[idx], xlab="Empirical p-value", ylab="Exponential model p-value")
abline(0,1)
dev.off()

message("Writing output table ", outfile)
write.table(dat, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)

