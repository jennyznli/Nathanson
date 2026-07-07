#!/usr/bin/env Rscript
library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)

# script to create manhattan plot of either a single chromosome of the whole genome
# this script greatly increases the speed of plotting by subsampling a small fraction
# of the non-significant snps for snptest 2.5.6

# INPUT: INFILE (text file) SNPTEST 2.5.6 output
# OUTPUT: OUTFILE (string), the name of the pdf to output

# ------------------------------ get.GIF ------------------------------ #
get.GIF <- function( p.vec )
    # compute genomic inflation factor
    # input: p.vec (numeric), vector of p-values
    # output: GIF
{
    median( qchisq(1- p.vec,1))/qchisq(0.5,1)
}
# --------------------------------------------------------------------- #

# -------------------------- gg.qq ------------------------------------ #
gg.qq <- function( p.vec, x.pos , y.pos  )
    # recreation of qq from qqman, but using ggplot
    # input:
    #   p.vec (numeric), vector of p-values
    #.  x.pos, y.pos (numeric), x/y coordinates for the position of the GIF text
    # output:
    #   p1, a ggplot object
{

    if( (!is.numeric(p.vec)))
        stop("input must be numeric")

    p.vec <- p.vec[!is.na(p.vec) & !is.nan(p.vec) & !is.null(p.vec) & is.finite(p.vec) &
                       p.vec < 1 & p.vec > 0]
    o = -log10(sort(p.vec, decreasing = FALSE))
    e = -log10(ppoints(length(p.vec)))

    txt = paste("lambda", round(get.GIF(p.vec),2), sep = " == ")
    p1 <- ggplot(data = data.frame(o = o, e = e), aes(x = e, y = o)) + geom_point() +
        theme_minimal() +
        xlab("Expected") +
        ylab("Observed") +
        geom_abline(slope=1, intercept=0, color="red") +
        annotate("text", x = x.pos, y = y.pos,
                 label = txt,
                 parse = T,
                 size = unit(10, "pt"))
    return(p1)
}
# --------------------------------------------------------------------- #

args = commandArgs(trailingOnly=TRUE)
if( length(args) < 2 )
{
    print("need to provide 2 arguments: INFILE OUTFILE")
    print("INFILE should be a SNPTEST 2.5.6 file")
    stop()
}

INFILE  = args[1]
OUTFILE = args[2]

# PARAMS
# the threshold for p-values; sub-sample any values with p > p.thresh
p.thresh = 0.01

# proportion of snps to keep
k.thresh = 0.2

MANOUTFILE <- paste(OUTFILE, "_man.png", sep="")
QQOUTFILE <- paste(OUTFILE, "_qq.png", sep = "")

# ---

print("reading association testing data...")

dat <- fread(INFILE, header=TRUE)

dat <- dat[,c("rsid", "frequentist_add_pvalue", "chromosome", "position", "all_maf")]

if( any(dat$chromosome == "X"))
{
    dat$chromosome[which(dat$chromosome == "X")] <- 23
}

colnames(dat) <- c("rsid", "p", "chr", "bp", "maf")
dat$chr <- as.numeric(dat$chr)
chrlabs <- as.character(seq(1:length(unique(dat$chr))))
dat$p <- as.numeric(dat$p)

# remove any rsids that snuck in
dat <- dat[!is.na(dat$chr),]
dat <- dat[!is.na(dat$p),]

print("done!")

# ---

print("subsampling association data...")

# s is the sub sample of nonsignificant
x <- dim(dat[which(dat$p > p.thresh),])[1]
s <- sample(dat$rsid[which(dat$p > p.thresh)], round(k.thresh * x), replace=F)

print("done!")

# ---

print("creating manhattan plot...")

png(MANOUTFILE, height=600, width=1200)
manhattan(dat[(dat$rsid %in% dat$rsid[match(s, dat$rsid)]) | (dat$p <= p.thresh),],
          chr="chr", bp="bp", p="p", snp="rsid", col=c("gray", "black"),
          logp=TRUE, chrlabs = chrlabs)
dev.off()

print("done!")

print("creating qq plot...")

png(QQOUTFILE)
print(gg.qq(dat$p,1,3))
dev.off()
print("done!")



INFILE  = here("breast", "data",  "ALL_snptest_chr10_maf05")
x <- as.data.frame(fread(INFILE))
# x <- x %>% sort(x$frequentist_add_pvalue)

x_sorted <- x[order(x$frequentist_add_pvalue, decreasing = FALSE), ]

# x_sorted <- x[, c("frequentist_add_pvalue")]

OUTFILE = args[2]



