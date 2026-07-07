#!/usr/bin/env Rscript
library(data.table)
library(qqman)
library(ggplot2)

# script to create manhattan plot of each chromosome from REGENIE output
# this script greatly increases the speed of plotting by subsampling a small fraction
# of the non-significant snps

# INPUT: INFILE .regenie file
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

# args = commandArgs(trailingOnly=TRUE)
# if( length(args) < 2 )
# {
#     print("need to provide 2 arguments: INFILE OUTFILE")
#     print("INFILE should be a REGENIE format text file")
#     print("OUTFILE should be the output prefix (without extension)")
#     stop()
# }

INFILE  = here("breast", "data", "breast_AFR_genome.regenie")
OUTFILE = "breast_genome"

print(paste("Input file:", INFILE))
print(paste("Output prefix:", OUTFILE))

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
dat <- as.data.frame(dat)

# Handle X chromosome conversion
if( any(dat$CHROM == "X"))
{
    dat$CHROM[which(dat$CHROM == "X")] <- 23
}

dat$ID <- as.character(dat$ID)
dat$CHROM <- as.integer(dat$CHROM)
dat$GENPOS <- as.integer(dat$GENPOS)
dat$CHISQ <- as.numeric(dat$CHISQ)
dat$LOG10P <- as.numeric(dat$LOG10P)
dat$INFO <- as.numeric(dat$INFO)
dat$PVAL <- as.numeric(10^(-dat$LOG10P))

chrlabs <- as.character(seq(1:length(unique(dat$CHROM))))

print("done!")

# ---

print("subsampling association data...")

# s is the sub sample of nonsignificant
x <- dim(dat[which(dat$PVAL > p.thresh),])[1]
s <- sample(dat$ID[which(dat$PVAL > p.thresh)], round(k.thresh * x), replace=F)

print("done!")

# ---

print("creating manhattan plot...")

png(MANOUTFILE, height=600, width=1200)
manhattan(dat[(dat$ID %in% dat$ID[match(s, dat$ID)]) | (dat$PVAL <= p.thresh),],
          chr="CHROM", bp="GENPOS", p="PVAL", snp="ID", col=c("gray", "black"),
          logp=TRUE, chrlabs = chrlabs)
dev.off()

print("done!")

print("creating qq plot...")

png(QQOUTFILE)
print(gg.qq(dat$PVAL,1,3))
dev.off()
print("done!")

print(paste("Manhattan plot saved as:", MANOUTFILE))
print(paste("QQ plot saved as:", QQOUTFILE))


