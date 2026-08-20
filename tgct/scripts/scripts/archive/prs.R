
# read in file 
x <- read.table("GWAS2026_PRS_SNPs.txt", header = TRUE)
x$Z <- x$Effect / x$StdErr

write.table(x, "GWAS2026_PRS_SNPs2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
