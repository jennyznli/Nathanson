# ========================
# Packages
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "geneplotter", "plotly"
)
library(data.table)
library(qqman)
library(ggplot2)

purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

DATE <- format(Sys.Date(), "%Y%m%d")

# ========================================
# LOAD AND PREPARE DATA
# ========================================
chr1 <- as.data.frame(read.table(here("breast", "EUR_step2_chr1_Y1.regenie"), header = TRUE))
chr2 <- as.data.frame(read.table(here("breast", "EUR_step2_chr2_Y1.regenie"), header = TRUE))
chromosome <- "2"
ancestry_group <- "EUR"

# ========================================
# SUMMARY STATISTICS
# ========================================

# make a job that can taken in the file, extract the genomic inflation factor, make a plot to identify cany chunks misisng
# make manhatten plot and then a qq plot

# [1] "CHROM"   "GENPOS"  "ID"      "ALLELE0" "ALLELE1" "A1FREQ"  "INFO"    "N"       "TEST"    "BETA"    "SE"
# [12] "CHISQ"   "LOG10P"  "EXTRA"

# Basic stats
cat("Variants tested:", nrow(chr1), "\n")
cat("Mean LOG10P:", round(mean(chr1$LOG10P, na.rm = TRUE), 3), "\n")
cat("Median LOG10P:", round(median(chr1$LOG10P, na.rm = TRUE), 3), "\n")
cat("Max LOG10P:", round(max(chr1$LOG10P, na.rm = TRUE), 3), "\n")
cat("Min p-value:", format(min(chr1$PVAL, na.rm = TRUE), scientific = TRUE), "\n")
cat("Mean INFO:", round(mean(chr1$INFO, na.rm = TRUE), 3), "\n")
cat("Mean A1FREQ:", round(mean(chr1$A1FREQ, na.rm = TRUE), 3), "\n")

# Genomic inflation factor (lambda) for this chromosome
if ("CHISQ" %in% names(chr1)) {
    chisq_median <- median(chr1$CHISQ, na.rm = TRUE)
    lambda <- chisq_median / qchisq(0.5, df = 1)
    cat("Chr", chromosome, "inflation factor (lambda):", round(lambda, 4), "\n")
} else {
    chisq_calc <- qchisq(1 - chr1$PVAL, df = 1)
    lambda <- median(chisq_calc, na.rm = TRUE) / qchisq(0.5, df = 1)
    cat("Chr", chromosome, "inflation factor (lambda):", round(lambda, 4), "\n")
}
# Chr 1 inflation factor (lambda): 1.0841

# Significance counts
genome_wide_sig <- filter(chr2, LOG10P >= 7.3)
suggestive <- filter(chr2, LOG10P >= 5.0)
nominal <- filter(chr1, chr1$LOG10P >= 2.0, na.rm = TRUE)

cat("Genome-wide significant (p < 5e-7.3):", dim(genome_wide_sig), "\n")
cat("Suggestive (p < 1e-5):", dim(suggestive), "\n")
cat("Nominally significant (p < 0.01):", dim(nominal), "\n")

# Position range
pos_range <- range(chr1$GENPOS, na.rm = TRUE)
cat("Position range:", format(pos_range[1], big.mark = ","), "-",
    format(pos_range[2], big.mark = ","), "bp\n")

# ========================================
# CREATE PLOTS
# ========================================

### GENOMIC INFLATION FACTOR ###
pval <- chr2[,"LOG10P"]
chisq <- chr2[,"CHISQ"]

gif <- median(chisq)/qchisq(0.5,1)
# 1.084077 chr1
# 1.0838 chr2

### MANHATTEN PLOT ###
p <- manhattan(chr2, chr="CHROM", bp="GENPOS", snp="ID", p="LOG10P" )
ggsave(here("breast", "20250626_bc_chr2_man.png"), plot = p, width = 8, height = 6, dpi = 300)

## BIG ASS GAP??
pos <- chr1$GENPOS
pos <- sort(pos)
gaps <- diff(pos)
largest_gap <- max(gaps)

# 4. Find the position where the largest gap occurs (i.e., between which positions)
index_of_largest_gap <- which(gaps == largest_gap)

# 5. Get the positions involved in the largest gap
start_position <- pos[index_of_largest_gap]
end_position <- pos[index_of_largest_gap + 1]

# Print the results
cat("The largest gap is between positions", start_position, "and", end_position,
    "with a gap size of", largest_gap, "base pairs.\n")

# The largest gap is between positions 121609248 and 144423246 with a gap size of 22813998 base pairs.


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
    print("INFILE should be a ggman format text file, with two columns: MarkerName (chrpos or rsid) and P-value")
    #	print("y.max, the maximum value for y on the -log10 scale")
    stop()
}

chr1$CHRPOS <- paste(chr1$CHROM, chr1$GENPOS, sep = ":")

INFILE  = chr1[,c("CHRPOS", "LOG10P")]
OUTFILE = "bc_chr1"

#y.max = args[3]

# PARAMS
# the threshold for p-values; sub-sample any values with p > p.thresh
p.thresh = 0.01

# proportion of snps to keep
k.thresh = 0.2


MANOUTFILE <- paste(OUTFILE, "_man.png", sep="")
QQOUTFILE <- paste(OUTFILE, "_qq.png", sep = "")
# ---

print("reading association testing data...")

# dat <- fread(INFILE, header=TRUE)
# dat <- as.data.frame(dat)
# dat <- dat[1:46033,]

# check input validity
# important for consistency with locuszoom
if( colnames(dat)[1] != "MarkerName")
{
    stop("Colnames of INFILE should be 'MarkerName' and 'P-value' ")
}

dat$MarkerName <- as.character(dat$MarkerName)
dat$CHR <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[1]))
dat$BP <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[2]))
dat$BP <- as.integer(dat$BP)
dat$`P-value` <- as.numeric(dat$`P-value`)
if( any(dat$CHR == "X"))
{
    dat$CHR[which(dat$CHR == "X")] <- 23
}

dat$CHR <- as.numeric(dat$CHR)
colnames(dat) <- c("rsid", "p", "chr", "bp")
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











# Manhattan Plot
p_manhattan <- ggplot(chr1, aes(x = GENPOS, y = LOG10P)) +
    geom_point(alpha = 0.7, size = 0.8, color = "darkblue") +
    geom_hline(yintercept = 7.3, color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 5.0, color = "orange", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(labels = function(x) paste0(round(x/1e6, 1), "M")) +
    labs(
        x = paste("Position on Chromosome", chromosome, "(Mb)"),
        y = "-log10(p)",
        title = paste("Manhattan Plot - Chromosome", chromosome, "-", ancestry_group, "ancestry"),
        subtitle = "Red line: genome-wide significance (5×10⁻⁸), Orange line: suggestive (1×10⁻⁵)"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        panel.grid.minor.x = element_blank()
    )

x <- chr1[chr1$LOG10P >= 7.3,]
if (genome_wide_sig > 0) {
    p_manhattan <- p_manhattan +
        geom_point(data = x,
                   aes(x = GENPOS, y = as.numeric(chr1$LOG10P)),
                   color = "red", size = 1.5, alpha = 0.8)
}

y <- chr1[chr1$LOG10P >= 5.0,]
if (suggestive > 0 & genome_wide_sig == 0) {
    p_manhattan <- p_manhattan +
        geom_point(data = y,
                   aes(x = GENPOS, y = as.numeric(chr1$LOG10P)),
                   color = "orange", size = 1.2, alpha = 0.8)
}

ggsave(here("breast", paste0("20250624_", "chr", chromosome, "_manhattan.png")),
       p_manhattan, width = 12, height = 6, dpi = 300)

# 2. QQ Plot for this chromosome
cat("Creating chromosome", chromosome, "QQ plot...\n")

observed <- sort(chr1$LOG10P)
expected <- -log10(ppoints(length(observed)))

qq_data <- data.frame(
    expected = expected,
    observed = observed
)

p_qq <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.7, size = 1, color = "darkblue") +
    geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
    labs(
        x = "Expected -log10(p)",
        y = "Observed -log10(p)",
        title = paste("QQ Plot - Chromosome", chromosome, "-", ancestry_group, "ancestry"),
        subtitle = paste("λ =", round(lambda, 3))
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
    )

ggsave(file.path(plot_dir, paste0("chr", chromosome, "_qq_plot.png")),
       p_qq, width = 8, height = 8, dpi = 300)

# 3. Effect Size Distribution
cat("Creating effect size distribution...\n")
p_beta <- ggplot(chr1, aes(x = BETA)) +
    geom_histogram(bins = 30, fill = "skyblue", alpha = 0.7, color = "black") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    labs(
        x = "Beta (Effect Size)",
        y = "Frequency",
        title = paste("Effect Size Distribution - Chr", chromosome)
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(plot_dir, paste0("chr", chromosome, "_beta_dist.png")),
       p_beta, width = 8, height = 6, dpi = 300)

# 4. INFO Score vs P-value
cat("Creating INFO vs p-value plot...\n")
p_info_p <- ggplot(chr1, aes(x = INFO, y = LOG10P)) +
    geom_point(alpha = 0.5, size = 0.8, color = "darkgreen") +
    geom_smooth(method = "loess", color = "red", se = TRUE) +
    geom_vline(xintercept = 0.3, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 5.0, color = "orange", linetype = "dashed") +
    labs(
        x = "INFO Score",
        y = "-log10(p)",
        title = paste("INFO Score vs P-value - Chr", chromosome)
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(plot_dir, paste0("chr", chromosome, "_info_pval.png")),
       p_info_p, width = 8, height = 6, dpi = 300)

# ========================================
# TOP ASSOCIATIONS
# ========================================

cat("Creating top associations table...\n")

# Top 20 associations for this chromosome
top_hits <- chr1 %>%
    arrange(desc(LOG10P)) %>%
    head(20) %>%
    select(CHROM, GENPOS, ID, A1, A0, A1FREQ, INFO, N, BETA, SE, LOG10P, PVAL) %>%
    mutate(
        PVAL = format(PVAL, scientific = TRUE, digits = 3),
        BETA = round(BETA, 4),
        SE = round(SE, 4),
        LOG10P = round(LOG10P, 2),
        A1FREQ = round(A1FREQ, 3),
        INFO = round(INFO, 3),
        GENPOS_MB = round(GENPOS / 1e6, 2)
    ) %>%
    select(-GENPOS) %>%
    select(CHROM, GENPOS_MB, everything())

# Save top hits
write.table(top_hits, file.path(plot_dir, paste0("chr", chromosome, "_top20.txt")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Display top 10
cat("\nTop 10 associations on chromosome", chromosome, ":\n")
print(head(top_hits, 10))

# ========================================
# SUMMARY REPORT
# ========================================

# Create summary report
report_file <- file.path(plot_dir, paste0("chr", chromosome, "_summary.txt"))

cat("Creating summary report...\n")
writeLines(c(
    paste("=== CHROMOSOME", chromosome, "ANALYSIS SUMMARY ==="),
    paste("Date:", Sys.Date()),
    paste("Ancestry:", ancestry_group),
    paste("Input file:", chr_file),
    "",
    "STATISTICS:",
    paste("- Total variants tested:", nrow(chr1)),
    paste("- Genome-wide significant (p < 5e-8):", genome_wide_sig),
    paste("- Suggestive (p < 1e-5):", suggestive),
    paste("- Nominally significant (p < 0.01):", nominal),
    paste("- Inflation factor (λ):", round(lambda, 4)),
    paste("- Position range:", format(pos_range[1], big.mark = ","), "-",
          format(pos_range[2], big.mark = ",")),
    "",
    "QUALITY METRICS:",
    paste("- Mean INFO score:", round(mean(chr1$INFO, na.rm = TRUE), 3)),
    paste("- Mean MAF:", round(mean(chr1$A1FREQ, na.rm = TRUE), 3)),
    paste("- Min p-value:", format(min(chr1$PVAL, na.rm = TRUE), scientific = TRUE)),
    "",
    "FILES CREATED:",
    paste("- chr", chromosome, "_manhattan.png", sep = ""),
    paste("- chr", chromosome, "_qq_plot.png", sep = ""),
    paste("- chr", chromosome, "_beta_dist.png", sep = ""),
    paste("- chr", chromosome, "_info_pval.png", sep = ""),
    paste("- chr", chromosome, "_top20.txt", sep = "")
), report_file)

# ========================================
# FINAL SUMMARY
# ========================================

cat("\n========================================\n")
cat("CHROMOSOME", chromosome, "ANALYSIS COMPLETE!\n")
cat("========================================\n")
cat("chr1 saved in:", plot_dir, "\n")
cat("Files created:\n")
cat("- Manhattan plot:", paste0("chr", chromosome, "_manhattan.png"), "\n")
cat("- QQ plot:", paste0("chr", chromosome, "_qq_plot.png"), "\n")
cat("- Effect size distribution:", paste0("chr", chromosome, "_beta_dist.png"), "\n")
cat("- INFO vs p-value:", paste0("chr", chromosome, "_info_pval.png"), "\n")
cat("- Top 20 associations:", paste0("chr", chromosome, "_top20.txt"), "\n")
cat("- Summary report:", paste0("chr", chromosome, "_summary.txt"), "\n")
cat("\nKey statistics:\n")
cat("- Variants tested:", nrow(chr1), "\n")
cat("- Inflation factor (λ):", round(lambda, 3), "\n")
cat("- Genome-wide significant:", genome_wide_sig, "\n")
cat("- Suggestive:", suggestive, "\n")

if (genome_wide_sig > 0) {
    cat("\n🎉 Found", genome_wide_sig, "genome-wide significant association(s)!\n")
} else if (suggestive > 0) {
    cat("\n📊 Found", suggestive, "suggestive association(s)\n")
} else {
    cat("\nNo genome-wide or suggestive associations found on this chromosome\n")
}
