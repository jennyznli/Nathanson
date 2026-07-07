library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(knitr)
library(kableExtra)
library(ggrepel)
library(qqman)
library(here)


single_variant_file <- here("simplexo", "data", "ALL_step2_genome.regenie")
gene_based_file <- here("simplexo", "data", "ALL_step2_gene_genome.regenie")

# Load single variant data
single_variant <- read_delim(single_variant_file, delim = " ", col_names = TRUE, trim_ws = TRUE) %>%
    mutate(across(c(CHROM, GENPOS, N, N_CASES, N_CONTROLS), as.integer),
           across(c(A1FREQ, A1FREQ_CASES, A1FREQ_CONTROLS, BETA, SE, CHISQ, LOG10P), as.numeric)) %>%
    mutate(CHROM = as.factor(CHROM)) %>%
    arrange(CHROM, GENPOS) %>%
    mutate(PVAL = 10^(-LOG10P),
           PVAL_ADJ = p.adjust(PVAL, method = "BH"),
           LOG10P_ADJ = -log10(PVAL_ADJ)) %>%
    mutate(POS_INDEX = row_number())

# Load gene-based data
gene_based <- read_delim(gene_based_file, delim = " ", col_names = TRUE, trim_ws = TRUE, skip = 1) %>%
    mutate(across(c(CHROM, GENPOS, N), as.integer),
           across(c(A1FREQ, BETA, SE, CHISQ, LOG10P), as.numeric)) %>%
    mutate(CHROM = as.factor(CHROM),
           PVAL = 10^(-LOG10P),
           PVAL_ADJ = p.adjust(PVAL, method = "BH"),
           LOG10P_ADJ = -log10(PVAL_ADJ),
           OR = exp(BETA)) %>%
    # Extract gene and mask information from ID column
    mutate(
        GENE = sub("\\..*", "", ID),
        MASK = sub(".*\\.(M[0-9]+).*", "\\1", ID),
        THRESH = sub(".*\\.(\\d+\\.\\d+)$", "\\1", ID)
    ) %>%
    # Filter to ADD test only for BETA-based analysis
    filter(TEST == "ADD")


sorted_genes <- gene_based %>%
    arrange(PVAL_ADJ, desc(PVAL_ADJ))
# write.csv(sorted_genes, file.path("data", "breast_regenie_gene_sorted.csv"))

### SEE IF ALL GENES ARE IN RESULTS ###
ref_genes <- readLines(file.path(DATA_DIR, "ref_gene.txt"))
hit_genes <- unique(as.factor(gene_based$GENE))

sum(hit_genes %in% ref_genes)
setdiff(ref_genes, hit_genes)
#  [1] "CDR1" - fine     "RAD51D"   "PFKFB1"   "PTEN"     "AR*"      "C1orf168" - fine "BRCC3" - fine
# [8] "KIAA1919" - fine "DIRC2" - fine   "GPR146"   "CDK9"     "PPM1D"    "USP7"

# things to find: GPR146, CDK9, PPM1D, RAD51D, PTEN
# *: AR only has one pathogenic variant (c.2395C>G_p.Gln799Glu) in Simplexo and has freq >0.001 , so calculate OR from this variant for AR only (from gnomAD).

# 96 out of 109 total genes to watch out for..





