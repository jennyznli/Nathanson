# ============================================================
# Gene CHIP analysis — Firth logistic in female, cancer-free
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(data.table)
library(MatchIt)
library(logistf)
library(splines)
library(tidyverse)
library(survival)
library(broom)
library(sandwich)
library(lmtest)
# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "chip_cov_minad4.xlsx"))

# join carrier status onto variant-level data
vars_cov <- vars %>%
    left_join(
        cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                              BRCA2_Case, Carrier, Sample_age,
                              Batch, Smoke_History, starts_with("PC")),
        by = c("Sample.ID" = "person_id")
    ) %>%
    filter(!is.na(BRCA12_Case))

cat("Cohort: N =", nrow(cov), "| CHIP+ =", sum(cov$CHIP_Binary),
    "| Variants =", nrow(vars), "\n")
# Cohort: N = 3004 | CHIP+ = 193 | Variants = 217

# ============================================================
# CONFIG
# ============================================================
MIN_CARRIERS  <- 3

exposure_labels <- c(
    "BRCA12_Case" = "BRCA1/2",
    "BRCA1_Case"  = "BRCA1",
    "BRCA2_Case"  = "BRCA2"
)

# shared ggplot theme
theme_ch <- theme_classic(base_size = 12) +
    theme(
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 10, color = "grey40"),
        axis.text        = element_text(color = "black"),
        legend.position  = "bottom"
    )

# ============================================================
# GENE FREQUENCY TABLE
# ============================================================
gene_person <- vars_cov %>%
    dplyr::select(Sample.ID, Gene, BRCA12_Case) %>%
    distinct()

n_carrier_total    <- sum(cov$BRCA12_Case == 1)
n_noncarrier_total <- sum(cov$BRCA12_Case == 0)

gene_freq <- gene_person %>%
    group_by(Gene, BRCA12_Case) %>%
    summarise(n_mutated = n(), .groups = "drop") %>%
    pivot_wider(names_from  = BRCA12_Case,
                values_from = n_mutated,
                values_fill = 0,
                names_prefix = "n_") %>%
    rename(n_noncarrier = n_0, n_carrier = n_1) %>%
    mutate(
        pct_noncarrier = n_noncarrier / n_noncarrier_total * 100,
        pct_carrier    = n_carrier    / n_carrier_total    * 100,
        n_total        = n_carrier + n_noncarrier
    ) %>%
    arrange(desc(n_total))

dim(gene_freq)
# 46  6

# ============================================================
# FIRTH LOGISTIC PER GENE
# ============================================================
source(here("ch", "scripts", "gene_sensitivity.R"))

# full cohort
res_all <- run_gene_sensitivity(cov, vars, cohort_label = "all")

# females only
res_f <- run_gene_sensitivity(cov, vars, females_only = TRUE)

write_xlsx(res_all, file.path("ch", "data", "ch_all_gene_firth_sensitivity.xlsx"))
write_xlsx(res_f, file.path("ch", "data", "ch_female_gene_firth_sensitivity.xlsx"))

# counts of those significant enough
gene_freq %>% filter(Gene %in% c("GATA1", "ASXL1", "DNMT3A", "TET2", "NF1", "KMT2A"))
# Gene   n_noncarrier n_carrier pct_noncarrier pct_carrier n_total
# <chr>         <int>     <int>          <dbl>       <dbl>   <int>
#     1 DNMT3A           40        17         1.72         2.51       57
# 2 TET2             23        13         0.988        1.92       36
# 3 ASXL1            12         4         0.516        0.591      16
# 4 NF1               6         3         0.258        0.443       9
# 5 GATA1             3         3         0.129        0.443       6
# 6 KMT2A             1         3         0.0430       0.443       4




