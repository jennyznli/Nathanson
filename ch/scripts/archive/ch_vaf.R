# ============================================================
# VAF & clonal characteristics — carrier vs non-carrier
# (restricted to CHIP-positive individuals)
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(coin)

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path(OUT_DIR, "pmbb_brca12_cov_chip_df.xlsx")) %>% filter(Strata == 1)
dim(cov)

vars <- read_excel(file.path(OUT_DIR, "ch_seq_wl_art_minad4_vars.xlsx"))

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

cat("Overall cohort: N =", nrow(cov),
    "| CHIP+ =", sum(cov$CHIP_Binary),
    sprintf("(%.1f%%)\n", 100 * mean(cov$CHIP_Binary)))
# Overall cohort: N = 3004 | CHIP+ = 193 (6.4%)

VAF_COL <- "Sample.AltFrac"   # <-- update if needed

# ============================================================
# RESTRICT TO CHIP+ INDIVIDUALS
# ============================================================
chip_pos <- cov %>% filter(CHIP_Binary == TRUE)
cat("CHIP+ individuals:", nrow(chip_pos), "\n")
cat("  Carriers:     ", sum(chip_pos$BRCA12_Case == 1), "\n")
cat("  Non-carriers: ", sum(chip_pos$BRCA12_Case == 0), "\n")

# Join variant data to CHIP+ individuals
vars_chip <- vars %>%
    inner_join(cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                                     BRCA2_Case, Carrier),
               by = c("Sample.ID" = "person_id"))

# Per-person summary: max VAF and total mutation count
person_summary <- vars_chip %>%
    group_by(Sample.ID, BRCA12_Case, BRCA1_Case, BRCA2_Case, Carrier, Sample_age) %>%
    summarise(
        max_vaf  = max(.data[[VAF_COL]], na.rm = TRUE),
        n_muts   = n(),
        .groups  = "drop"
    )

# ============================================================
# V1: MAX VAF COMPARISON (dominant clone size)
# ============================================================
cat("\n=== V1: Max VAF by carrier status ===\n")
vaf_summary <- person_summary %>%
    group_by(BRCA12_Case) %>%
    summarise(
        n         = n(),
        median_vaf = median(max_vaf),
        q25_vaf   = quantile(max_vaf, 0.25),
        q75_vaf   = quantile(max_vaf, 0.75),
        mean_vaf  = mean(max_vaf),
        .groups   = "drop"
    )
print(vaf_summary)

# Wilcoxon test
wt_vaf <- wilcox.test(max_vaf ~ BRCA12_Case, data = person_summary, exact = FALSE)
cat("Wilcoxon p-value (max VAF, carrier vs non-carrier):", wt_vaf$p.value, "\n")

# Also BRCA1 vs BRCA2 vs non-carrier (3-group Kruskal-Wallis)
kw_vaf <- kruskal.test(max_vaf ~ Carrier, data = person_summary)
cat("Kruskal-Wallis p-value (max VAF, by carrier type):", kw_vaf$p.value, "\n")

# Pairwise Wilcoxon if KW significant
if (kw_vaf$p.value < 0.05) {
    cat("\nPairwise Wilcoxon (Bonferroni corrected):\n")
    pw <- pairwise.wilcox.test(person_summary$max_vaf,
                               person_summary$Carrier,
                               p.adjust.method = "bonferroni")
    print(pw)
}

# ============================================================
# V2: VAF DISTRIBUTION — small vs large clone
# categorize into <2%, 2-10%, >10% (common CHIP thresholds)
# ============================================================
vars_chip <- vars_chip %>%
    mutate(vaf_cat = case_when(
        .data[[VAF_COL]] < 0.02 ~ "<2%",
        .data[[VAF_COL]] < 0.10 ~ "2-10%",
        TRUE                    ~ ">10%"
    ) %>% factor(levels = c("<2%", "2-10%", ">10%")))

cat("\n=== V2: VAF category distribution by carrier status ===\n")
vaf_cat_tab <- table(vars_chip$vaf_cat, vars_chip$BRCA12_Case)
print(vaf_cat_tab)
cat("Chi-square test:\n")
print(chisq.test(vaf_cat_tab))

# ============================================================
# V3: MUTATION BURDEN (n mutations per CHIP+ person)
# ============================================================
cat("\n=== V3: Mutation burden among CHIP+ ===\n")
burden_summary <- person_summary %>%
    group_by(BRCA12_Case) %>%
    summarise(
        n            = n(),
        median_count = median(n_muts),
        q25          = quantile(n_muts, 0.25),
        q75          = quantile(n_muts, 0.75),
        pct_multi    = mean(n_muts > 1) * 100,
        .groups      = "drop"
    )
print(burden_summary)

wt_burden <- wilcox.test(n_muts ~ BRCA12_Case, data = person_summary, exact = FALSE)
cat("Wilcoxon p-value (mutation burden):", wt_burden$p.value, "\n")

# ============================================================
# V4: LINEAR REGRESSION — VAF ~ carrier status + age
# (among CHIP+ individuals, controls for age)
# ============================================================
cat("\n=== V4: Linear regression — max VAF ~ BRCA12_Case + age ===\n")
lm_vaf <- lm(max_vaf ~ BRCA12_Case + Sample_age, data = person_summary)
print(summary(lm_vaf))

lm_vaf_b1 <- lm(max_vaf ~ BRCA1_Case + Sample_age, data = person_summary)
lm_vaf_b2 <- lm(max_vaf ~ BRCA2_Case + Sample_age, data = person_summary)
cat("\nBRCA1 coefficient (p-value):",
    coef(summary(lm_vaf_b1))["BRCA1_Case", "Pr(>|t|)"], "\n")
cat("BRCA2 coefficient (p-value):",
    coef(summary(lm_vaf_b2))["BRCA2_Case", "Pr(>|t|)"], "\n")

# ============================================================
# SAVE
# ============================================================
write_xlsx(
    list(
        person_summary = person_summary,
        vaf_by_group   = vaf_summary,
        vaf_categories = as.data.frame(vaf_cat_tab)
    ),
    file.path("ch", "data", "ch_vaf_clonal_results.xlsx")
)

cat("\nDone: 04_vaf_clonal.R\n")

