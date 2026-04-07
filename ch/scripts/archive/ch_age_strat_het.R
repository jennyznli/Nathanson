# ============================================================
# Age-stratified ORs + BRCA1 vs BRCA2 heterogeneity test
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(tidyverse)
library(writexl)
library(splines)
library(broom)
library(meta)       # for heterogeneity test across strata

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx")) %>%
    filter(Sequenced_gender == "Female", Strata == 1)

vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

base_covs <- "Batch + Smoke_History + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ============================================================
# A1: AGE INTERACTION TEST
# tests whether the BRCA12_Case OR changes with age
# ============================================================
cat("\n=== A1: Age x BRCA12 interaction test ===\n")

fit_interact <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case * Sample_age +", base_covs)),
    data = cov, family = binomial()
)
cat("Interaction term (BRCA12_Case:Sample_age):\n")
print(tidy(fit_interact, conf.int = TRUE, exponentiate = TRUE) %>%
          filter(grepl("BRCA12_Case", term)))

# Likelihood ratio test for interaction
fit_no_interact <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age +", base_covs)),
    data = cov, family = binomial()
)
lrt <- anova(fit_no_interact, fit_interact, test = "LRT")
cat("\nLRT p-value for age interaction:", lrt$`Pr(>Chi)`[2], "\n")

# ============================================================
# A2: STRATIFIED ORs BY AGE DECADE
# ============================================================
cat("\n=== A2: Stratified ORs by age decade ===\n")

cov <- cov %>%
    mutate(age_decade = cut(Sample_age,
                            breaks = c(0, 50, 60, 70, Inf),
                            labels = c("<50", "50-59", "60-69", "70+"),
                            right  = FALSE))

cat("N per age stratum:\n")
print(table(cov$age_decade, cov$BRCA12_Case))

fit_age_stratum <- function(stratum, df) {
    df_s <- df %>% filter(age_decade == stratum)
    if (sum(df_s$BRCA12_Case == 1) < 5 | sum(df_s$CHIP_Binary) < 5) {
        return(data.frame(age_decade = stratum, OR = NA, CI_lo = NA,
                          CI_hi = NA, p = NA, n = nrow(df_s)))
    }
    fit <- tryCatch(
        glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs)),
            data = df_s, family = binomial()),
        error = function(e) NULL
    )
    if (is.null(fit)) return(data.frame(age_decade = stratum, OR = NA,
                                        CI_lo = NA, CI_hi = NA, p = NA,
                                        n = nrow(df_s)))
    res <- tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term == "BRCA12_Case")
    data.frame(
        age_decade = stratum,
        n          = nrow(df_s),
        n_chip     = sum(df_s$CHIP_Binary),
        OR         = res$estimate,
        CI_lo      = res$conf.low,
        CI_hi      = res$conf.high,
        p          = res$p.value,
        OR_fmt     = sprintf("%.2f (%.2f\u2013%.2f)", res$estimate,
                             res$conf.low, res$conf.high)
    )
}

age_strata_results <- map_dfr(levels(cov$age_decade), fit_age_stratum, df = cov)
cat("\nAge-stratified ORs:\n")
print(age_strata_results)

# ============================================================
# A3: BRCA1 vs BRCA2 HETEROGENEITY TEST
# LRT comparing separate-gene model vs combined model
# ============================================================
cat("\n=== A3: BRCA1 vs BRCA2 heterogeneity test ===\n")

# Model with separate BRCA1 and BRCA2 terms (allows different ORs)
fit_separate <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case + BRCA2_Case + Sample_age +", base_covs)),
    data = cov, family = binomial()
)

# Model constraining BRCA1 = BRCA2 (combined term)
fit_combined <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age +", base_covs)),
    data = cov, family = binomial()
)

lrt_het <- anova(fit_combined, fit_separate, test = "LRT")
cat("LRT p-value (BRCA1 vs BRCA2 heterogeneity):", lrt_het$`Pr(>Chi)`[2], "\n")
cat("  p < 0.05 → BRCA1 and BRCA2 have significantly different ORs\n")
cat("  p >= 0.05 → No evidence of heterogeneity; combined model preferred\n\n")

# Individual ORs from separate model
cat("BRCA1 OR (from separate model):\n")
print(tidy(fit_separate, conf.int = TRUE, exponentiate = TRUE) %>%
          filter(term == "BRCA1_Case") %>%
          dplyr::select(term, estimate, conf.low, conf.high, p.value))
cat("BRCA2 OR (from separate model):\n")
print(tidy(fit_separate, conf.int = TRUE, exponentiate = TRUE) %>%
          filter(term == "BRCA2_Case") %>%
          dplyr::select(term, estimate, conf.low, conf.high, p.value))

# ============================================================
# A4: CANCER-NAIVE SENSITIVITY ANALYSIS
# restrict to individuals with blood draw >= 1 year before
# any cancer diagnosis (or no cancer ever)
# ============================================================
cat("\n=== A4: Cancer-naive sensitivity analysis ===\n")

# NOTE: adjust column names to match your data
# Expected: Cancer_Dx_Age (age at cancer diagnosis), Sample_age (age at blood draw)
# If no cancer: Cancer_Dx_Age should be NA

if ("Cancer_Dx_Age" %in% colnames(cov)) {
    cov_naive <- cov %>%
        filter(is.na(Cancer_Dx_Age) |
                   (Cancer_Dx_Age - Sample_age) >= 1)

    cat("Cancer-naive N:", nrow(cov_naive), "\n")
    cat("  Carriers:    ", sum(cov_naive$BRCA12_Case == 1), "\n")
    cat("  Non-carriers:", sum(cov_naive$BRCA12_Case == 0), "\n")

    fit_naive <- glm(
        as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age +", base_covs)),
        data = cov_naive, family = binomial()
    )
    cat("\nCancer-naive OR:\n")
    print(tidy(fit_naive, conf.int = TRUE, exponentiate = TRUE) %>%
              filter(term == "BRCA12_Case") %>%
              dplyr::select(term, estimate, conf.low, conf.high, p.value))
} else {
    cat("Cancer_Dx_Age column not found — skip or rename as needed\n")
}

# ============================================================
# SAVE
# ============================================================
write_xlsx(
    list(
        age_interaction  = tidy(fit_interact, conf.int = TRUE, exponentiate = TRUE),
        age_strata       = age_strata_results,
        het_test_lrt     = as.data.frame(lrt_het),
        brca_separate    = tidy(fit_separate, conf.int = TRUE, exponentiate = TRUE)
    ),
    file.path("ch", "data", "ch_age_stratified_het_results.xlsx")
)

cat("\nDone: 05_age_stratified_het.R\n")
