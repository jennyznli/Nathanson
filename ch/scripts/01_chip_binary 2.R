# ============================================================
# SCRIPT 1: CHIP BINARY OUTCOME — LOGISTIC MODELS
# BRCA1/2 Carrier CHIP Association Analysis - only for females...
# ============================================================
# Models implemented:
#   M1  — Weighted quasi-binomial, linear age
#   M1s — Weighted quasi-binomial, natural spline age (df=4)  [primary]
#   M2  — Conditional logistic within matched sets (cluster SE)
#   M3  — Firth logistic (separation/rare-event sensitivity)
#   M4  — Interaction model: Carrier × Age, Carrier × Sex
#   M5  — Multinomial (BRCA1 / BRCA2 / control)
#
# Weights: ATT weights from MatchIt match.data()
# Cluster SE: on subclass (matched set) via sandwich/lmtest
# ============================================================

library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(data.table)
library(MatchIt)
library(survival)     # clogit
library(logistf)      # Firth
library(nnet)         # multinom
library(splines)      # ns()
library(sandwich)     # vcovCL
library(lmtest)       # coeftest
library(pROC)
library(tidyverse)
library(writexl)

# ============================================================
# READ DATA
# ============================================================
m.out4   <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
m.data   <- match.data(m.out4)   # includes 'weights' and 'subclass'

cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx")) %>%
    left_join(
        m.data %>% dplyr::select(person_id, distance, weights, subclass),
        by = "person_id"
    ) %>%
    filter(!is.na(weights), Strata == 1, Sequenced_gender == "Female")
cat("Matched N:", nrow(cov), "\n")
# 1924
cat("Cases:", sum(cov$BRCA12_Case), "| Controls:", sum(cov$BRCA12_Case == 0), "\n")
# 405, 1519
cat("CHIP positive:", sum(cov$CHIP_Binary), "\n")
# 149

# ============================================================
# FORMULA COMPONENTS
# ============================================================
# NOTE: distance (propensity score) dropped from regression covariates —
# it was used for matching; adjusting for it post-hoc is redundant and
# can introduce bias (double-adjustment). We keep Sex/Batch/Smoking/PCs
# as they were exactly-matched on Sex/Batch/Smoking and PS-matched on PCs.
base_covs_linear <- "Sample_age + Sequenced_gender + Batch +
                      Smoke_History +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# natural spline for the variable Sample_age
base_covs_spline <- "ns(Sample_age, df = 4) + Sequenced_gender + Batch +
                      Smoke_History +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ============================================================
# HELPER: CLUSTER-ROBUST SE FOR WEIGHTED GLM
# ============================================================
# Uses sandwich::vcovCL clustered on matched set (subclass).
# This properly accounts for within-set correlation from 1:k matching.
extract_or_robust <- function(fit, term, data) {
    vcov_cl <- vcovCL(fit, cluster = ~subclass, data = data)
    ct      <- coeftest(fit, vcov = vcov_cl)
    ci      <- coefci(fit, vcov = vcov_cl)

    auc <- tryCatch(as.numeric(pROC::auc(fit$y, fitted(fit))), error = function(e) NA_real_)

    data.frame(
        term   = term,
        OR     = exp(ct[term, "Estimate"]),
        CI_lo  = exp(ci[term, 1]),
        CI_hi  = exp(ci[term, 2]),
        p      = ct[term, "Pr(>|z|)"],
        AUC    = auc
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

extract_or_firth <- function(fit, term) {
    idx <- which(names(coef(fit)) == term)
    data.frame(
        term   = term,
        OR     = exp(coef(fit)[idx]),
        CI_lo  = exp(fit$ci.lower[idx]),
        CI_hi  = exp(fit$ci.upper[idx]),
        p      = fit$prob[idx],
        AUC    = tryCatch(as.numeric(pROC::auc(fit$y, fit$predict)), error = function(e) NA_real_)
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

# ============================================================
# M1: WEIGHTED QUASI-BINOMIAL — LINEAR AGE (baseline)
# ============================================================
cat("\n--- M1: Linear age ---\n")
fit_m1_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_linear)),
                 data = cov, weights = weights, family = quasibinomial())
fit_m1_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +", base_covs_linear)),
                 data = cov, weights = weights, family = quasibinomial())
fit_m1_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +", base_covs_linear)),
                 data = cov, weights = weights, family = quasibinomial())

# ============================================================
# M1s: WEIGHTED QUASI-BINOMIAL — SPLINE AGE (primary model)
# ============================================================
# Preferred: spline captures nonlinear age-CHIP relationship
cat("\n--- M1s: Spline age (df=4) [PRIMARY] ---\n")
fit_m1s_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_spline)),
                  data = cov, weights = weights, family = quasibinomial())
fit_m1s_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +", base_covs_spline)),
                  data = cov, weights = weights, family = quasibinomial())
fit_m1s_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +", base_covs_spline)),
                  data = cov, weights = weights, family = quasibinomial())

# LRT: does spline improve over linear? (use F-test for quasibinomial)
cat("\n=== LRT: linear vs spline (BRCA1/2) ===\n")
print(anova(fit_m1_12, fit_m1s_12, test = "F"))

# ============================================================
# M2: CONDITIONAL LOGISTIC — within matched sets
# ============================================================
# Conditions on matched set (subclass), giving conditional OR.
# No weights needed; the conditioning handles the matched design.
# For 1:k unequal sets use clogit with strata(subclass).
cat("\n--- M2: Conditional logistic ---\n")

# clogit needs numeric subclass
cov <- cov %>% mutate(subclass_num = as.numeric(as.factor(subclass)))

fit_m2_12 <- clogit(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +",
                     "strata(subclass_num)")),
    data = cov, method = "efron"
)
fit_m2_1 <- clogit(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case + Sample_age +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +",
                     "strata(subclass_num)")),
    data = cov, method = "efron"
)
fit_m2_2 <- clogit(
    as.formula(paste("CHIP_Binary ~ BRCA2_Case + Sample_age +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +",
                     "strata(subclass_num)")),
    data = cov, method = "efron"
)

extract_or_clogit <- function(fit, term) {
    s  <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(
        term   = term,
        OR     = exp(s[term, "coef"]),
        CI_lo  = exp(ci[term, 1]),
        CI_hi  = exp(ci[term, 2]),
        p      = s[term, "Pr(>|z|)"],
        AUC    = NA_real_
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

# ============================================================
# M3: FIRTH LOGISTIC — separation/rare-event sensitivity
# ============================================================
cat("\n--- M3: Firth logistic ---\n")
# NOTE: logistf does not support weights — run on matched sample unweighted.
# This is a sensitivity check; results should broadly agree with M1s.
firth_12 <- logistf(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_spline)), data = cov)
firth_1  <- logistf(as.formula(paste("CHIP_Binary ~ BRCA1_Case +",  base_covs_spline)), data = cov)
firth_2  <- logistf(as.formula(paste("CHIP_Binary ~ BRCA2_Case +",  base_covs_spline)), data = cov)

# ============================================================
# M4: INTERACTION MODEL — Carrier × Age, Carrier × Sex
# ============================================================
cat("\n--- M4: Interaction models ---\n")
fit_int_age <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case * ns(Sample_age, df = 4) +",
                     "Sequenced_gender + Batch + Smoke_History +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6")),
    data = cov, weights = weights, family = quasibinomial()
)
fit_int_sex <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case * Sequenced_gender +",
                     base_covs_spline)),
    data = cov, weights = weights, family = quasibinomial()
)

# Joint test of interaction terms (F-test)
cat("\n=== Interaction: Carrier × Age ===\n")
print(anova(fit_m1s_12, fit_int_age, test = "F"))
# Analysis of Deviance Table

# Model 1: CHIP_Binary ~ BRCA12_Case + ns(Sample_age, df = 4) + Sequenced_gender +
#     Batch + Smoke_History + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
# Model 2: CHIP_Binary ~ BRCA12_Case * ns(Sample_age, df = 4) + Sequenced_gender +
#     Batch + Smoke_History + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
# Resid. Df Resid. Dev Df Deviance      F Pr(>F)
# 1      2989     1519.0
# 2      2985     1513.7  4   5.2836 1.3185 0.2606

cat("\n=== Interaction: Carrier × Sex ===\n")
print(anova(fit_m1s_12, fit_int_sex, test = "F"))

# Model 1: CHIP_Binary ~ BRCA12_Case + ns(Sample_age, df = 4) + Sequenced_gender +
#     Batch + Smoke_History + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
# Model 2: CHIP_Binary ~ BRCA12_Case * Sequenced_gender + ns(Sample_age,
#                                                            df = 4) + Sequenced_gender + Batch + Smoke_History + PC1 +
#     PC2 + PC3 + PC4 + PC5 + PC6
# Resid. Df Resid. Dev Df Deviance      F Pr(>F)
# 1      2989     1519.0
# 2      2988     1517.8  1   1.2059 1.2113 0.2712

# ============================================================
# M5: MULTINOMIAL — BRCA1 vs BRCA2 vs Control
# ============================================================
cat("\n--- M5: Multinomial logistic ---\n")
cov <- cov %>%
    mutate(Carrier3 = case_when(
        BRCA1_Case == 1 & BRCA2_Case == 0 ~ "BRCA1",
        BRCA2_Case == 1 & BRCA1_Case == 0 ~ "BRCA2",
        TRUE ~ "Control"
    )) %>%
    mutate(Carrier3 = relevel(factor(Carrier3), ref = "Control"))

fit_m5 <- multinom(
    as.formula(paste("CHIP_Binary ~ Carrier3 +", base_covs_spline)),
    data = cov, weights = weights, trace = FALSE
)
summary(fit_m5)

# ============================================================
# MODEL COMPARISON TABLE
# ============================================================
cat("\n=== MODEL COMPARISON ===\n")

results_binary <- bind_rows(
    # M1 linear
    extract_or_robust(fit_m1_12,  "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1_linear"),
    extract_or_robust(fit_m1_1,   "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1_linear"),
    extract_or_robust(fit_m1_2,   "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1_linear"),
    # M1s spline (PRIMARY)
    extract_or_robust(fit_m1s_12, "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1s_spline_PRIMARY"),
    extract_or_robust(fit_m1s_1,  "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1s_spline_PRIMARY"),
    extract_or_robust(fit_m1s_2,  "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1s_spline_PRIMARY"),
    # M2 conditional
    extract_or_clogit(fit_m2_12, "BRCA12_Case") %>% mutate(model = "BRCA1/2", spec = "M2_conditional"),
    extract_or_clogit(fit_m2_1,  "BRCA1_Case")  %>% mutate(model = "BRCA1",   spec = "M2_conditional"),
    extract_or_clogit(fit_m2_2,  "BRCA2_Case")  %>% mutate(model = "BRCA2",   spec = "M2_conditional"),
    # M3 Firth
    extract_or_firth(firth_12, "BRCA12_Case") %>% mutate(model = "BRCA1/2", spec = "M3_firth"),
    extract_or_firth(firth_1,  "BRCA1_Case")  %>% mutate(model = "BRCA1",   spec = "M3_firth"),
    extract_or_firth(firth_2,  "BRCA2_Case")  %>% mutate(model = "BRCA2",   spec = "M3_firth")
) %>%
    dplyr::select(spec, model, OR_fmt, OR, CI_lo, CI_hi, p, AUC) %>%
    arrange(model, spec)

print(results_binary)

# ============================================================
# MODEL dplyr::selectION GUIDANCE
# ============================================================
# Primary = M1s (weighted quasi-binomial + spline + cluster-robust SE):
#   - Accounts for matching weights (ATT estimand)
#   - Handles nonlinear age-CHIP relationship
#   - Cluster SE on subclass corrects for within-set correlation
# M2 (clogit) as supporting evidence for conditional OR
# M3 (Firth) as sensitivity for separation
# If M1s ≈ M2 ≈ M3 → results are robust

# AUC comparison
cat("\n=== AUC comparison (M1 linear vs M1s spline) ===\n")
roc1  <- roc(cov$CHIP_Binary, fitted(fit_m1_12),  quiet = TRUE)
roc1s <- roc(cov$CHIP_Binary, fitted(fit_m1s_12), quiet = TRUE)
print(roc.test(roc1, roc1s))

# ============================================================
# SAVE
# ============================================================
write_xlsx(results_binary, file.path("ch", "data", "ch_chip_binary_results.xlsx"))
saveRDS(list(
    m1_linear = list(b12 = fit_m1_12, b1 = fit_m1_1, b2 = fit_m1_2),
    m1s_spline = list(b12 = fit_m1s_12, b1 = fit_m1s_1, b2 = fit_m1s_2),
    m2_clogit  = list(b12 = fit_m2_12, b1 = fit_m2_1, b2 = fit_m2_2),
    m3_firth   = list(b12 = firth_12, b1 = firth_1, b2 = firth_2),
    m4_int     = list(age = fit_int_age, sex = fit_int_sex),
    m5_multinom = fit_m5,
    data = cov
), file.path("ch", "data", "ch_binary_model_fits.rds"))

cat("\nDone: 01_chip_binary.R\n")
