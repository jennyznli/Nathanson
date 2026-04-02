# ============================================================
# SCRIPT 2: CHIP COUNT OUTCOME — POISSON / NEGATIVE BINOMIAL
# BRCA1/2 Carrier CHIP Association Analysis
# ============================================================
# Outcome: number of CHIP variants per person (0, 1, 2, ...)
# Models:
#   C1 — Weighted Poisson (test overdispersion)
#   C2 — Negative binomial (preferred if overdispersed)
#   C3 — Zero-inflated NB (if excess zeros)
#
# NOTE: glm.nb() does not support frequency weights.
# Use unweighted NB on matched sample — this is acceptable because
# the matched sample already has near-equal covariate distributions.
# Sensitivity: re-run with IPTW weights using pscl::zeroinfl or VGAM.
# ============================================================

library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(data.table)
library(MatchIt)
library(MASS)         # glm.nb
library(pscl)         # zeroinfl, vuong
library(splines)
library(sandwich)
library(lmtest)
library(tidyverse)
library(writexl)

# ============================================================
# READ DATA
# ============================================================
m.out4 <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
m.data <- match.data(m.out4)

cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx")) %>%
    left_join(
        m.data %>% select(person_id, distance, weights, subclass),
        by = "person_id"
    ) %>%
    filter(!is.na(weights))   # keep only matched individuals

cat("Matched N:", nrow(cov), "\n")
cat("CHIP count distribution:\n")
print(table(cov$CHIP_Count))
cat("Zero fraction:", mean(cov$CHIP_Count == 0), "\n")
cat("Max count:", max(cov$CHIP_Count), "\n")

# ============================================================
# FORMULA
# ============================================================
base_covs_spline <- "ns(Sample_age, df = 4) + Sequenced_gender + Batch +
                      Smoke_History +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ============================================================
# C1: WEIGHTED POISSON — test for overdispersion
# ============================================================
cat("\n--- C1: Weighted Poisson ---\n")
fit_pois_12 <- glm(as.formula(paste("CHIP_Count ~ BRCA12_Case +", base_covs_spline)),
                   data = cov, weights = weights, family = poisson())
fit_pois_1  <- glm(as.formula(paste("CHIP_Count ~ BRCA1_Case +",  base_covs_spline)),
                   data = cov, weights = weights, family = poisson())
fit_pois_2  <- glm(as.formula(paste("CHIP_Count ~ BRCA2_Case +",  base_covs_spline)),
                   data = cov, weights = weights, family = poisson())

# Overdispersion: residual deviance / df should be ~1
cat("Overdispersion ratio (BRCA1/2 Poisson):",
    round(fit_pois_12$deviance / fit_pois_12$df.residual, 2), "\n")
# > 1 → use NB; ~1 → Poisson fine

# ============================================================
# C2: NEGATIVE BINOMIAL — unweighted on matched sample
# ============================================================
cat("\n--- C2: Negative binomial ---\n")
fit_nb_12 <- glm.nb(as.formula(paste("CHIP_Count ~ BRCA12_Case +", base_covs_spline)), data = cov)
fit_nb_1  <- glm.nb(as.formula(paste("CHIP_Count ~ BRCA1_Case +",  base_covs_spline)), data = cov)
fit_nb_2  <- glm.nb(as.formula(paste("CHIP_Count ~ BRCA2_Case +",  base_covs_spline)), data = cov)

# LRT: Poisson vs NB (NB nests Poisson as theta → ∞)
# Use unweighted Poisson for fair comparison with NB
fit_pois_unw <- glm(as.formula(paste("CHIP_Count ~ BRCA12_Case +", base_covs_spline)),
                    data = cov, family = poisson())
lrt_p <- pchisq(2 * (logLik(fit_nb_12) - logLik(fit_pois_unw)), df = 1, lower.tail = FALSE)
cat("Poisson vs NB LRT p-value:", lrt_p, "\n")
# Significant → prefer NB

# ============================================================
# C3: ZERO-INFLATED NB — if zero fraction is high (>60%)
# ============================================================
zero_frac <- mean(cov$CHIP_Count == 0)
cat("\nZero fraction:", round(zero_frac, 3), "\n")

if (zero_frac > 0.60) {
    cat("High zero fraction — fitting ZINB\n")
    fit_zinb_12 <- zeroinfl(
        as.formula(paste("CHIP_Count ~ BRCA12_Case +", base_covs_spline,
                         "| BRCA12_Case + Sample_age + Batch")),
        data = cov, dist = "negbin"
    )
    # Vuong test: ZINB vs NB
    cat("\n=== Vuong test: ZINB vs NB ===\n")
    print(vuong(fit_zinb_12, fit_nb_12))
} else {
    cat("Zero fraction OK — NB preferred\n")
    fit_zinb_12 <- NULL
}

# ============================================================
# EXTRACT IRR (incidence rate ratio)
# ============================================================
extract_irr <- function(fit, term) {
    s  <- coef(summary(fit))
    ci <- confint(fit)   # profile CI for NB
    data.frame(
        term    = term,
        IRR     = exp(s[term, "Estimate"]),
        CI_lo   = exp(ci[term, 1]),
        CI_hi   = exp(ci[term, 2]),
        p       = s[term, grep("Pr", colnames(s))],
        IRR_fmt = sprintf("%.2f (%.2f\u2013%.2f)",
                          exp(s[term, "Estimate"]),
                          exp(ci[term, 1]),
                          exp(ci[term, 2]))
    )
}

# ============================================================
# RESULTS TABLE
# ============================================================
count_results <- bind_rows(
    extract_irr(fit_nb_12, "BRCA12_Case") %>% mutate(model = "BRCA1/2", spec = "NB"),
    extract_irr(fit_nb_1,  "BRCA1_Case")  %>% mutate(model = "BRCA1",   spec = "NB"),
    extract_irr(fit_nb_2,  "BRCA2_Case")  %>% mutate(model = "BRCA2",   spec = "NB")
) %>%
    dplyr::select(spec, model, IRR_fmt, IRR, CI_lo, CI_hi, p)

cat("\n=== COUNT MODEL RESULTS (IRR) ===\n")
print(count_results)

# NB dispersion parameter (theta): lower = more overdispersion
cat("\nNB theta (BRCA1/2):", fit_nb_12$theta, "\n")

# ============================================================
# PREDICTED COUNT PLOT DATA (for figures)
# ============================================================
# Marginal predicted count by age and carrier status
age_grid <- seq(min(cov$Sample_age), max(cov$Sample_age), length.out = 100)

pred_grid <- expand.grid(
    Sample_age       = age_grid,
    BRCA12_Case      = c(0, 1),
    Sequenced_gender = "Female",
    Batch            = factor(2, levels = levels(factor(cov$Batch))),
    Smoke_History    = 0
) %>%
    dplyr::mutate(across(paste0("PC", 1:6), ~ 0))

pred_grid$pred_count <- predict(fit_nb_12, newdata = pred_grid, type = "response")

write.csv(pred_grid, file.path("ch", "data", "ch_nb_pred_grid.csv"), row.names = FALSE)

# ============================================================
# SAVE
# ============================================================
write_xlsx(count_results, file.path("ch", "data", "ch_chip_count_results.xlsx"))
saveRDS(list(
    pois  = fit_pois_12,
    nb    = list(b12 = fit_nb_12, b1 = fit_nb_1, b2 = fit_nb_2),
    zinb  = fit_zinb_12,
    data  = cov
), file.path("ch", "data", "ch_count_model_fits.rds"))

cat("\nDone: 02_chip_count.R\n")
