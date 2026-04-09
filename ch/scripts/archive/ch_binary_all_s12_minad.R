# ============================================================
# CHIP BINARY ASSOCIATION - MINAD SENSITIVITY ANALYSIS
# ALL SEXES
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table)
library(MatchIt)
library(survival)
library(splines)
library(sandwich)
library(lmtest)
library(pROC)
library(logistf)

# ============================================================
# READ DATA
# ============================================================
m.out4 <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
m.data <- match.data(m.out4)
cov_base <- m.data %>%
    mutate(
        Sample_age   = as.numeric(Sample_age),
        Sample_age2  = Sample_age^2,
        subclass_num = as.numeric(as.factor(subclass))
    )
# 677 groups

# If you saved your matchit object:
summary(m.out4)         # shows match ratio
m.out4$call           # tells you: "nearest", "full", "optimal", etc.
# matchit(formula = match_formula, data = unaff, method = "nearest",
#         distance = "glm", exact = ~Sequenced_gender + Batch + Smoke_History,
#         replace = FALSE, caliper = c(0.2, Sample_age = 5), std.caliper = c(TRUE,
#                                                                            FALSE), ratio = 4)

# Look at the matched data
table(m.data$subclass)  # group sizes tell you the structure

# Count controls per treated unit
m.data %>%
    group_by(subclass, BRCA12_Case) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = BRCA12_Case, values_from = n)
# subclass   `0`   `1`
# <fct>    <int> <int>
#     1 1            4     1
# 2 2            4     1
# 3 3            3     1
# 4 4            3     1
# 5 5            4     1
# 6 6            3     1
# 7 7            4     1
# 8 8            4     1
# 9 9            3     1
# 10 10           2     1
# def variable ratio?

# Count controls per treated unit
matched_data %>%
    group_by(subclass, treatment) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = treatment, values_from = n)

# ============================================================
# FORMULAS
# ============================================================
base_covs_linear <- "Sample_age + Batch + Smoke_History + Sequenced_gender +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

base_covs_spline <- "ns(Sample_age, df = 4) + Batch + Smoke_History + Sequenced_gender +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

base_covs_squared <- "Sample_age + Sample_age2 + Batch + Smoke_History + Sequenced_gender +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ============================================================
# HELPERS
# ============================================================
extract_or_robust <- function(fit, term, data) {
    vcov_cl <- vcovCL(fit, cluster = ~subclass, data = data)
    ct      <- coeftest(fit, vcov = vcov_cl)
    ci      <- coefci(fit, vcov = vcov_cl)
    auc     <- tryCatch(as.numeric(pROC::auc(fit$y, fitted(fit))),
                        error = function(e) NA_real_)

    p_col <- ifelse("Pr(>|z|)" %in% colnames(ct), "Pr(>|z|)", "Pr(>|t|)")  # <-- fix

    data.frame(
        term  = term,
        OR    = exp(ct[term, "Estimate"]),
        CI_lo = exp(ci[term, 1]),
        CI_hi = exp(ci[term, 2]),
        p     = ct[term, p_col],                                             # <-- fix
        AUC   = auc
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

extract_or_clogit <- function(fit, term) {
    s  <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(
        term  = term,
        OR    = exp(s[term, "coef"]),
        CI_lo = exp(ci[term, 1]),
        CI_hi = exp(ci[term, 2]),
        p     = s[term, "Pr(>|z|)"],
        AUC   = NA_real_
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

extract_or_firth <- function(fit, term) {
    ci <- confint(fit)
    idx <- which(names(fit$coefficients) == term)
    data.frame(
        term  = term,
        OR    = exp(fit$coefficients[idx]),
        CI_lo = exp(ci[idx, 1]),
        CI_hi = exp(ci[idx, 2]),
        p     = fit$prob[idx],
        AUC   = NA_real_
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

# ============================================================
# MAIN LOOP
# ============================================================
minad_thresholds <- c(3, 4, 5)

all_results <- list()

for (minad in minad_thresholds) {

    cat("\n============================================================\n")
    cat("Running minAD =", minad, "\n")
    cat("============================================================\n")

    # --- Load vars ---
    vars_file <- file.path("ch", "data",
                           sprintf("ch_seq_wl_art_minad%d_vars.xlsx", minad))
    vars <- read_excel(vars_file)

    # --- Build cov ---
    cov <- cov_base %>%
        mutate(
            CHIP_Binary = person_id %in% vars$Sample.ID,
            CHIP_Count  = {
                chip_counts <- vars %>% count(Sample.ID, name = "CHIP_Count")
                idx <- match(person_id, chip_counts$Sample.ID)
                ifelse(is.na(idx), 0L, chip_counts$CHIP_Count[idx])
            }
        )

    cat("N:", nrow(cov),
        "| Cases:", sum(cov$BRCA12_Case),
        "| Controls:", sum(cov$BRCA12_Case == 0),
        "| CHIP+:", sum(cov$CHIP_Binary), "\n")

    # --- M1: Unweighted quasi-binomial ---
    fit_m1u_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_linear)),
                      data = cov, family = quasibinomial())
    fit_m1u_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +",  base_covs_linear)),
                      data = cov, family = quasibinomial())
    fit_m1u_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +",  base_covs_linear)),
                      data = cov,  family = quasibinomial())

    # --- M1: Weighted quasi-binomial ---
    fit_m1_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_linear)),
                     data = cov, weights = weights, family = quasibinomial())
    fit_m1_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +",  base_covs_linear)),
                     data = cov, weights = weights, family = quasibinomial())
    fit_m1_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +",  base_covs_linear)),
                     data = cov, weights = weights, family = quasibinomial())

    # # --- M1s: Weighted quasi-binomial spline ---
    # fit_m1s_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_spline)),
    #                   data = cov, weights = weights, family = quasibinomial())
    # fit_m1s_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +",  base_covs_spline)),
    #                   data = cov, weights = weights, family = quasibinomial())
    # fit_m1s_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +",  base_covs_spline)),
    #                   data = cov, weights = weights, family = quasibinomial())
    #
    # # --- M1sq: Weighted quasi-binomial squared ---
    # fit_m1sq_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_squared)),
    #                    data = cov, weights = weights, family = quasibinomial())
    # fit_m1sq_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +",  base_covs_squared)),
    #                    data = cov, weights = weights, family = quasibinomial())
    # fit_m1sq_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +",  base_covs_squared)),
    #                    data = cov, weights = weights, family = quasibinomial())

    # age spec comparisons
    # cat("Linear vs spline:\n");  print(anova(fit_m1_12, fit_m1s_12,  test = "F"))
    # cat("Linear vs squared:\n"); print(anova(fit_m1_12, fit_m1sq_12, test = "F"))

    # --- M2: Conditional logistic ---
    fit_m2_12 <- clogit(
        as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age + Sequenced_gender +",
                         "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch +",
                         "strata(subclass_num)")),
        data = cov, weights = weights
    )
    fit_m2_1 <- clogit(
        as.formula(paste("CHIP_Binary ~ BRCA1_Case + Sample_age + Sequenced_gender +",
                         "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch +",
                         "strata(subclass)")),
        data = cov, weights = weights
    )
    fit_m2_2 <- clogit(
        as.formula(paste("CHIP_Binary ~ BRCA2_Case + Sample_age + Sequenced_gender +",
                         "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch +",
                         "strata(subclass_num)")),
        data = cov, weights = weights
    )

    # --- M3: Firth penalized logistic ---
    fit_m3_12 <- logistf(
        as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_linear)),
        data = cov, weights = weights
    )
    fit_m3_1 <- logistf(
        as.formula(paste("CHIP_Binary ~ BRCA1_Case +", base_covs_linear)),
        data = cov, weights = weights
    )
    fit_m3_2 <- logistf(
        as.formula(paste("CHIP_Binary ~ BRCA2_Case +", base_covs_linear)),
        data = cov, weights = weights
    )


    # --- Compile results ---
    results <- bind_rows(
        extract_or_robust(fit_m1u_12, "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1_unweighted_linear"),
        extract_or_robust(fit_m1u_1,  "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1_unweighted_linear"),
        extract_or_robust(fit_m1u_2,  "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1_unweighted_linear"),
        extract_or_robust(fit_m1_12, "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1_linear"),
        extract_or_robust(fit_m1_1,  "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1_linear"),
        extract_or_robust(fit_m1_2,  "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1_linear"),
        # extract_or_robust(fit_m1s_12,  "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1s_spline"),
        # extract_or_robust(fit_m1s_1,   "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1s_spline"),
        # extract_or_robust(fit_m1s_2,   "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1s_spline"),
        # extract_or_robust(fit_m1sq_12, "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1sq_squared"),
        # extract_or_robust(fit_m1sq_1,  "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1sq_squared"),
        # extract_or_robust(fit_m1sq_2,  "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1sq_squared"),
        extract_or_clogit(fit_m2_12, "BRCA12_Case") %>% mutate(model = "BRCA1/2", spec = "M2_conditional"),
        extract_or_clogit(fit_m2_1,  "BRCA1_Case")  %>% mutate(model = "BRCA1",   spec = "M2_conditional"),
        extract_or_clogit(fit_m2_2,  "BRCA2_Case")  %>% mutate(model = "BRCA2",   spec = "M2_conditional"),
        extract_or_firth(fit_m3_12, "BRCA12_Case") %>% mutate(model = "BRCA1/2", spec = "M3_firth"),
        extract_or_firth(fit_m3_1,  "BRCA1_Case")  %>% mutate(model = "BRCA1",   spec = "M3_firth"),
        extract_or_firth(fit_m3_2,  "BRCA2_Case")  %>% mutate(model = "BRCA2",   spec = "M3_firth")
    ) %>%
        mutate(minAD = minad) %>%
        dplyr::select(minAD, spec, model, OR_fmt, OR, CI_lo, CI_hi, p, AUC)

    all_results[[paste0("minAD", minad)]] <- results
}

# ============================================================
# SAVE ALL RESULTS
# ============================================================
write_xlsx(all_results, file.path("ch", "data", "chip_binary_all_s12_minad_sensitivity.xlsx"))

