# ============================================================
# CHIP binary association analysis
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
library(tidyverse)

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx")) %>%
    filter(Sequenced_gender == "Female", Strata == 1)
dim(cov)
# 2203

vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx")) %>%
    filter(Sample.ID %in% cov$person_id)
dim(vars)
# 149

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

m.out4 <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
m.data <- match.data(m.out4)

# ============================================================
# JOIN COVs
# ============================================================
cov <- cov %>%
    left_join(
        m.data %>% dplyr::select(person_id, distance, weights, subclass),
        by = "person_id"
    ) %>%
    filter(!is.na(weights))

cov$Sample_age <- as.numeric(cov$Sample_age)
cov <- cov %>% mutate(Sample_age2 = Sample_age^2)

cat("Matched N:", nrow(cov), "\n")
# 2203
cat("Cases:", sum(cov$BRCA12_Case), "| Controls:", sum(cov$BRCA12_Case == 0), "\n")
# 516, 1687
cat("CHIP positive:", sum(cov$CHIP_Binary), "\n")
# 132

# ============================================================
# FORMULAS
# ============================================================
# linear age
base_covs_linear <- "Sample_age + Batch + Smoke_History +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# natural spline
base_covs_spline <- "ns(Sample_age, df = 4) + Batch + Smoke_History +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# age squared
base_covs_squared <- "Sample_age + Sample_age2 + Batch + Smoke_History +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ============================================================
# HELPER: CLUSTER-ROBUST SE FOR WEIGHTED GLM
# ============================================================
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

# ============================================================
# M1: WEIGHTED QUASI-BINOMIAL — LINEAR AGE
# ============================================================
fit_m1_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_linear)),
                 data = cov, weights = weights, family = quasibinomial())
fit_m1_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +", base_covs_linear)),
                 data = cov, weights = weights, family = quasibinomial())
fit_m1_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +", base_covs_linear)),
                 data = cov, weights = weights, family = quasibinomial())

# ============================================================
# M1s: WEIGHTED QUASI-BINOMIAL — SPLINE AGE
# ============================================================
fit_m1s_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_spline)),
                  data = cov, weights = weights, family = quasibinomial())
fit_m1s_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +", base_covs_spline)),
                  data = cov, weights = weights, family = quasibinomial())
fit_m1s_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +", base_covs_spline)),
                  data = cov, weights = weights, family = quasibinomial())

# ============================================================
# M1sq: WEIGHTED QUASI-BINOMIAL — SQUARED AGE
# ============================================================
fit_m1sq_12 <- glm(as.formula(paste("CHIP_Binary ~ BRCA12_Case +", base_covs_squared)),
                  data = cov, weights = weights, family = quasibinomial())
fit_m1sq_1  <- glm(as.formula(paste("CHIP_Binary ~ BRCA1_Case +", base_covs_squared)),
                  data = cov, weights = weights, family = quasibinomial())
fit_m1sq_2  <- glm(as.formula(paste("CHIP_Binary ~ BRCA2_Case +", base_covs_squared)),
                  data = cov, weights = weights, family = quasibinomial())

# ============================================================
# M1 - AGE METHO COMPARISONS
# ============================================================
# linear vs spline
print(anova(fit_m1_12, fit_m1s_12, test = "F"))
# 0.782

# linear vs squared
print(anova(fit_m1_12, fit_m1sq_12, test = "F"))
# 0.311

# spline vs squared
print(anova(fit_m1s_12, fit_m1sq_12, test = "F"))
# 0.97

# so linear is fine..

# ============================================================
# M2: CONDITIONAL LOGISTIC
# ============================================================
cov <- cov %>% mutate(subclass_num = as.numeric(as.factor(subclass)))

fit_m2_12 <- clogit(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch + ",
                     "strata(subclass_num)")),
    data = cov, method = "efron"
)
fit_m2_1 <- clogit(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case + Sample_age +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch + ",
                     "strata(subclass_num)")),
    data = cov, method = "efron"
)
fit_m2_2 <- clogit(
    as.formula(paste("CHIP_Binary ~ BRCA2_Case + Sample_age +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch +",
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
# M3: INTERACTION MODEL - linear
# ============================================================
fit_int_12 <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case * Sample_age +",
                     "Batch + Smoke_History +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6")),
    data = cov, weights = weights, family = quasibinomial()
)

fit_int_1 <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case * Sample_age +",
                     "Batch + Smoke_History +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6")),
    data = cov, weights = weights, family = quasibinomial()
)

fit_int_2 <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA2_Case * Sample_age +",
                     "Batch + Smoke_History +",
                     "PC1 + PC2 + PC3 + PC4 + PC5 + PC6")),
    data = cov, weights = weights, family = quasibinomial()
)

### COMPARE WITH PREV MODELS
print(anova(fit_m1_12, fit_int_12, test = "F"))
# 0.361

print(anova(fit_m1s_12, fit_int_12, test = "F"))
# 0.878

print(anova(fit_m1sq_12, fit_int_12, test = "F"))
# -0.219

# ============================================================
# MODEL COMPARISON TABLE
# ============================================================
results_binary <- bind_rows(
    # M1 linear
    extract_or_robust(fit_m1_12,  "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1_linear"),
    extract_or_robust(fit_m1_1,   "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1_linear"),
    extract_or_robust(fit_m1_2,   "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1_linear"),
    # M1s spline
    extract_or_robust(fit_m1s_12, "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1s_spline"),
    extract_or_robust(fit_m1s_1,  "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1s_spline"),
    extract_or_robust(fit_m1s_2,  "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1s_spline"),
    # M1s squared
    extract_or_robust(fit_m1sq_12, "BRCA12_Case", cov) %>% mutate(model = "BRCA1/2", spec = "M1s_squared"),
    extract_or_robust(fit_m1sq_1,  "BRCA1_Case",  cov) %>% mutate(model = "BRCA1",   spec = "M1s_squared"),
    extract_or_robust(fit_m1sq_2,  "BRCA2_Case",  cov) %>% mutate(model = "BRCA2",   spec = "M1s_squared"),
    # M2 conditional
    extract_or_clogit(fit_m2_12, "BRCA12_Case") %>% mutate(model = "BRCA1/2", spec = "M2_conditional"),
    extract_or_clogit(fit_m2_1,  "BRCA1_Case")  %>% mutate(model = "BRCA1",   spec = "M2_conditional"),
    extract_or_clogit(fit_m2_2,  "BRCA2_Case")  %>% mutate(model = "BRCA2",   spec = "M2_conditional"),
    # M3 age * carrier interaction
    extract_or_robust(fit_int_12,  "BRCA12_Case", cov)  %>% mutate(model = "BRCA1/2",   spec = "M1_interaction"),
    extract_or_robust(fit_int_1,  "BRCA1_Case", cov)  %>% mutate(model = "BRCA1",   spec = "M1_interaction"),
    extract_or_robust(fit_int_2,  "BRCA2_Case", cov)  %>% mutate(model = "BRCA2",   spec = "M1_interaction")
) %>%
    dplyr::select(spec, model, OR_fmt, OR, CI_lo, CI_hi, p, AUC) %>%
    arrange(model, spec)

write_xlsx(results_binary, file.path("ch", "data", "ch_chip_binary_s1_results.xlsx"))

# ============================================================
# FOREST PLOT
# ============================================================
# ── Data: pull directly from primary model results ──
forest_df <- results_binary %>%
    filter(spec == "M2_conditional") %>%
    mutate(
        label = factor(model, levels = rev(c("BRCA1/2", "BRCA1", "BRCA2"))),
        p_fmt = case_when(
            p < 0.001 ~ "p<0.001",
            TRUE      ~ sprintf("p=%.3f", p)
        ),
        annotation = sprintf("%.2f (%.2f, %.2f), %s", OR, CI_lo, CI_hi, p_fmt)
    )

x_max    <-  4     # right edge of the plot panel
x_label  <- 5     # where the annotation text starts (in log10 space, beyond x_max)

fig_forest <- ggplot(forest_df,
                     aes(x = OR, y = label, xmin = CI_lo, xmax = CI_hi)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_errorbarh(height = 0.15, linewidth = 0.7, color = "#2C3E50") +
    geom_point(size = 3.5, color = "#C2185B") +
    # annotation text sitting outside the panel to the right
    geom_text(aes(x = x_label, label = annotation),
              hjust = 0, size = 3.3, color = "grey30") +
    scale_x_log10(
        breaks = c(0.5, 1, 2, 4),
        labels = c("0.5", "1", "2", "4"),
        limits = c(0.4, x_label)          # extend limit to make room
    ) +
    coord_cartesian(clip = "off") +        # allow text to render beyond panel edge
    labs(
        title = "Association between germline BRCA1/2 status and CHIP",
        x     = "Log odds ratio",
        y     = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.margin        = margin(5, 160, 5, 5),   # right margin room for labels
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        axis.text.y        = element_text(size = 11, face = "bold"),
        axis.text.x        = element_text(size = 10)
    )

ggsave(file.path("ch", "figures", "binary_association_forest_s1.pdf"), fig_forest, width = 8, height = 1.75)


