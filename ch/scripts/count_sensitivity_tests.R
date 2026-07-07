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
library(lmtest)
library(pROC)
library(lmtest)
library(sandwich)
library(brglm2)  # for penalized Poisson (Firth-type)
library(MASS)        # glm.nb
library(pscl)        # zeroinfl, vuong
library(AER)         # dispersiontest

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "chip_cov_minad4.xlsx"))
cov_f <- cov %>% filter(Sequenced_gender == "Female")

# ============================================================
# TEST DISPERSON - DETERMINE MODEL
# ============================================================

# ── 1. basic count summary ────────────────────────────────────
cat("Mean:", mean(cov$CHIP_Count), "\n")
# 0.07
cat("Variance:", var(cov$CHIP_Count), "\n")
# 0.098
cat("Variance/Mean ratio:", var(cov$CHIP_Count) / mean(cov$CHIP_Count), "\n")
# Poisson assumes var ≈ mean (ratio ~1). If ratio >> 1, overdispersed → NB.
# 1.35 - overdisperssed

table(cov$CHIP_Count)
# 0    1    2    3    6
# 2811  179   10    2    2

prop.table(table(cov$CHIP_Count))  # how many zeros?
# 0           1           2           3           6
# 0.935752330 0.059587217 0.003328895 0.000665779 0.000665779

# ── 2. fit both and formally test overdispersion ──────────────
fit_pois <- glm(CHIP_Count ~ BRCA12_Case + Sample_age + Sample_age2 +
                    Batch + Smoke_History + Sequenced_gender +
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                data = cov, family = poisson())

sum(residuals(fit_pois, "pearson")^2) / df.residual(fit_pois)
# 1.18

fit_nb <- glm.nb(CHIP_Count ~ BRCA12_Case + Sample_age + Sample_age2 +
                     Batch + Smoke_History + Sequenced_gender +
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                 data = cov)



# formal overdispersion test on Poisson fit
dispersiontest(fit_pois)
# p < 0.05 → overdispersed → use NB
# Overdispersion test
#
# data:  fit_pois
# z = 1.726, p-value = 0.04218
# alternative hypothesis: true dispersion is greater than 1
# sample estimates:
#     dispersion
# 1.147629

# likelihood ratio test Poisson vs NB (NB reduces to Poisson when theta → Inf)
pchisq(2 * (logLik(fit_nb) - logLik(fit_pois)), df = 1, lower.tail = FALSE)
# p < 0.05 → NB significantly better fit
# 'log Lik.' 1.806845e-07 (df=14)

# ============================================================
# RUN LOOPS - POISSON
# ============================================================
source(here("ch", "scripts", "count_sensitivity_models.R"))
res <- run_chip_qp(cov_base = cov, minad_thresholds = 4)
write_xlsx(res, file.path("ch", "data", "chip_count_all_s12_sensitivity.xlsx"))

# res_all_s12 <- run_chip_count(
#     cov_base      = cov,
#     cohort_label  = "all_s12",
#     females_only = FALSE
# )

# res_f_s12 <- run_chip_count(
#     cov_base      = cov_f,
#     cohort_label  = "female_s12",
#     females_only = TRUE
# )

write_xlsx(res_f_s12, file.path("ch", "data", "chip_count_female_s12_sensitivity.xlsx"))

