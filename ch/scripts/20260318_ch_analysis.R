library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)

# ========================
# READ IN DATA
# ========================
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x

# 3004
#### MATCHED DATA
matching <- read.csv(file.path("ch", "data", "ch_psm_matched4_data.csv"))
m.out4 <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
weights <- m.out4$weights
m.data <- match_data(m.out4)

### VARS
vars <- read_excel(file.path("ch", "data", "ch_wl_art.xlsx"))
cat("Variants:", nrow(vars), "\n")
cat("Unique carriers:", length(unique(vars$Sample.ID)), "\n")

cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"))
cov <- cov %>% filter(person_id %in% all_ids, Strata == 1)

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

cat("\nCHIP prevalence:\n")
print(table(cov$CHIP_Binary))
cat("\nCHIP by case/control:\n")
tab <- table(cov$BRCA12_Case, cov$CHIP_Binary)
print(tab)
# FALSE TRUE
# 0  2114  213
# 1   600   77

prop.table(tab, margin = 1)
# FALSE       TRUE
# 0 0.94198539 0.05801461
# 1 0.90398818 0.09601182

tab2 <- table(cov$Carrier, cov$CHIP_Binary)
prop.table(tab2, margin = 1)
# FALSE       TRUE
# BRCA1       0.90804598 0.09195402
# BRCA1+BRCA2 1.00000000 0.00000000
# BRCA2       0.92363636 0.07636364
# Non-carrier 0.94146105 0.05853895

# ========================
# CHIP PREVALENCE BY DECADE - ALL
# ========================
prev_by_decade <- cov %>%
    mutate(
        age_group = cut(Sample_age,
                        breaks = c(-Inf, 30, 40, 50, 60, 70, Inf),
                        labels = c("≤30", "30-40", "40-50", "50-60",
                                   "60-70", "≥70"),
                        right  = TRUE)
    ) %>%
    group_by(age_group) %>%
    summarise(
        n_total  = n(),
        n_chip   = sum(CHIP_Binary),
        prev_pct = 100 * mean(CHIP_Binary),
        .groups  = "drop"
    ) %>%
    mutate(x_num = as.numeric(age_group))

p_prev <- ggplot(prev_by_decade, aes(x = x_num, y = prev_pct)) +
    geom_smooth(method = "loess", se = FALSE, color = "#C0392B", linewidth = 1) +
    geom_point(color = "#C0392B", size = 3) +
    scale_x_continuous(
        breaks = prev_by_decade$x_num,
        labels = prev_by_decade$age_group
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    labs(
        title = "CHIP prevalence by decade",
        x     = "Age group (years)",
        y     = "CHIP prevalence (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold", hjust = 0.5)
    )

ggsave(file.path("ch", "figures", "chip_prevalence_by_decade.pdf"),
       p_prev, width = 7, height = 5)

# ========================
# CHIP PREVALENCE BY DECADE - BRCA1/2 CARRIERS
# ========================
prev_by_decade <- cov %>%
    filter(BRCA12_Case == 0) %>%
    mutate(
        age_group = cut(Sample_age,
                        breaks = c(-Inf, 30, 40, 50, 60, 70, Inf),
                        labels = c("≤30", "30-40", "40-50", "50-60",
                                   "60-70", "≥70"),
                        right  = TRUE)
    ) %>%
    group_by(age_group) %>%
    summarise(
        n_total  = n(),
        n_chip   = sum(CHIP_Binary),
        prev_pct = 100 * mean(CHIP_Binary),
        .groups  = "drop"
    ) %>%
    mutate(x_num = as.numeric(age_group))

p_prev <- ggplot(prev_by_decade, aes(x = x_num, y = prev_pct)) +
    geom_smooth(method = "loess", se = FALSE, color = "#C0392B", linewidth = 1) +
    geom_point(color = "#C0392B", size = 3) +
    scale_x_continuous(
        breaks = prev_by_decade$x_num,
        labels = prev_by_decade$age_group
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    labs(
        title = "CHIP prevalence by decade",
        x     = "Age group (years)",
        y     = "CHIP prevalence (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold", hjust = 0.5)
    )

ggsave(file.path("ch", "figures", "chip_prevalence_by_decade_no_brca12.pdf"),
       p_prev, width = 7, height = 5)

# ========================
# CHIP PREVALENCE BY DECADE - BRCA1/2 CARRIERS
# ========================
prev_by_decade <- cov %>%
    filter(BRCA12_Case == 1) %>%
    mutate(
        age_group = cut(Sample_age,
                        breaks = c(-Inf, 30, 40, 50, 60, 70),
                        labels = c("≤30", "30-40", "40-50", "50-60",
                                   "60-70"),
                        right  = TRUE)
    ) %>%
    group_by(age_group) %>%
    summarise(
        n_total  = n(),
        n_chip   = sum(CHIP_Binary),
        prev_pct = 100 * mean(CHIP_Binary),
        .groups  = "drop"
    ) %>%
    mutate(x_num = as.numeric(age_group))

p_prev <- ggplot(prev_by_decade, aes(x = x_num, y = prev_pct)) +
    geom_smooth(method = "loess", se = FALSE, color = "#C0392B", linewidth = 1) +
    geom_point(color = "#C0392B", size = 3) +
    scale_x_continuous(
        breaks = prev_by_decade$x_num,
        labels = prev_by_decade$age_group
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    labs(
        title = "CHIP prevalence by decade",
        x     = "Age group (years)",
        y     = "CHIP prevalence (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold", hjust = 0.5)
    )

ggsave(file.path("ch", "figures", "chip_prevalence_by_decade_brca12.pdf"),
       p_prev, width = 7, height = 5)

# ========================
# HISTOGRAM OF AGE
# ========================
p <- ggplot(cov, aes(x = Sample_age)) +
    geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Age distribution", x = "Age (years)", y = "Frequency") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path("ch", "figures", "strata1_age_hist.png"),
       p, width = 7, height = 5)

# ========================
# CHIP LOGISTIC REGRESSION — MODEL COMPARISON
# ========================
library(splines)
library(lmtest)
library(sandwich)
library(logistf)
library(pROC)

# ========================
# HELPER: extract OR + cluster-robust SE + AUC
# ========================
extract_or <- function(fit, term) {
    s  <- coef(summary(fit))
    ci <- confint.default(fit)
    p_col <- colnames(s)[4]   # always the 4th col regardless of t vs z
    data.frame(
        OR    = exp(coef(fit)),
        CI_lo = exp(ci[, 1]),
        CI_hi = exp(ci[, 2]),
        p     = s[, p_col]
    ) %>%
        mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi)) %>%
        tibble::rownames_to_column("term") %>%
        filter(term == !!term)
}


# Base covariates (no age term — added per model)
base_covs <- "distance + Batch + Smoke_History + Sequenced_gender +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ========================
# MODEL 1: Linear age
# ========================
fit12_linear <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

fit1_linear <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case + Sample_age +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

fit2_linear <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA2_Case + Sample_age +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

# ========================
# MODEL 2: Age + age^2
# ========================
cov <- cov %>% mutate(Sample_age2 = Sample_age^2)

fit12_age2 <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age + Sample_age2 +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

fit1_age2 <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case + Sample_age + Sample_age2 +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

fit2_age2 <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA2_Case + Sample_age + Sample_age2 +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

# ========================
# MODEL 3: Natural spline on age (df=4)
# ========================
fit12_spline <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + ns(Sample_age, df=4) +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

fit1_spline <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case + ns(Sample_age, df=4) +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

fit2_spline <- glm(
    as.formula(paste("CHIP_Binary ~ BRCA2_Case + ns(Sample_age, df=4) +", base_covs)),
    data = cov, weights = weights, family = quasibinomial()
)

# ========================
# MODEL 4: Firth logistic (unweighted, sensitivity)
# ========================
firth12 <- logistf(
    as.formula(paste("CHIP_Binary ~ BRCA12_Case + Sample_age +", base_covs)),
    data = cov
)

firth1 <- logistf(
    as.formula(paste("CHIP_Binary ~ BRCA1_Case + Sample_age +", base_covs)),
    data = cov
)

firth2 <- logistf(
    as.formula(paste("CHIP_Binary ~ BRCA2_Case + Sample_age +", base_covs)),
    data = cov
)

# Extract Firth results separately (no cluster SE support)
extract_or_firth <- function(fit, term) {
    idx   <- which(names(coef(fit)) == term)
    or    <- exp(coef(fit)[idx])
    ci_lo <- exp(fit$ci.lower[idx])
    ci_hi <- exp(fit$ci.upper[idx])
    p     <- fit$prob[idx]
    data.frame(
        term   = term,
        OR     = or,
        CI_lo  = ci_lo,
        CI_hi  = ci_hi,
        OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", or, ci_lo, ci_hi),
        p      = p,
        AUC    = NA_real_
    )
}

# ========================
# LRT: test if spline improves over linear age
# ========================
cat("\n=== LRT: linear vs age^2 ===\n")
print(anova(fit12_linear, fit12_age2,  test = "F"))

cat("\n=== LRT: linear vs spline (df=4) ===\n")
print(anova(fit12_linear, fit12_spline, test = "F"))

# ========================
# RESULTS TABLE
# ========================
results_all <- bind_rows(
    # Linear age
    extract_or(fit12_linear, "BRCA12_Case") %>% mutate(model = "BRCA1/2", age_term = "linear"),
    extract_or(fit1_linear,  "BRCA1_Case") %>% mutate(model = "BRCA1",   age_term = "linear"),
    extract_or(fit2_linear,  "BRCA2_Case") %>% mutate(model = "BRCA2",   age_term = "linear"),
    # Age^2
    extract_or(fit12_age2,   "BRCA12_Case") %>% mutate(model = "BRCA1/2", age_term = "age+age2"),
    extract_or(fit1_age2,    "BRCA1_Case") %>% mutate(model = "BRCA1",   age_term = "age+age2"),
    extract_or(fit2_age2,    "BRCA2_Case") %>% mutate(model = "BRCA2",   age_term = "age+age2"),
    # Spline
    extract_or(fit12_spline, "BRCA12_Case") %>% mutate(model = "BRCA1/2", age_term = "spline_df4"),
    extract_or(fit1_spline,  "BRCA1_Case") %>% mutate(model = "BRCA1",   age_term = "spline_df4"),
    extract_or(fit2_spline,  "BRCA2_Case") %>% mutate(model = "BRCA2",   age_term = "spline_df4"),
    # Firth
    extract_or_firth(firth12, "BRCA12_Case") %>% mutate(model = "BRCA1/2", age_term = "firth_linear"),
    extract_or_firth(firth1,  "BRCA1_Case")  %>% mutate(model = "BRCA1",   age_term = "firth_linear"),
    extract_or_firth(firth2,  "BRCA2_Case")  %>% mutate(model = "BRCA2",   age_term = "firth_linear")
) %>%
    select(model, age_term, OR_fmt, p, AUC)

# model     age_term           OR_fmt          p AUC
# ...1        BRCA1/2       linear 1.23 (0.93–1.63) 0.13923883  NA
# ...2          BRCA1       linear 1.38 (0.97–1.96) 0.07521396  NA
# ...3          BRCA2       linear 1.03 (0.71–1.50) 0.86846070  NA
# ...4        BRCA1/2     age+age2 1.23 (0.94–1.63) 0.13707572  NA
# ...5          BRCA1     age+age2 1.38 (0.97–1.97) 0.07413423  NA
# ...6          BRCA2     age+age2 1.03 (0.71–1.50) 0.86535557  NA
# ...7        BRCA1/2   spline_df4 1.24 (0.94–1.63) 0.13691206  NA
# ...8          BRCA1   spline_df4 1.38 (0.97–1.97) 0.07472395  NA
# ...9          BRCA2   spline_df4 1.03 (0.71–1.51) 0.86244711  NA
# BRCA12_Case BRCA1/2 firth_linear 1.24 (0.93–1.63) 0.14274633  NA
# BRCA1_Case    BRCA1 firth_linear 1.38 (0.96–1.95) 0.07909735  NA
# BRCA2_Case    BRCA2 firth_linear 1.04 (0.71–1.49) 0.83272575  NA

print(results_all)
write.csv(results_all, file.path("ch", "data", "ch_glm_results.csv"), row.names = FALSE)
