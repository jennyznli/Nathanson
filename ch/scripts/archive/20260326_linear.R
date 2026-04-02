library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)
library(MatchIt)
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
vars <- read_excel(file.path("ch", "data", "ch_wl_art_checked.xlsx")) %>% filter(Keep == 1)
cat("Variants:", nrow(vars), "\n")
cat("Unique carriers:", length(unique(vars$Sample.ID)), "\n")
# 220

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
# 0  1933  134
# 1   483   55

prop.table(tab, margin = 1)
# FALSE       TRUE
# 0 0.93517175 0.06482825
# 1 0.89776952 0.10223048

tab2 <- table(cov$Carrier, cov$CHIP_Binary)
prop.table(tab2, margin = 1)
# FALSE       TRUE
# BRCA1       0.88888889 0.11111111
# BRCA1+BRCA2 1.00000000 0.00000000
# BRCA2       0.90545455 0.09454545
# Non-carrier 0.93517175 0.06482825

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
    p_col <- colnames(s)[4]

    auc <- tryCatch({
        as.numeric(pROC::auc(fit$y, fitted(fit)))
    }, error = function(e) NA_real_)

    data.frame(
        OR    = exp(coef(fit)),
        CI_lo = exp(ci[, 1]),
        CI_hi = exp(ci[, 2]),
        p     = s[, p_col],
        AUC   = auc
    ) %>%
        mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi)) %>%
        tibble::rownames_to_column("term") %>%
        filter(term == !!term)
}

extract_or_firth <- function(fit, term) {
    idx   <- which(names(coef(fit)) == term)
    or    <- exp(coef(fit)[idx])
    ci_lo <- exp(fit$ci.lower[idx])
    ci_hi <- exp(fit$ci.upper[idx])
    p     <- fit$prob[idx]

    auc <- tryCatch({
        as.numeric(pROC::auc(fit$y, fit$predict))
    }, error = function(e) NA_real_)

    data.frame(
        term   = term,
        OR     = or,
        CI_lo  = ci_lo,
        CI_hi  = ci_hi,
        OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", or, ci_lo, ci_hi),
        p      = p,
        AUC    = auc
    )
}

base_covs <- "distance + Batch + Smoke_History + Sequenced_gender +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ========================
# MODEL 1: Linear age
# ========================
distance_df <- data.frame(
    person_id = m.data$person_id,
    distance  = m.data$distance
)

cov <- cov %>%
    left_join(distance_df, by = "person_id")

# verify
cat("NAs in distance:", sum(is.na(cov$distance)), "\n")

# get weights aligned to your filtered cov
weights <- m.data %>%
    filter(person_id %in% cov$person_id) %>%
    pull(weights)

# sanity check
cat("nrow(cov):", nrow(cov), "\n")
cat("length(weights):", length(weights), "\n")

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

# model     age_term           OR_fmt           p       AUC
# ...1        BRCA1/2       linear 1.60 (1.15–2.24) 0.005625915 0.6640628
# ...2          BRCA1       linear 1.64 (1.07–2.51) 0.021997076 0.6539669
# ...3          BRCA2       linear 1.37 (0.88–2.13) 0.166552629 0.6569672
# ...4        BRCA1/2     age+age2 1.60 (1.15–2.23) 0.005903771 0.6614917
# ...5          BRCA1     age+age2 1.64 (1.07–2.51) 0.022974463 0.6524799
# ...6          BRCA2     age+age2 1.37 (0.88–2.13) 0.168344249 0.6545823
# ...7        BRCA1/2   spline_df4 1.62 (1.16–2.26) 0.004882221 0.6661367
# ...8          BRCA1   spline_df4 1.64 (1.07–2.51) 0.021986346 0.6532442
# ...9          BRCA2   spline_df4 1.38 (0.89–2.16) 0.149497771 0.6571950
# BRCA12_Case BRCA1/2 firth_linear 1.59 (1.13–2.20) 0.008145357 0.6651162
# BRCA1_Case    BRCA1 firth_linear 1.63 (1.05–2.45) 0.030481463 0.6556510
# BRCA2_Case    BRCA2 firth_linear 1.37 (0.87–2.09) 0.166748514 0.6589973

results_all
