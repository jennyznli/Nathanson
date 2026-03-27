# ========================
# CHIP - LOGISTIC REGRESSION
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))


library(data.table, quietly = T)
library(MatchIt)
library(MASS)
library(logistf)
library(p.adjust)
library(splines)
library(lmtest)
library(sandwich)
library(logistf)
library(pROC)

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

cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx"))

# ========================
# FUNCTIONS
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

# ========================
# CHIP COUNT: POISSON vs NEGATIVE BINOMIAL
# ========================
library(MASS)

# first test for overdispersion with Poisson
fit_pois <- glm(
    as.formula(paste("CHIP_Count ~ BRCA12_Case + ns(Sample_age, df=4) +", base_covs)),
    data = cov, weights = weights, family = poisson()
)

# overdispersion test: ratio of residual deviance to df should be ~1
overdisp_ratio <- fit_pois$deviance / fit_pois$df.residual
cat("Overdispersion ratio:", overdisp_ratio, "\n")
# > 1 = overdispersed → use NB

# negative binomial (unweighted — glm.nb doesn't support weights well)
fit_nb12 <- glm.nb(
    as.formula(paste("CHIP_Count ~ BRCA12_Case + ns(Sample_age, df=4) +", base_covs)),
    data = cov
)

fit_nb1 <- glm.nb(
    as.formula(paste("CHIP_Count ~ BRCA1_Case + ns(Sample_age, df=4) +", base_covs)),
    data = cov
)

fit_nb2 <- glm.nb(
    as.formula(paste("CHIP_Count ~ BRCA2_Case + ns(Sample_age, df=4) +", base_covs)),
    data = cov
)

# compare Poisson vs NB with LRT
cat("\n=== Poisson vs NB LRT ===\n")
pchisq(2 * (logLik(fit_nb12) - logLik(fit_pois)), df = 1, lower.tail = FALSE)

# extract IRR (incidence rate ratio) instead of OR for count models
extract_irr <- function(fit, term) {
    s  <- coef(summary(fit))
    ci <- confint.default(fit)
    data.frame(
        IRR   = exp(coef(fit)[term]),
        CI_lo = exp(ci[term, 1]),
        CI_hi = exp(ci[term, 2]),
        p     = s[term, 4],
        IRR_fmt = sprintf("%.2f (%.2f\u2013%.2f)",
                          exp(coef(fit)[term]),
                          exp(ci[term, 1]),
                          exp(ci[term, 2]))
    )
}

count_results <- bind_rows(
    extract_irr(fit_nb12, "BRCA12_Case") %>% mutate(model = "BRCA1/2"),
    extract_irr(fit_nb1,  "BRCA1_Case")  %>% mutate(model = "BRCA1"),
    extract_irr(fit_nb2,  "BRCA2_Case")  %>% mutate(model = "BRCA2")
) %>% select(model, IRR_fmt, p)

print(count_results)
write.csv(count_results, file.path("ch", "data", "ch_nb_results.csv"), row.names = FALSE)

# ========================
# GENE-SPECIFIC CHIP: FIRTH LOGISTIC PER GENE
# ========================

# get genes with enough carriers to fit
gene_counts <- vars %>%
    count(Gene, name = "n_carriers") %>%
    filter(n_carriers >= 5) %>%   # minimum 5 carriers to attempt fit
    arrange(desc(n_carriers))

cat("Genes with >= 5 carriers:", nrow(gene_counts), "\n")
print(gene_counts)

# add per-gene binary carrier flags to cov
run_gene_firth <- function(gene, term = "BRCA12_Case") {
    gene_carriers <- vars %>%
        filter(Gene == gene) %>%
        distinct(Sample.ID) %>%
        pull(Sample.ID)

    cov_gene <- cov %>%
        mutate(gene_chip = person_id %in% gene_carriers)

    # skip if no variance
    if (sum(cov_gene$gene_chip) < 3) return(NULL)

    tryCatch({
        fit <- logistf(
            as.formula(paste("gene_chip ~", term, "+ ns(Sample_age, df=4) +", base_covs)),
            data = cov_gene
        )
        idx <- which(names(coef(fit)) == term)
        data.frame(
            Gene   = gene,
            term   = term,
            n      = sum(cov_gene$gene_chip),
            OR     = exp(coef(fit)[idx]),
            CI_lo  = exp(fit$ci.lower[idx]),
            CI_hi  = exp(fit$ci.upper[idx]),
            p      = fit$prob[idx]
        ) %>%
            mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
    }, error = function(e) {
        cat("Failed for gene:", gene, "\n")
        NULL
    })
}

# run for BRCA1/2 combined, BRCA1, BRCA2
gene_results12 <- gene_counts$Gene %>%
    lapply(run_gene_firth, term = "BRCA12_Case") %>%
    bind_rows() %>%
    mutate(
        p_fdr        = p.adjust(p, method = "BH"),
        p_bonferroni = p.adjust(p, method = "bonferroni"),
        model        = "BRCA1/2"
    )

gene_results1 <- gene_counts$Gene %>%
    lapply(run_gene_firth, term = "BRCA1_Case") %>%
    bind_rows() %>%
    mutate(
        p_fdr        = p.adjust(p, method = "BH"),
        p_bonferroni = p.adjust(p, method = "bonferroni"),
        model        = "BRCA1"
    )

gene_results2 <- gene_counts$Gene %>%
    lapply(run_gene_firth, term = "BRCA2_Case") %>%
    bind_rows() %>%
    mutate(
        p_fdr        = p.adjust(p, method = "BH"),
        p_bonferroni = p.adjust(p, method = "bonferroni"),
        model        = "BRCA2"
    )

gene_results_all <- bind_rows(gene_results12, gene_results1, gene_results2) %>%
    arrange(Gene)

print(gene_results_all %>% arrange((gene_results_all$Gene)))
write.csv(gene_results_all, file.path("ch", "data", "ch_gene_firth_results.csv"), row.names = FALSE)




