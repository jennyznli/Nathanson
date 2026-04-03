# ============================================================
# SCRIPT 3: GENE-SPECIFIC CHIP — FIRTH LOGISTIC PER GENE
# BRCA1/2 Carrier CHIP Association Analysis
# ============================================================
# For each CHIP gene with ≥ 5 variant carriers:
#   - Firth penalized logistic (handles separation from rare genes)
#   - Exposure: BRCA12_Case (combined), BRCA1_Case, BRCA2_Case
#   - Outcome: binary CHIP in that gene (1 = any variant in gene)
#   - Multiple testing: BH-FDR and Bonferroni across genes
#
# NOTE: We run on the matched sample (no weights in logistf).
# This is the recommended approach for gene-level rare outcomes.
#
# Optional: weighted Firth via Firth's method is not available;
# instead, run standard weighted logistic as a sensitivity check
# for genes with sufficient events (n ≥ 20).
# ============================================================

library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(data.table)
library(MatchIt)
library(logistf)
library(splines)
library(tidyverse)
library(writexl)

# ============================================================
# READ DATA
# ============================================================
m.out4 <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
m.data <- match.data(m.out4)

cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx")) %>%
    left_join(
        m.data %>% dplyr::select(person_id, distance, weights, subclass),
        by = "person_id"
    ) %>%
    filter(!is.na(weights), Strata == 1)

vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))

# ============================================================
# FORMULA BASE
# ============================================================
base_covs_spline <- "Sample_age + Sequenced_gender + Batch +
                      Smoke_History +
                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ============================================================
# GENE UNIVERSE — minimum carrier threshold
# ============================================================
MIN_CARRIERS <- 4   # minimum CHIP+ in gene to attempt fit

gene_counts <- vars %>%
    count(Gene, name = "n_chip_carriers") %>%
    arrange(desc(n_chip_carriers))

cat("\nAll genes with variant counts:\n")
print(gene_counts)

genes_to_test <- gene_counts %>%
    filter(n_chip_carriers >= MIN_CARRIERS) %>%
    pull(Gene)

cat("\nGenes passing threshold (n ≥", MIN_CARRIERS, "):", length(genes_to_test), "\n")

# ============================================================
# CORE FUNCTION: FIRTH PER GENE
# ============================================================
run_gene_firth <- function(gene, term, data = cov, vars_df = vars) {

    gene_carrier_ids <- vars_df %>%
        filter(Gene == gene) %>%
        distinct(Sample.ID) %>%
        pull(Sample.ID)

    dat <- data %>%
        mutate(gene_chip = as.integer(person_id %in% gene_carrier_ids))

    n_events   <- sum(dat$gene_chip)
    n_cases    <- sum(dat[[gsub("_Case", "_Case", term)]])   # n BRCA carriers
    n_co_events <- dat %>% filter(!!sym(term) == 1) %>% pull(gene_chip) %>% sum()

    # Skip if term variable absent or no co-events at all
    if (!(term %in% names(dat))) return(NULL)
    if (n_events < 3)            return(NULL)

    result_base <- data.frame(
        Gene         = gene,
        exposure     = term,
        n_chip       = n_events,
        n_chip_cases = n_co_events
    )

    # --- Firth logistic ---
    fit <- tryCatch({
        logistf(as.formula(paste("gene_chip ~", term, "+", base_covs_spline)), data = dat)
    }, error = function(e) {
        message("Firth failed for gene: ", gene, " | term: ", term, " | ", conditionMessage(e))
        return(NULL)
    })

    if (is.null(fit)) return(result_base %>% mutate(OR = NA, CI_lo = NA, CI_hi = NA, p = NA, OR_fmt = NA, method = "firth_failed"))

    idx <- which(names(coef(fit)) == term)

    result_base %>% mutate(
        OR     = exp(coef(fit)[idx]),
        CI_lo  = exp(fit$ci.lower[idx]),
        CI_hi  = exp(fit$ci.upper[idx]),
        p      = fit$prob[idx],
        OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi),
        method = "firth"
    )
}

# ============================================================
# SENSITIVITY: WEIGHTED GLM FOR HIGHER-COUNT GENES
# ============================================================
run_gene_weighted_glm <- function(gene, term, data = cov, vars_df = vars,
                                  min_events = 20) {
    gene_carrier_ids <- vars_df %>%
        filter(Gene == gene) %>%
        distinct(Sample.ID) %>%
        pull(Sample.ID)

    dat <- data %>%
        mutate(gene_chip = as.integer(person_id %in% gene_carrier_ids))

    if (sum(dat$gene_chip) < min_events) return(NULL)

    fit <- tryCatch({
        glm(as.formula(paste("gene_chip ~", term, "+", base_covs_spline)),
            data = dat, weights = dat$weights, family = quasibinomial())
    }, error = function(e) NULL)

    if (is.null(fit)) return(NULL)

    vcov_cl <- tryCatch(sandwich::vcovCL(fit, cluster = ~subclass, data = dat), error = function(e) NULL)
    if (is.null(vcov_cl)) return(NULL)

    ct <- lmtest::coeftest(fit, vcov = vcov_cl)
    ci <- lmtest::coefci(fit, vcov = vcov_cl)

    data.frame(
        Gene     = gene,
        exposure = term,
        n_chip   = sum(dat$gene_chip),
        OR       = exp(ct[term, "Estimate"]),
        CI_lo    = exp(ci[term, 1]),
        CI_hi    = exp(ci[term, 2]),
        p        = ct[term, "Pr(>|z|)"],
        method   = "weighted_glm"
    ) %>%
        mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

# ============================================================
# RUN ALL GENES × ALL EXPOSURES
# ============================================================
exposures <- c("BRCA12_Case", "BRCA1_Case", "BRCA2_Case")
exposure_labels <- c("BRCA12_Case" = "BRCA1/2", "BRCA1_Case" = "BRCA1", "BRCA2_Case" = "BRCA2")
#
# exposures <- c("BRCA12_Case", "BRCA1_Case", "BRCA2_Case")
# exposure_labels <- c("BRCA12_Case" = "BRCA1/2", "BRCA1_Case" = "BRCA1", "BRCA2_Case" = "BRCA2")


cat("\nRunning Firth per gene...\n")

gene_results_list <- list()
for (exp in exposures) {
    cat("  Exposure:", exp, "\n")
    res <- lapply(genes_to_test, run_gene_firth, term = exp) %>%
        bind_rows() %>%
        filter(!is.na(OR)) %>%
        mutate(
            model        = exposure_labels[exp],
            p_fdr        = p.adjust(p, method = "BH"),
            p_bonferroni = p.adjust(p, method = "bonferroni")
        )
    gene_results_list[[exp]] <- res
    cat("    Done:", nrow(res), "genes fitted\n")
}

gene_results_all <- bind_rows(gene_results_list) %>%
    arrange(model, p_fdr)

# ============================================================
# SENSITIVITY: WEIGHTED GLM FOR HIGH-COUNT GENES
# ============================================================
cat("\nRunning weighted GLM sensitivity for high-count genes...\n")

gene_sens_list <- list()
for (exp in exposures) {
    res_sens <- lapply(genes_to_test, run_gene_weighted_glm, term = exp) %>%
        bind_rows() %>%
        mutate(model = exposure_labels[exp])
    gene_sens_list[[exp]] <- res_sens
}
gene_sens_all <- bind_rows(gene_sens_list) %>%
    arrange(model, p)

# ============================================================
# SUMMARY
# ============================================================
cat("\n=== GENE RESULTS — BRCA1/2 (FDR-sorted) ===\n")
print(gene_results_all %>% filter(model == "BRCA1/2") )

cat("\n=== GENE RESULTS — BRCA1 ===\n")
print(gene_results_all %>% filter(model == "BRCA1") )

cat("\n=== GENE RESULTS — BRCA2 ===\n")
print(gene_results_all %>% filter(model == "BRCA2"))

# Significant hits (FDR < 0.05)
sig_genes <- gene_results_all %>% filter(p_fdr < 0.05)
cat("\nFDR < 0.05 hits:\n")
print(sig_genes)

# Nominally significant (p < 0.05, not FDR corrected)
nom_genes <- gene_results_all %>% filter(p < 0.05)
cat("\nNominal p < 0.05:\n")
print(nom_genes )

# ============================================================
# SAVE
# ============================================================
write_xlsx(
    list(
        all_results      = gene_results_all,
        fdr_significant  = sig_genes,
        nominal_p05      = nom_genes,
        sensitivity_wglm = gene_sens_all
    ),
    file.path("ch", "data", "ch_gene_firth_results.xlsx")
)

saveRDS(gene_results_all, file.path("ch", "data", "ch_gene_results.rds"))

cat("\nDone: 03_gene_specific.R\n")
