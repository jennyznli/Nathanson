# ============================================================
# Gene CHIP analysis — Firth logistic in female, cancer-free
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(data.table)
library(MatchIt)
library(logistf)
library(splines)
library(tidyverse)
library(broom)
library(sandwich)
library(lmtest)

# ============================================================
# CONFIG
# ============================================================
MIN_CARRIERS  <- 3
OUT_DIR <- file.path("ch", "data")
FIG_DIR <- file.path("ch", "figures")

exposure_labels <- c(
    "BRCA12_Case" = "BRCA1/2",
    "BRCA1_Case"  = "BRCA1",
    "BRCA2_Case"  = "BRCA2"
)

# shared ggplot theme
theme_ch <- theme_classic(base_size = 12) +
    theme(
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 10, color = "grey40"),
        axis.text        = element_text(color = "black"),
        legend.position  = "bottom"
    )

# ============================================================
# READ DATA
# ============================================================
m.out4 <- readRDS(file.path(OUT_DIR, "ch_psm_matched4.rds"))
m.data <- match.data(m.out4)

cov <- read_excel(file.path(OUT_DIR, "pmbb_brca12_cov_chip_df.xlsx")) %>%
    filter(Sequenced_gender == "Female", Strata == 1) %>%
    left_join(
        m.data %>% dplyr::select(person_id, distance, weights, subclass),
        by = "person_id"
    ) %>%
    filter(!is.na(weights))
dim(cov)

vars <- read_excel(file.path(OUT_DIR, "ch_seq_wl_art_minad4_vars.xlsx")) %>%
    filter(Sample.ID %in% cov$person_id)

# add CHIP flags to cov
cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

# join carrier status onto variant-level data
vars_cov <- vars %>%
    left_join(
        cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                              BRCA2_Case, Carrier, Sample_age,
                              Batch, Smoke_History, starts_with("PC")),
        by = c("Sample.ID" = "person_id")
    ) %>%
    filter(!is.na(BRCA12_Case))

cat("Cohort: N =", nrow(cov), "| CHIP+ =", sum(cov$CHIP_Binary),
    "| Variants =", nrow(vars), "\n")
# Cohort: N = 1924 | CHIP+ = 113 | Variants = 120

# ============================================================
# GENE FREQUENCY TABLE
# ============================================================
gene_person <- vars_cov %>%
    dplyr::select(Sample.ID, Gene, BRCA12_Case) %>%
    distinct()

n_carrier_total    <- sum(cov$BRCA12_Case == 1)
n_noncarrier_total <- sum(cov$BRCA12_Case == 0)

gene_freq <- gene_person %>%
    group_by(Gene, BRCA12_Case) %>%
    summarise(n_mutated = n(), .groups = "drop") %>%
    pivot_wider(names_from  = BRCA12_Case,
                values_from = n_mutated,
                values_fill = 0,
                names_prefix = "n_") %>%
    rename(n_noncarrier = n_0, n_carrier = n_1) %>%
    mutate(
        pct_noncarrier = n_noncarrier / n_noncarrier_total * 100,
        pct_carrier    = n_carrier    / n_carrier_total    * 100,
        n_total        = n_carrier + n_noncarrier
    ) %>%
    arrange(desc(n_total))

dim(gene_freq)
# 34  6
cat("\nGene frequency table (top 15):\n")
print(head(gene_freq, 15))

# ============================================================
# FIRTH LOGISTIC PER GENE
# ============================================================
base_covs <- "Sample_age + Batch + Smoke_History"
pc_covs <- "Sample_age + Batch + Smoke_History + PC1 + PC2 + PC3 + PC4 + PC5 +PC6"

genes_to_test <- gene_freq %>%
    filter(n_carrier >= MIN_CARRIERS) %>%
    pull(Gene)
cat("Genes:", paste(genes_to_test, collapse = ", "), "\n")

run_gene_firth <- function(gene, term, data = cov, vars_df = vars) {
    gene_ids <- vars_df %>% filter(Gene == gene) %>% distinct(Sample.ID) %>% pull(Sample.ID)
    dat      <- data %>% mutate(gene_chip = as.integer(person_id %in% gene_ids))

    n_events    <- sum(dat$gene_chip)
    n_co_events <- dat %>% filter(!!sym(term) == 1) %>% pull(gene_chip) %>% sum()

    if (!(term %in% names(dat))) return(NULL)
    if (n_events < 3)            return(NULL)
    if (n_co_events < 2) return(NULL)        # need >=2 carrier CHIP events
    if ((n_events - n_co_events) < 2) return(NULL)  # need >=2 non-carrier CHIP events

    fit <- tryCatch(
        logistf(as.formula(paste("gene_chip ~", term, "+", base_covs)), data = dat),
        error = function(e) {
            message("Firth failed: ", gene, " / ", term, " — ", conditionMessage(e))
            return(NULL)
        }
    )
    if (is.null(fit)) return(NULL)

    idx <- which(names(coef(fit)) == term)

    data.frame(
        Gene         = gene,
        exposure     = term,
        n_chip       = n_events,
        n_chip_cases = n_co_events,
        OR           = exp(coef(fit)[idx]),
        CI_lo        = exp(fit$ci.lower[idx]),
        CI_hi        = exp(fit$ci.upper[idx]),
        p            = fit$prob[idx],
        method       = "firth"
    ) %>%
        mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

exposures <- c("BRCA12_Case", "BRCA1_Case", "BRCA2_Case")

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
    cat("    Fitted:", nrow(res), "genes\n")
}
gene_results_all <- bind_rows(gene_results_list) %>% arrange(model, p_fdr)

# ============================================================
# SUMMARY
# ============================================================
for (mod in c("BRCA1/2", "BRCA1", "BRCA2")) {
    cat("\n=== Gene results —", mod, "===\n")
    print(gene_results_all %>%
              filter(model == mod) %>%
              dplyr::select(Gene, n_chip, n_chip_cases, OR_fmt, p, p_fdr))
}

sig_genes <- gene_results_all %>% filter(p_fdr  < 0.05)
nom_genes <- gene_results_all %>% filter(p < 0.05)
cat("\nFDR < 0.05 hits:\n");  print(sig_genes)
cat("\nNominal p < 0.05:\n"); print(nom_genes)

write_xlsx(gene_results_all, file.path("ch", "data", "ch_gene_firth_results.xlsx"))

write_xlsx(
    list(
        gene_frequency   = gene_freq,
        gene_fisher      = gene_fisher,
        firth_all        = gene_results_all,
        firth_fdr_sig    = sig_genes,
        firth_nominal    = nom_genes
    ),
    file.path(OUT_DIR, "ch_gene_results.xlsx")
)

saveRDS(gene_results_all, file.path(OUT_DIR, "ch_gene_results.rds"))







