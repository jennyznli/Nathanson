# ============================================================
# CHIP prevalence + age association + gene distribution
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
cov <- read_excel(file.path(OUT_DIR, "pmbb_brca12_cov_chip_df.xlsx"))

vars <- read_excel(file.path(OUT_DIR, "ch_seq_wl_art_minad4_vars.xlsx"))

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
# Cohort: N = 3004 | CHIP+ = 193 | Variants = 217

# ============================================================
# G1: GENE FREQUENCY TABLE
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
# G2: FISHER'S EXACT TEST PER GENE
# ============================================================
gene_fisher <- gene_freq %>%
    rowwise() %>%
    mutate(
        fisher_p = fisher.test(matrix(
            c(n_carrier,
              n_carrier_total    - n_carrier,
              n_noncarrier,
              n_noncarrier_total - n_noncarrier),
            nrow = 2
        ))$p.value,
        OR_crude = (n_carrier / (n_carrier_total - n_carrier)) /
            (n_noncarrier / (n_noncarrier_total - n_noncarrier))
    ) %>%
    ungroup() %>%
    arrange(fisher_p)

# Separate reportable genes BEFORE computing FDR
gene_fisher_reportable <- gene_fisher %>%
    filter(is.finite(OR_crude), n_noncarrier >= 1) %>%
    mutate(fdr = p.adjust(fisher_p, method = "BH"))

# Keep full table with FDR for saving, but compute separately
gene_fisher <- gene_fisher %>%
    mutate(fdr = p.adjust(fisher_p, method = "BH"))

cat("\nFisher results (top 10):\n")
print(head(gene_fisher_reportable %>% dplyr::select(Gene, n_carrier, n_noncarrier,
                                                    pct_carrier, pct_noncarrier,
                                                    OR_crude, fisher_p, fdr), 10))

# ============================================================
# FIGURES
# ============================================================
dot_data <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    pivot_longer(cols = c(pct_carrier, pct_noncarrier),
                 names_to  = "group",
                 values_to = "pct") %>%
    mutate(
        group = recode(group,
                       pct_carrier    = "gBRCA1/2 Carrier",
                       pct_noncarrier = "Non-carrier"),
        Gene  = factor(Gene, levels = top_genes)  # highest total on left
    )

# add overall prevalence per gene as a third layer
dot_overall <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    mutate(
        pct   = (n_carrier + n_noncarrier) /
                (n_carrier_total + n_noncarrier_total) * 100,
        group = "Overall",
        Gene  = factor(Gene, levels = top_genes)
    )

fig1 <- ggplot() +
    # vertical reference lines per gene
    geom_vline(data = dot_data %>% filter(group == "gBRCA1/2 Carrier"),
               aes(xintercept = as.numeric(Gene)),
               color = "grey88", linewidth = 0.4) +
    # group dots
    geom_point(data = dot_data,
               aes(x = Gene, y = pct, color = group, shape = group),
               size = 3, alpha = 0.9) +
    # overall dot on top
    geom_point(data = dot_overall,
               aes(x = Gene, y = pct, color = group, shape = group),
               size = 3, alpha = 0.9) +
    scale_color_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B",
                   "Non-carrier"       = "#546E7A",
                   "Overall"           = "black"),
        name = NULL
    ) +
    scale_shape_manual(
        values = c("gBRCA1/2 Carrier" = 16,
                   "Non-carrier"       = 16,
                   "Overall"           = 18),  # diamond for overall, matches QC fig
        name = NULL
    ) +
    scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        limits = c(0, NA),
        expand = expansion(mult = c(0.02, 0.1))
    ) +
    labs(
        title    = "CHIP gene prevalence by carrier status",
        subtitle = sprintf("gBRCA1/2 Carrier n = %d | Non-carrier n = %d | Overall n = %d",
                           n_carrier_total, n_noncarrier_total,
                           n_carrier_total + n_noncarrier_total),
        x        = NULL,
        y        = "Carrier frequency"
    ) +
    theme_ch +
    theme(
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        legend.position  = "bottom"
    )

ggsave(file.path(FIG_DIR, "fig_gene_prevalence_dot.pdf"),
       fig1, width = 8, height = 5)
ggsave(file.path(FIG_DIR, "fig_gene_prevalence_dot.png"),
       fig1, width = 8, height = 5, dpi = 300)
cat("Saved: fig_gene_prevalence_dot\n")
