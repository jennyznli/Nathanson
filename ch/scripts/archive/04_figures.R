# ============================================================
# SCRIPT 4: THESIS FIGURES
# BRCA1/2 Carrier CHIP Association Analysis
# ============================================================
# Figures:
#   Fig 1 — Study design CONSORT/flowchart (text-based, done elsewhere)
#   Fig 2 — Covariate balance: Love plot (SMD before/after matching)
#   Fig 3 — CHIP prevalence by carrier status (bar + CI)
#   Fig 4 — CHIP prevalence by age × carrier status (spline smooth)
#   Fig 5 — Forest plot: OR across all binary models
#   Fig 6 — Forest plot: gene-specific OR (BRCA1/2 combined)
#   Fig 7 — Volcano plot: gene-specific results
#   Fig 8 — CHIP count distribution (violin/boxplot by carrier)
#   Fig 9 — Model comparison: OR stability across specifications
# ============================================================

library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(MatchIt)
library(cobalt)
library(ggplot2)
library(ggtext)
library(patchwork)
library(tidyverse)
library(scales)
library(splines)
library(forcats)

SAVE_DIR <- file.path("ch", "figures")
THEME_BASE <- theme_minimal(base_size = 12) +
    theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        axis.title    = element_text(size = 11),
        legend.title  = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text    = element_text(face = "bold")
    )

COL_B12  <- "#7B2D8B"   # purple for BRCA1/2 combined
COL_B1   <- "#C2185B"   # magenta for BRCA1
COL_B2   <- "#0277BD"   # blue for BRCA2
COL_CTRL <- "#546E7A"   # grey for controls

# ============================================================
# LOAD SAVED FITS AND DATA
# ============================================================
m.out4     <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
binary_fits <- readRDS(file.path("ch", "data", "ch_binary_model_fits.rds"))
cov         <- binary_fits$data
gene_res    <- readRDS(file.path("ch", "data", "ch_gene_results.rds"))

# ============================================================
# FIG 2: LOVE PLOT — covariate balance
# ============================================================
# cobalt::love.plot gives publication-ready SMD before/after
fig2_love <- love.plot(
    m.out4,
    stats      = "mean.diffs",
    threshold  = 0.1,
    abs        = TRUE,
    var.order  = "unadjusted",
    colors     = c("Unadjusted" = COL_B1, "Adjusted" = COL_B2),
    shapes     = c("Unadjusted" = 17, "Adjusted" = 16),
    size       = 3,
    title      = "Covariate Balance Before and After Propensity Score Matching",
    var.names  = c(
        distance         = "Propensity Score",
        Sample_age       = "Age at Sample Collection",
        Sequenced_genderFemale = "Female Sex",
        Smoke_History0   = "Non-smoker",
        Batch1           = "Sequencing Batch (Freeze 2)",
        PC1 = "PC1", PC2 = "PC2", PC3 = "PC3",
        PC4 = "PC4", PC5 = "PC5", PC6 = "PC6"
    )
) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(SAVE_DIR, "fig2_love_plot.pdf"), fig2_love, width = 8, height = 5, dpi = 300)
cat("Saved: fig2_love_plot.pdf\n")

# ============================================================
# FIG 3: CHIP PREVALENCE BY CARRIER STATUS (bar + 95% CI)
# ============================================================
# Wilson 95% CI for proportions
prev_df <- cov %>%
    mutate(
        Group = case_when(
            BRCA1_Case == 1 & BRCA2_Case == 0 ~ "BRCA1",
            BRCA2_Case == 1 & BRCA1_Case == 0 ~ "BRCA2",
            BRCA12_Case == 1 ~ "BRCA1+BRCA2",
            TRUE ~ "Non-carrier"
        ),
        Group = factor(Group, levels = c("BRCA1", "BRCA2", "BRCA1+BRCA2", "Non-carrier"))
    ) %>%
    group_by(Group) %>%
    summarise(
        n     = n(),
        n_pos = sum(CHIP_Binary),
        prev  = n_pos / n,
        # Wilson CI
        ci_lo = (prev + qnorm(0.975)^2 / (2 * n) -
                     qnorm(0.975) * sqrt(prev * (1 - prev) / n + qnorm(0.975)^2 / (4 * n^2))) /
                (1 + qnorm(0.975)^2 / n),
        ci_hi = (prev + qnorm(0.975)^2 / (2 * n) +
                     qnorm(0.975) * sqrt(prev * (1 - prev) / n + qnorm(0.975)^2 / (4 * n^2))) /
                (1 + qnorm(0.975)^2 / n),
        .groups = "drop"
    )

fig3_prev <- ggplot(prev_df, aes(x = Group, y = prev * 100, fill = Group)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(aes(ymin = ci_lo * 100, ymax = ci_hi * 100), width = 0.2, linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", prev * 100, n_pos),
                  y = ci_hi * 100 + 0.5), size = 3.2, vjust = 0) +
    scale_fill_manual(values = c(
        "BRCA1"       = COL_B1,
        "BRCA2"       = COL_B2,
        "BRCA1+BRCA2" = COL_B12,
        "Non-carrier" = COL_CTRL
    )) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.15))) +
    labs(
        title    = "CHIP Prevalence by Germline Carrier Status",
        subtitle = "Propensity score-matched cohort | 95% Wilson CI",
        x        = NULL,
        y        = "CHIP Prevalence (%)"
    ) +
    THEME_BASE +
    theme(legend.position = "none")

ggsave(file.path(SAVE_DIR, "fig3_chip_prevalence_by_carrier.pdf"), fig3_prev,
       width = 6, height = 5, dpi = 300)
cat("Saved: fig3_chip_prevalence_by_carrier.pdf\n")

# ============================================================
# FIG 4: CHIP PREVALENCE BY AGE × CARRIER STATUS (spline smooth)
# ============================================================
# Loess/GAM smooth of CHIP_Binary ~ Age, stratified by carrier status
fig4_data <- cov %>%
    mutate(
        Carrier = case_when(
            BRCA12_Case == 1 ~ "BRCA1/2 Carrier",
            TRUE ~ "Non-carrier"
        )
    )

fig4_age <- ggplot(fig4_data, aes(x = Sample_age, y = CHIP_Binary, color = Carrier, fill = Carrier)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                se = TRUE, alpha = 0.15, linewidth = 1.2) +
    geom_rug(data = fig4_data %>% filter(CHIP_Binary == 1),
             aes(x = Sample_age), sides = "t", alpha = 0.2, length = unit(0.03, "npc")) +
    geom_rug(data = fig4_data %>% filter(CHIP_Binary == 0),
             aes(x = Sample_age), sides = "b", alpha = 0.05, length = unit(0.03, "npc")) +
    scale_color_manual(values = c("BRCA1/2 Carrier" = COL_B12, "Non-carrier" = COL_CTRL)) +
    scale_fill_manual(values  = c("BRCA1/2 Carrier" = COL_B12, "Non-carrier" = COL_CTRL)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    labs(
        title    = "CHIP Prevalence by Age and Carrier Status",
        subtitle = "GAM smooth with 95% CI | rug = CHIP+ individuals",
        x        = "Age at Sample Collection (years)",
        y        = "CHIP Probability"
    ) +
    THEME_BASE +
    theme(legend.position = "bottom")

ggsave(file.path(SAVE_DIR, "fig4_chip_by_age_carrier.pdf"), fig4_age,
       width = 7, height = 5, dpi = 300)
cat("Saved: fig4_chip_by_age_carrier.pdf\n")

# Stratified by BRCA1 vs BRCA2
fig4_data2 <- cov %>%
    mutate(Carrier = case_when(
        BRCA1_Case == 1 & BRCA2_Case == 0 ~ "BRCA1",
        BRCA2_Case == 1 & BRCA1_Case == 0 ~ "BRCA2",
        TRUE ~ "Non-carrier"
    ))

fig4b_age <- ggplot(fig4_data2, aes(x = Sample_age, y = CHIP_Binary, color = Carrier, fill = Carrier)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                se = TRUE, alpha = 0.15, linewidth = 1.2) +
    scale_color_manual(values = c("BRCA1" = COL_B1, "BRCA2" = COL_B2, "Non-carrier" = COL_CTRL)) +
    scale_fill_manual(values  = c("BRCA1" = COL_B1, "BRCA2" = COL_B2, "Non-carrier" = COL_CTRL)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    labs(
        title    = "CHIP Prevalence by Age: BRCA1 vs BRCA2",
        subtitle = "GAM smooth with 95% CI",
        x        = "Age at Sample Collection (years)",
        y        = "CHIP Probability"
    ) +
    THEME_BASE +
    theme(legend.position = "bottom")

ggsave(file.path(SAVE_DIR, "fig4b_chip_by_age_brca1v2.pdf"), fig4b_age,
       width = 7, height = 5, dpi = 300)
cat("Saved: fig4b_chip_by_age_brca1v2.pdf\n")

# ============================================================
# FIG 5: FOREST PLOT — OR across binary model specs
# ============================================================
binary_results <- read_xlsx(file.path("ch", "data", "ch_chip_binary_results.xlsx"))

# Clean labels for figure
spec_labels <- c(
    "M1_linear"          = "M1: Weighted logistic\n(linear age)",
    "M1s_spline_PRIMARY" = "M1s: Weighted logistic\n(spline age) [PRIMARY]",
    "M2_conditional"     = "M2: Conditional logistic\n(within matched sets)",
    "M3_firth"           = "M3: Firth logistic\n(separation sensitivity)"
)

forest_df <- binary_results %>%
    filter(model %in% c("BRCA1/2", "BRCA1", "BRCA2")) %>%
    mutate(
        spec_lab = factor(spec_labels[spec], levels = rev(spec_labels)),
        model    = factor(model, levels = c("BRCA1/2", "BRCA1", "BRCA2")),
        sig      = p < 0.05
    )

fig5_forest <- ggplot(forest_df, aes(x = OR, xmin = CI_lo, xmax = CI_hi,
                                      y = spec_lab, color = model, shape = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(height = 0), linewidth = 0.8, position = position_dodge(width = 0.5)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = OR_fmt, x = max(CI_hi, na.rm = TRUE) + 0.3),
              size = 2.8, hjust = 0, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("BRCA1/2" = COL_B12, "BRCA1" = COL_B1, "BRCA2" = COL_B2)) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
                  labels = c("0.5", "0.75", "1.0", "1.5", "2.0", "3.0")) +
    facet_wrap(~ model, ncol = 3) +
    labs(
        title    = "CHIP Association with BRCA1/2 Carrier Status",
        subtitle = "Odds Ratios across model specifications | filled = p < 0.05",
        x        = "Odds Ratio (log scale)",
        y        = NULL
    ) +
    THEME_BASE +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(SAVE_DIR, "fig5_forest_binary_models.pdf"), fig5_forest,
       width = 12, height = 5, dpi = 300)
cat("Saved: fig5_forest_binary_models.pdf\n")

# ============================================================
# FIG 6: GENE-SPECIFIC FOREST PLOT (BRCA1/2 combined)
# ============================================================
gene_b12 <- gene_res %>%
    filter(model == "BRCA1/2", !is.na(OR)) %>%
    arrange(OR) %>%
    mutate(
        Gene    = fct_reorder(Gene, OR),
        sig_fdr = p_fdr < 0.05,
        sig_nom = p < 0.05 & !sig_fdr,
        sig_lab = case_when(
            sig_fdr ~ "FDR < 0.05",
            sig_nom ~ "Nominal p < 0.05",
            TRUE    ~ "NS"
        ),
        sig_lab = factor(sig_lab, levels = c("FDR < 0.05", "Nominal p < 0.05", "NS"))
    )

# Clip extreme CIs for display
CI_CLIP <- 10
gene_b12 <- gene_b12 %>%
    mutate(
        CI_lo_plot = pmax(CI_lo, 1 / CI_CLIP),
        CI_hi_plot = pmin(CI_hi, CI_CLIP),
        clipped    = CI_hi > CI_CLIP | CI_lo < 1 / CI_CLIP
    )

fig6_genes <- ggplot(gene_b12,
    aes(x = OR, xmin = CI_lo_plot, xmax = CI_hi_plot,
        y = Gene, color = sig_lab)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(height = 0), linewidth = 0.7) +
    geom_point(aes(size = n_chip), alpha = 0.85) +
    geom_text(aes(label = sprintf("n=%d", n_chip), x = CI_CLIP + 0.3),
              size = 2.5, hjust = 0, color = "grey40") +
    scale_color_manual(values = c(
        "FDR < 0.05"       = "#B71C1C",
        "Nominal p < 0.05" = "#E65100",
        "NS"               = "grey60"
    )) +
    scale_size_continuous(name = "N CHIP carriers", range = c(1.5, 5)) +
    scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5, 10),
                  labels = c("0.2", "0.5", "1", "2", "5", "10"),
                  limits = c(1/CI_CLIP, CI_CLIP + 1)) +
    labs(
        title    = "Gene-Specific CHIP Association with BRCA1/2 Carrier Status",
        subtitle = "Firth logistic regression | point size ∝ CHIP carrier count",
        x        = "Odds Ratio (log scale)",
        y        = "CHIP Gene",
        color    = "Significance"
    ) +
    THEME_BASE +
    theme(
        legend.position = "right",
        axis.text.y     = element_text(size = 9)
    )

h <- max(5, nrow(gene_b12) * 0.35 + 1.5)
ggsave(file.path(SAVE_DIR, "fig6_forest_gene_brca12.pdf"), fig6_genes,
       width = 9, height = h, dpi = 300)
cat("Saved: fig6_forest_gene_brca12.pdf\n")

# ============================================================
# FIG 7: VOLCANO PLOT — gene-specific results
# ============================================================
gene_volcano <- gene_res %>%
    filter(!is.na(OR), !is.na(p)) %>%
    mutate(
        log2OR  = log2(OR),
        neg_log10_p = -log10(p),
        sig_lab = case_when(
            p_fdr < 0.05    ~ "FDR < 0.05",
            p < 0.05        ~ "p < 0.05",
            TRUE            ~ "NS"
        ),
        sig_lab = factor(sig_lab, levels = c("FDR < 0.05", "p < 0.05", "NS")),
        label   = ifelse(p < 0.05, Gene, NA)
    )

fig7_volcano <- ggplot(gene_volcano, aes(x = log2OR, y = neg_log10_p,
                                          color = sig_lab, size = n_chip)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                              box.padding = 0.5, segment.color = "grey50") +
    scale_color_manual(values = c(
        "FDR < 0.05" = "#B71C1C",
        "p < 0.05"   = "#E65100",
        "NS"         = "grey70"
    )) +
    scale_size_continuous(name = "N CHIP carriers", range = c(1, 5)) +
    scale_x_continuous(breaks = seq(-4, 4, 1)) +
    facet_wrap(~ model, ncol = 3) +
    labs(
        title    = "Volcano Plot: Gene-Specific CHIP Association",
        subtitle = "x-axis: log2(OR) vs BRCA carrier | dashed line: p = 0.05",
        x        = "log2(Odds Ratio)",
        y        = "-log10(p-value)",
        color    = "Significance"
    ) +
    THEME_BASE +
    theme(legend.position = "bottom")

ggsave(file.path(SAVE_DIR, "fig7_volcano_gene.pdf"), fig7_volcano,
       width = 12, height = 5, dpi = 300)
cat("Saved: fig7_volcano_gene.pdf\n")

# ============================================================
# FIG 8: CHIP COUNT DISTRIBUTION (violin + boxplot)
# ============================================================
fig8_data <- cov %>%
    mutate(
        Carrier = case_when(
            BRCA1_Case == 1 & BRCA2_Case == 0 ~ "BRCA1",
            BRCA2_Case == 1 & BRCA1_Case == 0 ~ "BRCA2",
            TRUE ~ "Non-carrier"
        ),
        Carrier = factor(Carrier, levels = c("BRCA1", "BRCA2", "Non-carrier"))
    ) %>%
    filter(CHIP_Count > 0)   # show only CHIP+ individuals for count distribution

fig8_count <- ggplot(fig8_data, aes(x = Carrier, y = CHIP_Count, fill = Carrier, color = Carrier)) +
    geom_violin(alpha = 0.3, trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.25, size = 1) +
    scale_fill_manual(values  = c("BRCA1" = COL_B1, "BRCA2" = COL_B2, "Non-carrier" = COL_CTRL)) +
    scale_color_manual(values = c("BRCA1" = COL_B1, "BRCA2" = COL_B2, "Non-carrier" = COL_CTRL)) +
    scale_y_continuous(breaks = 1:max(cov$CHIP_Count, na.rm = TRUE)) +
    labs(
        title    = "CHIP Variant Count Among CHIP-Positive Individuals",
        subtitle = "CHIP+ only | violin + box + jitter",
        x        = NULL,
        y        = "Number of CHIP Variants"
    ) +
    THEME_BASE +
    theme(legend.position = "none")

ggsave(file.path(SAVE_DIR, "fig8_chip_count_distribution.pdf"), fig8_count,
       width = 6, height = 5, dpi = 300)
cat("Saved: fig8_chip_count_distribution.pdf\n")

# ============================================================
# FIG 9: MODEL COMPARISON — OR stability panel
# ============================================================
# Shows that all model specifications give consistent ORs
# (reassures readers the result isn't artifact of one model)
fig9_data <- binary_results %>%
    filter(!grepl("firth_failed", spec)) %>%
    mutate(
        spec_lab = factor(spec_labels[spec], levels = spec_labels),
        model    = factor(model, levels = c("BRCA1/2", "BRCA1", "BRCA2"))
    )

fig9_stability <- ggplot(fig9_data, aes(x = spec_lab, y = OR, ymin = CI_lo, ymax = CI_hi,
                                         color = model, group = model)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.6, alpha = 0.5, position = position_dodge(width = 0.3)) +
    geom_errorbar(width = 0.15, linewidth = 0.8, position = position_dodge(width = 0.3)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    scale_color_manual(values = c("BRCA1/2" = COL_B12, "BRCA1" = COL_B1, "BRCA2" = COL_B2)) +
    scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
    labs(
        title    = "Robustness Check: OR Stability Across Model Specifications",
        subtitle = "All specifications give consistent estimates",
        x        = "Model Specification",
        y        = "Odds Ratio (log scale)",
        color    = "Carrier Group"
    ) +
    THEME_BASE +
    theme(
        axis.text.x  = element_text(angle = 30, hjust = 1, size = 9),
        legend.position = "right"
    )

ggsave(file.path(SAVE_DIR, "fig9_or_stability.pdf"), fig9_stability,
       width = 9, height = 5, dpi = 300)
cat("Saved: fig9_or_stability.pdf\n")

# ============================================================
# SUPPLEMENTARY: GENE-SPECIFIC BRCA1 vs BRCA2 SIDE-BY-SIDE
# ============================================================
gene_b1b2 <- gene_res %>%
    filter(model %in% c("BRCA1", "BRCA2"), !is.na(OR)) %>%
    mutate(Gene = fct_reorder(Gene, OR))

figS_b1b2 <- ggplot(gene_b1b2,
    aes(x = OR, xmin = pmax(CI_lo, 0.05), xmax = pmin(CI_hi, 20),
        y = Gene, color = model)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(height = 0), linewidth = 0.7, position = position_dodge(width = 0.5)) +
    geom_point(aes(size = n_chip), alpha = 0.85, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("BRCA1" = COL_B1, "BRCA2" = COL_B2)) +
    scale_size_continuous(range = c(1.5, 5), name = "N CHIP carriers") +
    scale_x_log10() +
    labs(
        title    = "Gene-Specific CHIP: BRCA1 vs BRCA2",
        subtitle = "Firth logistic | CI clipped at [0.05, 20]",
        x        = "Odds Ratio (log scale)",
        y        = "CHIP Gene",
        color    = "Exposure"
    ) +
    THEME_BASE +
    theme(axis.text.y = element_text(size = 8))

h2 <- max(5, length(unique(gene_b1b2$Gene)) * 0.35 + 1.5)
ggsave(file.path(SAVE_DIR, "figS_gene_brca1_vs_brca2.pdf"), figS_b1b2,
       width = 9, height = h2, dpi = 300)
cat("Saved: figS_gene_brca1_vs_brca2.pdf\n")

cat("\nDone: 04_figures.R\n")
cat("All figures saved to:", SAVE_DIR, "\n")

