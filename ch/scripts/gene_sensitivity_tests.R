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
library(survival)
library(broom)
library(sandwich)
library(lmtest)
# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "chip_cov_minad4.xlsx"))
vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx")) %>%
    filter(Sample.ID %in% cov$person_id)

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
# CONFIG
# ============================================================
MIN_CARRIERS  <- 3

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
# GENE FREQUENCY TABLE
# ============================================================
gene_person <- vars_cov %>%
    dplyr::select(Sample.ID, Gene, BRCA12_Case) %>%
    distinct()

n_carrier_total    <- sum(cov$BRCA12_Case == 1)
n_noncarrier_total <- sum(cov$BRCA12_Case == 0)

gene_freq_strat <- gene_person %>%
    left_join(
        cov %>% dplyr::select(person_id, BRCA1_Case, BRCA2_Case),
        by = c("Sample.ID" = "person_id")
    ) %>%
    group_by(Gene) %>%
    summarise(
        n_noncarrier   = sum(BRCA12_Case == 0),
        n_brca12       = sum(BRCA12_Case == 1),
        n_brca1        = sum(BRCA1_Case  == 1),
        n_brca2        = sum(BRCA2_Case  == 1),
        pct_noncarrier = n_noncarrier / n_noncarrier_total * 100,
        pct_brca12     = n_brca12     / n_carrier_total    * 100,
        pct_brca1      = n_brca1      / sum(cov$BRCA1_Case == 1) * 100,
        pct_brca2      = n_brca2      / sum(cov$BRCA2_Case == 1) * 100,
        n_total        = n_noncarrier + n_brca12,
        .groups        = "drop"
    ) %>%
    arrange(desc(n_total))

# ============================================================
# FIRTH LOGISTIC PER GENE
# ============================================================
source(here("ch", "scripts", "gene_sensitivity_models.R"))

### FINAL VERSION
res_all <- run_gene_final(cov, vars, cohort_label = "all")

res_all %>%
    filter(cov_spec == "full", n_co >= 3) %>%
    mutate(
        exposure_label = factor(exposure_label, levels = c("BRCA1/2", "BRCA1", "BRCA2")),
        OR_fmt = sprintf("%.2f (%.2f–%.2f)", OR, CI_lo, CI_hi)
    ) %>%
    dplyr::select(
        Gene,
        Exposure      = exposure_label,
        CHIP_Count    = n_chip,
        CHIP_Cases    = n_co,
        OR,
        P             = p,
        `OR (95% CI)` = OR_fmt,
        P_FDR         = p_fdr
    ) %>%
    arrange(Exposure, P_FDR) %>%
    writexl::write_xlsx(file.path("ch", "data", "ch_gene_firth_final_table.xlsx"))

write_xlsx(res_all, file.path("ch", "data", "ch_all_gene_firth_final.xlsx"))

plot_final_brca12 <- function(results, fig_dir, cohort_label) {

    df <- results %>%
        filter(cov_spec == "full", exposure_label == "BRCA1/2") %>%
        mutate(
            sig   = p_fdr < 0.05,
            shape = if_else(sig, "FDR < 0.05", "n.s."),
            annot = sprintf("%.2f (%.2f–%.2f)  %s",
                            OR, CI_lo, CI_hi,
                            ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)))
        )

    gene_order <- df %>%
        arrange(OR) %>%
        pull(Gene)

    df <- df %>%
        mutate(Gene = factor(Gene, levels = gene_order))

    n_genes     <- length(gene_order)
    plot_height <- max(3, n_genes * 0.2)
    x_max       <- max(df$CI_hi, na.rm = TRUE)

    p <- ggplot(df, aes(x = OR, y = Gene,
                        xmin = CI_lo, xmax = CI_hi,
                        shape = shape, color = shape)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
        geom_errorbarh(height = 0.35, linewidth = 0.6) +
        geom_point(size = 3) +
        geom_text(
            aes(x = x_max * 1.12, y = Gene, label = annot),
            hjust = 0, size = 2.4, color = "grey25", inherit.aes = FALSE,
            data = df
        ) +
        scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                      labels = c("0.25", "0.5", "1", "2", "4", "8")) +
        scale_shape_manual(values = c("FDR < 0.05" = 17, "n.s." = 16), name = NULL) +
        scale_color_manual(values = c("FDR < 0.05" = "#b2182b", "n.s." = "#2166ac"), name = NULL) +
        coord_cartesian(clip = "off") +
        labs(x        = "Odds Ratio (log scale)",
            y        = NULL
        ) +
        theme_classic(base_size = 11) +
        theme(
            plot.title    = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            plot.margin   = margin(5, 200, 5, 5),
            legend.position = "bottom",
            axis.text.y   = element_text(size = 9),
            axis.text.x   = element_text(size = 9)
        )

    out_path <- file.path(fig_dir,
                          sprintf("fig_gene_brca12_full_%s.pdf", cohort_label))
    ggsave(out_path, p,
           width  = 8,
           height = plot_height,
           device = cairo_pdf)

    message("Saved: ", out_path)
    invisible(p)
}
plot_final_brca12(res_all, fig_dir = file.path("ch", "figures"), cohort_label = "all")

# counts of those significant enough
gene_freq %>% filter(Gene %in% c("GATA1", "ASXL1", "DNMT3A", "TET2", "NF1", "KMT2A"))
# Gene   n_noncarrier n_carrier pct_noncarrier pct_carrier n_total
# <chr>         <int>     <int>          <dbl>       <dbl>   <int>
#     1 DNMT3A           40        17         1.72         2.51       57
# 2 TET2             23        13         0.988        1.92       36
# 3 ASXL1            12         4         0.516        0.591      16
# 4 NF1               6         3         0.258        0.443       9
# 5 GATA1             3         3         0.129        0.443       6
# 6 KMT2A             1         3         0.0430       0.443       4



