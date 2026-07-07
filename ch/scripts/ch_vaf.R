# ============================================================
# VAF & clonal characteristics — carrier vs non-carrier
# (restricted to CHIP-positive individuals)
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(coin)

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "chip_cov_minad4.xlsx"))
cov_f <- cov %>% filter(Sequenced_gender == "Female")
vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))
vars_chip <- vars %>% inner_join(cov %>% dplyr::select(-Sample_age, -Batch), by = c("Sample.ID" = "person_id"))

chip_pos <- cov %>% filter(CHIP_Binary == TRUE)
cat("CHIP+ individuals:", nrow(chip_pos), "\n")
# 193
cat("  Carriers:     ", sum(chip_pos$BRCA12_Case == 1), "\n")
# 55
cat("  Non-carriers: ", sum(chip_pos$BRCA12_Case == 0), "\n")
# 138

person_summary <- vars_chip %>%
    group_by(Sample.ID, BRCA12_Case, BRCA1_Case, BRCA2_Case, Carrier, Sample_age, Sample_age2,
             PC1, PC2, PC3, PC4, PC5, PC6, Sequenced_gender, Smoke_History, Batch) %>%
    summarise(
        max_vaf  = max(.data[["Sample.AltFrac"]], na.rm = TRUE),
        n_muts   = n(),
        .groups  = "drop"
    )

# ============================================================
# VAF 1
# ============================================================
vaf_summary <- person_summary %>%
    group_by(BRCA12_Case) %>%
    summarise(
        n         = n(),
        median_vaf = median(max_vaf),
        q25_vaf   = quantile(max_vaf, 0.25),
        q75_vaf   = quantile(max_vaf, 0.75),
        mean_vaf  = mean(max_vaf),
        .groups   = "drop"
    )
print(vaf_summary)
# BRCA12_Case     n median_vaf q25_vaf q75_vaf mean_vaf
# <dbl> <int>      <dbl>   <dbl>   <dbl>    <dbl>
#     1           0   138      0.085  0.0662  0.153    0.134
# 2           1    55      0.074  0.057   0.0935   0.0863

# Wilcoxon test
wt_vaf <- wilcox.test(max_vaf ~ BRCA12_Case, data = person_summary, exact = FALSE)
cat("Wilcoxon p-value (max VAF, carrier vs non-carrier):", wt_vaf$p.value, "\n")
# Wilcoxon p-value (max VAF, carrier vs non-carrier): 0.009976294

# Also BRCA1 vs BRCA2 vs non-carrier (3-group Kruskal-Wallis)
kw_vaf <- kruskal.test(max_vaf ~ Carrier, data = person_summary)
cat("Kruskal-Wallis p-value (max VAF, by carrier type):", kw_vaf$p.value, "\n")
# Kruskal-Wallis p-value (max VAF, by carrier type): 0.03050661

# Pairwise Wilcoxon if KW significant
if (kw_vaf$p.value < 0.05) {
    cat("\nPairwise Wilcoxon (Bonferroni corrected):\n")
    pw <- pairwise.wilcox.test(person_summary$max_vaf,
                               person_summary$Carrier,
                               p.adjust.method = "bonferroni")
    print(pw)
}

### PLOT
hist_df <- person_summary %>%
    mutate(group = if_else(BRCA12_Case == 1, "Carrier", "Non-carrier"))

fig <- ggplot(hist_df, aes(x = max_vaf, fill = group)) +
    geom_histogram(bins = 30, alpha = 0.5, position = "identity", color = "white") +
    scale_fill_manual(values = c("Carrier" = "#C2185B", "Non-carrier" = "#3e77c1"),
                      name = NULL) +
    labs(title = "Max VAF distribution by BRCA1/2 carrier status",
         # subtitle = sprintf("CHIP+ individuals | Carriers: %d, Non-carriers: %d | Wilcoxon p = %.3f",
         #                    sum(hist_df$BRCA12_Case == 1),
         #                    sum(hist_df$BRCA12_Case == 0),
         #                    wt_vaf$p.value),
         x = "Max VAF", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position  = "bottom",
          plot.subtitle     = element_text(size = 9, color = "grey40"),
          plot.title        = element_text(face = "bold"))

ggsave(file.path("ch", "figures", "fig_vaf_hist_carrier_status.pdf"),
       fig, width = 7, height = 5)


# ============================================================
# V2: VAF DISTRIBUTION — small vs large clone
# categorize into <2%, 2-10%, >10% (common CHIP thresholds)
# ============================================================
vars_chip <- vars_chip %>%
    mutate(vaf_cat = case_when(
        .data[["Sample.AltFrac"]] < 0.02 ~ "<2%",
        .data[["Sample.AltFrac"]] < 0.10 ~ "2-10%",
        TRUE                    ~ ">10%"
    ) %>% factor(levels = c("<2%", "2-10%", ">10%")))

vaf_cat_tab <- table(vars_chip$vaf_cat, vars_chip$BRCA12_Case)
vaf_cat_tab <- vaf_cat_tab[rowSums(vaf_cat_tab) > 0, ]  # drop empty rows

print(vaf_cat_tab)
# 0  1
# 2-10% 84 62
# >10%  61 10
print(chisq.test(vaf_cat_tab))
# X-squared = 16.099, df = 1, p-value = 6.012e-05

# ============================================================
# V4: LINEAR REGRESSION — VAF ~ carrier status + age
# (among CHIP+ individuals, controls for age)
# ============================================================
full_covs <- "Sample_age + Sample_age2 + Batch + Smoke_History + Sequenced_gender +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

lm_vaf    <- lm(as.formula(paste("max_vaf ~ BRCA12_Case +", full_covs)), data = person_summary)
lm_vaf_b1 <- lm(as.formula(paste("max_vaf ~ BRCA1_Case  +", full_covs)), data = person_summary)
lm_vaf_b2 <- lm(as.formula(paste("max_vaf ~ BRCA2_Case  +", full_covs)), data = person_summary)

cat("\nBRCA1/2 coefficient (p-value):",
    coef(summary(lm_vaf))["BRCA12_Case", "Pr(>|t|)"], "\n")
cat("BRCA1 coefficient (p-value):",
    coef(summary(lm_vaf_b1))["BRCA1_Case", "Pr(>|t|)"], "\n")
cat("BRCA2 coefficient (p-value):",
    coef(summary(lm_vaf_b2))["BRCA2_Case", "Pr(>|t|)"], "\n")


# ── extract lm coefficients ───────────────────────────────────
extract_lm <- function(fit, term) {
    ct <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
    ci <- coefci(fit, vcov = vcovHC(fit, type = "HC3"))
    data.frame(
        est   = ct[term, "Estimate"],
        CI_lo = ci[term, 1],
        CI_hi = ci[term, 2],
        p     = ct[term, "Pr(>|t|)"]
    )
}

lm_results <- bind_rows(
    extract_lm(lm_vaf,    "BRCA12_Case") %>% mutate(exposure = "BRCA1/2"),
    extract_lm(lm_vaf_b1, "BRCA1_Case")  %>% mutate(exposure = "BRCA1"),
    extract_lm(lm_vaf_b2, "BRCA2_Case")  %>% mutate(exposure = "BRCA2")
) %>%
    mutate(
        sig    = p < 0.05,
        shape  = if_else(sig, "p<0.05", "ns"),
        p_fmt  = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
        annot  = sprintf("%.3f (%.3f–%.3f) %s", est, CI_lo, CI_hi, p_fmt),
        exposure = factor(exposure, levels = c("BRCA2", "BRCA1", "BRCA1/2"))
    )

fig_lm <- ggplot(lm_results,
                 aes(x = est, y = exposure, xmin = CI_lo, xmax = CI_hi,
                     shape = shape)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbarh(height = 0.2, linewidth = 0.6, color = "#2166ac") +
    geom_point(size = 3, color = "#2166ac") +
    geom_text(aes(x = max(lm_results$CI_hi) * 1.15 + 0.01, label = annot),
              hjust = 0, size = 2.8, color = "grey30") +
    scale_shape_manual(values = c("p<0.05" = 17, "ns" = 16), name = NULL) +
    coord_cartesian(clip = "off") +
    labs(title = "Max VAF ~ carrier status (linear regression, CHIP+ only)",
         x = "Coefficient (change in max VAF)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(face = "bold"),
          plot.margin     = margin(5, 210, 5, 5),
          legend.position = "bottom")

ggsave(file.path("ch", "figures", "fig_lm_vaf_carrier.pdf"),
       fig_lm, width = 8, height = 3.5)


# ============================================================
# PLOT VIOLIN
# ============================================================
library(ggbeeswarm)  # for nicer jitter

# ── shared carrier label ──────────────────────────────────────
person_summary <- person_summary %>%
    mutate(carrier_label = case_when(
        BRCA1_Case == 1 ~ "BRCA1",
        BRCA2_Case == 1 ~ "BRCA2",
        TRUE            ~ "Non-carrier"
    ) %>% factor(levels = c("Non-carrier", "BRCA1", "BRCA2")),
    carrier_binary = ifelse(BRCA12_Case == 1, "BRCA1/2", "Non-carrier") %>%
        factor(levels = c("Non-carrier", "BRCA1/2"))
    )

vars_chip <- vars_chip %>%
    mutate(carrier_label = case_when(
        BRCA1_Case == 1 ~ "BRCA1",
        BRCA2_Case == 1 ~ "BRCA2",
        TRUE            ~ "Non-carrier"
    ) %>% factor(levels = c("Non-carrier", "BRCA1", "BRCA2")),
    carrier_binary = ifelse(BRCA12_Case == 1, "BRCA1/2", "Non-carrier") %>%
        factor(levels = c("Non-carrier", "BRCA1/2"))
    )

group_colors  <- c("Non-carrier" = "#3e77c1", "BRCA1/2" = "#C2185B")
group_colors3 <- c("Non-carrier" = "#3e77c1", "BRCA1" = "#C2185B", "BRCA2" = "#C2185B")

# ── P1: violin — max VAF, BRCA1/2 vs non-carrier ─────────────
wt_p  <- format.pval(wt_vaf$p.value, digits = 2)

p1 <- ggplot(person_summary,
             aes(x = carrier_binary, y = max_vaf, fill = carrier_binary)) +
    geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 fill = "white", linewidth = 0.5) +
    scale_fill_manual(values  = group_colors, guide = "none") +
    scale_color_manual(values = group_colors, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    annotate("text", x = 1.5, y = max(person_summary$max_vaf) * 1.05,
             label = paste0("Wilcoxon p = ", wt_p), size = 3.2, color = "grey30") +
    labs(x = NULL, y = "Max VAF") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text     = element_text(color = "black"))

# ── P2: violin — max VAF, BRCA1 / BRCA2 / non-carrier ────────
kw_p <- format.pval(kw_vaf$p.value, digits = 2)

p2 <- ggplot(person_summary,
             aes(x = carrier_label, y = max_vaf, fill = carrier_label)) +
    geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 fill = "white", linewidth = 0.5) +
    scale_fill_manual(values  = group_colors3, guide = "none") +
    scale_color_manual(values = group_colors3, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    annotate("text", x = 2, y = max(person_summary$max_vaf) * 1.05,
             label = paste0("Kruskal-Wallis p = ", kw_p), size = 3.2, color = "grey30") +
    labs(x = NULL, y = "Max VAF") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text     = element_text(color = "black"))

# ── P3: stacked bar — VAF categories, per mutation ───────────
# vaf_bar <- vars_chip %>%
#     count(carrier_binary, vaf_cat) %>%
#     group_by(carrier_binary) %>%
#     mutate(pct = n / sum(n) * 100) %>%
#     ungroup()
#
# vaf_bar3 <- vars_chip %>%
#     count(carrier_label, vaf_cat) %>%
#     group_by(carrier_label) %>%
#     mutate(pct = n / sum(n) * 100) %>%
#     ungroup()
#
# bar_colors <- c("<2%" = "#f7f7f7", "2-10%" = "#92c5de", ">10%" = "#2166ac")
#
# p3 <- ggplot(vaf_bar, aes(x = carrier_binary, y = pct, fill = vaf_cat)) +
#     geom_col(width = 0.5, color = "white", linewidth = 0.3) +
#     geom_text(aes(label = sprintf("%.0f%%", pct)),
#               position = position_stack(vjust = 0.5),
#               size = 3, color = "grey20") +
#     scale_fill_manual(values = bar_colors, name = "VAF") +
#     scale_y_continuous(labels = scales::percent_format(scale = 1)) +
#     labs(title    = "VAF category distribution per mutation",
#          subtitle = "Per CHIP variant (not per person)",
#          x = NULL, y = "% of mutations") +
#     theme_classic(base_size = 12) +
#     theme(plot.title    = element_text(face = "bold"),
#           plot.subtitle = element_text(size = 9, color = "grey40"),
#           axis.text     = element_text(color = "black"))
#
# # same but split by BRCA1/BRCA2
# p4 <- ggplot(vaf_bar3, aes(x = carrier_label, y = pct, fill = vaf_cat)) +
#     geom_col(width = 0.5, color = "white", linewidth = 0.3) +
#     geom_text(aes(label = sprintf("%.0f%%", pct)),
#               position = position_stack(vjust = 0.5),
#               size = 3, color = "grey20") +
#     scale_fill_manual(values = bar_colors, name = "VAF") +
#     scale_y_continuous(labels = scales::percent_format(scale = 1)) +
#     labs(title    = "VAF category distribution per mutation",
#          subtitle = "Per CHIP variant (not per person)",
#          x = NULL, y = "% of mutations") +
#     theme_classic(base_size = 12) +
#     theme(plot.title    = element_text(face = "bold"),
#           plot.subtitle = element_text(size = 9, color = "grey40"),
#           axis.text     = element_text(color = "black"))

# ── save ─────────────────────────────────────────────────────
ggsave(file.path("ch", "figures", "fig_vaf_violin_binary.pdf"),  p1, width = 5, height = 5)
ggsave(file.path("ch", "figures", "fig_vaf_violin_3group.pdf"),  p2, width = 6, height = 5)
ggsave(file.path("ch", "figures", "fig_vaf_bar_binary.pdf"),     p3, width = 5, height = 5)
ggsave(file.path("ch", "figures", "fig_vaf_bar_3group.pdf"),     p4, width = 6, height = 5)


# ============================================================
# JUST VAF BY GENE
# ============================================================
top20_genes <- vars_chip %>%
    count(Gene) %>%
    slice_max(n, n = 20) %>%
    pull(Gene)

gene_order_top20 <- vars_chip %>%
    filter(Gene %in% top20_genes) %>%
    group_by(Gene) %>%
    summarise(med = median(Sample.AltFrac), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(Gene)

p_gene_top20 <- vars_chip %>%
    filter(Gene %in% top20_genes) %>%
    mutate(Gene = factor(Gene, levels = gene_order_top20)) %>%
    ggplot(aes(x = Gene, y = Sample.AltFrac)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey90",
                 linewidth = 0.4, color = "grey30") +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "#636363") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Gene", y = "VAF") +
    theme_classic(base_size = 11) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"))

ggsave(file.path("ch", "figures", "fig_vaf_gene_top20.pdf"),
       p_gene_top20, width = 10, height = 5)



# ── prep gene-level VAF data ──────────────────────────────────
# per person: max VAF per gene (since one person can have multiple genes)
person_gene <- vars_chip %>%
    group_by(Sample.ID, Gene, BRCA12_Case, BRCA1_Case, BRCA2_Case,
             carrier_label, carrier_binary) %>%
    summarise(max_vaf = max(Sample.AltFrac, na.rm = TRUE),
              n_muts  = n(), .groups = "drop")

# filter genes with >= 5 CHIP+ individuals total
gene_keep <- person_gene %>%
    count(Gene) %>%
    filter(n >= 5) %>%
    pull(Gene)

person_gene <- person_gene %>% filter(Gene %in% gene_keep)
cat("Genes with >= 5 CHIP+ individuals:", length(gene_keep), "\n")
cat("Genes:", paste(sort(gene_keep), collapse = ", "), "\n")

# order genes by median VAF (non-carriers) for consistent x-axis
gene_order <- person_gene %>%
    filter(BRCA12_Case == 0) %>%
    group_by(Gene) %>%
    summarise(med = median(max_vaf), .groups = "drop") %>%
    arrange(med) %>%
    pull(Gene)

person_gene <- person_gene %>%
    mutate(Gene = factor(Gene, levels = gene_order))

# # ── P5: violin by gene — binary carrier ──────────────────────
# p5 <- ggplot(person_gene,
#              aes(x = Gene, y = max_vaf, fill = carrier_binary)) +
#     geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.3,
#                 position = position_dodge(0.8)) +
#     geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
#                  linewidth = 0.4, position = position_dodge(0.8)) +
#     geom_beeswarm(aes(color = carrier_binary), size = 1.3, alpha = 0.7,
#                   dodge.width = 0.8) +
#     scale_fill_manual(values  = group_colors,  name = NULL) +
#     scale_color_manual(values = group_colors,  name = NULL) +
#     scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#     labs(title    = "Max VAF by gene and carrier status",
#          subtitle = "CHIP+ individuals only | genes with ≥5 CHIP+",
#          x = NULL, y = "Max VAF") +
#     theme_classic(base_size = 11) +
#     theme(plot.title    = element_text(face = "bold"),
#           plot.subtitle = element_text(size = 9, color = "grey40"),
#           axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"),
#           axis.text.y   = element_text(color = "black"),
#           legend.position = "bottom")
#
# # ── P6: violin by gene — BRCA1 / BRCA2 / non-carrier ─────────
# p6 <- ggplot(person_gene,
#              aes(x = Gene, y = max_vaf, fill = carrier_label)) +
#     geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.3,
#                 position = position_dodge(0.9)) +
#     geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
#                  linewidth = 0.4, position = position_dodge(0.9)) +
#     geom_beeswarm(aes(color = carrier_label), size = 1.3, alpha = 0.7,
#                   dodge.width = 0.9) +
#     scale_fill_manual(values  = group_colors3, name = NULL) +
#     scale_color_manual(values = group_colors3, name = NULL) +
#     scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#     labs(title    = "Max VAF by gene and carrier type",
#          subtitle = "CHIP+ individuals only | genes with ≥5 CHIP+",
#          x = NULL, y = "Max VAF") +
#     theme_classic(base_size = 11) +
#     theme(plot.title    = element_text(face = "bold"),
#           plot.subtitle = element_text(size = 9, color = "grey40"),
#           axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"),
#           axis.text.y   = element_text(color = "black"),
#           legend.position = "bottom")
#
# # ── save ──────────────────────────────────────────────────────
# n_genes <- length(gene_keep)
#
# ggsave(file.path("ch", "figures", "fig_vaf_gene_violin_binary.pdf"),
#        p5, width = max(8, n_genes * 0.7), height = 5)
# ggsave(file.path("ch", "figures", "fig_vaf_gene_violin_3group.pdf"),
#        p6, width = max(8, n_genes * 0.8), height = 5)

library(patchwork)

# overall gene violin (no carrier stratification)
gene_order_overall <- person_gene %>%
    group_by(Gene) %>%
    filter(n() >= 3) %>%
    summarise(med = median(max_vaf), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(Gene)

person_gene_overall <- person_gene %>%
    filter(Gene %in% gene_order_overall) %>%
    mutate(Gene = factor(Gene, levels = gene_order_overall))
p_gene_overall <- ggplot(person_gene_overall, aes(x = Gene, y = max_vaf)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey90",
                 linewidth = 0.4, color = "grey30") +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "#636363") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title    = "Max VAF by gene",
         # subtitle = "CHIP+ individuals only | genes with ≥5 CHIP+",
         x = NULL, y = "Max VAF") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y   = element_text(color = "black"))

ggsave(file.path("ch", "figures", "fig_vaf_gene_violin_overall.pdf"),
       p_gene_overall, width = max(8, n_genes * 0.7), height = 5)



