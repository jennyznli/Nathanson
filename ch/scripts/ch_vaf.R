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
vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))

# ============================================================
# RESTRICT TO CHIP+ INDIVIDUALS
# ============================================================
chip_pos <- cov %>% filter(CHIP_Binary == TRUE)
cat("CHIP+ individuals:", nrow(chip_pos), "\n")
cat("  Carriers:     ", sum(chip_pos$BRCA12_Case == 1), "\n")
cat("  Non-carriers: ", sum(chip_pos$BRCA12_Case == 0), "\n")

# Join variant data to CHIP+ individuals
vars_chip <- vars %>%
    inner_join(cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                                     BRCA2_Case, Carrier),
               by = c("Sample.ID" = "person_id"))

# Per-person summary: max VAF and total mutation count
person_summary <- vars_chip %>%
    group_by(Sample.ID, BRCA12_Case, BRCA1_Case, BRCA2_Case, Carrier, Sample_age) %>%
    summarise(
        max_vaf  = max(.data[["Sample.AltFrac"]], na.rm = TRUE),
        n_muts   = n(),
        .groups  = "drop"
    )

# ============================================================
# VAF 1
# ============================================================
cat("\n=== V1: Max VAF by carrier status ===\n")
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

cat("\n=== V2: VAF category distribution by carrier status ===\n")
vaf_cat_tab <- table(vars_chip$vaf_cat, vars_chip$BRCA12_Case)
print(vaf_cat_tab)
# 0  1
# <2%    0  0
# 2-10% 84 62
# >10%  61 10

cat("Chi-square test:\n")
print(chisq.test(vaf_cat_tab))

# ============================================================
# V3: MUTATION BURDEN (n mutations per CHIP+ person)
# ============================================================
cat("\n=== V3: Mutation burden among CHIP+ ===\n")
burden_summary <- person_summary %>%
    group_by(BRCA12_Case) %>%
    summarise(
        n            = n(),
        median_count = median(n_muts),
        q25          = quantile(n_muts, 0.25),
        q75          = quantile(n_muts, 0.75),
        pct_multi    = mean(n_muts > 1) * 100,
        .groups      = "drop"
    )
print(burden_summary)

wt_burden <- wilcox.test(n_muts ~ BRCA12_Case, data = person_summary, exact = FALSE)
cat("Wilcoxon p-value (mutation burden):", wt_burden$p.value, "\n")

# BRCA12_Case     n median_count   q25   q75 pct_multi
# <dbl> <int>        <dbl> <dbl> <dbl>     <dbl>
#     1           0   138            1     1     1      4.35
# 2           1    55            1     1     1     14.5

# ============================================================
# V4: LINEAR REGRESSION — VAF ~ carrier status + age
# (among CHIP+ individuals, controls for age)
# ============================================================
cat("\n=== V4: Linear regression — max VAF ~ BRCA12_Case + age ===\n")
lm_vaf <- lm(max_vaf ~ BRCA12_Case + Sample_age, data = person_summary)
print(summary(lm_vaf))

lm_vaf_b1 <- lm(max_vaf ~ BRCA1_Case + Sample_age, data = person_summary)
lm_vaf_b2 <- lm(max_vaf ~ BRCA2_Case + Sample_age, data = person_summary)
cat("\nBRCA1 coefficient (p-value):",
    coef(summary(lm_vaf_b1))["BRCA1_Case", "Pr(>|t|)"], "\n")
cat("BRCA2 coefficient (p-value):",
    coef(summary(lm_vaf_b2))["BRCA2_Case", "Pr(>|t|)"], "\n")

write_xlsx(
    list(
        person_summary = person_summary,
        vaf_by_group   = vaf_summary,
        vaf_categories = as.data.frame(vaf_cat_tab)
    ),
    file.path("ch", "data", "ch_vaf_clonal_results.xlsx")
)

cat("\nDone: 04_vaf_clonal.R\n")

# ============================================================
# PLOT VIOLIN?
# ============================================================
library(ggplot2)
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
vaf_bar <- vars_chip %>%
    count(carrier_binary, vaf_cat) %>%
    group_by(carrier_binary) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

vaf_bar3 <- vars_chip %>%
    count(carrier_label, vaf_cat) %>%
    group_by(carrier_label) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

bar_colors <- c("<2%" = "#f7f7f7", "2-10%" = "#92c5de", ">10%" = "#2166ac")

p3 <- ggplot(vaf_bar, aes(x = carrier_binary, y = pct, fill = vaf_cat)) +
    geom_col(width = 0.5, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f%%", pct)),
              position = position_stack(vjust = 0.5),
              size = 3, color = "grey20") +
    scale_fill_manual(values = bar_colors, name = "VAF") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(title    = "VAF category distribution per mutation",
         subtitle = "Per CHIP variant (not per person)",
         x = NULL, y = "% of mutations") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text     = element_text(color = "black"))

# same but split by BRCA1/BRCA2
p4 <- ggplot(vaf_bar3, aes(x = carrier_label, y = pct, fill = vaf_cat)) +
    geom_col(width = 0.5, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f%%", pct)),
              position = position_stack(vjust = 0.5),
              size = 3, color = "grey20") +
    scale_fill_manual(values = bar_colors, name = "VAF") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(title    = "VAF category distribution per mutation",
         subtitle = "Per CHIP variant (not per person)",
         x = NULL, y = "% of mutations") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text     = element_text(color = "black"))

# ── save ─────────────────────────────────────────────────────
ggsave(file.path("ch", "figures", "fig_vaf_violin_binary.pdf"),  p1, width = 5, height = 5)
ggsave(file.path("ch", "figures", "fig_vaf_violin_3group.pdf"),  p2, width = 6, height = 5)
ggsave(file.path("ch", "figures", "fig_vaf_bar_binary.pdf"),     p3, width = 5, height = 5)
ggsave(file.path("ch", "figures", "fig_vaf_bar_3group.pdf"),     p4, width = 6, height = 5)


# ============================================================
# BY GENE
# ============================================================


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

# ── P5: violin by gene — binary carrier ──────────────────────
p5 <- ggplot(person_gene,
             aes(x = Gene, y = max_vaf, fill = carrier_binary)) +
    geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.3,
                position = position_dodge(0.8)) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
                 linewidth = 0.4, position = position_dodge(0.8)) +
    geom_beeswarm(aes(color = carrier_binary), size = 1.3, alpha = 0.7,
                  dodge.width = 0.8) +
    scale_fill_manual(values  = group_colors,  name = NULL) +
    scale_color_manual(values = group_colors,  name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title    = "Max VAF by gene and carrier status",
         subtitle = "CHIP+ individuals only | genes with ≥5 CHIP+",
         x = NULL, y = "Max VAF") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y   = element_text(color = "black"),
          legend.position = "bottom")

# ── P6: violin by gene — BRCA1 / BRCA2 / non-carrier ─────────
p6 <- ggplot(person_gene,
             aes(x = Gene, y = max_vaf, fill = carrier_label)) +
    geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.3,
                position = position_dodge(0.9)) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
                 linewidth = 0.4, position = position_dodge(0.9)) +
    geom_beeswarm(aes(color = carrier_label), size = 1.3, alpha = 0.7,
                  dodge.width = 0.9) +
    scale_fill_manual(values  = group_colors3, name = NULL) +
    scale_color_manual(values = group_colors3, name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title    = "Max VAF by gene and carrier type",
         subtitle = "CHIP+ individuals only | genes with ≥5 CHIP+",
         x = NULL, y = "Max VAF") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y   = element_text(color = "black"),
          legend.position = "bottom")

# ── save ──────────────────────────────────────────────────────
n_genes <- length(gene_keep)

ggsave(file.path("ch", "figures", "fig_vaf_gene_violin_binary.pdf"),
       p5, width = max(8, n_genes * 0.7), height = 5)
ggsave(file.path("ch", "figures", "fig_vaf_gene_violin_3group.pdf"),
       p6, width = max(8, n_genes * 0.8), height = 5)

library(patchwork)

# overall gene violin (no carrier stratification)
gene_order_overall <- person_gene %>%
    group_by(Gene) %>%
    summarise(med = median(max_vaf), .groups = "drop") %>%
    arrange(med) %>%
    pull(Gene)

person_gene_overall <- person_gene %>%
    mutate(Gene = factor(Gene, levels = gene_order_overall))

p_gene_overall <- ggplot(person_gene_overall,
                         aes(x = Gene, y = max_vaf)) +
    geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.3,
                fill = "#636363") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 linewidth = 0.4) +
    geom_beeswarm(size = 1.3, alpha = 0.5, color = "#636363") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title    = "Max VAF by gene",
         subtitle = "CHIP+ individuals only | genes with ≥5 CHIP+",
         x = NULL, y = "Max VAF") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y   = element_text(color = "black"))

ggsave(file.path("ch", "figures", "fig_vaf_gene_violin_overall.pdf"),
       p_gene_overall, width = max(8, n_genes * 0.7), height = 5)



