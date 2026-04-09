# ============================================================
# CHIP PREVALENCE - UNAFFECTED MATCHED COHORT
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(writexl)

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx"))
vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx")) %>%
    filter(Sample.ID %in% cov$person_id)

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

n_total         <- nrow(cov)
n_carrier_total    <- sum(cov$BRCA12_Case == 1)
n_noncarrier_total <- sum(cov$BRCA12_Case == 0)

# ============================================================
# SHARED THEME
# ============================================================
theme_ch <- theme_minimal(base_size = 12) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(hjust = 0.5, color = "grey40", size = 10),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom"
    )

# ============================================================
# LOG - PREVALENCE SUMMARY
# ============================================================
sink(file.path("ch", "data", "chip_unaff_prevalence_summary.log"), split = TRUE)
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("=== OVERALL ===\n")
cat("N =", n_total,
    "| CHIP+ =", sum(cov$CHIP_Binary),
    sprintf("(%.1f%%)\n\n", 100 * mean(cov$CHIP_Binary)))

cat("=== BY CARRIER STATUS ===\n")
cov %>%
    group_by(BRCA12_Case) %>%
    summarise(
        N      = n(),
        n_chip = sum(CHIP_Binary),
        pct    = round(100 * mean(CHIP_Binary), 1)
    ) %>% print()
cat("\n")

cov %>%
    group_by(Carrier) %>%
    summarise(
        N      = n(),
        n_chip = sum(CHIP_Binary),
        pct    = round(100 * mean(CHIP_Binary), 1)
    ) %>% print()
cat("\n")

cat("=== PREVALENCE BY AGE DECADE ===\n")
age_breaks <- c(-Inf, 30, 40, 50, 60, Inf)
age_labels <- c("≤30", "30–40", "40–50", "50–60", "≥60")

prev_overall <- cov %>%
    mutate(age_group = cut(Sample_age, breaks = age_breaks,
                           labels = age_labels, right = TRUE)) %>%
    group_by(age_group) %>%
    summarise(n_total  = n(),
              n_chip   = sum(CHIP_Binary),
              prev_pct = round(100 * mean(CHIP_Binary), 1),
              .groups  = "drop") %>%
    mutate(x_num = as.numeric(age_group))

print(prev_overall)
cat("\n")

cat("=== PREVALENCE BY AGE DECADE AND CARRIER STATUS ===\n")
prev_by_carrier <- cov %>%
    mutate(age_group = cut(Sample_age, breaks = age_breaks,
                           labels = age_labels, right = TRUE),
           Carrier_label = ifelse(BRCA12_Case == 1,
                                  "gBRCA1/2 Carrier", "Non-carrier")) %>%
    group_by(age_group, Carrier_label) %>%
    summarise(n_total  = n(),
              n_chip   = sum(CHIP_Binary),
              prev_pct = round(100 * mean(CHIP_Binary), 1),
              .groups  = "drop") %>%
    mutate(x_num = as.numeric(age_group))

print(prev_by_carrier)
cat("\n")

cat("=== CHIP COUNT PER INDIVIDUAL (CHIP+ ONLY) ===\n")
chip_count_dist <- cov %>%
    filter(CHIP_Binary) %>%
    count(CHIP_Count) %>%
    mutate(
        pct   = round(100 * n / n_total, 1),
        label = paste0(n, " (", pct, "%)")
    )
print(chip_count_dist)
cat("\n")

cat("=== PREVALENCE BY COVARIATE ===\n")

# age
age_wilcox <- wilcox.test(Sample_age ~ BRCA12_Case, data = cov)
age_stats <- cov %>%
    group_by(BRCA12_Case) %>%
    summarise(median = median(Sample_age), q1 = quantile(Sample_age, 0.25),
              q3 = quantile(Sample_age, 0.75), .groups = "drop")
cat("Age (Wilcoxon p =", round(age_wilcox$p.value, 4), ")\n")
print(age_stats)
cat("\n")

# categorical
format_p <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p == 0)   return("< 2.2e-16")
    if (p < 0.001) return(formatC(p, format = "e", digits = 2))
    as.character(round(p, 3))
}

cat_vars <- c("BRCA12_Case", "Sequenced_gender", "Smoke_History",
              "Batch", "Strata", "Class")

prev_by_cov <- lapply(cat_vars, function(var) {
    df_sub <- cov[!is.na(cov[[var]]), ]
    tab    <- table(df_sub[[var]], df_sub$CHIP_Binary)
    expected <- suppressWarnings(chisq.test(tab)$expected)

    if (any(expected < 5) || nrow(tab) > 2) {
        test   <- fisher.test(tab, simulate.p.value = TRUE)
        method <- "Fisher"
    } else {
        test   <- chisq.test(tab)
        method <- "Chi-squared"
    }

    counts <- df_sub %>%
        group_by(.data[[var]]) %>%
        summarise(N = n(), n_chip = sum(CHIP_Binary),
                  prev_pct = round(100 * mean(CHIP_Binary), 1),
                  .groups = "drop") %>%
        mutate(Variable = var, P_value = format_p(test$p.value), Method = method)

    counts
}) %>% bind_rows()

print(prev_by_cov, n = 50)
cat("\n")

sink()

# ============================================================
# GENE FREQUENCY TABLE
# ============================================================
vars_cov <- vars %>%
    left_join(
        cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                              BRCA2_Case, Sample_age, Batch, Smoke_History),
        by = c("Sample.ID" = "person_id")
    ) %>%
    filter(!is.na(BRCA12_Case))

gene_person <- vars_cov %>%
    dplyr::select(Sample.ID, Gene, BRCA12_Case) %>%
    distinct()

gene_freq <- gene_person %>%
    group_by(Gene, BRCA12_Case) %>%
    summarise(n_mutated = n(), .groups = "drop") %>%
    pivot_wider(names_from  = BRCA12_Case,
                values_from = n_mutated,
                values_fill = 0,
                names_prefix = "n_") %>%
    rename(n_noncarrier = n_0, n_carrier = n_1) %>%
    mutate(
        n_total        = n_carrier + n_noncarrier,
        pct_overall    = round(n_total    / (n_carrier_total + n_noncarrier_total) * 100, 2),
        pct_carrier    = round(n_carrier    / n_carrier_total    * 100, 2),
        pct_noncarrier = round(n_noncarrier / n_noncarrier_total * 100, 2)
    ) %>%
    dplyr::select(Gene,
                  n_total, pct_overall,
                  n_carrier, pct_carrier,
                  n_noncarrier, pct_noncarrier) %>%
    arrange(desc(n_total))

top_genes <- gene_freq %>%
    slice_max(order_by = n_total, n = 15, with_ties = FALSE) %>%
    pull(Gene)

write_xlsx(gene_freq, file.path("ch", "data", "gene_frequency_table.xlsx"))

# ============================================================
# PREVALENCE PLOTS
# ============================================================
scale_y_prev <- scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.1))
)
scale_x_prev <- scale_x_continuous(
    breaks = prev_overall$x_num,
    labels = prev_overall$age_group
)

# ── Overall ──
fig_prev_overall <- ggplot(prev_overall, aes(x = x_num, y = prev_pct)) +
    geom_smooth(method = "loess", span = 1.2, se = FALSE,
                color = "#546E7A", linewidth = 0.9) +
    geom_point(size = 3, color = "#546E7A") +
    scale_x_prev + scale_y_prev +
    labs(x = "Age group (years)", y = "CHIP prevalence (%)") +
    theme_ch

ggsave(file.path("ch", "figures", "fig_chip_prevalence_by_decade.pdf"),
       fig_prev_overall, width = 6, height = 5)

# ── Stratified by carrier (with legend) ──
fig_prev_carrier <- ggplot() +
    geom_smooth(data = prev_by_carrier,
                aes(x = x_num, y = prev_pct, color = Carrier_label),
                method = "loess", span = 1.2, se = FALSE, linewidth = 0.9) +
    geom_point(data = prev_by_carrier,
               aes(x = x_num, y = prev_pct, color = Carrier_label), size = 3) +
    geom_smooth(data = prev_overall,
                aes(x = x_num, y = prev_pct, group = 1),
                method = "loess", span = 1.2, se = FALSE,
                color = "darkgray", linewidth = 0.8) +
    geom_point(data = prev_overall,
               aes(x = x_num, y = prev_pct),
               color = "darkgray", size = 3, shape = 18) +
    scale_color_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B", "Non-carrier" = "#3e77c1"),
        name = NULL
    ) +
    scale_x_prev + scale_y_prev +
    labs(x = "Age group (years)", y = "CHIP prevalence (%)") +
    theme_ch

ggsave(file.path("ch", "figures", "fig_chip_prevalence_by_decade_carrier.pdf"),
       fig_prev_carrier, width = 6, height = 5)

ggsave(file.path("ch", "figures", "fig_chip_prevalence_by_decade_carrier_noleg.pdf"),
       fig_prev_carrier + theme(legend.position = "none"), width = 6, height = 5)

# ============================================================
# GENE PREVALENCE PLOTS
# ============================================================
scale_y_gene <- scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
)

theme_gene <- list(
    theme_ch,
    theme(
        axis.text.x        = element_text(angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
)

bar_data_strat <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    pivot_longer(cols = c(pct_carrier, pct_noncarrier),
                 names_to = "group", values_to = "pct") %>%
    mutate(
        group = recode(group,
                       pct_carrier    = "gBRCA1/2 Carrier",
                       pct_noncarrier = "Non-carrier"),
        Gene  = factor(Gene, levels = top_genes)
    )

fig_gene_strat <- ggplot(bar_data_strat, aes(x = Gene, y = pct, fill = group)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.9) +
    scale_fill_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B", "Non-carrier" = "#546E7A"),
        name = NULL
    ) +
    scale_y_gene +
    labs(x = NULL, y = "Prevalence (%)") +
    theme_gene

ggsave(file.path("ch", "figures", "fig_gene_prevalence_bar_stratified.pdf"),
       fig_gene_strat + theme(legend.position = "bottom"), width = 6, height = 5)

ggsave(file.path("ch", "figures", "fig_gene_prevalence_bar_stratified_noleg.pdf"),
       fig_gene_strat + theme(legend.position = "none"), width = 6, height = 5)

# ============================================================
# CHIP COUNT PER INDIVIDUAL PLOT
# ============================================================
fig_chip_count <- ggplot(chip_count_dist, aes(x = factor(CHIP_Count), y = n)) +
    geom_col(fill = "#546E7A", width = 0.7, alpha = 0.9) +
    geom_text(aes(label = label), vjust = -0.4, size = 3.5, color = "grey30") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Number of CHIP mutations", y = "Number of individuals") +
    theme_ch

ggsave(file.path("ch", "figures", "fig_chip_count_dist.pdf"),
       fig_chip_count, width = 6, height = 5)


