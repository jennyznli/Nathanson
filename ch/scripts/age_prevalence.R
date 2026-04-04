# ============================================================
# CHIP prevalence + age analysis
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(binom)

OUT_DIR <- file.path("ch", "data")
FIG_DIR <- file.path("ch", "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path(OUT_DIR, "pmbb_brca12_cov_chip_df.xlsx"))

vars <- read_excel(file.path(OUT_DIR, "ch_seq_wl_art_minad4_vars.xlsx"))

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

cat("Overall cohort: N =", nrow(cov),
    "| CHIP+ =", sum(cov$CHIP_Binary),
    sprintf("(%.1f%%)\n", 100 * mean(cov$CHIP_Binary)))
# Overall cohort: N = 3004 | CHIP+ = 193 (6.4%)

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
# QUICK SUMMARY STATS
# ============================================================
sink(file.path(OUT_DIR, "ch_overall_summary_stats.log"), split = TRUE)
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("N =", nrow(cov), "\n")

basic_crosstab <- function(df, var) {
    tab  <- table(df[[var]], useNA = "ifany")
    prop <- prop.table(tab) * 100
    out  <- data.frame(
        Level = names(tab),
        N     = as.integer(tab),
        Pct   = sprintf("%.1f%%", prop)
    )
    cat("\n", var, ":\n")
    print(out, row.names = FALSE)
    invisible(out)
}

chip_crosstab <- function(df, strat_var) {
    tab_n    <- table(df[[strat_var]], df$CHIP_Binary)
    tab_prop <- prop.table(tab_n, margin = 1) * 100
    combined <- matrix(
        sprintf("%d (%.1f%%)", tab_n, tab_prop),
        nrow = nrow(tab_n),
        dimnames = list(rownames(tab_n), c("No CHIP", "CHIP"))
    )
    total_n    <- colSums(tab_n)
    total_prop <- total_n / sum(total_n) * 100
    total_row  <- matrix(
        sprintf("%d (%.1f%%)", total_n, total_prop),
        nrow = 1,
        dimnames = list("Total", c("No CHIP", "CHIP"))
    )
    out <- rbind(combined, total_row)
    cat("\nCHIP by", strat_var, ":\n")
    print(out, quote = FALSE, right = FALSE)
    invisible(out)
}

cat("\n--- Basic counts ---\n")
for (sv in c("Sequenced_gender", "Carrier", "Batch")) {
    basic_crosstab(cov, sv)
}

cat("\nOverall CHIP prevalence:\n")
ov_n    <- table(cov$CHIP_Binary)
ov_prop <- prop.table(ov_n) * 100
cat(sprintf("  No CHIP : %d (%.1f%%)\n", ov_n["FALSE"], ov_prop["FALSE"]))
cat(sprintf("  CHIP    : %d (%.1f%%)\n", ov_n["TRUE"],  ov_prop["TRUE"]))

for (sv in c("BRCA12_Case", "Carrier", "Sequenced_gender", "Batch")) {
    chip_crosstab(cov, sv)
}

cat("\nAge summary:\n")
cat(sprintf("  Median (IQR): %.1f (%.1f-%.1f)\n",
            median(cov$Sample_age, na.rm = TRUE),
            quantile(cov$Sample_age, 0.25, na.rm = TRUE),
            quantile(cov$Sample_age, 0.75, na.rm = TRUE)))

sink()

# ============================================================
# AGE DISTRIBUTION HISTOGRAM
# ============================================================
fig_age <- ggplot(cov, aes(x = Sample_age,
                           fill = as.factor(BRCA12_Case))) +
    geom_histogram(binwidth = 2, color = "white", linewidth = 0.3,
                   alpha = 0.75, position = "identity") +
    geom_density(mapping = aes(x = Sample_age, y = after_stat(count) * 2),
                 color = "darkgray", linewidth = 0.4,
                 inherit.aes = FALSE) +
    scale_fill_manual(
        values = c("0" = "darkgray", "1" = "#C2185B"),
        labels = c("0" = "Non-carrier", "1" = "gBRCA1/2 Carrier"),
        name   = NULL
    ) +
    scale_x_continuous(breaks = seq(20, 90, by = 10), limits = c(18, 92)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
        title    = "Age distribution at sample collection",
        x = "Age at sample collection (years)",
        y = "Count"
    ) +
    theme_ch

ggsave(file.path(FIG_DIR, "fig_age_hist.pdf"), fig_age, width = 5.5, height = 4.5)

# ============================================================
# CHIP PREVALENCE BY AGE
# ============================================================
age_breaks <- c(-Inf, 30, 40, 50, 60, Inf)
age_labels <- c("≤30", "30–40", "40–50", "50–60", "≥60")

prev_overall <- cov %>%
    mutate(age_group = cut(Sample_age, breaks = age_breaks,
                           labels = age_labels, right = TRUE)) %>%
    group_by(age_group) %>%
    summarise(n_total  = n(),
              n_chip   = sum(CHIP_Binary),
              prev_pct = 100 * mean(CHIP_Binary),
              .groups  = "drop") %>%
    mutate(x_num = as.numeric(age_group))

fig_prev <- ggplot(prev_overall, aes(x = x_num, y = prev_pct)) +
    geom_smooth(method = "loess", span = 1.2, se = FALSE,
                color = "#546E7A", linewidth = 0.9) +
    geom_point(size = 3, color = "#546E7A") +
    scale_x_continuous(breaks = prev_overall$x_num,
                       labels = prev_overall$age_group) +
    scale_y_continuous(limits = c(0, NA),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0.02, 0.1))) +
    labs(
        title   = "CHIP prevalence by age decade",
        x       = "Age group (years)",
        y       = "CHIP prevalence (%)"
    ) +
    theme_ch +
    theme(plot.caption = element_text(size = 8, color = "grey50", hjust = 0))

ggsave(file.path(FIG_DIR, "fig_chip_prevalence_by_decade.pdf"),
       fig_prev, width = 5, height = 4)

# shared y-axis and x-axis scale for both figures
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
    labs(title = "CHIP prevalence by age decade",
         x = "Age group (years)", y = "CHIP prevalence (%)") +
    theme_ch

ggsave(file.path(FIG_DIR, "fig_chip_prevalence_by_decade.pdf"),
       fig_prev_overall, width = 5, height = 4)

# ── Stratified by carrier ──
fig_prev_carrier <- ggplot() +
    geom_smooth(data = prev_by_carrier,
                aes(x = x_num, y = prev_pct, color = Carrier_label),
                method = "loess", span = 1.2, se = FALSE, linewidth = 0.9) +
    geom_point(data = prev_by_carrier,
               aes(x = x_num, y = prev_pct, color = Carrier_label),
               size = 3) +
    geom_smooth(data = prev_overall,
                aes(x = x_num, y = prev_pct, group = 1),
                method = "loess", span = 1.2, se = FALSE,
                color = "darkgray", linewidth = 0.8, linetype = "dashed") +
    geom_point(data = prev_overall,
               aes(x = x_num, y = prev_pct),
               color = "darkgray", size = 3, shape = 18) +
    scale_color_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B", "Non-carrier" = "#546E7A"),
        name = NULL
    ) +
    scale_x_prev + scale_y_prev +
    labs(title = "CHIP prevalence by age decade",
         x = "Age group (years)", y = "CHIP prevalence (%)") +
    theme_ch +
    theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig_chip_prevalence_by_decade_carrier_nokey.pdf"),
       fig_prev_carrier, width = 5, height = 4)

fig_prev_carrier <- ggplot() +
    geom_smooth(data = prev_by_carrier,
                aes(x = x_num, y = prev_pct, color = Carrier_label),
                method = "loess", span = 1.2, se = FALSE, linewidth = 0.9) +
    geom_point(data = prev_by_carrier,
               aes(x = x_num, y = prev_pct, color = Carrier_label),
               size = 3) +
    geom_smooth(data = prev_overall,
                aes(x = x_num, y = prev_pct, group = 1),
                method = "loess", span = 1.2, se = FALSE,
                color = "darkgray", linewidth = 0.8, linetype = "dashed") +
    geom_point(data = prev_overall,
               aes(x = x_num, y = prev_pct),
               color = "darkgray", size = 3, shape = 18) +
    scale_color_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B", "Non-carrier" = "#546E7A"),
        name = NULL
    ) +
    scale_x_prev + scale_y_prev +
    labs(title = "CHIP prevalence by age decade",
         x = "Age group (years)", y = "CHIP prevalence (%)") +
    theme_ch

ggsave(file.path(FIG_DIR, "fig_chip_prevalence_by_decade_carrier.pdf"),
       fig_prev_carrier, width = 5, height = 4)


# ============================================================
# GENE FREQUENCY ANALYSIS
# ============================================================
vars_cov <- vars %>%
    left_join(
        cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                              BRCA2_Case, Carrier, Sample_age,
                              Batch, Smoke_History, starts_with("PC")),
        by = c("Sample.ID" = "person_id")
    ) %>%
    filter(!is.na(BRCA12_Case))

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
# 46 6
cat("\nGene frequency table (top 15):\n")
print(head(gene_freq, 15))

# ============================================================
# FISHER'S EXACT TEST PER GENE
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
# GENE PREVALENCE PLOT - STRATIFIED
# ============================================================
# shared y scale for both gene figures
scale_y_gene <- scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
)

# shared theme additions for both
theme_gene <- list(
    theme_ch,
    theme(
        axis.text.x        = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
)

# ── Overall ──
bar_data_overall <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    mutate(
        pct_total = n_total / (n_carrier_total + n_noncarrier_total) * 100,
        Gene      = factor(Gene, levels = top_genes)
    )

fig_gene_overall <- ggplot(bar_data_overall, aes(x = Gene, y = pct_total)) +
    geom_col(fill = "#546E7A", width = 0.7, alpha = 0.9) +
    scale_y_gene +
    labs(title = "CHIP gene prevalence",
         x = NULL, y = "Prevalence (%)") +
    theme_gene

ggsave(file.path(FIG_DIR, "fig_gene_prevalence_bar.pdf"),
       fig_gene_overall, width = 6, height = 5)

# ── Stratified ──
bar_data_strat <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    pivot_longer(cols = c(pct_carrier, pct_noncarrier),
                 names_to  = "group",
                 values_to = "pct") %>%
    mutate(
        group = recode(group,
                       pct_carrier    = "gBRCA1/2 Carrier",
                       pct_noncarrier = "Non-carrier"),
        Gene  = factor(Gene, levels = top_genes)
    )

fig_gene_strat <- ggplot(bar_data_strat,
                         aes(x = Gene, y = pct, fill = group)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.9) +
    scale_fill_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B",
                   "Non-carrier"       = "#546E7A"),
        name = NULL
    ) +
    scale_y_gene +
    labs(title = "CHIP gene prevalence by carrier status",
         x = NULL, y = "Prevalence (%)") +
    theme_gene +
    theme(legend.position = "bottom") +
    theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig_gene_prevalence_bar_stratified_noleg.pdf"),
       fig_gene_strat, width = 6, height = 5)


bar_data_strat <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    pivot_longer(cols = c(pct_carrier, pct_noncarrier),
                 names_to  = "group",
                 values_to = "pct") %>%
    mutate(
        group = recode(group,
                       pct_carrier    = "gBRCA1/2 Carrier",
                       pct_noncarrier = "Non-carrier"),
        Gene  = factor(Gene, levels = top_genes)
    )

fig_gene_strat <- ggplot(bar_data_strat,
                         aes(x = Gene, y = pct, fill = group)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.9) +
    scale_fill_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B",
                   "Non-carrier"       = "#546E7A"),
        name = NULL
    ) +
    scale_y_gene +
    labs(title = "CHIP gene prevalence by carrier status",
         x = NULL, y = "Prevalence (%)") +
    theme_gene +
    theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig_gene_prevalence_bar_stratified.pdf"),
       fig_gene_strat, width = 6, height = 5)


