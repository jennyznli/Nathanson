# ========================
# Packages
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "data.table"
)
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

DATE <- format(Sys.Date(), "%Y%m%d")
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

setwd(here("breast"))

### OLD VS. NEW GRAFANC OUTPUT ###
# old <- fread(here("breast", "data", "grafanc_output"))
new <- fread(here("breast", "data", "grafanc_pmbb_results.txt"))

### PLOT NEW ANCESTRIES V. HET ###
het_file <- here("simplexo", "data", "f2_breast.het")
imiss_file <- here("simplexo", "data", "f2_breast.imiss")
ancestry_file <- here("simplexo", "data", "ancestry_breast.txt")

het_data <- fread(het_file)
setnames(het_data,
         old = c("O(HOM)", "E(HOM)", "N(NM)"),
         new = c("O_HOM", "E_HOM", "N_NM"))
imiss_data <- fread(imiss_file)
ancestry_data <- fread(ancestry_file)

het_data[, obs_het_rate := (N_NM - O_HOM) / N_NM]

# # Check for any issues with calculation
# cat("Heterozygosity rate statistics:\n")
# cat("  Min:", min(het_data$obs_het_rate, na.rm = TRUE), "\n")
# cat("  Max:", max(het_data$obs_het_rate, na.rm = TRUE), "\n")
# cat("  Mean:", round(mean(het_data$obs_het_rate, na.rm = TRUE), 4), "\n")
# cat("  SD:", round(sd(het_data$obs_het_rate, na.rm = TRUE), 4), "\n")
#
# # Check missingness rate statistics
# cat("\nMissingness rate statistics:\n")
# cat("  Min:", min(imiss_data$F_MISS, na.rm = TRUE), "\n")
# cat("  Max:", max(imiss_data$F_MISS, na.rm = TRUE), "\n")
# cat("  Mean:", round(mean(imiss_data$F_MISS, na.rm = TRUE), 4), "\n")
# cat("  SD:", round(sd(imiss_data$F_MISS, na.rm = TRUE), 4), "\n\n")

# Merge all data together
data1 <- merge(
    het_data[, .(FID, IID, obs_het_rate, N_NM, O_HOM)],
    imiss_data[, .(FID, IID, F_MISS, N_MISS, N_GENO)],
    by = c("FID", "IID")
)

combined_data <- merge(
    data1,
    ancestry_data,
    by.x = "IID",
    by.y = "Sample"
)

combined_data <- combined_data[CovAnc %in% c("UNKNOWN1", "UNKNOWN2", "UNKNOWN0"), CovAnc := "UNKNOWN"]

# missing_het <- sum(is.na(combined_data$obs_het_rate))
# missing_miss <- sum(is.na(combined_data$F_MISS))
# cat("Missing heterozygosity rates:", missing_het, "\n")
# cat("Missing missingness rates:", missing_miss, "\n\n")

cat("Ancestry distribution:\n")
print(table(combined_data$CovAnc, useNA = "ifany"))

het_stats <- list(
    mean = mean(combined_data$obs_het_rate, na.rm = TRUE),
    median = median(combined_data$obs_het_rate, na.rm = TRUE),
    sd = sd(combined_data$obs_het_rate, na.rm = TRUE),
    min = min(combined_data$obs_het_rate, na.rm = TRUE),
    max = max(combined_data$obs_het_rate, na.rm = TRUE),
    q25 = quantile(combined_data$obs_het_rate, 0.25, na.rm = TRUE),
    q75 = quantile(combined_data$obs_het_rate, 0.75, na.rm = TRUE)
)

print(het_stats)

### 0.30 TO 0.38
filtered <- combined_data %>% filter(obs_het_rate > 0.30)



# Add this code after the het_stats calculation and before the Grafanc histogram

# Calculate 3 SD cutoffs
het_mean <- het_stats$mean
het_sd <- het_stats$sd
lower_cutoff <- het_mean - 3 * het_sd
upper_cutoff <- het_mean + 3 * het_sd

# Print cutoff information
cat("\n=== HETEROZYGOSITY QC CUTOFFS (Mean ± 3 SD) ===\n")
cat("Mean heterozygosity rate:", round(het_mean, 4), "\n")
cat("Standard deviation:", round(het_sd, 4), "\n")
cat("Lower cutoff (Mean - 3*SD):", round(lower_cutoff, 4), "\n")
cat("Upper cutoff (Mean + 3*SD):", round(upper_cutoff, 4), "\n")

# Mean heterozygosity rate: 0.1764
# Standard deviation: 0.0055
# Lower cutoff (Mean - 3*SD): 0.1599

# Count samples outside cutoffs
samples_below <- sum(combined_data$obs_het_rate < lower_cutoff, na.rm = TRUE)
samples_above <- sum(combined_data$obs_het_rate > upper_cutoff, na.rm = TRUE)
total_excluded <- samples_below + samples_above
total_samples <- nrow(combined_data[!is.na(obs_het_rate)])

cat("\n=== EXCLUSION SUMMARY ===\n")
cat("Samples below lower cutoff:", samples_below, "\n")
cat("Samples above upper cutoff:", samples_above, "\n")
cat("Total samples excluded:", total_excluded, "\n")
cat("Total samples with het data:", total_samples, "\n")
cat("Exclusion rate:", round(100 * total_excluded / total_samples, 2), "%\n")
# Samples below lower cutoff: 11
# > cat("Samples above upper cutoff:", samples_above, "\n")
# Samples above upper cutoff: 209
# > cat("Total samples excluded:", total_excluded, "\n")
# Total samples excluded: 220

# Create summary by ancestry
exclusion_summary <- combined_data[!is.na(obs_het_rate), .(
    total = .N,
    below_cutoff = sum(obs_het_rate < lower_cutoff),
    above_cutoff = sum(obs_het_rate > upper_cutoff),
    excluded = sum(obs_het_rate < lower_cutoff | obs_het_rate > upper_cutoff)
), by = GrafancCont][, exclusion_rate := round(100 * excluded / total, 2)]

cat("\n=== EXCLUSION BY GRAFANC ANCESTRY ===\n")
print(exclusion_summary)

# Modified Grafanc stacked histogram with cutoff lines
p_hist_grafanc_qc <- ggplot(combined_data, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    geom_vline(
        xintercept = lower_cutoff,
        color = "red",
        linetype = "solid",
        size = 1.2
    ) +
    geom_vline(
        xintercept = upper_cutoff,
        color = "red",
        linetype = "solid",
        size = 1.2
    ) +
    annotate("text",
             x = lower_cutoff,
             y = Inf,
             label = paste0("Lower\n", round(lower_cutoff, 3)),
             vjust = 1.1,
             hjust = -0.1,
             color = "red",
             size = 3.5,
             fontface = "bold") +
    annotate("text",
             x = upper_cutoff,
             y = Inf,
             label = paste0("Upper\n", round(upper_cutoff, 3)),
             vjust = 1.1,
             hjust = 1.1,
             color = "red",
             size = 3.5,
             fontface = "bold") +
    labs(
        title = "Heterozygosity Rate Distribution by Grafanc Ancestry with QC Cutoffs",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), "), Red = Mean ± 3SD cutoffs, n = ", total_excluded, " excluded (", round(100 * total_excluded / total_samples, 1), "%)"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_grafanc_qc)
ggsave(here("simplexo", "plots", "het_histogram_grafanc_qc.png"), plot = p_hist_grafanc_qc, width = 12, height = 8, dpi = 300)

# Optional: Create a list of samples to exclude for QC
samples_to_exclude <- combined_data[obs_het_rate < lower_cutoff | obs_het_rate > upper_cutoff, .(FID, IID, obs_het_rate, CovAnc, GrafancCont)]
cat("\nFirst 10 samples flagged for exclusion:\n")
print(head(samples_to_exclude, 10))

# Save exclusion list if you want to use it later
fwrite(samples_to_exclude, here("simplexo", "data", paste0("het_exclusions_3sd_", DATE, ".txt")), sep = "\t")
cat("\nExclusion list saved to: het_exclusions_3sd_", DATE, ".txt\n")





### BOX PLOTS OVERALL
p_box <- ggplot(combined_data, aes(x = 1, y = obs_het_rate)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave(here("simplexo", "plots", "het_box_overall_pmbb.png"), plot = p_box, width = 10, height = 8, dpi = 300)

### BY PMBB ANCESTRY ###

### BOX PLOTS
p_box <- ggplot(combined_data, aes(x = CovAnc, y = obs_het_rate, fill = CovAnc)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("n = ", nrow(combined_data), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave(here("simplexo", "plots", "het_box_pmbb.png"), plot = p_box, width = 10, height = 8, dpi = 300)

### BOX PLOTS
p_box <- ggplot(combined_data, aes(x = CovAnc, y = obs_het_rate, fill = CovAnc)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("n = ", nrow(combined_data), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave(here("simplexo", "plots", "het_box_pmbb.png"), plot = p_box, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist <- ggplot(combined_data, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist)
ggsave(here("simplexo", "plots", "het_histogram_pmbb.png"), plot = p_hist, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(combined_data, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave(here("simplexo", "plots", "het_histogram_facet_pmbb.png"), plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### BY GRAFANC ANCESTRY ###
### BOX PLOTS
p_box <- ggplot(combined_data, aes(x = GrafancCont, y = obs_het_rate, fill = GrafancCont)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("n = ", nrow(combined_data), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave(here("simplexo", "plots", "het_box_grafanc.png"), plot = p_box, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist <- ggplot(combined_data, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist)
ggsave(here("simplexo", "plots", "het_histogram_grafanc.png"), plot = p_hist, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(combined_data, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave(here("simplexo", "plots", "het_histogram_facet_grafanc.png"), plot = p_hist_facet, width = 10, height = 8, dpi = 300)



### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(combined_data, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB/Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        # legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_pmbb_grafanc.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(combined_data, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB/Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        # legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_grafanc_pmbb.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)





# Calculate statistics by ancestry
ancestry_stats <- combined_data %>%
    group_by(CovAnc) %>%
    summarise(
        n = n(),
        mean_het = mean(obs_het_rate, na.rm = TRUE),
        median_het = median(obs_het_rate, na.rm = TRUE),
        sd_het = sd(obs_het_rate, na.rm = TRUE),
        min_het = min(obs_het_rate, na.rm = TRUE),
        max_het = max(obs_het_rate, na.rm = TRUE),
        .groups = "drop"
    )

cat("Heterozygosity statistics by ancestry:\n")
print(ancestry_stats)
cat("\n")

# Save plots
ggsave("heterozygosity_boxplot_by_ancestry.png", plot = p_box, width = 10, height = 8, dpi = 300)
ggsave("heterozygosity_histogram_stacked.png", plot = p_hist, width = 10, height = 8, dpi = 300)
ggsave("heterozygosity_histogram_faceted.png", plot = p_hist_facet, width = 12, height = 8, dpi = 300)








### HET OF FREEZE 3 GENOTYPES ###

het3 <- fread(here("breast", "data", "all_pmbb_3.0.het"))
anc3 <- fread(here("breast", "data", "ancestry_final.txt"))

setnames(het3,
         old = c("O(HOM)", "E(HOM)", "N(NM)"),
         new = c("O_HOM", "E_HOM", "N_NM"))
het3[, obs_het_rate := (N_NM - O_HOM) / N_NM]

cat("Heterozygosity rate statistics:\n")
cat("  Min:", min(het3$obs_het_rate, na.rm = TRUE), "\n")
cat("  Max:", max(het3$obs_het_rate, na.rm = TRUE), "\n")
cat("  Mean:", round(mean(het3$obs_het_rate, na.rm = TRUE), 4), "\n")
cat("  SD:", round(sd(het3$obs_het_rate, na.rm = TRUE), 4), "\n")

comb3 <- merge(
    het3,
    anc3,
    by.x = "IID",
    by.y = "Sample"
)

comb3 <- comb3[CovAnc %in% c("UNKNOWN1", "UNKNOWN2", "UNKNOWN0"), CovAnc := "UNKNOWN"]

# missing_het <- sum(is.na(comb3$obs_het_rate))
# missing_miss <- sum(is.na(comb3$F_MISS))
# cat("Missing heterozygosity rates:", missing_het, "\n")
# cat("Missing missingness rates:", missing_miss, "\n\n")

cat("Ancestry distribution:\n")
print(table(comb3$CovAnc, useNA = "ifany"))

het_stats <- list(
    mean = mean(comb3$obs_het_rate, na.rm = TRUE),
    median = median(comb3$obs_het_rate, na.rm = TRUE),
    sd = sd(comb3$obs_het_rate, na.rm = TRUE),
    min = min(comb3$obs_het_rate, na.rm = TRUE),
    max = max(comb3$obs_het_rate, na.rm = TRUE),
    q25 = quantile(comb3$obs_het_rate, 0.25, na.rm = TRUE),
    q75 = quantile(comb3$obs_het_rate, 0.75, na.rm = TRUE)
)

### BY PMBB ANCESTRY ###
### BOX PLOTS
p_box <- ggplot(comb3, aes(x = CovAnc, y = obs_het_rate, fill = CovAnc)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("n = ", nrow(comb3), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave("het_box_pmbb3.png", plot = p_box, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist <- ggplot(comb3, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist)
ggsave("het_histogram_pmbb3.png", plot = p_hist, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb3, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_pmbb3.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### BY GRAFANC ANCESTRY ###
### BOX PLOTS
p_box <- ggplot(comb3, aes(x = GrafancCont, y = obs_het_rate, fill = GrafancCont)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("n = ", nrow(comb3), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave("het_box_grafanc_pmbb3.png", plot = p_box, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist <- ggplot(comb3, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist)
ggsave("het_histogram_grafanc_pmbb3.png", plot = p_hist, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb3, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_grafanc_pmbb3.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb3, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB/Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        # legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_pmbb3_grafanc.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb3, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB/Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        # legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_grafanc_pmbb3.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)







### HET OF FREEZE 2 GENOTYPES

het2 <- fread(here("breast", "data", "all_pmbb_2.0.het"))
anc3 <- fread(here("breast", "data", "ancestry_final.txt"))

setnames(het2,
         old = c("O(HOM)", "E(HOM)", "N(NM)"),
         new = c("O_HOM", "E_HOM", "N_NM"))
het2[, obs_het_rate := (N_NM - O_HOM) / N_NM]

cat("Heterozygosity rate statistics:\n")
cat("  Min:", min(het2$obs_het_rate, na.rm = TRUE), "\n")
cat("  Max:", max(het2$obs_het_rate, na.rm = TRUE), "\n")
cat("  Mean:", round(mean(het2$obs_het_rate, na.rm = TRUE), 4), "\n")
cat("  SD:", round(sd(het2$obs_het_rate, na.rm = TRUE), 4), "\n")

comb2 <- merge(
    het2,
    anc3,
    by.x = "IID",
    by.y = "Sample"
)

comb2 <- comb2[CovAnc %in% c("UNKNOWN1", "UNKNOWN2", "UNKNOWN0"), CovAnc := "UNKNOWN"]

# missing_het <- sum(is.na(comb2$obs_het_rate))
# missing_miss <- sum(is.na(comb2$F_MISS))
# cat("Missing heterozygosity rates:", missing_het, "\n")
# cat("Missing missingness rates:", missing_miss, "\n\n")

cat("Ancestry distribution:\n")
print(table(comb2$CovAnc, useNA = "ifany"))

het_stats <- list(
    mean = mean(comb2$obs_het_rate, na.rm = TRUE),
    median = median(comb2$obs_het_rate, na.rm = TRUE),
    sd = sd(comb2$obs_het_rate, na.rm = TRUE),
    min = min(comb2$obs_het_rate, na.rm = TRUE),
    max = max(comb2$obs_het_rate, na.rm = TRUE),
    q25 = quantile(comb2$obs_het_rate, 0.25, na.rm = TRUE),
    q75 = quantile(comb2$obs_het_rate, 0.75, na.rm = TRUE)
)

### BY PMBB ANCESTRY ###
### BOX PLOTS
p_box <- ggplot(comb2, aes(x = CovAnc, y = obs_het_rate, fill = CovAnc)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("n = ", nrow(comb2), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave("het_box_pmbb2.png", plot = p_box, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist <- ggplot(comb2, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist)
ggsave("het_histogram_pmbb2.png", plot = p_hist, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb2, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_pmbb2.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### BY GRAFANC ANCESTRY ###
### BOX PLOTS
p_box <- ggplot(comb2, aes(x = GrafancCont, y = obs_het_rate, fill = GrafancCont)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("n = ", nrow(comb2), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave("het_box_grafanc_pmbb2.png", plot = p_box, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist <- ggplot(comb2, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist)
ggsave("het_histogram_grafanc_pmbb2.png", plot = p_hist, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb2, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_grafanc_pmbb2.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)



### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb2, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB/Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        # legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_pmbb2_grafanc.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb2, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB/Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        # legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave("het_histogram_facet_grafanc_pmbb2.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)












### HET OF FREEZE 2 GENOTYPES

het2 <- fread(here("breast", "data", "all_pmbb_2.0.het"))
anc3 <- fread(here("breast", "data", "ancestry_final.txt"))

setnames(het2,
         old = c("O(HOM)", "E(HOM)", "N(NM)"),
         new = c("O_HOM", "E_HOM", "N_NM"))
het2[, obs_het_rate := (N_NM - O_HOM) / N_NM]

cat("Heterozygosity rate statistics:\n")
cat("  Min:", min(het2$obs_het_rate, na.rm = TRUE), "\n")
cat("  Max:", max(het2$obs_het_rate, na.rm = TRUE), "\n")
cat("  Mean:", round(mean(het2$obs_het_rate, na.rm = TRUE), 4), "\n")
cat("  SD:", round(sd(het2$obs_het_rate, na.rm = TRUE), 4), "\n")

comb2 <- merge(
    het2,
    anc3,
    by.x = "IID",
    by.y = "Sample"
)

comb2 <- comb2[CovAnc %in% c("UNKNOWN1", "UNKNOWN2", "UNKNOWN0"), CovAnc := "UNKNOWN"]

# missing_het <- sum(is.na(comb2$obs_het_rate))
# missing_miss <- sum(is.na(comb2$F_MISS))
# cat("Missing heterozygosity rates:", missing_het, "\n")
# cat("Missing missingness rates:", missing_miss, "\n\n")

cat("Ancestry distribution:\n")
print(table(comb2$CovAnc, useNA = "ifany"))

het_stats <- list(
    mean = mean(comb2$obs_het_rate, na.rm = TRUE),
    median = median(comb2$obs_het_rate, na.rm = TRUE),
    sd = sd(comb2$obs_het_rate, na.rm = TRUE),
    min = min(comb2$obs_het_rate, na.rm = TRUE),
    max = max(comb2$obs_het_rate, na.rm = TRUE),
    q25 = quantile(comb2$obs_het_rate, 0.25, na.rm = TRUE),
    q75 = quantile(comb2$obs_het_rate, 0.75, na.rm = TRUE)
)

### BY PMBB ANCESTRY ###
### BOX PLOTS
p_box <- ggplot(comb2, aes(x = CovAnc, y = obs_het_rate, fill = CovAnc)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("n = ", nrow(comb2), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box)
ggsave(here("breast", "plot2", "2_het_box_pmbb2.png"), plot = p_box, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist <- ggplot(comb2, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        linewidth = 1  # Fixed: 'size' is deprecated, use 'linewidth'
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist)
ggsave(here("breast", "plot2", "2_het_histogram_pmbb2.png"), plot = p_hist, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet <- ggplot(comb2, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        linewidth = 1  # Fixed: 'size' is deprecated, use 'linewidth'
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet)
ggsave(here("breast", "plot2", "2_het_histogram_facet_pmbb2.png"), plot = p_hist_facet, width = 10, height = 8, dpi = 300)

### BY GRAFANC ANCESTRY ###
### BOX PLOTS
p_box_grafanc <- ggplot(comb2, aes(x = GrafancCont, y = obs_het_rate, fill = GrafancCont)) +
    geom_boxplot(
        alpha = 0.7,
        outlier.color = "blue",
        outlier.size = 1
    ) +
    geom_jitter(
        width = 0.2,
        alpha = 0.05,
        size = 0.5
    ) +
    labs(
        title = "Distribution of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("n = ", nrow(comb2), " individuals"),
        x = "Ancestry",
        y = "Heterozygosity Rate"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p_box_grafanc)
ggsave(here("breast", "plot2", "2_het_box_grafanc_pmbb2.png"), plot = p_box_grafanc, width = 10, height = 8, dpi = 300)

### STACKED BAR PLOT
p_hist_grafanc <- ggplot(comb2, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        linewidth = 1  # Fixed: 'size' is deprecated, use 'linewidth'
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by Grafanc Ancestry",
        subtitle = paste0("Blue = Median (", round(het_stats$median, 3), ")"),
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
print(p_hist_grafanc)
ggsave(here("breast", "plot2", "2_het_histogram_grafanc_pmbb2.png"), plot = p_hist_grafanc, width = 10, height = 8, dpi = 300)

### FACETED HISTOGRAM ###
p_hist_facet_grafanc <- ggplot(comb2, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        linewidth = 1  # Fixed: 'size' is deprecated, use 'linewidth'
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet_grafanc)
ggsave(here("breast", "plot2", "2_het_histogram_facet_grafanc_pmbb2.png"), plot = p_hist_facet_grafanc, width = 10, height = 8, dpi = 300)

### COMBINED ANCESTRY PLOTS ###
### FACETED HISTOGRAM (PMBB by Grafanc) ###
p_hist_facet_combined1 <- ggplot(comb2, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ GrafancCont, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        linewidth = 1  # Fixed: 'size' is deprecated, use 'linewidth'
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by PMBB/Grafanc Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "PMBB Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet_combined1)
ggsave(here("breast", "plot2", "3_het_histogram_facet_pmbb2_grafanc.png"), plot = p_hist_facet_combined1, width = 12, height = 8, dpi = 300)

### FACETED HISTOGRAM (Grafanc by PMBB) ###
p_hist_facet_combined2 <- ggplot(comb2, aes(x = obs_het_rate, fill = GrafancCont)) +
    geom_histogram(
        bins = 30,
        color = "black",
        alpha = 0.7
    ) +
    facet_wrap(~ CovAnc, scales = "free_y") +
    geom_vline(
        xintercept = het_stats$median,
        color = "blue",
        linetype = "dashed",
        linewidth = 1  # Fixed: 'size' is deprecated, use 'linewidth'
    ) +
    labs(
        title = "Heterozygosity Rate Distribution by Grafanc/PMBB Ancestry (Facet)",
        x = "Heterozygosity Rate",
        y = "Count",
        fill = "Grafanc Ancestry"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")

print(p_hist_facet_combined2)
ggsave(here("breast", "plot2", "3_het_histogram_facet_grafanc_pmbb2.png"), plot = p_hist_facet_combined2, width = 12, height = 8, dpi = 300)




