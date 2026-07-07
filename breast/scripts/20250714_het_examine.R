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
old <- fread(here("breast", "data", "grafanc_output"))
new <- fread(here("breast", "data", "grafanc_pmbb_results.txt"))

disagreements <- old$AncGroupID != new$AncGroupID
disagreement_df <- data.table(
    Sample_ID = old[disagreements, get(names(old)[1])],
    AncGroupID_old = old$AncGroupID[disagreements],
    AncGroupID_new = new$AncGroupID[disagreements],
    Row_Number = which(disagreements)
)

if(nrow(disagreement_df) > 0) {
    disagreement_summary <- disagreement_df[, .N, by = .(AncGroupID_old, AncGroupID_new)]
    disagreement_summary <- disagreement_summary[order(-N)]

    disagreement_df <- merge(disagreement_df,
                             disagreement_summary[, .(AncGroupID_old, AncGroupID_new, N)],
                             by = c("AncGroupID_old", "AncGroupID_new"))
    disagreement_df <- disagreement_df[order(-N, Sample_ID)]

    # Remove the N column from final output (it was just for sorting)
    disagreement_df[, N := NULL]
} else {
    cat("No disagreements found!\n")
}

fwrite(disagreement_df, here("breast", "data", "ancestry_disagreements.txt"), sep = "\t")

### PLOT NEW ANCESTRIES V. HET ###
het_file <- here("breast", "data", "common_snps_breast.het")
imiss_file <- here("breast", "data", "common_snps_breast.imiss")
ancestry_file <- here("breast", "data", "ancestry_breast.txt")

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
ggsave(here("breast", "plots", "het_box_overall_pmbb.png"), plot = p_box, width = 10, height = 8, dpi = 300)

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
ggsave(here("breast", "plots", "het_box_pmbb.png"), plot = p_box, width = 10, height = 8, dpi = 300)

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
ggsave(here("breast", "plots", "het_box_pmbb.png"), plot = p_box, width = 10, height = 8, dpi = 300)

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
ggsave(here("breast", "plots", "het_histogram_pmbb.png"), plot = p_hist, width = 10, height = 8, dpi = 300)

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
ggsave("het_histogram_facet_pmbb.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)

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
ggsave("het_box_grafanc.png", plot = p_box, width = 10, height = 8, dpi = 300)

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
ggsave("het_histogram_grafanc.png", plot = p_hist, width = 10, height = 8, dpi = 300)

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
ggsave("het_histogram_facet_grafanc.png", plot = p_hist_facet, width = 10, height = 8, dpi = 300)



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




