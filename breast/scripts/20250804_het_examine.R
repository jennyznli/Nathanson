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

old <- fread(here("breast", "data", "grafanc_output"))
new <- fread(here("breast", "data", "grafanc_pmbb_results.txt"))

### PLOT ALL PMBB ANCESTRIES V. HET ###
het_file <- here("breast", "data", "common_snps_sex.het")
imiss_file <- here("breast", "data", "common_snps_sex.imiss")
ancestry_file <- here("breast", "data", "ancestry_breast.txt")

het_data <- fread(het_file)
setnames(het_data,
         old = c("O(HOM)", "E(HOM)", "N(NM)"),
         new = c("O_HOM", "E_HOM", "N_NM"))
imiss_data <- fread(imiss_file)
ancestry_data <- fread(ancestry_file)

het_data[, obs_het_rate := (N_NM - O_HOM) / N_NM]

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
ggsave(here("breast", "plots", "het_histogram_all_pmbb.png"), plot = p_hist, width = 10, height = 8, dpi = 300)

### PLOT PMBB IMPUTED COMMON SNPS - ANCESTRIES V. HET ###
het_file <- here("breast", "data", "common_snps_sex.het")
imiss_file <- here("breast", "data", "common_snps_sex.imiss")
ancestry_file <- here("breast", "data", "ancestry_breast.txt")

het_data <- fread(het_file)
setnames(het_data,
         old = c("O(HOM)", "E(HOM)", "N(NM)"),
         new = c("O_HOM", "E_HOM", "N_NM"))
imiss_data <- fread(imiss_file)
ancestry_data <- fread(ancestry_file)

het_data[, obs_het_rate := (N_NM - O_HOM) / N_NM]

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
ggsave(here("breast", "plots", "het_histogram_all_pmbb.png"), plot = p_hist, width = 10, height = 8, dpi = 300)



### PLOT PMBB GENOTYPE 2.0 - ANCESTRIES V. HET ###
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

# Calculate mean and SD
het_mean <- mean(combined_data$obs_het_rate, na.rm = TRUE)
het_sd <- sd(combined_data$obs_het_rate, na.rm = TRUE)

# Compute thresholds
lower_thresh <- het_mean - 4 * het_sd
upper_thresh <- het_mean + 4 * het_sd

cat("Mean =", het_mean, "\n")
cat("SD =", het_sd, "\n")
cat("Lower threshold =", lower_thresh, "\n")
cat("Upper threshold =", upper_thresh, "\n")

# Identify individuals outside thresholds
het_outliers <- combined_data[
    combined_data$obs_het_rate < lower_thresh |
        combined_data$obs_het_rate > upper_thresh,
    .(FID, IID, obs_het_rate, CovAnc)
]
fwrite(het_outliers, file = here("breast", "qc", "4sd_het_outliers.txt"), sep = "\t")

p_hist <- ggplot(combined_data, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_mean,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    geom_vline(
        xintercept = c(lower_thresh, upper_thresh),
        color = "red",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0(
            "Blue = Mean (", round(het_mean, 3), "), ",
            "Red = ±3 SD (", round(lower_thresh, 3), ", ", round(upper_thresh, 3), ")"
        ),
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
ggsave(here("breast", "plots", "het_histogram_f2_4sd.png"),
       plot = p_hist, width = 10, height = 8, dpi = 300)



p_hist <- ggplot(combined_data, aes(x = obs_het_rate, fill = CovAnc)) +
    geom_histogram(
        bins = 50,
        color = "black",
        alpha = 0.7,
        position = "stack"
    ) +
    geom_vline(
        xintercept = het_mean,
        color = "blue",
        linetype = "dashed",
        size = 1
    ) +
    geom_vline(
        xintercept = c(lower_thresh, upper_thresh),
        color = "red",
        linetype = "dashed",
        size = 1
    ) +
    labs(
        title = "Histogram of Observed Heterozygosity Rates by PMBB Ancestry",
        subtitle = paste0(
            "Blue = Mean (", round(het_mean, 3), "), ",
            "Red = ±3 SD (", round(lower_thresh, 3), ", ", round(upper_thresh, 3), ")"
        ),
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
    scale_fill_brewer(type = "qual", palette = "Set2") +
    coord_cartesian(ylim = c(0, 250))


print(p_hist)
ggsave(here("breast", "plots", "het_histogram_f2_4sd.png"),
       plot = p_hist, width = 10, height = 8, dpi = 300)






