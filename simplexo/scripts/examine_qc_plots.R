library(tidyverse)
# ========================
# SIMPLEXO F4
# CASE CONTROL SELECTION
# ========================
library(here)
setwd(here())
source(here("R/config.R"))

pmbb4 <- here("PMBB", "4.0")

files <- list.files("examine", pattern = "^gene_qc_metrics\\..*\\.tsv$", full.names = TRUE)
# exclude a manually-combined gene_qc_metrics.tsv if one exists alongside the per-gene files
files <- files[basename(files) != "gene_qc_metrics.tsv"]

col_names <- c("CHROM", "POS", "ID", "SAMPLE", "GT", "DP", "AD", "GQ", "GENE")

qc <- map_dfr(files, ~ read_tsv(.x, col_names = col_names, na = ".", show_col_types = FALSE))

qc <- qc %>%
  filter(GT != "./.") %>%
  separate(AD, into = c("AD_REF", "AD_ALT"), sep = ",", convert = TRUE, fill = "right") %>%
  mutate(
    DP = as.numeric(DP),
    GQ = as.numeric(GQ),
    AB = AD_ALT / (AD_REF + AD_ALT)
  )

summary_stats <- qc %>%
  group_by(GENE) %>%
  summarise(
    N = n(),
    median_DP = median(DP, na.rm = TRUE),
    median_GQ = median(GQ, na.rm = TRUE),
    median_AB = median(AB, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_DP)

print(summary_stats)
write_tsv(summary_stats, "examine/gene_qc_summary.tsv")

# Faceted histogram for one metric, panels ordered by that metric's per-gene median
plot_metric <- function(data, metric, summary_tbl) {
  order_col <- paste0("median_", metric)
  gene_order <- summary_tbl %>% arrange(.data[[order_col]]) %>% pull(GENE)
  data %>%
    mutate(GENE = factor(GENE, levels = gene_order)) %>%
    ggplot(aes(x = .data[[metric]])) +
    geom_histogram(bins = 30) +
    facet_wrap(~GENE, scales = "free_y") +
    labs(title = paste(metric, "distribution by gene (ordered by median)"), x = metric) +
    theme_minimal()
}

ggsave("examine/plot_DP_by_gene.png", plot_metric(qc, "DP", summary_stats), width = 10, height = 8)
ggsave("examine/plot_GQ_by_gene.png", plot_metric(qc, "GQ", summary_stats), width = 10, height = 8)
ggsave("examine/plot_AB_by_gene.png", plot_metric(qc, "AB", summary_stats), width = 10, height = 8)
