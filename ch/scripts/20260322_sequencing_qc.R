# ========================
# Sequencing based QC + plotting
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)
library(gridExtra)
library(grid)

gList <- fread(file.path("ch", "data", "whitelist_filter_20230531", "Full_CHIP_gene_list_08262022.txt"))
# gList <- fread(file.path("ch", "data", "whitelist_filter_20230531", "NEJM_2017_genes_01262020.txt"))$Gene

# ========================
# READ IN SOMATIC CALLS
# ========================
vars3 <- read.csv(file.path("ch", "data", "f3_ch.csv"))
vars2 <- read.csv(file.path("ch", "data", "f2_ch.csv"))

vars2$Batch <- 1
vars3$Batch <- 2
vars <- rbind(vars2, vars3)
write.csv(vars, file.path("ch", "data", "ch_all_vars.csv"))

cat("Unique samples:", length(unique(vars$Sample.ID)), "\n")
# 3003
cat("Unique genes:",   length(unique(vars$Gene)), "\n")
# 74

genes_found <- unique(vars$Gene)
cat("Genes in list not found:", paste(gList$Gene[!(gList$Gene %in% genes_found)], collapse = ", "), "\n")
# SRCAP, ZBTB33, YLPM1 i don't think we need these?

# ========================
# INITIAL QC STATS - ALT DEPTH
# ========================
### TOGETHER
p_alt <- ggplot(vars, aes(x = Sample.AltDepth)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    coord_cartesian(xlim = c(0, 40)) +   # <-- add this
    labs(title = "Alt Depth distribution (pre-filter)",
         x     = "Alt Depth",
         y     = "Number of variants") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title        = element_text(face = "bold", hjust = 0.5))
ggsave(file.path("ch", "figures", "qc_altdepth_overall.pdf"),
       p_alt, width = 7, height = 5)

### BY BATCH
p_alt_batch <- ggplot(vars, aes(x = Sample.AltDepth, fill = as.factor(Batch))) +
    geom_histogram(binwidth = 1, color = "white", alpha = 0.4,
                   position = "identity", linewidth = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_fill_brewer(palette = "Set2", labels = c("Freeze 2", "Freeze 3")) +
    coord_cartesian(xlim = c(0, 40)) +   # <-- add this
    labs(title = "Alt Depth distribution by batch (pre-filter)",
         x     = "Alt Depth",
         y     = "Number of variants",
         fill  = "Batch") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title        = element_text(face = "bold", hjust = 0.5),
          legend.position   = "bottom")
ggsave(file.path("ch", "figures", "qc_altdepth_by_batch.pdf"),
       p_alt_batch, width = 7, height = 5)

### SUMMARY
vars %>%
    group_by(Batch) %>%
    summarise(n_variants = n(),
              mean_alt   = mean(Sample.AltDepth, na.rm = TRUE),
              median_alt = median(Sample.AltDepth, na.rm = TRUE),
              pct_lt3    = 100 * mean(Sample.AltDepth < 3, na.rm = TRUE),  # % very low alt
              pct_lt5    = 100 * mean(Sample.AltDepth < 5, na.rm = TRUE),
              .groups    = "drop") %>%
    print()

# Batch n_variants mean_alt median_alt pct_lt3 pct_lt5
# <dbl>      <int>    <dbl>      <dbl>   <dbl>   <dbl>
#     1     1      36398     9.84          4    36.1    54.4
# 2     2     105026     6.28          2    56.3    72.7

# ========================
# SEQUENCING QC
# ========================
no_zero_counts <- function(x) {
    sapply(x, function(v) {
        if (is.na(v) || v == ".") return(FALSE)
        vals <- as.numeric(strsplit(as.character(v), ",")[[1]])
        !any(is.na(vals)) && !any(vals == 0)
    })
}

process_ch <- function(data, dataset_label, AF_LOWER, AF_UPPER, DP, AD) {
    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_file       <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, ".txt"))
    summary_file   <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, "_summary.csv"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)

    # --- Upfront protein coding filter (not a flag, just remove) ---
    n_before <- nrow(data)
    data <- data %>%
        filter(
            Bio.type == "protein_coding",
            Sample.AltFrac  != ".",
            Sample.Depth    != ".",
            Sample.AltDepth != "."
        ) %>%
        dplyr::mutate(
            Sample.AltFrac  = as.numeric(Sample.AltFrac),
            Sample.Depth    = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth),
            Sample.TLOD     = as.numeric(TLOD)
        ) %>%
        filter(
            !is.na(Sample.AltFrac),
            !is.na(Sample.Depth),
            !is.na(Sample.AltDepth),
            Sample.AltFrac  != -1,
            Sample.Depth    != -1,
            Sample.AltDepth != -1
        )
    n_after_pc <- nrow(data)
    cat("Removed", n_before - n_after_pc, "non-protein-coding / missing variants\n")

    # --- Filter flags in new order: orientation -> alt depth -> total depth -> VAF ---
    data <- data %>%
        dplyr::mutate(
            pass_orientation = FALSE,
            pass_alt_depth   = FALSE,
            pass_depth       = FALSE,
            pass_vaf         = FALSE,
            qc_pass          = FALSE
        )

    data$pass_orientation <- no_zero_counts(data$Sample.F1R2) &
        no_zero_counts(data$Sample.F2R1)

    data$pass_alt_depth <- data$pass_orientation &
        (data$Sample.AltDepth >= AD)

    data$pass_depth <- data$pass_alt_depth &
        (data$Sample.Depth >= DP)

    data$pass_vaf <- data$pass_depth &
        (data$Sample.AltFrac >= AF_LOWER) &
        (data$Sample.AltFrac <= AF_UPPER)

    data$qc_pass <- data$pass_vaf

    summary_df <- data.frame(
        Step = c("After protein coding filter", "Pass orientation filter",
                 "Pass alt depth filter", "Pass total depth filter",
                 "Pass VAF filter"),
        Variants_Remaining = c(
            n_after_pc,
            sum(data$pass_orientation),
            sum(data$pass_alt_depth),
            sum(data$pass_depth),
            sum(data$pass_vaf)
        )
    )
    summary_df$Percent_of_Original <- round(100 * summary_df$Variants_Remaining / n_before, 1)

    sink(log_file)
    cat("=== POST-VEP GENE QC LOG ===\n")
    cat("Dataset:", dataset_label, "\n\n")
    print(summary_df, row.names = FALSE)
    sink()

    write.csv(summary_df, summary_file, row.names = FALSE)
    cat("Log saved to:", log_file, "\n")

    return(data)
}

vars <- process_ch(vars, "ch", 0.02, 0.45, 20, 3)

all_ch <- vars %>% filter(qc_pass)
cat("QC-pass variants:", nrow(all_ch), "\n")
# 18231

write.csv(all_ch %>% select(-TLOD, -pass_protein_coding, -pass_vaf, -pass_depth, -pass_orientation, -qc_pass),
          file.path("ch", "data", "ch_seq_vars.csv"))

# ========================
# DIAGNOSTIC PLOTS
# ========================
create_diagnostic_pdf(
    vars      = vars,
    steps     = c(
        "Initial"             = NA,
        "After Orientation"   = "pass_orientation",
        "After Alt Depth"     = "pass_alt_depth",
        "After Total Depth"   = "pass_depth",
        "Final (After VAF)"   = "pass_vaf"
    ),
    gene_name   = "CH",
    output_file = file.path("ch", "figures", "ch_seq_filtering.pdf")
)
