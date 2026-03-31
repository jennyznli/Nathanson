# ========================
# Sequencing based QC + plotting
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

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

# reaad in whitelist
gList <- fread(file.path("ch", "data", "whitelist_filter_20230531", "Full_CHIP_gene_list_08262022.txt"))
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
    geom_histogram(binwidth = 1, color = "white", alpha = 0.85,
                   position = "identity", linewidth = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_fill_brewer(palette = "Set2", labels = c("Freeze 2", "Freeze 3")) +
    coord_cartesian(xlim = c(0, 40)) +
    facet_wrap(~ Batch, ncol = 1,
               labeller = as_labeller(c("1" = "Freeze 2", "2" = "Freeze 3"))) +
    labs(title = "Alt Depth by Freeze (pre-filter)",
         x     = "Alt Depth",
         y     = "Number of variants",
         fill  = "Batch") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title        = element_text(face = "bold", hjust = 0.5),
          legend.position   = "none")
ggsave(file.path("ch", "figures", "qc_altdepth_by_batch.pdf"),
       p_alt_batch, width = 5, height = 6)

### SUMMARY
alt_summary <- vars %>%
    group_by(Batch) %>%
    summarise(n_variants = n(),
              mean_alt   = mean(Sample.AltDepth, na.rm = TRUE),
              median_alt = median(Sample.AltDepth, na.rm = TRUE),
              pct_lt3    = 100 * mean(Sample.AltDepth < 3, na.rm = TRUE),  # % very low alt
              pct_lt5    = 100 * mean(Sample.AltDepth < 5, na.rm = TRUE),
              .groups    = "drop") %>%
    print()
write_xlsx(alt_summary, file.path("ch", "data", "ch_seq_altdepth_summary.xlsx"))

# Batch n_variants mean_alt median_alt pct_lt3 pct_lt5
# <dbl>      <int>    <dbl>      <dbl>   <dbl>   <dbl>
#     1     1      36398     9.84          4    36.1    54.4
# 2     2     105026     6.28          2    56.3    72.7

# ========================
# SEQUENCE BASED QC
# ========================
no_zero_counts <- function(x) {
    sapply(x, function(v) {
        if (is.na(v) || v == ".") return(FALSE)
        vals <- as.numeric(strsplit(as.character(v), ",")[[1]])
        !any(is.na(vals)) && !any(vals == 0)
    })
}

process_ch <- function(data, dataset_label, AF_LOWER, AF_UPPER, DP, AD, VAF_MIN = 0.02) {
    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_file       <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, ".txt"))
    summary_file   <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, "_summary.csv"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)
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
        (data$Sample.AltFrac >= VAF_MIN)
    data$qc_pass <- data$pass_vaf
    summary_df <- data.frame(
        Step = c("After protein coding filter", "Pass orientation filter",
                 "Pass alt depth filter", "Pass total depth filter",
                 paste0("Pass VAF >= ", VAF_MIN, " filter")),
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

vars <- process_ch(vars, "ch", DP = 20, AD = 3, VAF_MIN = 0.02)
all_ch <- vars %>% filter(qc_pass)
cat("QC-pass variants:", nrow(all_ch), "\n")
# 29656

write.csv(all_ch %>% select(-Sample.TLOD, -pass_alt_depth, -pass_depth, -pass_orientation, -pass_vaf, -qc_pass),
          file.path("ch", "data", "ch_seq_vars.csv"))

# ========================
# PLOTS
# ========================
create_diagnostic_pdf <- function(vars, steps, gene_name, output_file) {
    batch_colors <- setNames(RColorBrewer::brewer.pal(3, "Set2")[1:2], c("1", "2"))
    batch_labels <- c("1" = "Freeze 2", "2" = "Freeze 3")

    metrics <- list(
        list(var = "Sample.Depth",    lab = "Total Depth (DP)",        bw = 5,    xlim = NULL),
        list(var = "Sample.AltDepth", lab = "Alt Depth (AD)",          bw = 2,    xlim = NULL),
        list(var = "Sample.AltFrac",  lab = "VAF",                     bw = 0.02, xlim = c(0,1)),
        list(var = "Sample.TLOD",     lab = "Tumor LOD Score",         bw = 5,    xlim = NULL)
    )

    get_data <- function(col) if (is.na(col)) vars else vars %>% filter(.data[[col]])

    make_hist <- function(col, m) {
        d <- get_data(col) %>% mutate(Batch = as.factor(Batch))
        p <- ggplot(d, aes(.data[[m$var]], fill = Batch)) +
            geom_histogram(binwidth = m$bw, color = "white", linewidth = 0.3,
                           alpha = 0.50, position = "identity") +
            scale_fill_manual(values = batch_colors, labels = batch_labels) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            labs(x = m$lab, y = "Count") +
            theme_minimal(base_size = 9) +
            theme(panel.grid.minor = element_blank(), legend.position = "none",
                  axis.title = element_text(size = 8), axis.text = element_text(size = 7))
        if (!is.null(m$xlim)) p <- p + coord_cartesian(xlim = m$xlim)
        p
    }

    void_label <- function(top, bot = NULL) {
        p <- ggplot() + theme_void() + xlim(0,1) + ylim(0,1) +
            annotate("text", x=0.5, y=0.65, label=top, fontface="bold", size=3.2, hjust=0.5)
        if (!is.null(bot))
            p <- p + annotate("text", x=0.5, y=0.28, label=bot, size=2.8, hjust=0.5, color="grey40")
        p
    }

    step_names <- names(steps); step_cols <- unname(steps)
    n_col <- length(metrics) + 1

    # header row
    plots <- c(
        list(ggplot() + theme_void()),
        lapply(metrics, \(m) void_label(m$lab))
    )

    # one row per step
    for (i in seq_along(steps)) {
        sc <- step_cols[[i]]
        n  <- length(unique(get_data(sc)$Sample.ID))
        plots <- c(plots,
                   list(void_label(step_names[[i]], paste0("n = ", n))),
                   lapply(metrics, \(m) make_hist(sc, m))
        )
    }

    # shared legend
    legend <- cowplot::get_legend(
        ggplot(vars %>% mutate(Batch = as.factor(Batch)), aes(Sample.Depth, fill = Batch)) +
            geom_histogram(binwidth = 5, alpha = 0.75) +
            scale_fill_manual(values = batch_colors, labels = batch_labels, name = "Batch") +
            theme(legend.position = "bottom")
    )

    grid <- wrap_plots(plots,
                       nrow    = length(steps) + 1,
                       ncol    = n_col,
                       widths  = c(0.6, rep(1, length(metrics))),
                       heights = c(0.3, rep(1, length(steps)))) +
        plot_annotation(title = paste(gene_name, "Filtering Pipeline"),
                        theme = theme(plot.title = element_text(face="bold", hjust=0.5, size=14)))

    final <- cowplot::plot_grid(grid, legend, ncol=1, rel_heights=c(1, 0.04))

    ggsave(output_file, final,
           width  = 3.5 * length(metrics) + 1,
           height = 2.8 * length(steps) + 1,
           device = "pdf")
    cat("Saved to:", output_file, "\n")
}

create_diagnostic_pdf(
    vars      = vars,
    steps     = c(
        "Initial"           = NA,
        "After Orientation" = "pass_orientation",
        "After Alt Depth"   = "pass_alt_depth",
        "After Depth"       = "pass_depth",
        "Final (VAF)"       = "pass_vaf"
    ),
    gene_name   = "CH",
    output_file = file.path("ch", "figures", "ch_seq_filtering.pdf")
)

### BY GENE - ALT DEPTHS
highlight_genes <- c("ASXL1", "PPM1D", "JAK2", "SF3B1", "TP53", "SRSF2")

gene_ad <- all_ch %>%
    filter(Gene %in% gList$Gene) %>%
    mutate(
        Gene      = factor(Gene, levels = names(sort(tapply(Sample.AltDepth, Gene, median), decreasing = TRUE))),
        highlight = Gene %in% highlight_genes
    )
ad_faces <- ifelse(levels(gene_ad$Gene) %in% highlight_genes, "bold", "plain")

p_ad_gene <- ggplot(gene_ad, aes(x = Gene, y = Sample.AltDepth, fill = highlight)) +
    geom_boxplot(aes(color = highlight), alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(values  = c("FALSE" = "grey70", "TRUE" = "steelblue")) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue4")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(ylim = c(0, 75)) +
    labs(title = "Alt Depth by Gene", x = NULL, y = "Alt Depth") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", hjust = 0.5),
          axis.text.x      = element_text(angle = 45, hjust = 1, size = 7, face = ad_faces),
          legend.position  = "none")
ggsave(file.path("ch", "figures", "qc_altdepth_by_gene.pdf"), p_ad_gene, width = 14, height = 5)

### BY GENE - TOTAL DEPTHS
gene_td <- all_ch %>%
    filter(Gene %in% gList$Gene) %>%
    mutate(
        Gene      = factor(Gene, levels = names(sort(tapply(Sample.Depth, Gene, median), decreasing = TRUE))),
        highlight = Gene %in% highlight_genes
    )
td_faces <- ifelse(levels(gene_td$Gene) %in% highlight_genes, "bold", "plain")

p_d_gene <- ggplot(gene_td, aes(x = Gene, y = Sample.Depth, fill = highlight)) +
    geom_boxplot(aes(color = highlight), alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(values  = c("FALSE" = "grey70", "TRUE" = "steelblue")) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue4")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(ylim = c(0, 175)) +
    labs(title = "Total Depth by Gene", x = NULL, y = "Total Depth") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", hjust = 0.5),
          axis.text.x      = element_text(angle = 45, hjust = 1, size = 7, face = td_faces),
          legend.position  = "none")
ggsave(file.path("ch", "figures", "qc_totaldepth_by_gene.pdf"), p_d_gene, width = 14, height = 5)
