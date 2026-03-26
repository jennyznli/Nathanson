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


bcorl1 <- vars %>% filter(Gene == "BCORL1")
table(bcorl1$MANE_SELECT)

table(bcorl1$CANONICAL)
table(bcorl1$MANE_PLUS_CLINICAL)


cat("Unique samples:", length(unique(vars$Sample.ID)), "\n")
# 3003
cat("Unique genes:",   length(unique(vars$Gene)), "\n")
# 74

genes_found <- unique(vars$Gene)
cat("Genes in list not found:", paste(gList$Gene[!(gList$Gene %in% genes_found)], collapse = ", "), "\n")
# SRCAP, ZBTB33, YLPM1 i don't think we need these?

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

    data <- data %>%
        dplyr::mutate(
            pass_protein_coding = FALSE,
            pass_vaf            = FALSE,
            pass_depth          = FALSE,
            pass_orientation    = FALSE,
            qc_pass             = FALSE
        )

    data$pass_protein_coding <- (data$Bio.type == "protein_coding") &
        (data$Sample.AltFrac  != ".") &
        (data$Sample.Depth    != ".") &
        (data$Sample.AltDepth != ".")

    data <- data %>%
        dplyr::mutate(
            Sample.AltFrac  = as.numeric(Sample.AltFrac),
            Sample.Depth    = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth),
            Sample.TLOD = as.numeric(TLOD)
        )

    data$pass_protein_coding <- data$pass_protein_coding &
        !is.na(data$Sample.AltFrac)  &
        !is.na(data$Sample.Depth)    &
        !is.na(data$Sample.AltDepth) &
        (data$Sample.AltFrac  != -1) &
        (data$Sample.Depth    != -1) &
        (data$Sample.AltDepth != -1)

    data$pass_vaf <- data$pass_protein_coding &
        (data$Sample.AltFrac >= AF_LOWER) &
        (data$Sample.AltFrac <= AF_UPPER)

    data$pass_depth <- data$pass_vaf &
        (data$Sample.Depth    >= DP) &
        (data$Sample.AltDepth >= AD)

    data$pass_orientation <- data$pass_depth &
        no_zero_counts(data$Sample.F1R2) &
        no_zero_counts(data$Sample.F2R1)

    data$qc_pass <- data$pass_orientation

    n_initial  <- nrow(data)
    summary_df <- data.frame(
        Step = c("Initial", "Pass protein coding", "Pass VAF filter",
                 "Pass depth filter", "Pass orientation filter"),
        Variants_Remaining = c(
            n_initial,
            sum(data$pass_protein_coding),
            sum(data$pass_vaf),
            sum(data$pass_depth),
            sum(data$pass_orientation)
        )
    )
    summary_df$Percent_of_Original <- round(100 * summary_df$Variants_Remaining / n_initial, 1)

    sink(log_file)
    cat("=== POST-VEP GENE QC LOG ===\n")
    cat("Dataset:", dataset_label, "\n\n")
    print(summary_df, row.names = FALSE)
    sink()

    write.csv(summary_df, summary_file, row.names = FALSE)
    cat("Log saved to:", log_file, "\n")

    return(data)
}

vars <- process_ch(vars, "ch", 0.02, 0.8, 20, 3)

all_ch <- vars %>% filter(qc_pass)
cat("QC-pass variants:", nrow(all_ch), "\n")
# 29604

write.csv(all_ch %>% select(-TLOD, -pass_protein_coding, -pass_vaf, -pass_depth, -pass_orientation, -qc_pass),
          file.path("ch", "data", "ch_seq_vars.csv"))

# ========================
# DIAGNOSTIC PLOTS
# ========================
create_diagnostic_pdf <- function(vars, steps, gene_name, output_file) {
    # steps: named character vector where names = display labels, values = filter column names
    # e.g. c("Initial" = NA, "After Depth/VAF" = "pass_depth", "Final" = "pass_orientation")

    plots <- list()
    metrics <- list(
        list(col = "Sample.Depth",    name = "Total Depth", xlab = "Total Depth (DP)",       color = "steelblue"),
        list(col = "Sample.AltDepth", name = "Alt Depth",   xlab = "Alt Depth (AD)",          color = "darkgreen"),
        list(col = "Sample.AltFrac",  name = "VAF",         xlab = "Variant Allele Fraction", color = "purple"),
        list(col = "Sample.TLOD",     name = "TLOD",        xlab = "Tumor LOD Score",         color = "darkorange")
    )

    # Build filtered data for each step by applying flags cumulatively
    step_list <- vector("list", length(steps))
    current_data <- vars
    for (s in seq_along(steps)) {
        flag_col <- steps[[s]]
        if (!is.na(flag_col)) {
            current_data <- current_data %>% filter(.data[[flag_col]])
        }
        step_list[[s]] <- list(
            data  = current_data,
            label = names(steps)[s]
        )
    }

    idx <- 1
    for (s in seq_along(step_list)) {
        step_data  <- step_list[[s]]$data
        step_label <- step_list[[s]]$label
        n          <- length(unique(step_data$Sample.ID))
        for (m in metrics) {
            if (m$col %in% names(step_data)) {
                plots[[idx]] <- ggplot(step_data, aes(x = .data[[m$col]])) +
                    geom_histogram(bins = 100, fill = m$color, alpha = 0.7, color = "white") +
                    labs(title    = paste0(step_label, ": ", m$name),
                         subtitle = paste("n =", n, "individuals"),
                         x = m$xlab, y = "Count") +
                    theme_minimal() +
                    theme(plot.title    = element_text(size = 10),
                          plot.subtitle = element_text(size = 8))
                idx <- idx + 1
            }
        }
    }

    n_steps   <- length(step_list)
    n_metrics <- length(metrics)

    # Build a flat list: [row_label, plot, plot, plot, plot,  row_label, ...]
    all_grobs <- list()
    idx <- 1
    for (s in seq_along(step_list)) {
        step_data  <- step_list[[s]]$data
        step_label <- step_list[[s]]$label
        n          <- length(unique(step_data$Sample.ID))

        # Left-side row label
        all_grobs[[idx]] <- textGrob(
            paste0(step_label, "\nn = ", n),
            rot = 90,
            gp  = gpar(fontsize = 10, fontface = "bold")
        )
        idx <- idx + 1

        for (m in metrics) {
            if (m$col %in% names(step_data)) {
                all_grobs[[idx]] <- ggplot(step_data, aes(x = .data[[m$col]])) +
                    geom_histogram(bins = 100, fill = m$color, alpha = 0.7, color = "white") +
                    labs(title = m$name, x = m$xlab, y = "Count") +
                    theme_minimal() +
                    theme(plot.title = element_text(size = 10))
                idx <- idx + 1
            }
        }
    }

    # 1 narrow label column + n_metrics plot columns
    col_widths <- c(0.3, rep(1, n_metrics))

    pdf(output_file, width = 1.5 + 5 * n_metrics, height = 4 * n_steps)
    grid.arrange(
        grobs  = all_grobs,
        ncol   = n_metrics + 1,
        nrow   = n_steps,
        widths = col_widths,
        top    = textGrob(paste(gene_name, "Filtering Pipeline"),
                          gp = gpar(fontsize = 16, fontface = "bold"))
    )
    dev.off()
}

create_diagnostic_pdf(
    vars      = vars,
    steps     = c(
        "Initial" = NA,
        "After Protein Coding" = "pass_protein_coding",
        "After VAF" = "pass_vaf",
        "After Depth" = "pass_depth",
        "After Orientation" = "pass_orientation",
        "Final" = "qc_pass"
    ),
    gene_name   = "CH",
    output_file = file.path("ch", "figures", "ch_seq_filtering.pdf")
)


